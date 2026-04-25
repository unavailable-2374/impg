//! Reproducer for the `impg depth` ↔ `impg query -x` sample-set divergence
//! reported against the VGP 10-sample subset.
//!
//! Gold standard: `impg query -x -m 0` on an anchor enumerates every sample
//! reachable through the implicit pangenome graph. `impg depth -r ANCHOR` must
//! report the same sample set under every supported `(transitive-mode, -m)`
//! combination — anything else means depth is either dropping or fabricating
//! transitive hits.
//!
//! All tests are `#[ignore]` because they need ~600 MB of cached PAFs + a
//! pre-built combined IMPG index that aren't checked in. Run manually:
//!
//! ```bash
//! cargo build --release --manifest-path /scratch/10779/shuocao2374/tool/impg/Cargo.toml
//! cargo test  --release --manifest-path /scratch/10779/shuocao2374/tool/impg/Cargo.toml \
//!     --test depth_chain -- --ignored --nocapture --test-threads=1
//! ```
//!
//! Override input paths via env vars (defaults wired to the VGP test bed):
//!   IMPG_DEPTH_DATADIR  → cwd to invoke impg in       (default `/scratch/10779/shuocao2374/VGP`)
//!   IMPG_DEPTH_ALIST    → `--alignment-list` argument (default `ab_test/subset10/filter_sub.txt`)
//!   IMPG_DEPTH_INDEX    → cached `-i` argument        (default `combined_5b4e05867bab4443.impg`)
//!
//! Anchor `ANCHOR` exercises a 2-hop chain
//!   GCA_016904835.1#0#CM029184.1[80805779-80806908]
//!     → GCF_017639655.2#0#NC_054055.1[59716759-59717909]
//!     → GCA_963264785.1#0#OY725422.1[71008826-71009933]
//! and `query -x -m 0` returns exactly four samples (anchor + 3 reachable).

use std::collections::BTreeSet;
use std::path::PathBuf;
use std::process::Command;

const DEFAULT_DATADIR: &str = "/scratch/10779/shuocao2374/VGP";
const DEFAULT_ALIST: &str = "ab_test/subset10/filter_sub.txt";
const DEFAULT_INDEX: &str = "combined_5b4e05867bab4443.impg";

const ANCHOR: &str = "GCA_016904835.1#0#CM029184.1:80805779-80806908";

/// Sample set returned by `impg query -r ANCHOR -x -m 0` on the VGP subset10
/// data — recorded 2026-04-25 against `dc5336a` + WIP IntervalSet BTreeMap diff.
fn expected_samples() -> BTreeSet<String> {
    [
        "GCA_016904835.1",
        "GCA_903797595.2",
        "GCA_963264785.1",
        "GCF_017639655.2",
    ]
    .into_iter()
    .map(str::to_string)
    .collect()
}

fn impg_binary() -> PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ] {
        if candidate.exists() {
            return candidate;
        }
    }
    PathBuf::from("impg")
}

fn datadir() -> PathBuf {
    PathBuf::from(std::env::var("IMPG_DEPTH_DATADIR").unwrap_or_else(|_| DEFAULT_DATADIR.into()))
}
fn alist() -> String {
    std::env::var("IMPG_DEPTH_ALIST").unwrap_or_else(|_| DEFAULT_ALIST.into())
}
fn index() -> String {
    std::env::var("IMPG_DEPTH_INDEX").unwrap_or_else(|_| DEFAULT_INDEX.into())
}

/// Panic with a readable message if the test bed is missing — these tests are
/// `#[ignore]`'d so this fires only on explicit `--ignored` runs.
fn require_data() {
    let dir = datadir();
    assert!(
        dir.is_dir(),
        "test bed not found: IMPG_DEPTH_DATADIR={} (override via env)",
        dir.display()
    );
    let alist_path = dir.join(alist());
    assert!(
        alist_path.exists(),
        "alignment list not found: {}",
        alist_path.display()
    );
    let index_path = dir.join(index());
    assert!(
        index_path.exists(),
        "cached index not found: {}",
        index_path.display()
    );
}

fn run_impg(args: &[&str]) -> String {
    let bin = impg_binary();
    let dir = datadir();
    let out = Command::new(&bin)
        .current_dir(&dir)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn {}: {e}", bin.display()));
    assert!(
        out.status.success(),
        "impg {:?} failed (exit={:?})\nstderr: {}",
        args,
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("non-utf8 stdout from impg")
}

/// Parse `impg query` BED-ish output and collect unique sample prefixes
/// (token before the first `#`).
fn samples_from_query(extra: &[&str]) -> BTreeSet<String> {
    let alist_str = alist();
    let index_str = index();
    let mut args: Vec<&str> = vec![
        "query",
        "--alignment-list",
        alist_str.as_str(),
        "-i",
        index_str.as_str(),
        "-r",
        ANCHOR,
        "-x",
    ];
    args.extend_from_slice(extra);
    let stdout = run_impg(&args);
    stdout
        .lines()
        .filter(|l| !l.is_empty())
        .filter_map(|l| l.split('\t').next())
        .filter_map(|seq| seq.split('#').next().map(str::to_string))
        .collect()
}

/// Parse `impg depth -r` table output (TSV with header `#ref_seq…depth + sample
/// columns`) and collect, **per row**, the samples with non-empty alignment
/// lists. Returns a `(reported_depth, contributing_samples)` per non-header row.
fn depth_rows(extra: &[&str]) -> Vec<(u32, BTreeSet<String>)> {
    let alist_str = alist();
    let index_str = index();
    let mut args: Vec<&str> = vec![
        "depth",
        "--alignment-list",
        alist_str.as_str(),
        "-i",
        index_str.as_str(),
        "-r",
        ANCHOR,
        "-x",
    ];
    args.extend_from_slice(extra);
    let stdout = run_impg(&args);

    let mut header_samples: Vec<String> = Vec::new();
    let mut rows: Vec<(u32, BTreeSet<String>)> = Vec::new();
    for (i, line) in stdout.lines().filter(|l| !l.is_empty()).enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        if i == 0 {
            // Header: cols 0..3 are ref_seq/ref_start/ref_end/depth, 4.. are sample names
            header_samples = cols.iter().skip(4).map(|s| s.to_string()).collect();
            continue;
        }
        assert!(
            cols.len() >= 4,
            "depth row has fewer than 4 columns: {:?}",
            cols
        );
        let depth: u32 = cols[3]
            .parse()
            .unwrap_or_else(|_| panic!("non-numeric depth column: {:?}", cols[3]));
        let mut contributing = BTreeSet::new();
        for (idx, sample) in header_samples.iter().enumerate() {
            let cell = cols.get(4 + idx).copied().unwrap_or("");
            if !cell.is_empty() {
                contributing.insert(sample.clone());
            }
        }
        rows.push((depth, contributing));
    }
    rows
}

/// Assert that every depth row reports the same sample set (and matching
/// `depth` count) as the gold-standard query.
fn assert_depth_matches_query(label: &str, depth_extra: &[&str]) {
    require_data();
    let gold = expected_samples();
    let rows = depth_rows(depth_extra);
    assert!(
        !rows.is_empty(),
        "{label}: depth produced no rows (expected at least one)"
    );

    let mut failures: Vec<String> = Vec::new();
    for (i, (depth, samples)) in rows.iter().enumerate() {
        if samples != &gold {
            let missing: Vec<&str> = gold.difference(samples).map(String::as_str).collect();
            let extra: Vec<&str> = samples.difference(&gold).map(String::as_str).collect();
            failures.push(format!(
                "row {i}: depth={depth} samples={:?} (missing={:?} extra={:?})",
                samples, missing, extra
            ));
        } else if (*depth as usize) != gold.len() {
            failures.push(format!(
                "row {i}: depth={depth} but sample set has {} entries (sets agree)",
                gold.len()
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "[{label}] {} of {} rows disagree with `query -x -m 0` gold ({:?}):\n  {}",
        failures.len(),
        rows.len(),
        gold,
        failures.join("\n  ")
    );
}

#[test]
#[ignore]
fn query_anchor_returns_expected_chain() {
    require_data();
    let got_unlimited = samples_from_query(&["-m", "0"]);
    assert_eq!(
        got_unlimited,
        expected_samples(),
        "query -x -m 0 baseline drifted"
    );
    let got_two_hop = samples_from_query(&["-m", "2"]);
    assert_eq!(
        got_two_hop,
        expected_samples(),
        "query -x -m 2 should reach all 4 samples in this 2-hop chain"
    );
}

#[test]
#[ignore]
fn depth_use_bfs_unlimited_matches_query() {
    assert_depth_matches_query("--use-BFS -m 0", &["--use-BFS", "-m", "0"]);
}

#[test]
#[ignore]
fn depth_use_bfs_default_max_depth_matches_query() {
    // Default -m is 2; this combination is the recommended production setting.
    assert_depth_matches_query("--use-BFS -m 2", &["--use-BFS", "-m", "2"]);
}

/// **Reproducer for the active P0 bug (raw-BFS over-reports at unlimited -m).**
///
/// On `dc5336a` + WIP this fails: every row reports depth=7 with 3 phantom
/// samples (GCA_009914755.4, GCF_009764565.3, GCF_028564815.1) that
/// `query -x -m 0` does not reach. Linear-interpolation BFS without a depth
/// cap drifts into unrelated chains. Once the §3 fix lands (or raw-BFS is
/// retired in favour of `--use-BFS`, plan §4 option A), this test must turn
/// green alongside the others.
#[test]
#[ignore]
fn depth_raw_bfs_unlimited_matches_query() {
    assert_depth_matches_query("raw-BFS -m 0", &["-m", "0"]);
}

#[test]
#[ignore]
fn depth_raw_bfs_default_max_depth_matches_query() {
    // Default -m=2 happens to coincide with the gold sample set on this anchor,
    // because the depth cap clips the linear-interpolation drift. Kept as a
    // regression guard so a future "fix" doesn't accidentally break the case
    // that already works.
    assert_depth_matches_query("raw-BFS -m 2", &["-m", "2"]);
}
