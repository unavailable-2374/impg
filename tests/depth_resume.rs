//! End-to-end test for `impg depth --resume` (CIGAR-precise BFS, global mode).
//!
//! Strategy
//! --------
//! 1. Build a tiny synthetic dataset (3 samples × 30 kb each, two PAFs forming
//!    an A↔B↔C chain) so transitive BFS through the implicit graph reaches
//!    every sample. CIGARs are included so PAF parsing succeeds.
//! 2. Run `impg depth -x --use-BFS -O <prefix>` once without `--resume` →
//!    baseline TSV. This is the "ground truth" output.
//! 3. Run again with `--resume` plus `IMPG_TEST_EXIT_AFTER_COMMIT=1` and
//!    `IMPG_CHECKPOINT_INTERVAL_OVERRIDE=1` → impg writes one ckpt commit then
//!    immediately exits with code 99, leaving partial `<prefix>.depth.tsv`,
//!    `<prefix>.depth.ckpt`, `<prefix>.depth.work.bin`.
//! 4. Re-run with `--resume` and clean env → impg loads the ckpt, replays the
//!    work-log, finishes the remaining chunks, deletes ckpt + work.bin on
//!    success.
//! 5. Compare the post-resume TSV to the baseline (after sorting and stripping
//!    the non-deterministic `#id` column).
//!
//! The "happy path" test runs `--resume` from scratch and verifies the run
//! still produces baseline-equivalent output and cleans up its artifacts.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::Command;

use impg::alignment_record::AlignmentRecord;
use impg::impg::Impg;
use impg::seqidx::SequenceIndex;
use tempfile::TempDir;

fn impg_binary() -> PathBuf {
    if let Ok(p) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(p);
    }
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for c in [
        manifest.join("target/release/impg"),
        manifest.join("target/debug/impg"),
    ] {
        if c.exists() {
            return c;
        }
    }
    PathBuf::from("impg")
}

fn write_paf(path: &Path, lines: &[String]) {
    let mut f = File::create(path).expect("create paf");
    for line in lines {
        writeln!(f, "{}", line).expect("write paf line");
    }
}

fn build_per_file_index(paf_path: &str, index_path: &str) {
    let mut seq_index = SequenceIndex::new();
    let file = File::open(paf_path).expect("open paf");
    let records: Vec<AlignmentRecord> = impg::paf::parse_paf_file(
        paf_path,
        file,
        NonZeroUsize::new(1).unwrap(),
        &mut seq_index,
    )
    .expect("parse paf");
    assert!(!records.is_empty(), "paf produced zero records");
    let by_file = vec![(records, paf_path.to_string())];
    let impg = Impg::from_multi_alignment_records(&by_file, seq_index, None, true)
        .expect("from_multi_alignment_records");
    let out = File::create(index_path).expect("create index");
    let mut w = BufWriter::new(out);
    impg.serialize_with_forest_map(&mut w)
        .expect("serialize_with_forest_map");
    w.flush().unwrap();
}

/// Build a 3-sample synthetic graph:
///   sampleA#0#chr1 ── B→A alignment ── sampleB#0#chr1
///   sampleB#0#chr1 ── C→B alignment ── sampleC#0#chr1
///
/// Returns (alignment_files, alist_path).
fn build_test_dataset(dir: &Path) -> (Vec<String>, PathBuf) {
    let len: i64 = 30_000;
    let pafs: [(&str, Vec<String>); 2] = [
        (
            "a_b.paf",
            vec![format!(
                "sampleB#0#chr1\t{len}\t0\t{len}\t+\tsampleA#0#chr1\t{len}\t0\t{len}\t{len}\t{len}\t60\tcg:Z:{len}="
            )],
        ),
        (
            "b_c.paf",
            vec![format!(
                "sampleC#0#chr1\t{len}\t0\t{len}\t+\tsampleB#0#chr1\t{len}\t0\t{len}\t{len}\t{len}\t60\tcg:Z:{len}="
            )],
        ),
    ];
    let mut paths = Vec::new();
    for (name, lines) in &pafs {
        let paf = dir.join(name);
        write_paf(&paf, lines);
        let idx_path = dir.join(format!("{}.impg", name));
        build_per_file_index(paf.to_str().unwrap(), idx_path.to_str().unwrap());
        paths.push(paf.to_string_lossy().into_owned());
    }
    let alist = dir.join("alist.txt");
    {
        let mut f = File::create(&alist).unwrap();
        for p in &paths {
            writeln!(f, "{}", p).unwrap();
        }
    }
    (paths, alist)
}

/// Run `impg depth` and return the exit status. stderr is captured and only
/// printed on assertion failure (helps debug CI breakage).
fn run_depth(
    alist: &Path,
    prefix: &Path,
    extra_args: &[&str],
    env_overrides: &[(&str, &str)],
) -> (std::process::ExitStatus, String) {
    let bin = impg_binary();
    let mut cmd = Command::new(&bin);
    cmd.arg("depth")
        .arg("--alignment-list")
        .arg(alist)
        .arg("-O")
        .arg(prefix)
        .arg("-x")
        .arg("--use-BFS")
        .arg("--ref")
        .arg("sampleA")
        .arg("--min-transitive-len")
        .arg("100");
    for a in extra_args {
        cmd.arg(a);
    }
    // Make sure no inherited env var disturbs the run.
    cmd.env_remove("IMPG_TEST_EXIT_AFTER_COMMIT");
    cmd.env_remove("IMPG_CHECKPOINT_INTERVAL_OVERRIDE");
    for (k, v) in env_overrides {
        cmd.env(k, v);
    }
    let out = cmd.output().expect("spawn impg");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (out.status, stderr)
}

/// Read a `<prefix>.depth.tsv`, drop the header, drop the per-row `#id`
/// column (non-deterministic across runs), and return sorted lines.
fn read_normalized_tsv(path: &Path) -> Vec<String> {
    let s = fs::read_to_string(path).expect("read tsv");
    let mut lines: Vec<String> = s
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            // Strip leading id field (everything up to and including the
            // first tab). The remaining columns — length, depth, positions —
            // are deterministic.
            match l.split_once('\t') {
                Some((_, rest)) => rest.to_string(),
                None => l.to_string(),
            }
        })
        .collect();
    lines.sort();
    lines
}

#[test]
fn resume_happy_path_completes_and_cleans_up() {
    let tmp = TempDir::new().unwrap();
    let (_paths, alist) = build_test_dataset(tmp.path());

    let baseline = tmp.path().join("baseline");
    let (st, err) = run_depth(&alist, &baseline, &[], &[]);
    assert!(st.success(), "baseline depth failed:\n{err}");
    let base_tsv =
        read_normalized_tsv(Path::new(&format!("{}.depth.tsv", baseline.display())));
    assert!(!base_tsv.is_empty(), "baseline produced empty TSV");

    let resumed = tmp.path().join("resumed");
    let (st, err) = run_depth(&alist, &resumed, &["--resume"], &[]);
    assert!(st.success(), "--resume depth failed:\n{err}");
    let resumed_tsv =
        read_normalized_tsv(Path::new(&format!("{}.depth.tsv", resumed.display())));

    assert_eq!(base_tsv, resumed_tsv, "TSV differs between baseline and --resume");

    // Cleanup: ckpt and work-log artifacts must not survive a successful run.
    assert!(
        !Path::new(&format!("{}.depth.ckpt", resumed.display())).exists(),
        "ckpt should be removed on success"
    );
    assert!(
        !Path::new(&format!("{}.depth.work.bin", resumed.display())).exists(),
        "work.bin should be removed on success"
    );
}

#[test]
fn resume_after_simulated_crash_matches_baseline() {
    let tmp = TempDir::new().unwrap();
    let (_paths, alist) = build_test_dataset(tmp.path());

    // Baseline (no --resume).
    let baseline = tmp.path().join("baseline");
    let (st, err) = run_depth(&alist, &baseline, &[], &[]);
    assert!(st.success(), "baseline depth failed:\n{err}");
    let base_tsv =
        read_normalized_tsv(Path::new(&format!("{}.depth.tsv", baseline.display())));

    // First --resume run with crash injection: every commit triggers
    // IMPG_TEST_EXIT_AFTER_COMMIT, so the process exits 99 right after the
    // first ckpt is durably saved.
    let resume_prefix = tmp.path().join("crash");
    let (st, err) = run_depth(
        &alist,
        &resume_prefix,
        &["--resume"],
        &[
            ("IMPG_CHECKPOINT_INTERVAL_OVERRIDE", "1"),
            ("IMPG_TEST_EXIT_AFTER_COMMIT", "1"),
        ],
    );
    assert_eq!(
        st.code(),
        Some(99),
        "expected forced exit 99 from test hook; stderr=\n{err}",
    );
    let ckpt_path = PathBuf::from(format!("{}.depth.ckpt", resume_prefix.display()));
    let work_path = PathBuf::from(format!("{}.depth.work.bin", resume_prefix.display()));
    assert!(ckpt_path.exists(), "ckpt missing after simulated crash");
    assert!(work_path.exists(), "work.bin missing after simulated crash");

    // Second --resume run with no crash hook → must complete and produce the
    // same TSV as the baseline (modulo the non-deterministic `#id` column,
    // which `read_normalized_tsv` strips).
    let (st, err) = run_depth(&alist, &resume_prefix, &["--resume"], &[]);
    assert!(st.success(), "second --resume run failed:\n{err}");
    let final_tsv =
        read_normalized_tsv(Path::new(&format!("{}.depth.tsv", resume_prefix.display())));

    assert_eq!(base_tsv, final_tsv, "post-resume TSV differs from baseline");

    // Cleanup happened on the successful second run.
    assert!(!ckpt_path.exists(), "ckpt should be cleaned up on success");
    assert!(
        !work_path.exists(),
        "work.bin should be cleaned up on success"
    );
}
