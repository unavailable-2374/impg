//! End-to-end checks for complete global homology depth. All generated PAFs
//! and sequences here are explicitly TEST fixtures, not experimental results.
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

fn fixture() -> (TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let files = [
        (
            "TEST_ab.paf",
            concat!(
            "TEST_B#0#chr1\t100\t0\t100\t+\tTEST_A#0#chr1\t100\t0\t100\t100\t100\t60\tcg:Z:100=\n",
            "TEST_B#0#chr1\t100\t0\t100\t+\tTEST_A#0#chr2\t100\t0\t100\t100\t100\t60\tcg:Z:100=\n"),
        ),
        (
            "TEST_ac.paf",
            "TEST_C#0#chr1\t100\t0\t100\t+\tTEST_A#0#chr2\t100\t0\t100\t0\t100\t60\tcg:Z:100X\n",
        ),
    ];
    let mut paths = Vec::new();
    for (name, contents) in files {
        let path = dir.path().join(name);
        fs::write(&path, contents).unwrap();
        paths.push(path);
    }
    let list = dir.path().join("TEST_alignments.list");
    fs::write(
        &list,
        paths
            .iter()
            .map(|p| format!("{}\n", p.display()))
            .collect::<String>(),
    )
    .unwrap();
    (dir, list)
}

fn run(list: &Path, prefix: &Path, args: &[&str], hook: Option<(&str, &str)>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_impg"));
    command
        .args(["depth", "--alignment-list"])
        .arg(list)
        .arg("-O")
        .arg(prefix)
        .args(args);
    command
        .env_remove("IMPG_TEST_EXIT_AFTER_COMMIT")
        .env_remove("IMPG_TEST_EXIT_AFTER_STAGE");
    if let Some((name, value)) = hook {
        command.env(name, value);
    }
    command.output().unwrap()
}
fn success(output: Output) -> String {
    assert!(
        output.status.success(),
        "TEST impg failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stderr).unwrap()
}
fn text(prefix: &Path) -> String {
    fs::read_to_string(format!("{}.depth.tsv", prefix.display())).unwrap()
}
fn body(text: &str) -> Vec<Vec<&str>> {
    text.lines()
        .filter(|s| !s.starts_with('#'))
        .map(|s| s.split('\t').collect())
        .collect()
}
fn parse_position(token: &str) -> (&str, u64, u64) {
    let (name, coordinates) = token.rsplit_once(':').unwrap();
    let (start, end) = coordinates.split_once('-').unwrap();
    (name, start.parse().unwrap(), end.parse().unwrap())
}
fn check_partition(text: &str, lengths: &[(&str, u64)]) {
    let mut positions: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
    for row in body(text) {
        let tokens: Vec<_> = row[3].split(';').collect();
        let samples: BTreeSet<_> = tokens
            .iter()
            .map(|t| parse_position(t).0.split('#').next().unwrap())
            .collect();
        assert_eq!(
            row[2].parse::<usize>().unwrap(),
            samples.len(),
            "TEST count samples, not copies"
        );
        for token in tokens {
            let (name, start, end) = parse_position(token);
            assert!(start < end);
            positions.entry(name).or_default().push((start, end));
        }
    }
    assert_eq!(positions.len(), lengths.len());
    for &(name, length) in lengths {
        let ranges = positions.get_mut(name).unwrap();
        ranges.sort_unstable();
        let mut end = 0;
        for &(start, stop) in ranges.iter() {
            assert_eq!(start, end, "TEST duplicate/gap on {name}: {ranges:?}");
            end = stop;
        }
        assert_eq!(end, length, "TEST incomplete sequence {name}");
    }
}
const LENGTHS: &[(&str, u64)] = &[
    ("TEST_A#0#chr1", 100),
    ("TEST_A#0#chr2", 100),
    ("TEST_B#0#chr1", 100),
    ("TEST_C#0#chr1", 100),
];

#[test]
fn global_depth_preserves_all_copies_and_counts_transitive_samples() {
    let (dir, list) = fixture();
    let prefix = dir.path().join("TEST_default");
    success(run(&list, &prefix, &["-t", "1"], None));
    let output = text(&prefix);
    let rows = body(&output);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].len(), 4);
    assert_eq!(rows[0][1], "100");
    assert_eq!(rows[0][2], "3");
    assert_eq!(rows[0][3].split(';').count(), 4);
    check_partition(&output, LENGTHS);

    let representatives = dir.path().join("TEST_representative");
    success(run(
        &list,
        &representatives,
        &["--copy-mode", "representative", "-t", "4"],
        None,
    ));
    let output = text(&representatives);
    let marked = body(&output);
    assert_eq!(marked[0].len(), 5);
    assert_eq!(&marked[0][..4], &rows[0][..]);
    assert_eq!(
        marked[0][4],
        "TEST_A#0#chr1:0-100;TEST_B#0#chr1:0-100;TEST_C#0#chr1:0-100"
    );
    check_partition(&output, LENGTHS);
}

#[test]
fn global_depth_is_invariant_to_threads_file_order_and_index_backend() {
    let (dir, list) = fixture();
    let baseline = dir.path().join("TEST_baseline");
    success(run(&list, &baseline, &["-t", "1"], None));
    let expected = text(&baseline);
    let reversed = dir.path().join("TEST_reversed.list");
    fs::write(
        &reversed,
        fs::read_to_string(&list)
            .unwrap()
            .lines()
            .rev()
            .map(|s| format!("{s}\n"))
            .collect::<String>(),
    )
    .unwrap();
    for (i, args) in [
        vec!["-t", "4", "--index-mode", "per-file"],
        vec!["-t", "2", "--ref-order", "TEST_A#0#chr2,TEST_A#0#chr1"],
        vec!["-t", "1", "-x", "-m", "0"],
    ]
    .iter()
    .enumerate()
    {
        let prefix = dir.path().join(format!("TEST_variant{i}"));
        success(run(&reversed, &prefix, args, None));
        assert_eq!(text(&prefix), expected);
    }
}

#[test]
fn global_depth_windows_fai_filters_reference_and_stats_are_consistent() {
    let (dir, list) = fixture();
    let fai = dir.path().join("TEST.fa.fai");
    fs::write(&fai, "TEST_Z#0#unaligned\t17\t0\t17\t18\n").unwrap();
    let fai_list = dir.path().join("TEST_fai.list");
    fs::write(&fai_list, format!("{}\n", fai.display())).unwrap();
    let prefix = dir.path().join("TEST_windows");
    success(run(
        &list,
        &prefix,
        &[
            "--fai-list",
            fai_list.to_str().unwrap(),
            "--window-size",
            "7",
            "--ref",
            "TEST_C",
            "--copy-mode",
            "representative",
        ],
        None,
    ));
    let output = text(&prefix);
    let mut lengths = LENGTHS.to_vec();
    lengths.push(("TEST_Z#0#unaligned", 17));
    check_partition(&output, &lengths);
    for row in body(&output) {
        assert!(row[1].parse::<u64>().unwrap() <= 7);
    }

    let filtered = dir.path().join("TEST_filtered");
    success(run(&list, &filtered, &["--samples", "TEST_A,TEST_C"], None));
    let output = text(&filtered);
    check_partition(
        &output,
        &LENGTHS
            .iter()
            .copied()
            .filter(|(name, _)| !name.starts_with("TEST_B"))
            .collect::<Vec<_>>(),
    );
    assert!(body(&output)
        .iter()
        .any(|row| row[2] == "1" && row[3] == "TEST_A#0#chr1:0-100"));

    let stats = dir.path().join("TEST_stats");
    success(run(&list, &stats, &["--stats"], None));
    assert!(
        fs::read_to_string(format!("{}.summary.txt", stats.display()))
            .unwrap()
            .contains("3\t100\t400")
    );
    assert_eq!(
        fs::read_to_string(format!("{}.depth3.bed", stats.display()))
            .unwrap()
            .lines()
            .count(),
        4
    );
}

#[test]
fn global_depth_approximation_is_explicit_and_accepts_missing_cigar() {
    let dir = tempfile::tempdir().unwrap();
    let paf = dir.path().join("TEST_without_cigar.paf");
    fs::write(
        &paf,
        "TEST_B#0#chr1\t20\t0\t20\t-\tTEST_A#0#chr1\t10\t0\t10\t10\t20\t60\n",
    )
    .unwrap();
    let list = dir.path().join("TEST.list");
    fs::write(&list, format!("{}\n", paf.display())).unwrap();
    let prefix = dir.path().join("TEST_approx");
    success(run(&list, &prefix, &["--coordinate-mode", "approx"], None));
    let output = text(&prefix);
    assert!(output.contains("##coordinate_mode=approx"));
    check_partition(&output, &[("TEST_A#0#chr1", 10), ("TEST_B#0#chr1", 20)]);
    let exact = dir.path().join("TEST_exact_error");
    let error = run(&list, &exact, &[], None);
    assert!(!error.status.success());
    assert!(String::from_utf8_lossy(&error.stderr).contains("CIGAR"));
    assert!(
        !Path::new(&format!("{}.depth.tsv", exact.display())).exists(),
        "TEST failed evidence must not produce downstream depth"
    );
    let conflicting = run(
        &list,
        &exact,
        &["--coordinate-mode", "approx", "--cigar-precise"],
        None,
    );
    assert!(!conflicting.status.success());
}

#[test]
fn global_depth_reuses_dumped_order_including_unmapped_test_fai_sequences() {
    let (dir, list) = fixture();
    let fai = dir.path().join("TEST_unmapped.fa.fai");
    fs::write(&fai, "TEST_Z#0#unaligned\t13\t0\t13\t14\n").unwrap();
    let fais = dir.path().join("TEST_fais.list");
    fs::write(&fais, format!("{}\n", fai.display())).unwrap();
    let order = dir.path().join("TEST_dumped_order.txt");
    let baseline = dir.path().join("TEST_order_baseline");
    success(run(
        &list,
        &baseline,
        &[
            "--fai-list",
            fais.to_str().unwrap(),
            "--dump-order",
            order.to_str().unwrap(),
        ],
        None,
    ));
    let contents = fs::read_to_string(&order).unwrap();
    assert!(contents.lines().any(|n| n == "TEST_Z#0#unaligned"));
    let csv = contents.lines().collect::<Vec<_>>().join(",");
    let prefix = dir.path().join("TEST_order_reused");
    success(run(
        &list,
        &prefix,
        &["--fai-list", fais.to_str().unwrap(), "--ref-order", &csv],
        None,
    ));
    assert_eq!(text(&prefix), text(&baseline));
    check_partition(
        &text(&prefix),
        &[
            ("TEST_A#0#chr1", 100),
            ("TEST_A#0#chr2", 100),
            ("TEST_B#0#chr1", 100),
            ("TEST_C#0#chr1", 100),
            ("TEST_Z#0#unaligned", 13),
        ],
    );
}

#[test]
fn global_depth_rejects_truncated_closure_and_inconsistent_lengths() {
    let (dir, list) = fixture();
    let prefix = dir.path().join("TEST_invalid");
    assert!(!run(&list, &prefix, &["--max-depth", "2"], None)
        .status
        .success());
    let fai = dir.path().join("TEST_bad.fai");
    fs::write(&fai, "TEST_A#0#chr1\t99\n").unwrap();
    let fai_list = dir.path().join("TEST_bad.list");
    fs::write(&fai_list, format!("{}\n", fai.display())).unwrap();
    let error = run(
        &list,
        &prefix,
        &["--fai-list", fai_list.to_str().unwrap()],
        None,
    );
    assert!(!error.status.success());
    assert!(String::from_utf8_lossy(&error.stderr).contains("length conflict"));
}

#[test]
fn global_depth_resume_after_source_commit_matches_complete_run() {
    let (dir, list) = fixture();
    let baseline = dir.path().join("TEST_complete");
    success(run(&list, &baseline, &[], None));
    let prefix = dir.path().join("TEST_resumed");
    let stopped = run(
        &list,
        &prefix,
        &["--resume", "--index-mode", "per-file"],
        Some(("IMPG_TEST_EXIT_AFTER_COMMIT", "1")),
    );
    assert_eq!(
        stopped.status.code(),
        Some(99),
        "{}",
        String::from_utf8_lossy(&stopped.stderr)
    );
    for suffix in [".depth.ckpt", ".depth.work.bin"] {
        assert!(Path::new(&format!("{}{suffix}", prefix.display())).exists());
    }
    let incompatible = run(
        &list,
        &prefix,
        &[
            "--resume",
            "--index-mode",
            "per-file",
            "--copy-mode",
            "representative",
        ],
        None,
    );
    assert!(!incompatible.status.success());
    success(run(
        &list,
        &prefix,
        &["--resume", "--index-mode", "per-file"],
        None,
    ));
    assert_eq!(text(&prefix), text(&baseline));
    for suffix in [".depth.ckpt", ".depth.work.bin"] {
        assert!(!Path::new(&format!("{}{suffix}", prefix.display())).exists());
    }
}

#[test]
fn global_depth_resumes_every_durable_stage_and_truncates_uncommitted_output() {
    use std::io::Write;
    let (dir, list) = fixture();
    let baseline = dir.path().join("TEST_stage_baseline");
    success(run(&list, &baseline, &[], None));
    for stage in [
        "sources",
        "boundaries",
        "refinement",
        "closure",
        "groups",
        "output",
    ] {
        let prefix = dir.path().join(format!("TEST_stage_{stage}"));
        let stopped = run(
            &list,
            &prefix,
            &["--resume"],
            Some(("IMPG_TEST_EXIT_AFTER_STAGE", stage)),
        );
        assert_eq!(
            stopped.status.code(),
            Some(99),
            "TEST {stage}: {}",
            String::from_utf8_lossy(&stopped.stderr)
        );
        if stage == "output" {
            fs::OpenOptions::new()
                .append(true)
                .open(format!("{}.depth.tsv", prefix.display()))
                .unwrap()
                .write_all(b"TEST_UNCOMMITTED_TRAILER\n")
                .unwrap();
        }
        success(run(&list, &prefix, &["--resume"], None));
        assert_eq!(text(&prefix), text(&baseline), "TEST resume after {stage}");
    }
    let prefix = dir.path().join("TEST_legacy_checkpoint");
    fs::write(
        format!("{}.depth.ckpt", prefix.display()),
        b"TEST_OLD_CHECKPOINT",
    )
    .unwrap();
    let error = run(&list, &prefix, &["--resume"], None);
    assert!(!error.status.success());
    assert!(String::from_utf8_lossy(&error.stderr).contains("legacy depth checkpoint"));
}

#[test]
fn global_depth_resumes_single_boundary_frontier_through_multiple_hops() {
    let dir = tempfile::tempdir().unwrap();
    let paf = dir.path().join("TEST_single_frontier.paf");
    fs::write(
        &paf,
        concat!(
            "TEST_B#0#x\t100\t0\t100\t+\tTEST_A#0#x\t100\t0\t100\t100\t100\t60\tcg:Z:100=\n",
            "TEST_C#0#x\t100\t0\t100\t+\tTEST_B#0#x\t100\t0\t100\t100\t100\t60\tcg:Z:100=\n",
            "TEST_D#0#x\t100\t0\t70\t+\tTEST_A#0#x\t100\t0\t70\t70\t70\t60\tcg:Z:70=\n"
        ),
    )
    .unwrap();
    let list = dir.path().join("TEST.list");
    fs::write(&list, format!("{}\n", paf.display())).unwrap();
    let baseline = dir.path().join("TEST_complete");
    success(run(&list, &baseline, &["-t", "4"], None));
    let prefix = dir.path().join("TEST_resumed");
    let stopped = run(
        &list,
        &prefix,
        &["--resume", "-t", "1"],
        Some(("IMPG_TEST_EXIT_AFTER_STAGE", "refinement")),
    );
    assert_eq!(
        stopped.status.code(),
        Some(99),
        "{}",
        String::from_utf8_lossy(&stopped.stderr)
    );
    success(run(&list, &prefix, &["--resume", "-t", "4"], None));
    assert_eq!(text(&prefix), text(&baseline));
    check_partition(
        &text(&prefix),
        &[
            ("TEST_A#0#x", 100),
            ("TEST_B#0#x", 100),
            ("TEST_C#0#x", 100),
            ("TEST_D#0#x", 100),
        ],
    );
}

#[test]
fn global_depth_streams_large_test_copy_group_and_resumes_exactly() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let paf = dir.path().join("TEST_large_copy_group.paf");
    let mut out = std::io::BufWriter::new(fs::File::create(&paf).unwrap());
    for i in 0..70000 {
        writeln!(
            out,
            "TEST_B#0#x\t1\t0\t1\t+\tTEST_A#0#x\t70000\t{i}\t{}\t1\t1\t60\tcg:Z:1=",
            i + 1
        )
        .unwrap();
    }
    out.flush().unwrap();
    drop(out);
    let list = dir.path().join("TEST_large.list");
    fs::write(&list, format!("{}\n", paf.display())).unwrap();
    let baseline = dir.path().join("TEST_large_complete");
    let args = [
        "--ref",
        "TEST_B",
        "--copy-mode",
        "representative",
        "-t",
        "4",
    ];
    success(run(&list, &baseline, &args, None));
    let contents = text(&baseline);
    let rows = body(&contents);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0][1], "1");
    assert_eq!(rows[0][2], "2");
    assert!(rows[0][3].starts_with("TEST_B#0#x:0-1;TEST_A#0#x:0-1;"));
    assert_eq!(rows[0][3].split(';').count(), 70001);
    assert_eq!(rows[0][4], "TEST_A#0#x:0-1;TEST_B#0#x:0-1");
    check_partition(&contents, &[("TEST_A#0#x", 70000), ("TEST_B#0#x", 1)]);
    let spilled = dir.path().join("TEST_large_spilled");
    success(run(
        &list,
        &spilled,
        &args,
        Some(("IMPG_DEPTH_MEMORY_BYTES", "4194304")),
    ));
    assert_eq!(
        text(&spilled),
        contents,
        "TEST external and memory tables must produce identical output"
    );
    let prefix = dir.path().join("TEST_large_resumed");
    let args = [
        "--ref",
        "TEST_B",
        "--copy-mode",
        "representative",
        "-t",
        "1",
        "--resume",
    ];
    let stopped = run(
        &list,
        &prefix,
        &args,
        Some(("IMPG_TEST_EXIT_AFTER_STAGE", "output")),
    );
    assert_eq!(
        stopped.status.code(),
        Some(99),
        "{}",
        String::from_utf8_lossy(&stopped.stderr)
    );
    success(run(&list, &prefix, &args, None));
    assert_eq!(text(&prefix), contents);
}

#[test]
fn global_depth_resumes_inside_reverse_oriented_windowed_group() {
    let dir = tempfile::tempdir().unwrap();
    let paf = dir.path().join("TEST_reverse.paf");
    fs::write(&paf, "TEST_B#0#chr1\t3000\t0\t3000\t-\tTEST_A#0#chr1\t3000\t0\t3000\t3000\t3000\t60\tcg:Z:3000M\n").unwrap();
    let list = dir.path().join("TEST.list");
    fs::write(&list, format!("{}\n", paf.display())).unwrap();
    let baseline = dir.path().join("TEST_reverse_complete");
    let args = ["--window-size", "1", "--ref", "TEST_B"];
    success(run(&list, &baseline, &args, None));
    let expected = text(&baseline);
    check_partition(
        &expected,
        &[("TEST_A#0#chr1", 3000), ("TEST_B#0#chr1", 3000)],
    );
    for row in body(&expected) {
        let tokens: Vec<_> = row[3].split(';').collect();
        let (b, s, e) = parse_position(tokens[0]);
        let (a, t, u) = parse_position(tokens[1]);
        assert_eq!((b, a), ("TEST_B#0#chr1", "TEST_A#0#chr1"));
        assert_eq!((t, u), (3000 - e, 3000 - s));
    }
    let prefix = dir.path().join("TEST_reverse_resumed");
    let args = ["--resume", "--window-size", "1", "--ref", "TEST_B"];
    let stopped = run(
        &list,
        &prefix,
        &args,
        Some(("IMPG_TEST_EXIT_AFTER_STAGE", "output")),
    );
    assert_eq!(stopped.status.code(), Some(99));
    assert_eq!(body(&text(&prefix)).len(), 1024);
    success(run(&list, &prefix, &args, None));
    assert_eq!(text(&prefix), expected);
}
