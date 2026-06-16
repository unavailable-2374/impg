//! Regression test for the bounded `sub_indices` cache in `MultiImpg`.
//!
//! The CIGAR-precise transitive depth path (`--cigar-precise` / `--use-BFS`)
//! keeps tree caching ON and only clears `sub_indices` at the *end* of Phase 1,
//! so at hundreds-of-thousands of per-file indices the cache grows monotonically
//! toward "every file" mid-phase and blows past `vm.max_map_count`. The fix
//! bounds residency: when a fresh `get_sub_index` miss would exceed
//! `IMPG_SUB_INDEX_CACHE_LIMIT`, every slot is evicted first.
//!
//! Cache eviction is a pure memoization detail — it only forces a reload, never
//! changes what `Impg::query` computes. This test pins that guarantee: it runs
//! the single-hop `query` AND the transitive BFS (`query_transitive_bfs`, the
//! path that actually walks many files) under an aggressive cap of 1 and under
//! no cap, and asserts the results are byte-identical (modulo order). It also
//! asserts the cap is actually enforced (residency never exceeds it).

use impg::alignment_record::AlignmentRecord;
use impg::impg::{AdjustedInterval, Impg};
use impg::impg_index::ImpgIndex;
use impg::multi_impg::MultiImpg;
use impg::seqidx::SequenceIndex;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use tempfile::TempDir;

fn write_paf(path: &PathBuf, lines: &[&str]) {
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

    let records_by_file = vec![(records, paf_path.to_string())];
    let impg = Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true)
        .expect("from_multi_alignment_records");

    let out = File::create(index_path).expect("create index file");
    let mut writer = BufWriter::new(out);
    impg.serialize_with_forest_map(&mut writer)
        .expect("serialize_with_forest_map");
    writer.flush().expect("flush");
}

/// Canonicalize an AdjustedInterval list (with CIGAR) for order-independent
/// comparison. Includes the CIGAR ops so a divergence in projected alignment
/// would be caught, not just coordinate drift.
fn canon(v: Vec<AdjustedInterval>) -> Vec<(i64, i64, u32, Vec<(char, i32)>, i64, i64, u32)> {
    let mut out: Vec<_> = v
        .into_iter()
        .map(|(q, cigar, t)| {
            let ops: Vec<(char, i32)> = cigar.iter().map(|o| (o.op(), o.len())).collect();
            (
                q.first, q.last, q.metadata, ops, t.first, t.last, t.metadata,
            )
        })
        .collect();
    out.sort();
    out
}

fn build_multi(tmp: &TempDir) -> (Vec<PathBuf>, Vec<String>) {
    // A transitive chain that spans files: querying seqA must hop
    // seqA -> seqB -> seqC -> seqD, pulling in a different per-file index at
    // each hop. With a cap of 1 this forces an eviction on every new file.
    let pafs: [(&str, &[&str]); 4] = [
        (
            "aln1.paf",
            &["seqA\t300\t10\t290\t+\tseqB\t300\t5\t285\t280\t280\t60\tcg:Z:280="],
        ),
        (
            "aln2.paf",
            &["seqB\t300\t20\t280\t+\tseqC\t300\t15\t275\t260\t260\t60\tcg:Z:260="],
        ),
        (
            "aln3.paf",
            &["seqC\t300\t30\t270\t+\tseqD\t300\t25\t265\t240\t240\t60\tcg:Z:240="],
        ),
        (
            "aln4.paf",
            &["seqD\t300\t0\t150\t-\tseqA\t300\t40\t190\t150\t150\t60\tcg:Z:150="],
        ),
    ];

    let mut index_paths = Vec::new();
    let mut alignment_files = Vec::new();
    for (name, lines) in &pafs {
        let paf_path = tmp.path().join(name);
        write_paf(&paf_path, lines);
        let idx_path = tmp.path().join(format!("{}.impg", name));
        build_per_file_index(paf_path.to_str().unwrap(), idx_path.to_str().unwrap());
        index_paths.push(idx_path);
        alignment_files.push(paf_path.to_string_lossy().into_owned());
    }
    (index_paths, alignment_files)
}

/// Collect every single-hop `query` + transitive BFS result over all targets.
fn collect_all(multi: &MultiImpg) -> Vec<Vec<(i64, i64, u32, Vec<(char, i32)>, i64, i64, u32)>> {
    let num_targets = multi.seq_index().len() as u32;
    let mut all = Vec::new();
    for tid in 0..num_targets {
        let len = multi.seq_index().get_len_from_id(tid).unwrap_or(0) as i64;
        if len == 0 {
            continue;
        }
        // Single-hop projected query (with CIGAR).
        let single = multi
            .query(tid, 0, len, true, None, None, false)
            .expect("query");
        all.push(canon(single));

        // Transitive BFS (with CIGAR) — the path that walks many files and is
        // the actual target of the residency bound. max_depth 0 = unlimited.
        let trans = multi
            .query_transitive_bfs(
                tid, 0, len, None, 0, 1, 0, None, true, None, None, false, None,
            )
            .expect("query_transitive_bfs");
        all.push(canon(trans));
    }
    all
}

#[test]
fn bounded_sub_index_cache_is_byte_identical_to_unbounded() {
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    // Unbounded reference (limit = 0).
    std::env::set_var("IMPG_SUB_INDEX_CACHE_LIMIT", "0");
    let unbounded =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load unbounded");
    let reference = collect_all(&unbounded);

    // Aggressive cap (limit = 1): every new file load evicts the previous, so
    // the transitive BFS reloads sub-indices repeatedly mid-walk.
    std::env::set_var("IMPG_SUB_INDEX_CACHE_LIMIT", "1");
    let bounded =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load bounded");
    let got = collect_all(&bounded);

    // Residency must never exceed the cap of 1 after a full sweep.
    assert!(
        bounded.loaded_sub_index_count() <= 1,
        "bounded cache exceeded its cap of 1: {} resident",
        bounded.loaded_sub_index_count()
    );

    assert_eq!(
        got.len(),
        reference.len(),
        "different number of result sets between bounded and unbounded"
    );
    for (i, (a, b)) in reference.iter().zip(got.iter()).enumerate() {
        assert_eq!(
            a, b,
            "bounded sub-index cache diverged from unbounded on result set {}",
            i
        );
    }

    // Restore default for any later tests sharing this process.
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
}
