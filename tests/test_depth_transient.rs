//! Regression test for the non-transitive depth path's transient query
//! variants (`query_raw_intervals_transient` / `query_raw_overlapping_transient`).
//!
//! These exist to bound peak memory on per-file indexing at ≥10⁴ files in the
//! `impg depth` non-transitive Phase 1/2 hot paths: they must load each
//! sub-index via `load_sub_index_transient`, walk its trees, and drop both
//! without writing `MultiImpg::sub_indices` or the `Impg::trees` cache.
//!
//! This test constructs a small synthetic multi-file index (3 PAF files with
//! overlapping target sequences), runs both the cached and transient query
//! variants, and asserts:
//!
//!   1. Result parity (byte-identical, modulo order) between cached and
//!      transient for every unified target and for every `query_raw_overlapping`
//!      range exercised here.
//!   2. The shared `sub_indices` vec stays empty after only-transient usage.
//!   3. The shared `sub_indices` vec is populated after cached usage (sanity
//!      check that we are actually measuring the right thing).

use impg::alignment_record::AlignmentRecord;
use impg::impg::Impg;
use impg::impg_index::{ImpgIndex, RawAlignmentInterval};
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

/// Build a per-file Impg index from a single PAF file and serialize it to disk.
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
    let impg = Impg::from_multi_alignment_records(
        &records_by_file,
        seq_index,
        None,
        true, // bidirectional (V2)
    )
    .expect("from_multi_alignment_records");

    let out = File::create(index_path).expect("create index file");
    let mut writer = BufWriter::new(out);
    impg.serialize_with_forest_map(&mut writer)
        .expect("serialize_with_forest_map");
    writer.flush().expect("flush");
}

/// Canonicalize a raw interval list so parity comparisons are order-independent.
fn sorted(mut v: Vec<RawAlignmentInterval>) -> Vec<(i64, i64, u32, i64, i64, bool)> {
    v.sort_by_key(|r| {
        (
            r.target_start,
            r.target_end,
            r.query_id,
            r.query_start,
            r.query_end,
            r.is_reverse,
        )
    });
    v.into_iter()
        .map(|r| {
            (
                r.target_start,
                r.target_end,
                r.query_id,
                r.query_start,
                r.query_end,
                r.is_reverse,
            )
        })
        .collect()
}

/// Load a fresh MultiImpg from the given per-file index paths. Returned
/// `MultiImpg` starts with an empty `sub_indices` cache.
fn fresh_multi(index_paths: &[PathBuf], alignment_files: &[String]) -> MultiImpg {
    MultiImpg::load_from_files(index_paths, alignment_files, None)
        .expect("MultiImpg::load_from_files")
}

#[test]
fn transient_raw_queries_match_cached_and_do_not_populate_cache() {
    let tmp = TempDir::new().expect("tempdir");
    let tmp_path = tmp.path().to_path_buf();

    // Three synthetic PAF files. Sequences deliberately overlap across files so
    // that the unified forest_map has multiple TreeLocation entries per target
    // (which is the scenario the transient path must dedupe by index_idx).
    // All alignments include CIGAR (cg:Z:) — required for PAF ingestion.
    let pafs: [(&str, &[&str]); 3] = [
        (
            "aln1.paf",
            &[
                "seqA\t200\t10\t190\t+\tseqB\t200\t5\t185\t180\t180\t60\tcg:Z:180=",
                "seqC\t150\t0\t150\t+\tseqB\t200\t20\t170\t150\t150\t60\tcg:Z:150=",
            ],
        ),
        (
            "aln2.paf",
            &[
                "seqD\t100\t0\t100\t+\tseqB\t200\t40\t140\t100\t100\t60\tcg:Z:100=",
                "seqD\t100\t0\t50\t-\tseqA\t200\t30\t80\t50\t50\t60\tcg:Z:50=",
            ],
        ),
        (
            "aln3.paf",
            &[
                "seqE\t120\t0\t120\t+\tseqA\t200\t60\t180\t120\t120\t60\tcg:Z:120=",
                "seqE\t120\t10\t90\t-\tseqC\t150\t20\t100\t80\t80\t60\tcg:Z:80=",
            ],
        ),
    ];

    let mut index_paths: Vec<PathBuf> = Vec::new();
    let mut alignment_files: Vec<String> = Vec::new();
    for (name, lines) in &pafs {
        let paf_path = tmp_path.join(name);
        write_paf(&paf_path, lines);
        let idx_path = tmp_path.join(format!("{}.impg", name));
        build_per_file_index(
            paf_path.to_str().unwrap(),
            idx_path.to_str().unwrap(),
        );
        index_paths.push(idx_path);
        alignment_files.push(paf_path.to_string_lossy().into_owned());
    }

    // =====================================================================
    // Parity check: cached vs transient must return identical intervals
    // for every unified target, modulo order.
    // =====================================================================
    let multi = fresh_multi(&index_paths, &alignment_files);
    let num_targets = multi.seq_index().len() as u32;
    assert!(num_targets > 0, "unified index has no sequences");

    for tid in 0..num_targets {
        let cached = multi.query_raw_intervals(tid);
        let transient = multi.query_raw_intervals_transient(tid);
        assert_eq!(
            sorted(cached.clone()),
            sorted(transient.clone()),
            "query_raw_intervals vs _transient mismatch for target {}",
            tid
        );

        // Exercise query_raw_overlapping_transient with a range that should
        // cover every interval, a middle-range, and an out-of-range window.
        for &(start, end) in &[(0i64, 1_000i64), (50, 150), (10_000, 20_000)] {
            let cached_range = multi.query_raw_overlapping(tid, start, end);
            let transient_range = multi.query_raw_overlapping_transient(tid, start, end);
            assert_eq!(
                sorted(cached_range.clone()),
                sorted(transient_range.clone()),
                "query_raw_overlapping vs _transient mismatch for target {} range [{}, {})",
                tid,
                start,
                end
            );
        }
    }

    // =====================================================================
    // Memory-bound check: transient calls must NOT populate sub_indices.
    // Use a FRESH MultiImpg so any pollution from the parity loop above
    // does not mask the result.
    // =====================================================================
    let multi_transient_only = fresh_multi(&index_paths, &alignment_files);
    assert_eq!(
        multi_transient_only.loaded_sub_index_count(),
        0,
        "pre-condition: fresh MultiImpg should have empty sub_indices"
    );
    for tid in 0..num_targets {
        let _ = multi_transient_only.query_raw_intervals_transient(tid);
        let _ = multi_transient_only.query_raw_overlapping_transient(tid, 0, 1_000);
    }
    assert_eq!(
        multi_transient_only.loaded_sub_index_count(),
        0,
        "transient query variants must not populate MultiImpg::sub_indices \
         (otherwise Phase 1 non-transitive depth will OOM on per-file scale)"
    );

    // =====================================================================
    // Sanity check: the cached variants DO populate sub_indices. If this
    // ever stops being true, the memory-bound assertion above becomes
    // vacuously true and silently stops testing anything, so fail loudly.
    // =====================================================================
    let multi_cached_only = fresh_multi(&index_paths, &alignment_files);
    assert_eq!(multi_cached_only.loaded_sub_index_count(), 0);
    for tid in 0..num_targets {
        let _ = multi_cached_only.query_raw_intervals(tid);
    }
    assert!(
        multi_cached_only.loaded_sub_index_count() > 0,
        "cached query_raw_intervals was expected to populate sub_indices — \
         if it no longer does, the transient-parity test above has nothing \
         to guard against and must be rewritten"
    );
}
