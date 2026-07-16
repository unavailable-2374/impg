//! Regression test for the bounded `sub_indices` cache in `MultiImpg`.
//!
//! The CIGAR-precise transitive depth path (`--cigar-precise` / `--use-BFS`)
//! keeps tree caching ON and only clears `sub_indices` at the *end* of Phase 1,
//! so at hundreds-of-thousands of per-file indices the cache grows monotonically
//! toward "every file" mid-phase and blows past `vm.max_map_count`. The fix
//! bounds residency: when a fresh `get_sub_index` miss would exceed
//! `IMPG_SUB_INDEX_CACHE_LIMIT` (or the byte budget), the oldest residents are
//! evicted one at a time (FIFO) until the new one fits — incrementally, so the
//! warm working set survives instead of being flushed wholesale.
//!
//! Cache eviction is a pure memoization detail — it only forces a reload, never
//! changes what `Impg::query` computes. This test pins that guarantee: it runs
//! the single-hop `query` AND the transitive BFS (`query_transitive_bfs`, the
//! path that actually walks many files) under an aggressive cap of 1 and under
//! no cap, and asserts the results are byte-identical (modulo order). It also
//! asserts the cap is actually enforced (residency never exceeds it).

use impg::alignment_record::AlignmentRecord;
use impg::impg::{AdjustedInterval, Impg};
use impg::impg_index::{ImpgIndex, RawAlignmentInterval};
use impg::multi_impg::MultiImpg;
use impg::seqidx::SequenceIndex;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Mutex;
use tempfile::TempDir;

/// Both tests mutate the process-global `IMPG_SUB_INDEX_CACHE_LIMIT` env var;
/// cargo runs tests in one binary on multiple threads, so serialize them.
static ENV_GUARD: Mutex<()> = Mutex::new(());

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

fn canon_raw(v: Vec<RawAlignmentInterval>) -> Vec<(i64, i64, u32, i64, i64, bool)> {
    let mut out: Vec<_> = v
        .into_iter()
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
            // Deliberately asymmetric CIGAR.  The V2 index stores the reverse
            // A/B entry at the same backing-file offset, but its CIGAR must be
            // inverted (I <-> D).  A cache key that ignores entry direction
            // can therefore make a multi-target batch reuse the wrong CIGAR;
            // a pure `280=` fixture cannot expose that bug.
            &["seqA\t300\t10\t290\t+\tseqB\t300\t5\t295\t270\t300\t60\tcg:Z:50=10I40=20D180="],
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

fn build_single(alignment_files: &[String]) -> Impg {
    let mut seq_index = SequenceIndex::new();
    let mut records_by_file = Vec::with_capacity(alignment_files.len());
    for path in alignment_files {
        let file = File::open(path).expect("open paf for combined index");
        let records =
            impg::paf::parse_paf_file(path, file, NonZeroUsize::new(1).unwrap(), &mut seq_index)
                .expect("parse paf for combined index");
        records_by_file.push((records, path.clone()));
    }
    Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true)
        .expect("build combined single index")
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
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
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

#[test]
fn multi_transitive_cigar_matches_single_index_composition() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load multi");
    let single = build_single(&alignment_files);

    let multi_tid = multi.seq_index().get_id("seqA").expect("multi seqA");
    let single_tid = single.seq_index.get_id("seqA").expect("single seqA");
    assert_eq!(multi_tid, single_tid, "fixture sequence IDs must agree");

    let multi_result = multi
        .query_transitive_bfs(
            multi_tid, 0, 300, None, 3, 1, 0, None, true, None, None, false, None,
        )
        .expect("multi transitive query");
    let single_result = single.query_transitive_bfs(
        single_tid, 0, 300, None, 3, 1, 0, None, true, None, None, false, None,
    );

    assert_eq!(
        canon(multi_result),
        canon(single_result),
        "multi-file hop>=1 CIGARs must be composed back to anchor coordinates"
    );

    let multi_dfs = multi
        .query_transitive_dfs(
            multi_tid, 0, 300, None, 3, 1, 0, None, true, None, None, false, None,
        )
        .expect("multi transitive DFS query");
    let single_dfs = single.query_transitive_dfs(
        single_tid, 0, 300, None, 3, 1, 0, None, true, None, None, false, None,
    );
    assert_eq!(
        canon(multi_dfs),
        canon(single_dfs),
        "multi-file DFS hop>=1 CIGARs must be composed back to anchor coordinates"
    );
}

/// `set_sub_index_cache_limit` is the adaptive lever the CIGAR-precise depth
/// path uses to shrink the cap to the thread count before its parallel region.
/// It must (1) only ever lower the cap, (2) actually bound residency to the new
/// cap, (3) ignore a `0` ("unbounded") request, and (4) honour an explicit
/// `IMPG_SUB_INDEX_CACHE_LIMIT` env override by refusing to change it.
#[test]
fn set_sub_index_cache_limit_only_lowers_and_bounds_residency() {
    use impg::impg_index::ImpgIndex;

    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    // Default (no env): cap = min(num_indices, 8192) = 4 here.
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load default");
    assert_eq!(multi.sub_index_cache_limit(), 4, "unexpected default cap");

    // Raising is a no-op (only lowers).
    multi.set_sub_index_cache_limit(100);
    assert_eq!(multi.sub_index_cache_limit(), 4, "cap must not be raised");

    // `0` (unbounded) is ignored — we never adaptively disable bounding.
    multi.set_sub_index_cache_limit(0);
    assert_eq!(multi.sub_index_cache_limit(), 4, "cap=0 must be ignored");

    // Lowering takes effect and is enforced: residency stays within the new cap
    // even after a transitive BFS that walks all four files.
    multi.set_sub_index_cache_limit(1);
    assert_eq!(multi.sub_index_cache_limit(), 1, "cap must be lowered to 1");
    let _ = collect_all(&multi);
    assert!(
        multi.loaded_sub_index_count() <= 1,
        "residency {} exceeded lowered cap of 1",
        multi.loaded_sub_index_count()
    );

    // An explicit env override is honoured exactly: adaptive lowering is a no-op.
    std::env::set_var("IMPG_SUB_INDEX_CACHE_LIMIT", "3");
    let pinned =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load pinned");
    assert_eq!(
        pinned.sub_index_cache_limit(),
        3,
        "explicit cap not applied"
    );
    pinned.set_sub_index_cache_limit(1);
    assert_eq!(
        pinned.sub_index_cache_limit(),
        3,
        "explicit cap must override adaptive lowering"
    );
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
}

/// `relax_sub_index_cache_count_cap` is the lever the CIGAR-precise depth path
/// pulls once a byte budget is in force, so the (VMA-motivated) slot-count cap
/// stops throttling the warm working set and the byte budget becomes the sole
/// governor. It must (1) *raise* the cap toward `num_indices` (the mirror of
/// `set_sub_index_cache_limit`'s lower-only rule), (2) let residency actually
/// grow to the relaxed cap, and (3) honour an explicit `IMPG_SUB_INDEX_CACHE_LIMIT`
/// by refusing to change it.
#[test]
fn relax_sub_index_cache_count_cap_raises_and_respects_explicit_env() {
    use impg::impg_index::ImpgIndex;

    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    // Default cap = min(num_indices, 8192) = 4 here. Lower it to 1 so the relax
    // has something to raise.
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load default");
    multi.set_sub_index_cache_limit(1);
    assert_eq!(multi.sub_index_cache_limit(), 1, "cap must be lowered to 1");

    // Relax: the cap is raised to num_indices (4), and residency may now grow to
    // hold every file after a transitive BFS that walks all four.
    multi.relax_sub_index_cache_count_cap();
    assert_eq!(
        multi.sub_index_cache_limit(),
        4,
        "relax must raise the cap to num_indices"
    );
    let _ = collect_all(&multi);
    assert!(
        multi.loaded_sub_index_count() > 1,
        "residency {} did not grow past the old cap of 1 after relax",
        multi.loaded_sub_index_count()
    );
    assert!(
        multi.loaded_sub_index_count() <= 4,
        "residency {} exceeded num_indices",
        multi.loaded_sub_index_count()
    );

    // An explicit env override is honoured exactly: relax is a no-op.
    std::env::set_var("IMPG_SUB_INDEX_CACHE_LIMIT", "2");
    let pinned =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load pinned");
    assert_eq!(
        pinned.sub_index_cache_limit(),
        2,
        "explicit cap not applied"
    );
    pinned.relax_sub_index_cache_count_cap();
    assert_eq!(
        pinned.sub_index_cache_limit(),
        2,
        "explicit cap must override relax"
    );
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
}

/// `query_overlapping_cigar_compressed` is the memory-bounding hop-0 CIGAR query:
/// it streams per file and compresses each alignment's CIGAR as it is built,
/// instead of holding every overlapping file's uncompressed CIGAR at once. With
/// an identity `compress` it must reproduce `query(.., store_cigar=true)` exactly
/// — same overlaps, same order, same CIGARs — so the depth hop-0 sweep that uses
/// it stays byte-identical to the plain-query path. This pins that equivalence.
#[test]
fn streaming_cigar_compressed_query_matches_plain_query() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    std::env::set_var("IMPG_SUB_INDEX_CACHE_LIMIT", "0");
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load multi");

    let num_targets = multi.seq_index().len() as u32;
    for tid in 0..num_targets {
        let len = multi.seq_index().get_len_from_id(tid).unwrap_or(0) as i64;
        if len == 0 {
            continue;
        }
        let plain = multi
            .query(tid, 0, len, true, None, None, false)
            .expect("query");
        // Identity compress: the result must equal the plain query verbatim.
        let streamed = multi.query_overlapping_cigar_compressed(tid, 0, len, &|ops| ops.to_vec());
        assert_eq!(
            canon(plain),
            canon(streamed),
            "streaming compressed query diverged from plain query for target {tid}"
        );

        let plain_again = multi
            .query(tid, 0, len, true, None, None, false)
            .expect("query");
        let raw = multi.query_raw_overlapping(tid, 0, len);
        let (combined, combined_raw) =
            multi.query_overlapping_cigar_compressed_with_raw(tid, 0, len, &|ops| ops.to_vec());
        assert_eq!(
            canon(plain_again),
            canon(combined),
            "combined compressed query diverged for target {tid}"
        );
        assert_eq!(
            canon_raw(raw),
            canon_raw(combined_raw),
            "combined raw extents diverged for target {tid}"
        );
    }
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
}

#[test]
fn batched_cigar_compressed_query_matches_individual_queries() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load multi");

    let mut queries = Vec::new();
    for tid in 0..multi.seq_index().len() as u32 {
        let len = multi.seq_index().get_len_from_id(tid).unwrap_or(0) as i64;
        if len > 0 {
            queries.push((tid, 0, len / 2));
            queries.push((tid, len / 3, len));
        }
    }

    let batched = multi
        .batch_query_overlapping_cigar_compressed_with_raw(&queries, &|ops| ops.to_vec())
        .expect("batch query");
    assert_eq!(batched.len(), queries.len());
    for ((tid, start, end), (batch_overlaps, batch_raw)) in
        queries.into_iter().zip(batched.into_iter())
    {
        let (single_overlaps, single_raw) =
            multi.query_overlapping_cigar_compressed_with_raw(tid, start, end, &|ops| ops.to_vec());
        assert_eq!(canon(batch_overlaps), canon(single_overlaps));
        assert_eq!(canon_raw(batch_raw), canon_raw(single_raw));
    }
}

#[test]
fn batched_cigar_query_propagates_alignment_read_failure() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load multi");

    // TEST corruption after index construction: metadata still points to the
    // original CIGAR offsets, but the backing alignment is now too short.
    std::fs::write(&alignment_files[0], b"TEST_CORRUPT\n").expect("corrupt TEST paf");
    let queries: Vec<_> = (0..multi.seq_index().len() as u32)
        .filter_map(|tid| {
            let len = multi.seq_index().get_len_from_id(tid)? as i64;
            (len > 0).then_some((tid, 0, len))
        })
        .collect();
    let err = multi
        .batch_query_overlapping_cigar_compressed_with_raw(&queries, &|ops| ops.to_vec())
        .expect_err("corrupt CIGAR backing file must abort the batch");
    assert!(
        err.to_string().contains("CIGAR"),
        "unexpected strict CIGAR error: {err}"
    );
}

#[test]
fn batched_raw_query_propagates_missing_sub_index() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load multi");
    std::fs::remove_file(&index_paths[0]).expect("remove TEST sub-index");
    let queries: Vec<_> = (0..multi.seq_index().len() as u32)
        .filter_map(|tid| {
            let len = multi.seq_index().get_len_from_id(tid)? as i64;
            (len > 0).then_some((tid, 0, len))
        })
        .collect();
    let err = multi
        .batch_query_raw_overlapping(&queries)
        .expect_err("missing sub-index must abort raw batch");
    assert!(
        err.to_string().contains("batch raw depth query failed"),
        "unexpected strict raw error: {err}"
    );
}

/// The byte-budget axis of the sub-index cache (added to bound RAM where a slot
/// *count* cap can't, because `.impg` sizes span orders of magnitude) must:
/// (1) keep resident bytes within the budget, (2) leave query results
/// byte-identical to the unbounded path, and (3) treat a single index larger
/// than the whole budget as "oversized" — never cached, but still queryable.
#[test]
fn byte_budgeted_sub_index_cache_is_byte_identical_and_bounded() {
    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    // Unbounded reference (no count cap, no byte budget).
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_LIMIT");
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_BYTES");
    let unbounded =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load unbounded");
    let reference = collect_all(&unbounded);

    // Per-file estimated resident = on-disk size × 125% (mirrors the cache
    // model's `estimated_resident_bytes` / SUB_INDEX_RESIDENT_EXPANSION_PCT).
    let max_est = index_paths
        .iter()
        .map(|p| std::fs::metadata(p).unwrap().len() * 125 / 100)
        .max()
        .unwrap();

    // Budget that fits exactly one (the largest) index but not two: loading a
    // second file forces a full flush. Bytes — not slot count — are the lever.
    std::env::set_var("IMPG_SUB_INDEX_CACHE_BYTES", max_est.to_string());
    let bounded =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load bounded");
    assert_eq!(bounded.sub_index_cache_byte_budget(), max_est);
    let got = collect_all(&bounded);
    assert!(
        bounded.sub_index_cache_bytes() <= max_est,
        "resident bytes {} exceeded budget {}",
        bounded.sub_index_cache_bytes(),
        max_est
    );
    assert_eq!(
        reference, got,
        "byte-budgeted cache diverged from unbounded"
    );

    // Oversized backstop: a budget below every file's estimate caches nothing,
    // yet every query still resolves identically (the Arc is returned uncached).
    std::env::set_var("IMPG_SUB_INDEX_CACHE_BYTES", "1");
    let oversized =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load oversized");
    let got2 = collect_all(&oversized);
    assert_eq!(
        oversized.loaded_sub_index_count(),
        0,
        "no index should cache when all exceed the byte budget"
    );
    assert_eq!(oversized.sub_index_cache_bytes(), 0);
    assert_eq!(
        reference, got2,
        "oversized-backstop path diverged from unbounded"
    );

    std::env::remove_var("IMPG_SUB_INDEX_CACHE_BYTES");
}

/// `set_sub_index_cache_byte_budget` is the adaptive lever the CIGAR-precise
/// depth path uses to shrink the cache to a RAM fraction. Like the count-cap
/// lever it must (1) only ever lower, (2) ignore a `0` ("unbounded") request,
/// and (3) honour an explicit `IMPG_SUB_INDEX_CACHE_BYTES` override (incl. its
/// `K`/`M`/`G` suffix parsing) by refusing to change it.
#[test]
fn set_sub_index_cache_byte_budget_only_lowers_and_respects_explicit_env() {
    use impg::impg_index::ImpgIndex;

    let _g = ENV_GUARD.lock().unwrap_or_else(|e| e.into_inner());
    let tmp = TempDir::new().expect("tempdir");
    let (index_paths, alignment_files) = build_multi(&tmp);

    // Default: no byte budget (0 = unbounded on this axis).
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_BYTES");
    let multi =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load default");
    assert_eq!(
        multi.sub_index_cache_byte_budget(),
        0,
        "default budget should be unbounded"
    );

    // Setting from unbounded takes effect; raising afterwards is a no-op.
    multi.set_sub_index_cache_byte_budget(1_000_000);
    assert_eq!(multi.sub_index_cache_byte_budget(), 1_000_000);
    multi.set_sub_index_cache_byte_budget(2_000_000);
    assert_eq!(
        multi.sub_index_cache_byte_budget(),
        1_000_000,
        "budget must not be raised"
    );

    // `0` (unbounded) is ignored — never adaptively disable bounding.
    multi.set_sub_index_cache_byte_budget(0);
    assert_eq!(
        multi.sub_index_cache_byte_budget(),
        1_000_000,
        "budget=0 must be ignored"
    );

    // Lowering takes effect.
    multi.set_sub_index_cache_byte_budget(500_000);
    assert_eq!(multi.sub_index_cache_byte_budget(), 500_000);

    // Explicit env override (with suffix) is honoured exactly; adaptive
    // lowering is a no-op against it.
    std::env::set_var("IMPG_SUB_INDEX_CACHE_BYTES", "4M");
    let pinned =
        MultiImpg::load_from_files(&index_paths, &alignment_files, None).expect("load pinned");
    assert_eq!(
        pinned.sub_index_cache_byte_budget(),
        4 * 1024 * 1024,
        "suffix parse / explicit budget not applied"
    );
    pinned.set_sub_index_cache_byte_budget(1000);
    assert_eq!(
        pinned.sub_index_cache_byte_budget(),
        4 * 1024 * 1024,
        "explicit budget must override adaptive lowering"
    );
    std::env::remove_var("IMPG_SUB_INDEX_CACHE_BYTES");
}
