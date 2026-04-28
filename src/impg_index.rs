// src/impg_index.rs
//! Trait abstraction for IMPG index operations.
//!
//! This trait allows both single-file `Impg` and multi-file `MultiImpg` to provide
//! the same interface to query commands, making the multi-index logic invisible
//! to callers.

use crate::impg::{AdjustedInterval, CigarOp, Impg, QueryMetadata, SortedRanges};
use crate::multi_impg::MultiImpg;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::BasicCOITree;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;

/// Raw alignment interval without CIGAR projection, for fast depth computation.
/// Contains only the coordinate metadata stored in the interval tree nodes.
#[derive(Debug, Clone)]
pub struct RawAlignmentInterval {
    pub target_start: i64,
    pub target_end: i64,
    pub query_id: u32,
    pub query_start: i64,
    pub query_end: i64,
    pub is_reverse: bool,
}

/// Trait for IMPG index operations.
///
/// Both `Impg` (single-file) and `MultiImpg` (multi-file) implement this trait,
/// allowing query commands to work transparently with either.
pub trait ImpgIndex: Send + Sync {
    /// Get a reference to the sequence index (unified for MultiImpg).
    fn seq_index(&self) -> &SequenceIndex;

    /// Query a single region without transitivity.
    fn query(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval>;

    /// Query with pre-populated CIGAR cache for efficiency.
    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval>;

    /// Populate CIGAR cache for a region.
    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    );

    /// Transitive query using depth-first search.
    fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i64,
        min_distance_between_ranges: i64,
        min_output_length: Option<i64>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval>;

    /// Transitive query using breadth-first search.
    fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i64,
        min_distance_between_ranges: i64,
        min_output_length: Option<i64>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval>;

    /// Get or load an interval tree for a target sequence.
    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>>;

    /// Get all target IDs that have interval trees (for iteration).
    fn target_ids(&self) -> Vec<u32>;

    /// Remove a cached tree from memory (for memory management).
    fn remove_cached_tree(&self, target_id: u32);

    /// Get number of targets with trees.
    fn num_targets(&self) -> usize {
        self.target_ids().len()
    }

    /// Get the sequence files (FASTA/AGC) associated with this index.
    fn sequence_files(&self) -> &[String];

    /// Get the alignment files (PAF/.1aln) associated with this index.
    /// Used for caching by paths (e.g., the depth-degrees sidecar cache).
    fn alignment_files(&self) -> &[String];

    /// Query alignments where the specified sequence is the QUERY (reverse direction).
    /// Returns: Vec of (our_q_start, our_q_end, other_t_start, other_t_end, other_seq_id)
    /// - our_q_start/our_q_end: coordinates on our sequence (appears as query in the alignment)
    /// - other_t_start/other_t_end: coordinates on the other sequence (appears as target)
    /// - other_seq_id: the other sequence's unified ID
    fn query_reverse_for_depth(&self, query_id: u32) -> Vec<(i64, i64, i64, i64, u32)>;

    /// Build a lightweight reverse index: query_id -> [target_ids that have alignments with this query]
    fn build_query_to_targets_map(&self) -> FxHashMap<u32, Vec<u32>>;

    /// Query reverse alignments using a pre-built query_to_targets map.
    fn query_reverse_for_depth_with_map(
        &self,
        query_id: u32,
        query_to_targets: &FxHashMap<u32, Vec<u32>>,
    ) -> Vec<(i64, i64, i64, i64, u32)>;

    /// Clear tree cache to free memory.
    fn clear_tree_cache(&self);

    /// Clear sub-index cache (MultiImpg only) to free memory.
    /// For single Impg, this is a no-op.
    /// Useful for depth computation with many alignment files to bound peak memory.
    fn clear_sub_index_cache(&self);

    /// Clear the transient per-file header cache used by `MultiImpg`'s
    /// chunked Phase 1/2 hot paths (`load_sub_index_transient`). For single
    /// `Impg` this is a no-op; the default implementation matches that.
    ///
    /// Call between distinct phases (e.g. after the global degree pre-scan)
    /// to release retained `Arc<Impg>` headers that are no longer needed,
    /// freeing the kernel mmap regions they hold. Independent of
    /// `clear_sub_index_cache`, which targets the BFS/transitive cache.
    fn clear_transient_header_cache(&self) {}

    /// Enable or disable tree caching.
    /// When disabled, trees loaded from disk are not stored in the cache,
    /// bounding peak memory for transitive queries.
    fn set_tree_cache_enabled(&self, enabled: bool);

    /// Check if this index was built with bidirectional mode.
    /// Bidirectional indices contain both A→B and B→A entries for each alignment,
    /// eliminating the need for a separate reverse index in depth calculations.
    fn is_bidirectional(&self) -> bool;

    /// Iterate all alignment intervals for a target sequence, returning raw coordinates
    /// without CIGAR projection. Much faster than query() for bulk operations like depth.
    /// For MultiImpg, this merges results from all sub-indices that contain the target.
    fn query_raw_intervals(&self, target_id: u32) -> Vec<RawAlignmentInterval>;

    /// Query only alignment intervals that overlap a specific range [start, end) on a target sequence.
    /// Returns raw coordinates without CIGAR projection, using coitrees range query for O(n+k) performance.
    /// Much more efficient than query_raw_intervals() when only a subset of intervals is needed (e.g., BFS).
    fn query_raw_overlapping(&self, target_id: u32, start: i64, end: i64) -> Vec<RawAlignmentInterval>;

    /// Transient variant of `query_raw_intervals` that MUST NOT populate any
    /// long-lived sub-index or tree cache. `MultiImpg` overrides this to load
    /// each sub-index via `load_sub_index_transient` and drop it immediately
    /// after walking its trees, so peak retained memory per call stays at one
    /// sub-index + one tree regardless of `forest_map` breadth. The default
    /// implementation delegates to `query_raw_intervals` (single `Impg` has
    /// bounded per-file state and does not need the transient path).
    ///
    /// Required for the depth command's non-transitive Phase 1/2 hot paths
    /// when running with `--index-mode per-file` and ≫ 10⁴ alignment files.
    ///
    /// Output order is **not** guaranteed to match `query_raw_intervals`: the
    /// `MultiImpg` override groups locations by sub-index file and iterates the
    /// group map, which for `FxHashMap` is non-deterministic. Callers that need
    /// a deterministic order must sort the result themselves (the depth command
    /// already does so via `sort_unstable_by_key`).
    fn query_raw_intervals_transient(&self, target_id: u32) -> Vec<RawAlignmentInterval> {
        self.query_raw_intervals(target_id)
    }

    /// Transient variant of `query_raw_overlapping`. Same rationale as
    /// `query_raw_intervals_transient`: `MultiImpg` overrides it to avoid
    /// writing the sub-index/tree caches. Default delegates to the cached path.
    ///
    /// Output order is not guaranteed (see `query_raw_intervals_transient`).
    fn query_raw_overlapping_transient(
        &self,
        target_id: u32,
        start: i64,
        end: i64,
    ) -> Vec<RawAlignmentInterval> {
        self.query_raw_overlapping(target_id, start, end)
    }

    /// Batch variant of `query_raw_overlapping_transient`.
    ///
    /// Given a slice of `(unified_target_id, start, end)` queries, returns a
    /// parallel Vec of results — `output[i]` holds the raw intervals for
    /// `queries[i]`. The key optimisation for `MultiImpg`: every query that
    /// touches the same alignment file is served by a single transient
    /// sub-index load, so each file is read from disk at most once per call
    /// regardless of how many queries reference it. Peak live memory is
    /// bounded to one sub-index at a time (sequential file processing).
    ///
    /// Default: calls `query_raw_overlapping_transient` individually for each
    /// query — correct for single-file `Impg` which has no file sharing to
    /// exploit.
    fn batch_query_raw_overlapping(
        &self,
        queries: &[(u32, i64, i64)],
    ) -> Vec<Vec<RawAlignmentInterval>> {
        queries
            .iter()
            .map(|&(target_id, start, end)| {
                self.query_raw_overlapping_transient(target_id, start, end)
            })
            .collect()
    }

    /// Pre-scan: for each unified target in `seq_included`, count unique OTHER samples
    /// directly aligned to it. Used by the depth command to auto-detect hub sequences.
    ///
    /// Parameters:
    /// - `seq_included`: length `num_unified`, true if the sequence participates in this run.
    /// - `seq_to_sample`: unified_seq_id -> sample_id mapping, length `num_unified`.
    ///
    /// Returns a `Vec<u16>` of length `num_unified` giving the degree (unique OTHER samples) per
    /// sequence; entries for excluded sequences are 0.
    ///
    /// The default implementation iterates targets in parallel via `query_raw_intervals`. This
    /// is fine for single-file `Impg` but explodes memory for `MultiImpg` with hundreds of
    /// thousands of per-file sub-indices, which is why `MultiImpg` overrides it with a
    /// file-parallel strategy that bounds retained sub-indices to `num_threads`.
    fn compute_sample_degrees(
        &self,
        seq_included: &[bool],
        seq_to_sample: &[u16],
    ) -> Vec<u16> {
        (0..seq_included.len() as u32)
            .into_par_iter()
            .map(|seq_id| {
                if !seq_included.get(seq_id as usize).copied().unwrap_or(false) {
                    return 0u16;
                }
                let raw_alns = self.query_raw_intervals(seq_id);
                if raw_alns.is_empty() {
                    return 0;
                }
                let self_sample = seq_to_sample.get(seq_id as usize).copied().unwrap_or(0);
                let mut samples: FxHashSet<u16> = FxHashSet::default();
                for aln in &raw_alns {
                    if !seq_included
                        .get(aln.query_id as usize)
                        .copied()
                        .unwrap_or(false)
                    {
                        continue;
                    }
                    let sample_id = seq_to_sample
                        .get(aln.query_id as usize)
                        .copied()
                        .unwrap_or(0);
                    if sample_id != self_sample {
                        samples.insert(sample_id);
                    }
                }
                samples.len().min(u16::MAX as usize) as u16
            })
            .collect()
    }
}

/// Enum wrapper that can hold either a single `Impg` or a `MultiImpg`.
///
/// This allows the CLI to work with either index type without changing
/// function signatures throughout the codebase.
pub enum ImpgWrapper {
    Single(Impg),
    Multi(MultiImpg),
}

impl ImpgWrapper {
    /// Create a wrapper from a single Impg
    pub fn from_single(impg: Impg) -> Self {
        ImpgWrapper::Single(impg)
    }

    /// Create a wrapper from a MultiImpg
    pub fn from_multi(multi: MultiImpg) -> Self {
        ImpgWrapper::Multi(multi)
    }

    /// Get a reference to the inner Impg if this is a Single variant.
    /// Returns None for Multi variant.
    pub fn as_single(&self) -> Option<&Impg> {
        match self {
            ImpgWrapper::Single(impg) => Some(impg),
            ImpgWrapper::Multi(_) => None,
        }
    }
}

impl ImpgIndex for ImpgWrapper {
    fn seq_index(&self) -> &SequenceIndex {
        match self {
            ImpgWrapper::Single(impg) => impg.seq_index(),
            ImpgWrapper::Multi(multi) => multi.seq_index(),
        }
    }

    fn query(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            ),
            ImpgWrapper::Multi(multi) => multi.query(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            ),
        }
    }

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_with_cache(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                cigar_cache,
            ),
            ImpgWrapper::Multi(multi) => multi.query_with_cache(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                cigar_cache,
            ),
        }
    }

    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
        match self {
            ImpgWrapper::Single(impg) => impg.populate_cigar_cache(
                target_id,
                range_start,
                range_end,
                min_gap_compressed_identity,
                sequence_index,
                cache,
            ),
            ImpgWrapper::Multi(multi) => multi.populate_cigar_cache(
                target_id,
                range_start,
                range_end,
                min_gap_compressed_identity,
                sequence_index,
                cache,
            ),
        }
    }

    fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i64,
        min_distance_between_ranges: i64,
        min_output_length: Option<i64>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_transitive_dfs(
                target_id,
                range_start,
                range_end,
                masked_regions,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            ),
            ImpgWrapper::Multi(multi) => multi.query_transitive_dfs(
                target_id,
                range_start,
                range_end,
                masked_regions,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            ),
        }
    }

    fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i64,
        min_distance_between_ranges: i64,
        min_output_length: Option<i64>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_transitive_bfs(
                target_id,
                range_start,
                range_end,
                masked_regions,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            ),
            ImpgWrapper::Multi(multi) => multi.query_transitive_bfs(
                target_id,
                range_start,
                range_end,
                masked_regions,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            ),
        }
    }

    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        match self {
            ImpgWrapper::Single(impg) => impg.get_or_load_tree(target_id),
            ImpgWrapper::Multi(multi) => multi.get_or_load_tree(target_id),
        }
    }

    fn target_ids(&self) -> Vec<u32> {
        match self {
            ImpgWrapper::Single(impg) => impg.target_ids(),
            ImpgWrapper::Multi(multi) => multi.target_ids(),
        }
    }

    fn remove_cached_tree(&self, target_id: u32) {
        match self {
            ImpgWrapper::Single(impg) => impg.remove_cached_tree(target_id),
            ImpgWrapper::Multi(multi) => multi.remove_cached_tree(target_id),
        }
    }

    fn num_targets(&self) -> usize {
        match self {
            ImpgWrapper::Single(impg) => impg.num_targets(),
            ImpgWrapper::Multi(multi) => multi.num_targets(),
        }
    }

    fn sequence_files(&self) -> &[String] {
        match self {
            ImpgWrapper::Single(impg) => impg.sequence_files(),
            ImpgWrapper::Multi(multi) => multi.sequence_files(),
        }
    }

    fn alignment_files(&self) -> &[String] {
        match self {
            ImpgWrapper::Single(impg) => impg.alignment_files(),
            ImpgWrapper::Multi(multi) => multi.alignment_files(),
        }
    }

    fn query_reverse_for_depth(&self, query_id: u32) -> Vec<(i64, i64, i64, i64, u32)> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_reverse_for_depth(query_id),
            ImpgWrapper::Multi(multi) => multi.query_reverse_for_depth(query_id),
        }
    }

    fn build_query_to_targets_map(&self) -> FxHashMap<u32, Vec<u32>> {
        match self {
            ImpgWrapper::Single(impg) => impg.build_query_to_targets_map(),
            ImpgWrapper::Multi(multi) => multi.build_query_to_targets_map(),
        }
    }

    fn query_reverse_for_depth_with_map(
        &self,
        query_id: u32,
        query_to_targets: &FxHashMap<u32, Vec<u32>>,
    ) -> Vec<(i64, i64, i64, i64, u32)> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_reverse_for_depth_with_map(query_id, query_to_targets),
            ImpgWrapper::Multi(multi) => multi.query_reverse_for_depth_with_map(query_id, query_to_targets),
        }
    }

    fn clear_tree_cache(&self) {
        match self {
            ImpgWrapper::Single(impg) => impg.clear_tree_cache(),
            ImpgWrapper::Multi(multi) => multi.clear_tree_cache(),
        }
    }

    fn clear_sub_index_cache(&self) {
        match self {
            ImpgWrapper::Single(impg) => impg.clear_sub_index_cache(),
            ImpgWrapper::Multi(multi) => multi.clear_sub_index_cache(),
        }
    }

    fn clear_transient_header_cache(&self) {
        // Without this override, calls on a `&dyn ImpgIndex` / `&impl ImpgIndex`
        // typed as `ImpgWrapper` would silently use the trait default no-op,
        // leaking cached `Arc<Impg>` headers across phase boundaries and
        // exhausting `vm.max_map_count` at large `--index-mode per-file`
        // alignment-file counts (depth `memory allocation of N bytes failed`
        // crash long before RSS approaches the host limit).
        match self {
            ImpgWrapper::Single(_) => {}
            ImpgWrapper::Multi(multi) => multi.clear_transient_header_cache(),
        }
    }

    fn set_tree_cache_enabled(&self, enabled: bool) {
        match self {
            ImpgWrapper::Single(impg) => impg.set_tree_cache_enabled(enabled),
            ImpgWrapper::Multi(multi) => multi.set_tree_cache_enabled(enabled),
        }
    }

    fn is_bidirectional(&self) -> bool {
        match self {
            ImpgWrapper::Single(impg) => impg.is_bidirectional(),
            ImpgWrapper::Multi(multi) => multi.is_bidirectional(),
        }
    }

    fn query_raw_intervals(&self, target_id: u32) -> Vec<RawAlignmentInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_raw_intervals(target_id),
            ImpgWrapper::Multi(multi) => multi.query_raw_intervals(target_id),
        }
    }

    fn query_raw_overlapping(&self, target_id: u32, start: i64, end: i64) -> Vec<RawAlignmentInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_raw_overlapping(target_id, start, end),
            ImpgWrapper::Multi(multi) => multi.query_raw_overlapping(target_id, start, end),
        }
    }

    fn query_raw_intervals_transient(&self, target_id: u32) -> Vec<RawAlignmentInterval> {
        match self {
            // Single Impg: no per-file sub-index cache, trait default (== cached path) is fine.
            ImpgWrapper::Single(impg) => impg.query_raw_intervals_transient(target_id),
            // MultiImpg: override that skips writing the sub-index/tree caches.
            ImpgWrapper::Multi(multi) => multi.query_raw_intervals_transient(target_id),
        }
    }

    fn query_raw_overlapping_transient(
        &self,
        target_id: u32,
        start: i64,
        end: i64,
    ) -> Vec<RawAlignmentInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_raw_overlapping_transient(target_id, start, end),
            ImpgWrapper::Multi(multi) => multi.query_raw_overlapping_transient(target_id, start, end),
        }
    }

    fn batch_query_raw_overlapping(
        &self,
        queries: &[(u32, i64, i64)],
    ) -> Vec<Vec<RawAlignmentInterval>> {
        match self {
            ImpgWrapper::Single(impg) => impg.batch_query_raw_overlapping(queries),
            ImpgWrapper::Multi(multi) => multi.batch_query_raw_overlapping(queries),
        }
    }

    fn compute_sample_degrees(
        &self,
        seq_included: &[bool],
        seq_to_sample: &[u16],
    ) -> Vec<u16> {
        match self {
            // Single Impg: default trait impl (parallel-by-target) is already optimal.
            ImpgWrapper::Single(impg) => impg.compute_sample_degrees(seq_included, seq_to_sample),
            // MultiImpg: file-parallel override to bound peak memory.
            ImpgWrapper::Multi(multi) => multi.compute_sample_degrees(seq_included, seq_to_sample),
        }
    }
}
