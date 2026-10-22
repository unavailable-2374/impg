// src/multi_impg.rs
//! Multi-file IMPG index implementation.
//!
//! `MultiImpg` coordinates queries across multiple per-file `.impg` indices,
//! presenting a unified view while internally managing ID translation.

use crate::alignment_record::Strand;
use crate::forest_map::ForestMap;
use crate::impg::{
    compose_hop, extend_frontier_from_hit, AdjustedInterval, BfsHit, CigarCacheKey, CigarOp, Impg,
    NextCarry, QueryMetadata, SortedRanges, TransitiveRange,
};
use crate::impg_index::{ImpgIndex, RawAlignmentInterval};
use crate::seqidx::SequenceIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use log::{debug, info, warn};
use parking_lot::{Condvar as PlCondvar, Mutex as PlMutex};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, Once, RwLock};
use std::time::SystemTime;

/// Location of a tree within a specific sub-index.
///
/// Packed as `(index_idx as u64) << 32 | local_target_id as u64`.
///
/// On 64-bit Linux the natural `{ usize, u32 }` layout is 16 B due to padding;
/// at CHM13 / 580-file scale the unified `forest_map` holds ~280 K targets ×
/// ~580 locations each, so 16 → 8 B per entry saves ~1.3 GB of constant
/// resident memory. `index_idx` is bounded by the number of input alignment
/// files (always < 2³² in practice) and `local_target_id` is the per-file
/// target id (already u32 in the on-disk format), so the pack is lossless.
#[derive(Debug, Clone, Copy)]
#[repr(transparent)]
struct TreeLocation(u64);

impl TreeLocation {
    #[inline]
    fn new(index_idx: usize, local_target_id: u32) -> Self {
        debug_assert!(index_idx <= u32::MAX as usize, "index_idx overflows u32");
        TreeLocation(((index_idx as u64) << 32) | local_target_id as u64)
    }

    #[inline]
    fn index_idx(&self) -> usize {
        (self.0 >> 32) as usize
    }

    #[inline]
    fn local_target_id(&self) -> u32 {
        self.0 as u32
    }
}

/// Serializable version of TreeLocation for cache.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct TreeLocationSer {
    index_idx: u32,
    local_target_id: u32,
}

impl From<&TreeLocation> for TreeLocationSer {
    fn from(loc: &TreeLocation) -> Self {
        TreeLocationSer {
            index_idx: loc.index_idx() as u32,
            local_target_id: loc.local_target_id(),
        }
    }
}

impl From<TreeLocationSer> for TreeLocation {
    fn from(ser: TreeLocationSer) -> Self {
        TreeLocation::new(ser.index_idx as usize, ser.local_target_id)
    }
}

/// File entry for staleness detection in the cache manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct FileEntry {
    /// Path to the index file (relative or absolute as stored in list)
    path: String,
    /// File size in bytes
    size: u64,
    /// Modification time as duration since UNIX_EPOCH
    mtime_secs: u64,
}

/// Cache for MultiImpg unified data.
///
/// This cache stores precomputed unified sequence index, forest map, and
/// local-to-unified translation tables to speed up repeated queries with
/// the same set of per-file indices.
#[derive(Debug, Serialize, Deserialize)]
pub struct MultiImpgCache {
    /// Magic bytes for identification
    magic: [u8; 10],
    /// Version for format evolution
    version: u32,
    /// Manifest of index files for staleness detection
    manifest: Vec<FileEntry>,
    /// The alignment list file path (for auto-detection)
    list_file_path: String,
    /// Unified sequence index
    unified_seq_index: SequenceIndex,
    /// Unified forest map: target_id → list of tree locations
    unified_forest_map: Vec<(u32, Vec<TreeLocationSer>)>,
    /// Local-to-unified translation tables per index
    local_to_unified: Vec<Vec<u32>>,
    /// Whether all sub-indices are bidirectional (V2 format)
    is_bidirectional: bool,
}

const CACHE_MAGIC: &[u8; 10] = b"MIMPGCACH1";
const CACHE_VERSION: u32 = 2; // bumped: added is_bidirectional field

/// Metadata loaded from a per-file index header (seq_index + forest_map only).
struct IndexHeader {
    /// Sequence index from this file
    seq_index: SequenceIndex,
    /// Forest map from this file
    forest_map: ForestMap,
    /// Whether this index is bidirectional (V2 format)
    is_bidirectional: bool,
}

/// Multi-file IMPG index.
///
/// Coordinates queries across multiple per-file indices while presenting
/// a unified interface via the `ImpgIndex` trait.
pub struct MultiImpg {
    // ============ UNIFIED VIEW (what callers see) ============
    /// Unified sequence index: name ↔ unified_id
    pub seq_index: SequenceIndex,

    /// Unified forest map: unified_target_id → Vec<TreeLocation>
    /// Multiple indices may have trees for the same sequence
    forest_map: FxHashMap<u32, Vec<TreeLocation>>,

    // ============ PER-INDEX DATA (internal only) ============
    /// Per-index file paths
    index_paths: Vec<PathBuf>,

    /// Per-index alignment file paths
    alignment_files: Vec<String>,

    /// Per-index sequence files (if any)
    sequence_files: Vec<String>,

    /// Per-index ID translation: local_id → unified_id
    /// Uses Vec for O(1) indexed access since local IDs are dense 0..n
    local_to_unified: Vec<Vec<u32>>,

    /// Lazily-loaded sub-indices (only loaded when tree data is needed)
    sub_indices: RwLock<Vec<Option<Arc<Impg>>>>,

    /// Per-index singleflight locks. Without these, concurrent misses all
    /// opened and deserialized the same large index before racing to populate
    /// one cache slot, multiplying disk traffic and allocation pressure.
    sub_index_load_locks: Vec<PlMutex<()>>,

    /// FIFO load-order of the currently-resident `sub_indices` slots. Enables
    /// *incremental* eviction (drop the oldest resident until back under budget)
    /// instead of flushing every slot when a soft bound is tripped. At all-vs-all
    /// scale (10^5 per-file indices, ~580 files touched per query × ~100 threads)
    /// a full flush wiped the entire warm working set on every miss past the cap,
    /// forcing constant re-decompression from disk — the dominant cost of
    /// CIGAR-precise depth. Only mutated under the `sub_indices` write lock, so
    /// its own mutex is effectively uncontended (present only for interior
    /// mutability of a `&self` method).
    sub_index_lru: PlMutex<VecDeque<usize>>,

    /// Number of populated (`Some`) slots in `sub_indices`.
    ///
    /// Mirrors `transient_cache_count` but guards the BFS/transitive cache.
    /// Maintained under the `sub_indices` write lock so it stays consistent
    /// with the actual slot population without an O(num_indices) scan per miss.
    sub_index_cache_count: AtomicUsize,

    /// Soft upper bound on `sub_index_cache_count`. When a fresh `get_sub_index`
    /// miss would push residency above this, the oldest residents are evicted
    /// one at a time (FIFO, via `sub_index_lru`) until the new one fits.
    ///
    /// Why this matters: the CIGAR-precise transitive depth path
    /// (`--cigar-precise` / `--use-BFS`) keeps tree caching ON for re-use across
    /// BFS hops (see `depth.rs` `set_tree_cache_enabled(true)`), and
    /// `clear_sub_index_cache()` only fires at the *end* of Phase 1. With
    /// hundreds of thousands of per-file indices that means `sub_indices` grows
    /// monotonically toward "every file" mid-phase — each tree-pinned `Arc<Impg>`
    /// backs several mmap'd allocations, so residency blows past the kernel
    /// `vm.max_map_count` limit (default 65530 on Linux) and the allocator
    /// aborts with `memory allocation of N bytes failed` long before RSS nears
    /// the host limit. Bounding residency trades a periodic header/tree reload
    /// for survival on these workloads. Eviction is transparent to query
    /// results (it only forces a reload), so output is byte-identical to the
    /// unbounded path.
    ///
    /// Tunable via `IMPG_SUB_INDEX_CACHE_LIMIT` (number of slots). Default:
    /// `min(num_indices, 8192)`. Set to `0` to disable bounding (recovers the
    /// prior unbounded behaviour). For typical single-region queries this is
    /// never reached, so their cached-tree re-use is unchanged.
    ///
    /// Stored as an atomic so a caller that knows the concurrency it is about
    /// to drive (e.g. the CIGAR-precise depth path, which fans out one chunk
    /// per rayon thread) can adaptively *lower* the cap before the parallel
    /// region via [`MultiImpg::set_sub_index_cache_limit`]. Peak residency in
    /// this cache is bounded by the limit (a full-flush fires once the count
    /// reaches it), so the limit — not the thread count — is the dominant
    /// memory lever for that path.
    sub_index_cache_limit: AtomicUsize,

    /// True when `sub_index_cache_limit` came from an explicit
    /// `IMPG_SUB_INDEX_CACHE_LIMIT` env override. When set, the adaptive
    /// lowering in [`MultiImpg::set_sub_index_cache_limit`] is suppressed so
    /// the user's chosen value is honoured exactly.
    sub_index_cache_limit_explicit: bool,

    /// On-disk size (bytes) of each per-file index, parallel to `index_paths`.
    /// Used as a proxy for the resident footprint of a loaded sub-index when
    /// byte-budgeting the `sub_indices` cache (see `sub_index_cache_byte_budget`).
    index_sizes: Vec<u64>,

    /// Running sum of estimated resident bytes for the populated `sub_indices`
    /// slots (`index_sizes[i] * SUB_INDEX_RESIDENT_EXPANSION`). Maintained under
    /// the `sub_indices` write lock alongside `sub_index_cache_count`.
    sub_index_cache_bytes: AtomicU64,

    /// Soft upper bound (estimated resident bytes) on the `sub_indices` cache.
    /// When a fresh miss would push `sub_index_cache_bytes` over this, every slot
    /// is evicted first — the same full-flush policy as the count cap, but
    /// bounding RAM instead of slot count. `0` = no byte bound (count cap only).
    ///
    /// Why a byte budget in addition to the count cap: `.impg` sizes span ~4
    /// orders of magnitude (KB to ~1 GB) on all-vs-all workloads, so a *count*
    /// cap cannot bound RAM. The CIGAR-precise depth path anchors Phase 1 on the
    /// highest-degree hub sequences, whose neighbour indices are exactly the
    /// large tail of the size distribution, so a few thousand cached slots can
    /// reach >100 GB and OOM. The byte budget bounds RAM directly; the count cap
    /// independently bounds mmap/VMA pressure (`vm.max_map_count`) for the
    /// opposite regime of very many tiny files. Both flush full on exceed.
    ///
    /// Tunable via `IMPG_SUB_INDEX_CACHE_BYTES` (accepts a raw byte count or a
    /// `K`/`M`/`G`/`T` suffix). Default `0`; the CIGAR-precise depth path lowers
    /// it to a fraction of system RAM via `set_sub_index_cache_byte_budget`.
    sub_index_cache_byte_budget: AtomicU64,

    /// True when `sub_index_cache_byte_budget` came from an explicit
    /// `IMPG_SUB_INDEX_CACHE_BYTES` override, suppressing the adaptive lowering
    /// in [`MultiImpg::set_sub_index_cache_byte_budget`].
    sub_index_cache_byte_budget_explicit: bool,

    /// Per-file lazy header cache for the **transient** query path.
    ///
    /// Each slot holds an `Arc<Impg>` parsed from one alignment file's index
    /// header (seq_index + forest_map only). Trees are NOT pinned: each cached
    /// `Impg` has `set_tree_cache_enabled(false)`, so subsequent
    /// `get_or_load_tree` calls fetch the COITree from disk and the
    /// `Arc<COITree>` drops as soon as the caller releases it.
    ///
    /// Why a separate cache from `sub_indices`:
    ///   - `sub_indices` services the BFS / transitive query path which
    ///     deliberately keeps tree caching ON for re-use across BFS hops; we
    ///     can't safely flip the cache flag once an `Arc<Impg>` has been
    ///     handed out from there.
    ///   - The transient path (`load_sub_index_transient`, used by chunked
    ///     non-transitive Phase 1/2) was previously a fresh `File::open` +
    ///     bincode-decode per call. Once Phase 1 was chunked at 5 MB the
    ///     same file gets revisited dozens of times per chromosome × every
    ///     hub chromosome, multiplying the O(N_files) header-parse cost by
    ///     two orders of magnitude.
    ///
    /// Per-slot `PlMutex` (vs a single `RwLock<Vec<...>>`) so the first miss
    /// for file A doesn't block the first miss for file B.
    transient_header_cache: Vec<PlMutex<Option<Arc<Impg>>>>,

    /// FIFO of populated transient-header slots.  Unlike the old full-cache
    /// flush, this lets a miss evict only as many entries as needed and avoids
    /// scanning/locking every one of 10^5+ slots at each capacity crossing.
    transient_header_fifo: PlMutex<VecDeque<usize>>,

    /// Number of populated slots in `transient_header_cache`.
    ///
    /// Tracked separately from the slot mutexes so we can decide on a cheap
    /// upper bound without scanning all `num_indices` slots. Updated under
    /// the same per-slot lock that flips a slot from `None` → `Some` (or
    /// vice-versa) so it stays consistent with the cache contents.
    transient_cache_count: AtomicUsize,

    /// Soft upper bound on `transient_cache_count`. When a fresh miss would
    /// push the cache above this, `load_sub_index_transient` evicts every
    /// slot before populating the new one.
    ///
    /// Why this matters: with hundreds of thousands of per-file indices,
    /// retaining one `Arc<Impg>` per file blows past the kernel
    /// `vm.max_map_count` limit (default 65530 on Linux) — each cached `Impg`
    /// holds several internal allocations and glibc spreads them across
    /// per-thread arenas, each of which costs VMA slots. Long before RSS
    /// approaches the host limit the allocator returns ENOMEM and Rust
    /// aborts with `memory allocation of N bytes failed`. A bounded cache
    /// trades a small amount of redundant header parsing for survival on
    /// these workloads.
    ///
    /// Tunable via `IMPG_TRANSIENT_HEADER_CACHE_LIMIT` (number of slots).
    /// Default: `min(num_indices, 8192)`. Set to `0` to disable bounding
    /// (recovers the prior unbounded behaviour).
    transient_cache_limit: usize,

    /// Whether all indices are bidirectional (V2 format)
    /// True only if ALL sub-indices are V2 format
    is_bidirectional: bool,

    /// Whether tree caching is enabled for sub-indices.
    /// Propagated to newly lazy-loaded sub-indices.
    tree_cache_enabled: std::sync::atomic::AtomicBool,
}

/// Resolve the transient-header-cache size limit from the environment, with a
/// safe default. Returns `0` to mean "unbounded" (legacy behaviour).
fn resolve_transient_cache_limit(num_indices: usize) -> usize {
    if let Ok(s) = std::env::var("IMPG_TRANSIENT_HEADER_CACHE_LIMIT") {
        if let Ok(v) = s.parse::<usize>() {
            return v;
        } else {
            warn!(
                "IMPG_TRANSIENT_HEADER_CACHE_LIMIT='{}' is not a non-negative integer; using default",
                s
            );
        }
    }
    // 8192 keeps the cache well under typical vm.max_map_count budgets even
    // when several internal allocations per cached Impg back into mmap, and
    // is large enough that file-locality re-hits dominate at chunked depth
    // workloads on per-file indices in the few-thousand-files regime.
    num_indices.min(8192)
}

/// Bound simultaneous transient sub-index/tree/CIGAR working sets.  Rayon
/// thread count is a CPU setting and can be 112+ on HPC nodes; using it as the
/// file-I/O concurrency also allowed 112 large per-file indices and CIGAR
/// caches to coexist.  Keep a conservative default and allow an explicit
/// workload-specific override.
fn transient_file_query_concurrency() -> usize {
    std::env::var("IMPG_FILE_QUERY_CONCURRENCY")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or_else(|| rayon::current_num_threads().min(16).max(1))
}

/// Resident-byte budget for concurrently processed file-first query tasks.
/// This is deliberately separate from the long-lived sub-index cache budget:
/// these tasks own transient trees and CIGAR buffers, not `sub_indices` entries.
/// An explicit `0` leaves the count limit as the sole governor; otherwise the
/// default is 25% of the smallest detected physical/cgroup/Slurm memory limit.
fn transient_file_query_byte_budget() -> u64 {
    if let Ok(raw) = std::env::var("IMPG_FILE_QUERY_BYTES") {
        if let Some(bytes) = parse_byte_size(&raw) {
            return bytes;
        }
        warn!(
            "IMPG_FILE_QUERY_BYTES='{}' is not a valid byte size (e.g. 16G); using adaptive default",
            raw
        );
    }
    detected_process_memory_limit_bytes()
        .map(|bytes| bytes / 4)
        .unwrap_or(0)
}

/// Best-effort memory available to this process. Shared HPC nodes frequently
/// expose the node's full `MemTotal` while Slurm grants the job only a subset,
/// so take the minimum of physical RAM, cgroup limits and common Slurm memory
/// variables. Values of zero and kernel "unlimited" sentinels are ignored.
fn detected_process_memory_limit_bytes() -> Option<u64> {
    fn take_limit(current: &mut Option<u64>, candidate: u64) {
        if candidate > 0 && candidate < (1u64 << 62) {
            *current = Some(current.map_or(candidate, |old| old.min(candidate)));
        }
    }

    let mut limit = std::fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|content| {
            content.lines().find_map(|line| {
                line.strip_prefix("MemTotal:")?
                    .trim()
                    .trim_end_matches("kB")
                    .trim()
                    .parse::<u64>()
                    .ok()
                    .map(|kb| kb.saturating_mul(1024))
            })
        });

    for path in [
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    ] {
        if let Ok(raw) = std::fs::read_to_string(path) {
            if let Ok(bytes) = raw.trim().parse::<u64>() {
                take_limit(&mut limit, bytes);
            }
        }
    }

    const MIB: u64 = 1024 * 1024;
    if let Ok(raw) = std::env::var("SLURM_MEM_PER_NODE") {
        if let Ok(mib) = raw.parse::<u64>() {
            take_limit(&mut limit, mib.saturating_mul(MIB));
        }
    }
    if let (Ok(mem_raw), Ok(cpus_raw)) = (
        std::env::var("SLURM_MEM_PER_CPU"),
        std::env::var("SLURM_CPUS_ON_NODE"),
    ) {
        if let (Ok(mib_per_cpu), Ok(cpus)) = (mem_raw.parse::<u64>(), cpus_raw.parse::<u64>()) {
            take_limit(
                &mut limit,
                mib_per_cpu.saturating_mul(cpus).saturating_mul(MIB),
            );
        }
    }

    limit
}

fn file_query_limiter() -> FileQueryLimiter {
    let max_active = transient_file_query_concurrency();
    let max_bytes = transient_file_query_byte_budget();
    static LOGGED: Once = Once::new();
    LOGGED.call_once(|| {
        info!(
            "Transient file-query scheduler: concurrency={}, resident-byte-budget={} ({})",
            max_active,
            max_bytes,
            if max_bytes == 0 {
                "unbounded".to_string()
            } else {
                format!("{:.1} GiB", max_bytes as f64 / (1_u64 << 30) as f64)
            }
        );
    });
    FileQueryLimiter::new(max_active, max_bytes)
}

#[derive(Default)]
struct FileQueryLimitState {
    active: usize,
    bytes: u64,
}

/// Count + byte weighted semaphore used inside rayon file-parallel queries.
/// It removes the fixed-size `.chunks(N)` barriers: when one large file is
/// slow, completed workers can immediately start later small files instead of
/// waiting for the entire cohort.
struct FileQueryLimiter {
    max_active: usize,
    max_bytes: u64,
    state: PlMutex<FileQueryLimitState>,
    changed: PlCondvar,
}

impl FileQueryLimiter {
    fn new(max_active: usize, max_bytes: u64) -> Self {
        Self {
            max_active: max_active.max(1),
            max_bytes,
            state: PlMutex::new(FileQueryLimitState::default()),
            changed: PlCondvar::new(),
        }
    }

    fn acquire(&self, estimated_bytes: u64) -> FileQueryPermit<'_> {
        let weight = estimated_bytes.max(1);
        let mut state = self.state.lock();
        while state.active >= self.max_active
            || (self.max_bytes > 0
                && state.active > 0
                && state.bytes.saturating_add(weight) > self.max_bytes)
        {
            self.changed.wait(&mut state);
        }
        state.active += 1;
        state.bytes = state.bytes.saturating_add(weight);
        FileQueryPermit {
            limiter: self,
            weight,
        }
    }
}

struct FileQueryPermit<'a> {
    limiter: &'a FileQueryLimiter,
    weight: u64,
}

impl Drop for FileQueryPermit<'_> {
    fn drop(&mut self) {
        let mut state = self.limiter.state.lock();
        state.active -= 1;
        state.bytes = state.bytes.saturating_sub(self.weight);
        self.limiter.changed.notify_all();
    }
}

/// Restores a transient sub-index to header-only mode even when a strict CIGAR
/// read returns early with an error.
struct TransientTreeCacheGuard<'a> {
    impg: &'a Impg,
}

impl<'a> TransientTreeCacheGuard<'a> {
    fn new(impg: &'a Impg) -> Self {
        impg.set_tree_cache_enabled(true);
        Self { impg }
    }
}

impl Drop for TransientTreeCacheGuard<'_> {
    fn drop(&mut self) {
        self.impg.clear_tree_cache();
        self.impg.set_tree_cache_enabled(false);
    }
}

/// Resolve the BFS/transitive sub-index cache size limit from the environment.
/// Returns `0` to mean "unbounded" (legacy behaviour). Mirrors
/// `resolve_transient_cache_limit` but for the `sub_indices` cache.
///
/// Returns `(limit, explicit)` where `explicit` is `true` only when the value
/// came from a successfully-parsed env override, so the caller can suppress the
/// adaptive lowering applied to the default.
fn resolve_sub_index_cache_limit(num_indices: usize) -> (usize, bool) {
    if let Ok(s) = std::env::var("IMPG_SUB_INDEX_CACHE_LIMIT") {
        if let Ok(v) = s.parse::<usize>() {
            return (v, true);
        } else {
            warn!(
                "IMPG_SUB_INDEX_CACHE_LIMIT='{}' is not a non-negative integer; using default",
                s
            );
        }
    }
    // Same default rationale as the transient header cache: bound resident
    // tree-pinned sub-indices to keep mmap pressure under vm.max_map_count at
    // hundreds-of-thousands-of-files scale, while staying high enough that
    // ordinary few-thousand-file workloads never evict.
    (num_indices.min(8192), false)
}

/// Estimated resident-bytes multiplier (percent) over the on-disk `.impg` size.
///
/// A `.impg` is essentially the serialized COITree forest: the CIGAR strings
/// live in the alignment (PAF) file and are read transiently per query (never
/// cached), and `load_from_file` decodes only `seq_index` + `forest_map` plus,
/// on query, the per-target COITrees that the tree cache pins. The
/// bincode→in-memory expansion of those trees was **measured at 1.06×** the
/// on-disk size (jemalloc, settled RSS), stable across 300- and 1160-file
/// samples of the VGP all-vs-all ref working set (`examples/measure_resident`):
/// the full 1160-file ref neighbour set was 9287 MB on disk → 9845 MB resident,
/// header/base a negligible ~25 KB/file. We use `125%` (≈ +18% over the measured
/// 1.06×) as a safety margin covering the minority of files cached in both
/// alignment directions (their reverse-direction trees load too) and allocator
/// noise. Over-estimating only costs extra cache reload churn; under-estimating
/// risks the OOM this budget prevents.
const SUB_INDEX_RESIDENT_EXPANSION_PCT: u64 = 125;

/// Estimated resident bytes for the sub-index loaded from a `.impg` of `size`.
#[inline]
fn estimated_resident_bytes(size: u64) -> u64 {
    // size ≤ ~1 GB in practice, so `size * 125` cannot overflow u64; the
    // saturating mul guards the pathological case regardless.
    size.saturating_mul(SUB_INDEX_RESIDENT_EXPANSION_PCT) / 100
}

/// Parse a byte size that optionally carries a `K`/`M`/`G`/`T` (× 1024) suffix.
/// `"0"` is valid and means "unbounded". Returns `None` on malformed input.
fn parse_byte_size(s: &str) -> Option<u64> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    let (num, mult) = match s.chars().last().unwrap().to_ascii_uppercase() {
        'K' => (&s[..s.len() - 1], 1u64 << 10),
        'M' => (&s[..s.len() - 1], 1u64 << 20),
        'G' => (&s[..s.len() - 1], 1u64 << 30),
        'T' => (&s[..s.len() - 1], 1u64 << 40),
        '0'..='9' => (s, 1u64),
        _ => return None,
    };
    num.trim()
        .parse::<u64>()
        .ok()
        .map(|v| v.saturating_mul(mult))
}

/// Resolve the BFS/transitive sub-index cache **byte** budget from the
/// environment. Returns `(budget_bytes, explicit)`; `0` means "no byte bound"
/// (the count cap still applies). The CIGAR-precise depth path narrows the
/// default to a RAM fraction via `set_sub_index_cache_byte_budget`.
fn resolve_sub_index_cache_byte_budget() -> (u64, bool) {
    if let Ok(s) = std::env::var("IMPG_SUB_INDEX_CACHE_BYTES") {
        if let Some(v) = parse_byte_size(&s) {
            return (v, true);
        } else {
            warn!(
                "IMPG_SUB_INDEX_CACHE_BYTES='{}' is not a valid byte size (e.g. 32G); ignoring",
                s
            );
        }
    }
    (0, false)
}

impl MultiImpg {
    /// Return one deterministic alignment-file connectivity label per input
    /// sequence. Only `seq_ids` participate: a Phase-1 hub that is already
    /// fully processed must not collapse otherwise independent Phase-2 leaf
    /// components merely because every leaf has a pairwise file with that hub.
    ///
    /// Union-find scans the existing compact `forest_map` locations in place;
    /// it allocates O(sequences + files), not another copy of the potentially
    /// hundreds-of-millions of sequence/file incidences.
    pub(crate) fn depth_locality_components(&self, seq_ids: &[u32]) -> Vec<u32> {
        if seq_ids.is_empty() {
            return Vec::new();
        }

        fn find(parent: &mut [usize], mut node: usize) -> usize {
            while parent[node] != node {
                parent[node] = parent[parent[node]];
                node = parent[node];
            }
            node
        }

        fn union(parent: &mut [usize], rank: &mut [u8], left: usize, right: usize) {
            let mut left_root = find(parent, left);
            let mut right_root = find(parent, right);
            if left_root == right_root {
                return;
            }
            if rank[left_root] < rank[right_root] {
                std::mem::swap(&mut left_root, &mut right_root);
            }
            parent[right_root] = left_root;
            if rank[left_root] == rank[right_root] {
                rank[left_root] = rank[left_root].saturating_add(1);
            }
        }

        let mut parent: Vec<usize> = (0..seq_ids.len()).collect();
        let mut rank = vec![0u8; seq_ids.len()];
        let mut first_seq_by_file = vec![usize::MAX; self.index_paths.len()];

        for (seq_pos, seq_id) in seq_ids.iter().copied().enumerate() {
            if let Some(locations) = self.forest_map.get(&seq_id) {
                for location in locations {
                    let file_idx = location.index_idx();
                    let first = first_seq_by_file[file_idx];
                    if first == usize::MAX {
                        first_seq_by_file[file_idx] = seq_pos;
                    } else {
                        union(&mut parent, &mut rank, seq_pos, first);
                    }
                }
            }
        }

        // Canonicalise labels to the smallest sequence ID in each component so
        // the result is independent of union-by-rank choices.
        let mut min_seq_by_root: FxHashMap<usize, u32> = FxHashMap::default();
        for (idx, seq_id) in seq_ids.iter().copied().enumerate() {
            let root = find(&mut parent, idx);
            min_seq_by_root
                .entry(root)
                .and_modify(|current| *current = (*current).min(seq_id))
                .or_insert(seq_id);
        }
        seq_ids
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let root = find(&mut parent, idx);
                min_seq_by_root[&root]
            })
            .collect()
    }

    /// Load headers from multiple per-file indices and build unified mappings.
    ///
    /// This loads ONLY the headers (seq_index + forest_map) from each file,
    /// NOT the tree data. Trees are loaded on demand.
    pub fn load_from_files(
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
    ) -> std::io::Result<Self> {
        let num_indices = index_paths.len();
        // Load headers in parallel
        let headers: Vec<IndexHeader> = index_paths
            .par_iter()
            .map(|path| {
                Self::load_header(path).map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load header from {:?}: {}", path, e),
                    )
                })
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Build unified sequence index
        let mut unified_seq_index = SequenceIndex::new();
        let mut local_to_unified: Vec<Vec<u32>> = Vec::with_capacity(num_indices);

        for header in &headers {
            // Pre-allocate Vec with capacity for all local IDs
            let mut l2u = Vec::with_capacity(header.seq_index.len());

            for local_id in 0..header.seq_index.len() as u32 {
                if let Some(name) = header.seq_index.get_name(local_id) {
                    let len = header.seq_index.get_len_from_id(local_id);
                    let unified_id = unified_seq_index.get_or_insert_id(name, len);
                    l2u.push(unified_id);
                } else {
                    // Should not happen for valid indices, but handle gracefully
                    l2u.push(u32::MAX);
                }
            }

            local_to_unified.push(l2u);
        }

        // Build unified forest map
        let mut unified_forest_map: FxHashMap<u32, Vec<TreeLocation>> = FxHashMap::default();

        for (index_idx, header) in headers.iter().enumerate() {
            let l2u = &local_to_unified[index_idx];

            for &local_target_id in header.forest_map.entries.keys() {
                let unified_id = l2u[local_target_id as usize];
                if unified_id != u32::MAX {
                    unified_forest_map
                        .entry(unified_id)
                        .or_default()
                        .push(TreeLocation::new(index_idx, local_target_id));
                }
            }
        }

        // Check if all indices are bidirectional
        let all_bidirectional = headers.iter().all(|h| h.is_bidirectional);
        if !all_bidirectional {
            let v1_count = headers.iter().filter(|h| !h.is_bidirectional).count();
            warn!(
                "{} of {} indices are V1 (unidirectional). Rebuild with default settings for full bidirectional support.",
                v1_count, num_indices
            );
        }

        info!(
            "Built unified index with {} sequences and {} targets (bidirectional: {})",
            unified_seq_index.len(),
            unified_forest_map.len(),
            all_bidirectional
        );

        let (sub_index_cache_limit, sub_index_cache_limit_explicit) =
            resolve_sub_index_cache_limit(num_indices);
        let (sub_index_cache_byte_budget, sub_index_cache_byte_budget_explicit) =
            resolve_sub_index_cache_byte_budget();
        // On-disk size per index, parallel to `index_paths`, for byte-budgeting
        // the sub-index cache. Cheap relative to the parallel header load above
        // (one `stat` per file vs a full bincode header decode); a stat failure
        // maps to 0 so a missing file never inflates the accounting.
        let index_sizes: Vec<u64> = index_paths
            .par_iter()
            .map(|p| fs::metadata(p).map(|m| m.len()).unwrap_or(0))
            .collect();
        Ok(Self {
            seq_index: unified_seq_index,
            forest_map: unified_forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            sub_index_load_locks: (0..num_indices).map(|_| PlMutex::new(())).collect(),
            sub_index_lru: PlMutex::new(VecDeque::new()),
            sub_index_cache_count: AtomicUsize::new(0),
            sub_index_cache_limit: AtomicUsize::new(sub_index_cache_limit),
            sub_index_cache_limit_explicit,
            index_sizes,
            sub_index_cache_bytes: AtomicU64::new(0),
            sub_index_cache_byte_budget: AtomicU64::new(sub_index_cache_byte_budget),
            sub_index_cache_byte_budget_explicit,
            transient_header_cache: (0..num_indices).map(|_| PlMutex::new(None)).collect(),
            transient_header_fifo: PlMutex::new(VecDeque::new()),
            transient_cache_count: AtomicUsize::new(0),
            transient_cache_limit: resolve_transient_cache_limit(num_indices),
            is_bidirectional: all_bidirectional,
            tree_cache_enabled: std::sync::atomic::AtomicBool::new(true),
        })
    }

    /// Load MultiImpg from a cache file if valid, or build from scratch.
    ///
    /// Auto-detects cache file as `{list_file}.multi_impg` and validates
    /// staleness before using. If cache is stale or missing, builds from
    /// scratch and saves a new cache.
    pub fn load_with_cache(
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
        list_file: &Path,
    ) -> std::io::Result<Self> {
        let list_str = list_file.to_string_lossy();
        let cache_path = if list_str.starts_with("/proc/") || list_str.starts_with("/dev/fd/") {
            // Process substitution: derive deterministic cache path from alignment file paths
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};
            let mut hasher = DefaultHasher::new();
            for f in alignment_files {
                f.hash(&mut hasher);
            }
            let hash = hasher.finish();
            let dir = Path::new(&alignment_files[0])
                .parent()
                .unwrap_or(Path::new("."));
            let cache_file = dir.join(format!("impg_cache_{:016x}.multi_impg", hash));
            debug!(
                "Alignment list is a process substitution ({}), using deterministic cache: {}",
                list_str,
                cache_file.display()
            );
            cache_file
        } else {
            list_file.with_extension("multi_impg")
        };

        // Try to load from cache
        if cache_path.exists() {
            match MultiImpgCache::load(&cache_path) {
                Ok(cache) => {
                    if cache.is_valid(index_paths, list_file)? {
                        return Self::from_cache(
                            cache,
                            index_paths,
                            alignment_files,
                            sequence_files,
                        );
                    } else {
                        info!("Cache stale, rebuilding...");
                    }
                }
                Err(e) => {
                    warn!(
                        "Failed to load cache {:?}: {}, rebuilding...",
                        cache_path, e
                    );
                }
            }
        }

        // Build from scratch
        let multi = Self::load_from_files(index_paths, alignment_files, sequence_files)?;

        // Save cache for next time
        if let Err(e) = multi.save_cache(&cache_path, list_file) {
            warn!("Failed to save cache {:?}: {}", cache_path, e);
        }

        Ok(multi)
    }

    /// Create MultiImpg from a validated cache.
    fn from_cache(
        cache: MultiImpgCache,
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
    ) -> std::io::Result<Self> {
        let num_indices = index_paths.len();

        // Convert serialized forest map back to FxHashMap<u32, Vec<TreeLocation>>
        let forest_map: FxHashMap<u32, Vec<TreeLocation>> = cache
            .unified_forest_map
            .into_iter()
            .map(|(target_id, locs)| {
                (
                    target_id,
                    locs.into_iter().map(TreeLocation::from).collect(),
                )
            })
            .collect();

        info!(
            "Loaded {} sequences and {} targets from cache",
            cache.unified_seq_index.len(),
            forest_map.len()
        );

        let (sub_index_cache_limit, sub_index_cache_limit_explicit) =
            resolve_sub_index_cache_limit(num_indices);
        let (sub_index_cache_byte_budget, sub_index_cache_byte_budget_explicit) =
            resolve_sub_index_cache_byte_budget();
        // Reuse the cache manifest's per-file sizes — `is_valid` (run before
        // `from_cache`) already confirmed `manifest[i]` corresponds to
        // `index_paths[i]` and matches its on-disk size, so no extra `stat`s are
        // needed. Fall back to a parallel stat only if the manifest is somehow
        // the wrong length (defensive; should not happen post-validation).
        let index_sizes: Vec<u64> = if cache.manifest.len() == num_indices {
            cache.manifest.iter().map(|e| e.size).collect()
        } else {
            index_paths
                .par_iter()
                .map(|p| fs::metadata(p).map(|m| m.len()).unwrap_or(0))
                .collect()
        };
        Ok(Self {
            seq_index: cache.unified_seq_index,
            forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified: cache.local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            sub_index_load_locks: (0..num_indices).map(|_| PlMutex::new(())).collect(),
            sub_index_lru: PlMutex::new(VecDeque::new()),
            sub_index_cache_count: AtomicUsize::new(0),
            sub_index_cache_limit: AtomicUsize::new(sub_index_cache_limit),
            sub_index_cache_limit_explicit,
            index_sizes,
            sub_index_cache_bytes: AtomicU64::new(0),
            sub_index_cache_byte_budget: AtomicU64::new(sub_index_cache_byte_budget),
            sub_index_cache_byte_budget_explicit,
            transient_header_cache: (0..num_indices).map(|_| PlMutex::new(None)).collect(),
            transient_header_fifo: PlMutex::new(VecDeque::new()),
            transient_cache_count: AtomicUsize::new(0),
            transient_cache_limit: resolve_transient_cache_limit(num_indices),
            is_bidirectional: cache.is_bidirectional,
            tree_cache_enabled: std::sync::atomic::AtomicBool::new(true),
        })
    }

    /// Save the unified index data to a cache file.
    pub fn save_cache(&self, cache_path: &Path, list_file: &Path) -> std::io::Result<()> {
        // Build manifest from index files
        let manifest: Vec<FileEntry> = self
            .index_paths
            .iter()
            .map(|path| {
                let metadata = fs::metadata(path)?;
                let mtime = metadata
                    .modified()?
                    .duration_since(SystemTime::UNIX_EPOCH)
                    .map_err(std::io::Error::other)?;
                Ok(FileEntry {
                    path: path.to_string_lossy().to_string(),
                    size: metadata.len(),
                    mtime_secs: mtime.as_secs(),
                })
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Convert forest_map to serializable format
        let unified_forest_map: Vec<(u32, Vec<TreeLocationSer>)> = self
            .forest_map
            .iter()
            .map(|(&target_id, locs)| (target_id, locs.iter().map(TreeLocationSer::from).collect()))
            .collect();

        let cache = MultiImpgCache {
            magic: *CACHE_MAGIC,
            version: CACHE_VERSION,
            manifest,
            list_file_path: list_file.to_string_lossy().to_string(),
            unified_seq_index: self.seq_index.clone(),
            unified_forest_map,
            local_to_unified: self.local_to_unified.clone(),
            is_bidirectional: self.is_bidirectional,
        };

        let file = File::create(cache_path)?;
        let mut writer = BufWriter::new(file);

        bincode::serde::encode_into_std_write(&cache, &mut writer, bincode::config::standard())
            .map_err(std::io::Error::other)?;

        writer.flush()?;
        Ok(())
    }

    /// Load only the header (seq_index + forest_map) from a single index file.
    fn load_header(path: &Path) -> std::io::Result<IndexHeader> {
        const MAGIC_V1: &[u8] = b"IMPGIDX1";
        const MAGIC_V2: &[u8] = b"IMPGIDX2";

        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read and verify magic bytes (support both V1 and V2)
        let mut magic_buf = [0u8; 8];
        reader.read_exact(&mut magic_buf)?;
        if magic_buf != MAGIC_V1 && magic_buf != MAGIC_V2 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Invalid magic bytes in {:?}", path),
            ));
        }

        // Read forest map offset
        let mut offset_buf = [0u8; 8];
        reader.read_exact(&mut offset_buf)?;
        let forest_map_offset = u64::from_le_bytes(offset_buf);

        // Read sequence index
        let seq_index: SequenceIndex =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load sequence index: {e}"),
                    )
                })?;

        // Seek to forest map and read it
        reader.seek(SeekFrom::Start(forest_map_offset))?;
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load forest map: {e}"),
                    )
                })?;

        let is_bidirectional = magic_buf == MAGIC_V2;
        Ok(IndexHeader {
            seq_index,
            forest_map,
            is_bidirectional,
        })
    }

    /// Get or load a sub-index.
    fn get_sub_index(&self, index_idx: usize) -> std::io::Result<Arc<Impg>> {
        // Fast path: check if already loaded
        {
            let indices = self.sub_indices.read().unwrap();
            if let Some(ref impg) = indices[index_idx] {
                return Ok(Arc::clone(impg));
            }
        }

        // Only one worker may perform the expensive load for a given slot.
        // Re-check after acquiring because another worker may have populated it
        // while this one waited.
        let _load_guard = self.sub_index_load_locks[index_idx].lock();
        {
            let indices = self.sub_indices.read().unwrap();
            if let Some(ref impg) = indices[index_idx] {
                return Ok(Arc::clone(impg));
            }
        }

        // Slow path: load the index
        let path = &self.index_paths[index_idx];
        let alignment_files = vec![self.alignment_files[index_idx].clone()];
        let seq_files = if self.sequence_files.is_empty() {
            None
        } else {
            Some(self.sequence_files.as_slice())
        };

        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let impg = Impg::load_from_file(
            reader,
            &alignment_files,
            path.to_string_lossy().to_string(),
            seq_files,
        )?;
        let impg = Arc::new(impg);

        // Propagate tree cache setting to the newly loaded sub-index
        let cache_enabled = self
            .tree_cache_enabled
            .load(std::sync::atomic::Ordering::Relaxed);
        impg.set_tree_cache_enabled(cache_enabled);

        // Store in sub-index cache, bounding residency on two independent axes:
        // a slot *count* cap (mmap/VMA pressure under vm.max_map_count, for the
        // many-tiny-files regime) and an estimated-resident *byte* budget (RAM,
        // for the few-large-files regime that a count cap cannot bound).
        let new_bytes = estimated_resident_bytes(self.index_sizes[index_idx]);
        let byte_budget = self
            .sub_index_cache_byte_budget
            .load(AtomicOrdering::Relaxed);

        // Backstop: never let a single oversized index into the cache. If this
        // one file's estimated resident already exceeds the whole byte budget,
        // caching it would either bust the budget outright or force a full flush
        // of every other slot on each access. Hand back the Arc (the caller
        // still needs it to answer the query) but leave it uncached.
        let oversized = byte_budget > 0 && new_bytes > byte_budget;

        // Evicted `Arc<Impg>`s are moved here and dropped only *after* the write
        // lock is released. Deallocating a large sub-index's COITrees can be
        // expensive; doing it inside the lock would serialize every other
        // `get_sub_index` behind a free() storm (the exact write-lock contention
        // this cache change targets). The slot is already `None` by then, so the
        // memory is logically evicted the instant we take() it — only the
        // physical free is deferred.
        let mut evicted: Vec<Arc<Impg>> = Vec::new();

        if !oversized {
            let mut indices = self.sub_indices.write().unwrap();

            // Count only genuine None -> Some transitions (a racing thread may
            // have populated this slot before we took the lock).
            if indices[index_idx].is_none() {
                // Before storing the fresh resident, evict the *oldest* residents
                // one at a time until adding this one stays within both soft
                // bounds. This replaces the previous full-flush (drop every slot
                // on the first miss past the cap), which wiped the entire warm
                // working set repeatedly: at all-vs-all scale the concurrent
                // working set (~580 files/query × thread count) far exceeds the
                // cap, so a full flush thrashed the cache to disk on nearly every
                // miss. FIFO incremental eviction keeps the hot majority resident.
                //
                // All mutation happens under the `sub_indices` write lock, so the
                // counters and the LRU deque stay exactly consistent with the slot
                // population. Eviction only drops cache references — any
                // `Arc<Impg>` already handed to an in-flight query stays alive
                // until that query drops it — so this is transparent to results.
                let cache_limit = self.sub_index_cache_limit.load(AtomicOrdering::Relaxed);
                {
                    let mut lru = self.sub_index_lru.lock();
                    loop {
                        let count = self.sub_index_cache_count.load(AtomicOrdering::Relaxed);
                        let bytes = self.sub_index_cache_bytes.load(AtomicOrdering::Relaxed);
                        let count_exceeds = cache_limit > 0 && count + 1 > cache_limit;
                        let bytes_exceed =
                            byte_budget > 0 && bytes.saturating_add(new_bytes) > byte_budget;
                        if !(count_exceeds || bytes_exceed) {
                            break;
                        }
                        // Pop the oldest resident. Skip stale deque entries whose
                        // slot is already empty (evicted via another path).
                        let victim = match lru.pop_front() {
                            Some(v) => v,
                            None => break, // nothing left to evict
                        };
                        if let Some(old) = indices[victim].take() {
                            let freed = estimated_resident_bytes(self.index_sizes[victim]);
                            self.sub_index_cache_count
                                .fetch_sub(1, AtomicOrdering::Relaxed);
                            self.sub_index_cache_bytes
                                .fetch_sub(freed, AtomicOrdering::Relaxed);
                            // Defer the physical free until the lock is released.
                            evicted.push(old);
                        }
                    }

                    self.sub_index_cache_count
                        .fetch_add(1, AtomicOrdering::Relaxed);
                    self.sub_index_cache_bytes
                        .fetch_add(new_bytes, AtomicOrdering::Relaxed);
                    indices[index_idx] = Some(Arc::clone(&impg));
                    lru.push_back(index_idx);
                }
            }
        }

        // Write lock (and `lru`) are now released; free evicted sub-indices
        // outside the critical section.
        drop(evicted);

        Ok(impg)
    }

    /// Load a sub-index WITHOUT storing it in `self.sub_indices`.
    ///
    /// Originally used by the file-parallel pre-scan (`compute_sample_degrees`)
    /// to bound peak memory; now also the workhorse for chunked Phase 1/2
    /// non-transitive depth.
    ///
    /// The returned `Arc<Impg>` has `tree_cache_enabled = false`, so any tree
    /// walked via `get_or_load_tree` is freed when the caller drops the
    /// `Arc<COITree>` it received — `self.trees` never accumulates.
    ///
    /// Header caching: the parsed `Impg` (seq_index + forest_map; trees not
    /// included) is stashed in `self.transient_header_cache[index_idx]` after
    /// the first call. Subsequent calls return a clone of the cached `Arc`,
    /// avoiding the `File::open` + bincode-decode round trip. With chunked
    /// Phase 1 hitting the same file once per 5 MB chunk per hub chromosome,
    /// the cache cuts header-parse work by 2–3 orders of magnitude on
    /// CHM13-scale workloads. Per-file `PlMutex` keeps misses for distinct
    /// files independent.
    fn load_sub_index_transient(&self, index_idx: usize) -> std::io::Result<Arc<Impg>> {
        // Cache hit fast path.
        {
            let slot = self.transient_header_cache[index_idx].lock();
            if let Some(ref impg) = *slot {
                return Ok(Arc::clone(impg));
            }
        }

        // Cache miss. Incrementally evict oldest populated slots until this
        // entry fits.  The previous implementation dropped the entire cache
        // and scanned all `num_indices` mutexes at every capacity crossing;
        // with 331k indices that destroyed the warm set and created repeated
        // O(num_indices) lock storms.
        if self.transient_cache_limit > 0 {
            while self.transient_cache_count.load(AtomicOrdering::Relaxed)
                >= self.transient_cache_limit
                && self.evict_one_transient_header()
            {}
        }

        let arc = self.load_sub_index_uncached(index_idx)?;

        // Stash for next time. Race-tolerant: if another thread populated the
        // slot first, we drop our copy and use theirs (functionally identical).
        let mut slot = self.transient_header_cache[index_idx].lock();
        let (result, inserted) = match *slot {
            Some(ref existing) => (Arc::clone(existing), false),
            None => {
                *slot = Some(Arc::clone(&arc));
                self.transient_cache_count
                    .fetch_add(1, AtomicOrdering::Relaxed);
                (arc, true)
            }
        };
        drop(slot);
        if inserted {
            self.transient_header_fifo.lock().push_back(index_idx);
            // Concurrent misses may overshoot the soft cap by at most the
            // in-flight loader count. Bring it back without a full scan. This
            // runs after releasing our slot lock so evicting this just-added
            // entry cannot self-deadlock.
            if self.transient_cache_limit > 0 {
                while self.transient_cache_count.load(AtomicOrdering::Relaxed)
                    > self.transient_cache_limit
                    && self.evict_one_transient_header()
                {}
            }
        }
        Ok(result)
    }

    /// Evict one populated transient header from the FIFO. Stale queue entries
    /// are possible after an explicit full clear and are skipped cheaply.
    fn evict_one_transient_header(&self) -> bool {
        loop {
            let Some(victim) = self.transient_header_fifo.lock().pop_front() else {
                return false;
            };
            let mut slot = self.transient_header_cache[victim].lock();
            if slot.take().is_some() {
                self.transient_cache_count
                    .fetch_sub(1, AtomicOrdering::Relaxed);
                return true;
            }
        }
    }

    /// Drop every cached header and reset the populated-slot counter.
    ///
    /// Concurrency: another thread may walk past the limit check and start
    /// loading a fresh sub-index while we're evicting. That's safe — its
    /// final write-back acquires the per-slot mutex *after* our eviction
    /// has released it, so the two operations serialize and the final
    /// `transient_cache_count` matches the populated-slot total.
    fn evict_transient_header_cache(&self) {
        let mut dropped: usize = 0;
        for slot in &self.transient_header_cache {
            let mut s = slot.lock();
            if s.is_some() {
                *s = None;
                dropped += 1;
            }
        }
        // Decrement by what we actually dropped, so concurrent populate-misses
        // don't drive the counter negative.
        if dropped > 0 {
            self.transient_cache_count
                .fetch_sub(dropped, AtomicOrdering::Relaxed);
        }
        self.transient_header_fifo.lock().clear();
    }

    /// Like `load_sub_index_transient` but never writes to
    /// `transient_header_cache`. Use for one-shot scans that visit each file
    /// exactly once (e.g. the degree pre-scan in `compute_sample_degrees`),
    /// where caching only ratchets up retained allocations / VMAs.
    ///
    /// Will still honour an existing cache hit: if another caller has already
    /// stashed this `index_idx`, we hand back that `Arc` to avoid a redundant
    /// disk parse. The new behaviour is purely "do not pollute the cache on a
    /// miss" — peak retained sub-indices for a parallel pre-scan stays bounded
    /// by the rayon worker count, regardless of `index_paths.len()`.
    ///
    /// Why this matters: with hundreds of thousands of per-file indices,
    /// caching every loaded `Impg` in `transient_header_cache` blows past
    /// `vm.max_map_count` (default 65530 on Linux) long before RSS gets close
    /// to the host limit, because each cached `Impg` retains several
    /// glibc-mmap'd allocations. The allocator then aborts the process with
    /// `memory allocation of N bytes failed` even though physical memory is
    /// nowhere near exhausted.
    fn load_sub_index_uncached(&self, index_idx: usize) -> std::io::Result<Arc<Impg>> {
        let path = &self.index_paths[index_idx];
        let alignment_files = vec![self.alignment_files[index_idx].clone()];
        let seq_files = if self.sequence_files.is_empty() {
            None
        } else {
            Some(self.sequence_files.as_slice())
        };
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let impg = Impg::load_from_file(
            reader,
            &alignment_files,
            path.to_string_lossy().to_string(),
            seq_files,
        )?;
        // Always disable tree caching on transient sub-indices so that trees
        // walked through the loaded header are released as soon as the caller
        // drops the Arc returned by get_or_load_tree.
        impg.set_tree_cache_enabled(false);
        Ok(Arc::new(impg))
    }

    /// Translate an AdjustedInterval from local IDs to unified IDs.
    fn translate_to_unified(
        &self,
        interval: AdjustedInterval,
        index_idx: usize,
    ) -> Option<AdjustedInterval> {
        let (query_interval, cigar, target_interval) = interval;
        let l2u = &self.local_to_unified[index_idx];

        // Vec indexing with bounds check
        let unified_query_id = *l2u.get(query_interval.metadata as usize)?;
        let unified_target_id = *l2u.get(target_interval.metadata as usize)?;

        // Check for invalid sentinel values
        if unified_query_id == u32::MAX || unified_target_id == u32::MAX {
            return None;
        }

        Some((
            Interval {
                first: query_interval.first,
                last: query_interval.last,
                metadata: unified_query_id,
            },
            cigar,
            Interval {
                first: target_interval.first,
                last: target_interval.last,
                metadata: unified_target_id,
            },
        ))
    }

    /// Query all sub-indices that have trees for the given target.
    fn query_all_indices(
        &self,
        unified_target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> io::Result<Vec<AdjustedInterval>> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => {
                return Ok(vec![self.make_self_interval(
                    unified_target_id,
                    range_start,
                    range_end,
                    store_cigar,
                )])
            }
        };

        // Query all relevant sub-indices in parallel
        let results: Vec<Vec<AdjustedInterval>> = locations
            .par_iter()
            .map(|loc| -> io::Result<Vec<AdjustedInterval>> {
                let impg = self.get_sub_index(loc.index_idx()).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "query failed to load sub-index '{}': {e}",
                            self.index_paths[loc.index_idx()].display()
                        ),
                    )
                })?;

                // Query using local target ID
                let local_results = impg.query_without_self(
                    loc.local_target_id(),
                    range_start,
                    range_end,
                    store_cigar,
                    min_gap_compressed_identity,
                    sequence_index,
                    approximate_mode,
                );

                // Translate results to unified IDs
                let unified_results: Vec<AdjustedInterval> = local_results
                    .into_iter()
                    .filter_map(|r| self.translate_to_unified(r, loc.index_idx()))
                    .collect();

                Ok(unified_results)
            })
            .collect::<io::Result<Vec<_>>>()?;

        // Merge all results with exactly one synthetic self-interval. Building
        // it up front lets us sort only the non-self tail; the previous
        // remove(0) + insert(0) pair shifted every result twice per query.
        let result_capacity = results.iter().map(Vec::len).sum::<usize>() + 1;
        let mut final_results = Vec::with_capacity(result_capacity);
        final_results.push(self.make_self_interval(
            unified_target_id,
            range_start,
            range_end,
            store_cigar,
        ));

        for result_set in results {
            final_results.extend(result_set);
        }

        // Sort results for deterministic ordering, leaving self first.
        // Sort by: query_id, query_start, query_end, target_start, target_end
        if final_results.len() > 2 {
            final_results[1..].sort_by(|a, b| {
                let a_key = (a.0.metadata, a.0.first, a.0.last, a.2.first, a.2.last);
                let b_key = (b.0.metadata, b.0.first, b.0.last, b.2.first, b.2.last);
                a_key.cmp(&b_key)
            });
        }

        Ok(final_results)
    }

    /// Per-file-streaming, CIGAR-compressing variant of `query_all_indices`.
    ///
    /// Identical merged/sorted output to `query_all_indices(.., store_cigar=true)`
    /// except every overlap's CIGAR is run through `compress` the moment its file
    /// finishes reconstructing — so the uncompressed per-base CIGAR of any single
    /// alignment file is freed before the next file's is built. The accumulated
    /// `results` therefore holds only compressed CIGARs (tiny) instead of all
    /// `degree` files' uncompressed CIGARs at once. This is what keeps the depth
    /// hop-0 CIGAR sweep from OOMing on all-vs-all at 10^5 per-file indices.
    fn query_all_indices_compressed_impl(
        &self,
        unified_target_id: u32,
        range_start: i64,
        range_end: i64,
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
        collect_raw: bool,
    ) -> (Vec<AdjustedInterval>, Vec<RawAlignmentInterval>) {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => {
                return (
                    vec![self.make_self_interval(unified_target_id, range_start, range_end, true)],
                    Vec::new(),
                )
            }
        };

        // Query all relevant sub-indices in parallel, compressing each file's
        // CIGARs before they leave the worker so only one file's uncompressed
        // CIGAR is ever resident per thread.
        let results: Vec<(Vec<AdjustedInterval>, Vec<RawAlignmentInterval>)> = locations
            .par_iter()
            .filter_map(|loc| {
                let impg = match self.get_sub_index(loc.index_idx()) {
                    Ok(i) => i,
                    Err(e) => {
                        warn!("Failed to load sub-index {}: {}", loc.index_idx(), e);
                        return None;
                    }
                };

                let (local_results, local_raw) = if collect_raw {
                    impg.query_without_self_with_raw(
                        loc.local_target_id(),
                        range_start,
                        range_end,
                        true, // store_cigar
                        None,
                        None,
                        false,
                    )
                } else {
                    (
                        impg.query_without_self(
                            loc.local_target_id(),
                            range_start,
                            range_end,
                            true,
                            None,
                            None,
                            false,
                        ),
                        Vec::new(),
                    )
                };

                let unified_results: Vec<AdjustedInterval> = local_results
                    .into_iter()
                    .filter_map(|r| self.translate_to_unified(r, loc.index_idx()))
                    .map(|mut r| {
                        // Compress and drop the uncompressed CIGAR immediately.
                        r.1 = compress(&r.1);
                        r
                    })
                    .collect();

                // Raw extents were collected in the same COITree traversal as
                // the projected CIGARs; only ID translation remains here.
                let unified_raw = if collect_raw {
                    let l2u = &self.local_to_unified[loc.index_idx()];
                    local_raw
                        .into_iter()
                        .filter_map(|mut raw| {
                            let unified_query_id = *l2u.get(raw.query_id as usize)?;
                            if unified_query_id == u32::MAX {
                                return None;
                            }
                            raw.query_id = unified_query_id;
                            Some(raw)
                        })
                        .collect()
                } else {
                    Vec::new()
                };

                Some((unified_results, unified_raw))
            })
            .collect();

        // Merge with one synthetic unified self interval + deterministic tail
        // sort — identical to `query_all_indices` so downstream `cigar_idx`
        // assignment matches.
        let result_capacity = results.iter().map(|(set, _)| set.len()).sum::<usize>() + 1;
        let raw_capacity = results.iter().map(|(_, raw)| raw.len()).sum();
        let mut final_results = Vec::with_capacity(result_capacity);
        final_results.push(self.make_self_interval(
            unified_target_id,
            range_start,
            range_end,
            true,
        ));
        let mut final_raw = Vec::with_capacity(raw_capacity);
        for (result_set, mut raw_set) in results {
            final_raw.append(&mut raw_set);
            final_results.extend(result_set);
        }
        if final_results.len() > 2 {
            final_results[1..].sort_by(|a, b| {
                let a_key = (a.0.metadata, a.0.first, a.0.last, a.2.first, a.2.last);
                let b_key = (b.0.metadata, b.0.first, b.0.last, b.2.first, b.2.last);
                a_key.cmp(&b_key)
            });
        }
        (final_results, final_raw)
    }

    fn query_all_indices_compressed(
        &self,
        unified_target_id: u32,
        range_start: i64,
        range_end: i64,
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
    ) -> Vec<AdjustedInterval> {
        self.query_all_indices_compressed_impl(
            unified_target_id,
            range_start,
            range_end,
            compress,
            false,
        )
        .0
    }

    fn query_all_indices_compressed_with_raw(
        &self,
        unified_target_id: u32,
        range_start: i64,
        range_end: i64,
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
    ) -> (Vec<AdjustedInterval>, Vec<RawAlignmentInterval>) {
        self.query_all_indices_compressed_impl(
            unified_target_id,
            range_start,
            range_end,
            compress,
            true,
        )
    }

    /// Create a self-referential interval for the query region.
    fn make_self_interval(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
    ) -> AdjustedInterval {
        (
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            if store_cigar {
                CigarOp::new_run(range_end - range_start, '=')
            } else {
                Vec::new()
            },
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
        )
    }
}

impl ImpgIndex for MultiImpg {
    fn seq_index(&self) -> &SequenceIndex {
        &self.seq_index
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        self.query_all_indices(
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
        )
    }

    fn query_overlapping_cigar_compressed(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
    ) -> Vec<AdjustedInterval> {
        self.query_all_indices_compressed(target_id, range_start, range_end, compress)
    }

    fn query_overlapping_cigar_compressed_with_raw(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
    ) -> (Vec<AdjustedInterval>, Vec<RawAlignmentInterval>) {
        self.query_all_indices_compressed_with_raw(target_id, range_start, range_end, compress)
    }

    fn batch_query_overlapping_cigar_compressed_with_raw(
        &self,
        queries: &[(u32, i64, i64)],
        compress: &(dyn Fn(&[CigarOp]) -> Vec<CigarOp> + Sync),
    ) -> io::Result<Vec<(Vec<AdjustedInterval>, Vec<RawAlignmentInterval>)>> {
        let n = queries.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        // Invert query→file fan-out once for the whole batch. This is the core
        // locality transformation needed by depth Phase 2: millions of small
        // gaps no longer deserialize the same pairwise index independently.
        let mut by_file: FxHashMap<usize, Vec<(usize, u32, i64, i64)>> = FxHashMap::default();
        for (qi, &(target_id, start, end)) in queries.iter().enumerate() {
            if let Some(locs) = self.forest_map.get(&target_id) {
                for loc in locs {
                    by_file.entry(loc.index_idx()).or_default().push((
                        qi,
                        loc.local_target_id(),
                        start,
                        end,
                    ));
                }
            }
        }

        let results: Vec<PlMutex<(Vec<AdjustedInterval>, Vec<RawAlignmentInterval>)>> = (0..n)
            .map(|_| PlMutex::new((Vec::new(), Vec::new())))
            .collect();
        let mut file_order: Vec<usize> = by_file.keys().copied().collect();
        // Largest-first scheduling reduces the long tail while the weighted
        // limiter bounds concurrent residency. Index id is a deterministic
        // tie-break; rayon execution order was already intentionally free.
        file_order.sort_unstable_by_key(|&i| (std::cmp::Reverse(self.index_sizes[i]), i));
        let limiter = file_query_limiter();

        file_order
            .par_iter()
            .try_for_each(|&file_idx| -> io::Result<()> {
                let _permit = limiter.acquire(estimated_resident_bytes(self.index_sizes[file_idx]));
                let impg = self.load_sub_index_transient(file_idx).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "load sub-index '{}' for batched CIGAR depth: {e}",
                            self.index_paths[file_idx].display()
                        ),
                    )
                })?;
                let l2u = &self.local_to_unified[file_idx];

                // Keep each referenced target tree resident for this file's
                // whole query batch, and read each distinct CIGAR once. The
                // transient loader normally disables tree caching because the
                // legacy caller executes only one query per load; batching
                // reverses that assumption.
                let _tree_cache_guard = TransientTreeCacheGuard::new(&impg);
                let mut cigar_cache: FxHashMap<CigarCacheKey, Vec<CigarOp>> = FxHashMap::default();
                let local_queries: Vec<(u32, i64, i64)> = by_file[&file_idx]
                    .iter()
                    .map(|&(_, local_target_id, start, end)| (local_target_id, start, end))
                    .collect();
                impg.populate_cigar_cache_strict_batch(
                    &local_queries,
                    None,
                    compress,
                    &mut cigar_cache,
                )?;

                for &(qi, local_target_id, start, end) in &by_file[&file_idx] {
                    let (local_overlaps, local_raw) = impg.query_with_cache_and_raw(
                        local_target_id,
                        start,
                        end,
                        true,
                        None,
                        None,
                        &cigar_cache,
                    );
                    let mut overlaps = Vec::with_capacity(local_overlaps.len());
                    for r in local_overlaps {
                        if let Some(mut translated) = self.translate_to_unified(r, file_idx) {
                            // The batch owns one synthetic self interval per query.
                            if translated.0.metadata == translated.2.metadata
                                && translated.0.first == start
                                && translated.0.last == end
                                && translated.2.metadata == queries[qi].0
                            {
                                continue;
                            }
                            // Cache values were compressed once above; the
                            // clipped CIGAR returned here is compressed again
                            // because clipping may split runs at range edges.
                            translated.1 = compress(&translated.1);
                            overlaps.push(translated);
                        }
                    }

                    let mut raw = Vec::with_capacity(local_raw.len());
                    for mut r in local_raw {
                        let Some(&unified_query_id) = l2u.get(r.query_id as usize) else {
                            continue;
                        };
                        if unified_query_id == u32::MAX {
                            continue;
                        }
                        r.query_id = unified_query_id;
                        raw.push(r);
                    }

                    let mut slot = results[qi].lock();
                    slot.0.append(&mut overlaps);
                    slot.1.append(&mut raw);
                }
                Ok(())
            })?;

        Ok(results
            .into_iter()
            .enumerate()
            .map(|(qi, slot)| {
                let (mut overlaps, raw) = slot.into_inner();
                overlaps.sort_unstable_by_key(|r| {
                    (r.0.metadata, r.0.first, r.0.last, r.2.first, r.2.last)
                });
                let (target_id, start, end) = queries[qi];
                overlaps.insert(0, self.make_self_interval(target_id, start, end, true));
                (overlaps, raw)
            })
            .collect())
    }

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        _cigar_cache: &FxHashMap<CigarCacheKey, Vec<CigarOp>>,
    ) -> io::Result<Vec<AdjustedInterval>> {
        // For MultiImpg, we don't use the shared cache since each sub-index
        // has its own file offsets. Just do a normal query.
        self.query(
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            false, // approximate_mode
        )
    }

    fn populate_cigar_cache(
        &self,
        _target_id: u32,
        _range_start: i64,
        _range_end: i64,
        _min_gap_compressed_identity: Option<f64>,
        _sequence_index: Option<&UnifiedSequenceIndex>,
        _cache: &mut FxHashMap<CigarCacheKey, Vec<CigarOp>>,
    ) {
        // For MultiImpg, CIGAR caching is not implemented since each sub-index
        // has different file offsets. This is a no-op.
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        // Transitive query implementation using DFS with deterministic ordering
        self.transitive_query_impl(
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
            true, // use_dfs
        )
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        // Transitive query implementation using BFS with deterministic ordering
        self.transitive_query_impl(
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
            false, // use_dfs
        )
    }

    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        // For MultiImpg, we can't return a single tree since multiple indices
        // may have data for the same target. Return None and let callers use query().
        // This is only used by stats and similarity commands which may need adaptation.

        // Try to get a tree from any sub-index that has this target
        let locations = self.forest_map.get(&target_id)?;

        // Get the first location and try to load its tree
        let loc = locations.first()?;
        let local_target_id = loc.local_target_id();

        let impg = self.get_sub_index(loc.index_idx()).ok()?;
        impg.get_or_load_tree(local_target_id)
    }

    fn target_ids(&self) -> Vec<u32> {
        self.forest_map.keys().copied().collect()
    }

    fn remove_cached_tree(&self, _target_id: u32) {
        // For MultiImpg, we don't cache trees at this level.
        // Sub-indices manage their own tree caches.
    }

    fn num_targets(&self) -> usize {
        self.forest_map.len()
    }

    fn sequence_files(&self) -> &[String] {
        &self.sequence_files
    }

    fn alignment_files(&self) -> &[String] {
        &self.alignment_files
    }

    fn query_reverse_for_depth(&self, query_id: u32) -> Vec<(i64, i64, i64, i64, u32)> {
        let mut results = Vec::new();

        // Iterate through all target_ids in the forest map
        for &target_id in self.forest_map.keys() {
            // Skip if querying self
            if target_id == query_id {
                continue;
            }

            if let Some(tree) = self.get_or_load_tree(target_id) {
                for interval in tree.iter() {
                    if interval.metadata.query_id() == query_id {
                        let query_start = interval.metadata.query_start();
                        let query_end = interval.metadata.query_end();
                        let target_start = interval.first as i64;
                        let target_end = interval.last as i64;
                        results.push((query_start, query_end, target_start, target_end, target_id));
                    }
                }
            }
        }

        results
    }

    fn build_query_to_targets_map(&self) -> FxHashMap<u32, Vec<u32>> {
        let mut query_to_targets: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

        for &target_id in self.forest_map.keys() {
            if let Some(tree) = self.get_or_load_tree(target_id) {
                let mut seen_queries: rustc_hash::FxHashSet<u32> = rustc_hash::FxHashSet::default();
                for interval in tree.iter() {
                    let qid = interval.metadata.query_id();
                    if qid != target_id && seen_queries.insert(qid) {
                        query_to_targets.entry(qid).or_default().push(target_id);
                    }
                }
            }
        }

        // Clear tree cache after building the map
        self.clear_tree_cache();

        query_to_targets
    }

    fn query_reverse_for_depth_with_map(
        &self,
        query_id: u32,
        query_to_targets: &FxHashMap<u32, Vec<u32>>,
    ) -> Vec<(i64, i64, i64, i64, u32)> {
        let mut results = Vec::new();

        if let Some(target_ids) = query_to_targets.get(&query_id) {
            for &target_id in target_ids {
                if let Some(tree) = self.get_or_load_tree(target_id) {
                    for interval in tree.iter() {
                        if interval.metadata.query_id() == query_id {
                            let query_start = interval.metadata.query_start();
                            let query_end = interval.metadata.query_end();
                            let target_start = interval.first as i64;
                            let target_end = interval.last as i64;
                            results.push((
                                query_start,
                                query_end,
                                target_start,
                                target_end,
                                target_id,
                            ));
                        }
                    }
                }
            }
        }

        results
    }

    fn clear_tree_cache(&self) {
        // For MultiImpg, clear caches of all loaded sub-indices
        let indices = self.sub_indices.read().unwrap();
        for sub_index in indices.iter().flatten() {
            sub_index.clear_tree_cache();
        }
    }

    fn clear_sub_index_cache(&self) {
        // Clear the sub-index cache to free memory
        let mut indices = self.sub_indices.write().unwrap();
        for slot in indices.iter_mut() {
            *slot = None;
        }
        self.sub_index_lru.lock().clear();
        self.sub_index_cache_count.store(0, AtomicOrdering::Relaxed);
        self.sub_index_cache_bytes.store(0, AtomicOrdering::Relaxed);
    }

    fn set_sub_index_cache_limit(&self, limit: usize) {
        // Adaptive *lowering* only: a caller that knows the concurrency it is
        // about to drive (the CIGAR-precise depth path runs the BFS / transitive
        // query path, which keeps tree caching ON and so pins each resident
        // sub-index's COITrees) can shrink the cap to bound peak residency. Peak
        // memory in this cache is `cap × per-sub-index-bytes`: once the populated
        // count reaches the cap, each fresh miss evicts the oldest resident
        // (FIFO, incremental — no full flush), so the cap is the dominant memory
        // lever for that path.
        //
        // An explicit `IMPG_SUB_INDEX_CACHE_LIMIT` override is honoured exactly
        // (never overridden here). We only ever move the cap DOWN: raising it
        // mid-run could let residency exceed the budget the caller picked, and
        // `clear_sub_index_cache` already handles the count, so a lower cap just
        // means the next miss past it triggers an earlier flush. A `limit` of 0
        // would mean "unbounded" elsewhere, so it is ignored here (we never
        // adaptively *disable* bounding).
        if self.sub_index_cache_limit_explicit || limit == 0 {
            return;
        }
        let current = self.sub_index_cache_limit.load(AtomicOrdering::Relaxed);
        if current == 0 || limit < current {
            self.sub_index_cache_limit
                .store(limit, AtomicOrdering::Relaxed);
        }
    }

    fn relax_sub_index_cache_count_cap(&self) {
        // Raise the slot-count cap to `num_indices` so the byte budget is the
        // sole residency governor. Safe because sub-indices are heap-loaded and
        // jemalloc arena-packs them (~0.1 VMAs/file measured), so even the full
        // byte-budget working set stays far under `vm.max_map_count`. An explicit
        // `IMPG_SUB_INDEX_CACHE_LIMIT` is honoured exactly (never overridden).
        if self.sub_index_cache_limit_explicit {
            return;
        }
        let num_indices = self.index_paths.len();
        let current = self.sub_index_cache_limit.load(AtomicOrdering::Relaxed);
        if current == 0 || num_indices > current {
            self.sub_index_cache_limit
                .store(num_indices, AtomicOrdering::Relaxed);
        }
    }

    fn set_sub_index_cache_byte_budget(&self, budget_bytes: u64) {
        // Adaptive *lowering* only, mirroring `set_sub_index_cache_limit` but for
        // the RAM (byte) axis. The CIGAR-precise depth path picks a fraction of
        // system RAM as the budget before its parallel region; this is the
        // dominant memory lever at all-vs-all scale, where index sizes span
        // orders of magnitude and the slot-count cap cannot bound RAM.
        //
        // An explicit `IMPG_SUB_INDEX_CACHE_BYTES` override is honoured exactly.
        // `budget_bytes == 0` means "unbounded" and is ignored here (we never
        // adaptively *disable* bounding). We only ever move the budget DOWN.
        if self.sub_index_cache_byte_budget_explicit || budget_bytes == 0 {
            return;
        }
        let current = self
            .sub_index_cache_byte_budget
            .load(AtomicOrdering::Relaxed);
        if current == 0 || budget_bytes < current {
            self.sub_index_cache_byte_budget
                .store(budget_bytes, AtomicOrdering::Relaxed);
        }
    }

    fn clear_transient_header_cache(&self) {
        // Drop every Arc<Impg> stashed by load_sub_index_transient so the
        // backing mmap regions are returned to the kernel. Independent from
        // clear_sub_index_cache, which targets the BFS/transitive cache.
        // Routes through evict_transient_header_cache so the populated-slot
        // counter that gates our soft cap stays in sync.
        self.evict_transient_header_cache();
    }

    fn set_tree_cache_enabled(&self, enabled: bool) {
        // Store at MultiImpg level so newly lazy-loaded sub-indices inherit the setting
        self.tree_cache_enabled
            .store(enabled, std::sync::atomic::Ordering::Relaxed);
        // Also propagate to already-loaded sub-indices
        let indices = self.sub_indices.read().unwrap();
        for sub_index in indices.iter().flatten() {
            sub_index.set_tree_cache_enabled(enabled);
        }
    }

    fn is_bidirectional(&self) -> bool {
        self.is_bidirectional
    }

    fn query_raw_intervals(&self, unified_target_id: u32) -> Vec<RawAlignmentInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => return Vec::new(),
        };

        let mut results = Vec::new();
        for loc in locations {
            let impg = match self.get_sub_index(loc.index_idx()) {
                Ok(i) => i,
                Err(_) => continue,
            };
            if let Some(tree) = impg.get_or_load_tree(loc.local_target_id()) {
                let l2u = &self.local_to_unified[loc.index_idx()];
                for interval in tree.iter() {
                    let m = &interval.metadata;
                    let unified_query_id = match l2u.get(m.query_id() as usize) {
                        Some(&id) if id != u32::MAX => id,
                        _ => continue,
                    };
                    results.push(RawAlignmentInterval {
                        target_start: interval.first,
                        target_end: interval.last,
                        query_id: unified_query_id,
                        query_start: m.query_start(),
                        query_end: m.query_end(),
                        is_reverse: m.is_reverse_strand(),
                    });
                }
            }
        }
        results
    }

    /// File-parallel pre-scan of unique-sample degrees.
    ///
    /// Overrides the default trait implementation (which iterates targets in parallel
    /// and loads every sub-index that touches each target). That default scales as
    /// O(num_files) retained sub-indices and OOMs with hundreds of thousands of
    /// per-file indices. This override iterates files in parallel: each worker loads
    /// ONE sub-index via `load_sub_index_uncached`, walks its trees, accumulates
    /// degree contributions into a shared per-target aggregator, then drops the
    /// sub-index. The no-cache loader bypasses `transient_header_cache` so peak
    /// retained sub-indices = rayon worker count, regardless of `num_files`.
    ///
    /// Correctness vs. default:
    /// - Both count "unique OTHER samples with direct alignments to this target".
    /// - The aggregator merges contributions from every file that contains a tree for
    ///   the same unified target (which is exactly what `query_raw_intervals` does
    ///   under the hood via the unified forest_map).
    /// - Excluded sequences (`seq_included[id] == false`) are skipped on both sides.
    fn compute_sample_degrees(
        &self,
        seq_included: &[bool],
        seq_to_sample: &[u16],
    ) -> io::Result<Vec<u16>> {
        use std::sync::atomic::{AtomicU64, Ordering};

        let num_unified = self.seq_index.len();

        // Aggregator is a flat lock-free bitset: one row of `chunks_per_row`
        // u64 chunks per unified target, one bit per sample.
        //
        // Why this beats the previous `Vec<PlMutex<FxHashSet<u16>>>`:
        //   1. Memory: one contiguous Vec<AtomicU64> instead of `num_unified`
        //      independent FxHashSet heap allocations + per-target Mutex.
        //      For 280K targets × ~600 samples that's ~21 MB vs. ≥160 MB at
        //      full occupancy, and — more importantly on hosts with low
        //      `vm.max_map_count` — collapses tens of thousands of small
        //      allocations into a single mmap-backed VMA.
        //   2. Concurrency: `fetch_or` on a u64 chunk is wait-free, so the
        //      per-file workers no longer queue on a per-target mutex when
        //      two files report alignments to the same target.
        //   3. Insertion is O(1) per (target, sample) pair instead of
        //      hash-and-rehash inside the local FxHashSet plus a mutex
        //      critical section to merge it.
        let max_sample = seq_to_sample.iter().copied().max().unwrap_or(0);
        let num_samples = max_sample as usize + 1;
        let chunks_per_row = num_samples.div_ceil(64);
        let total_chunks = num_unified.checked_mul(chunks_per_row).unwrap_or(0);
        let mut bitset: Vec<AtomicU64> = Vec::with_capacity(total_chunks);
        bitset.resize_with(total_chunks, || AtomicU64::new(0));

        (0..self.index_paths.len()).into_par_iter().try_for_each(
            |index_idx| -> io::Result<()> {
                // Pre-scan visits each file exactly once. Use the no-cache
                // loader so the transient header cache doesn't accumulate
                // hundreds of thousands of `Arc<Impg>` headers (each with its
                // own glibc-mmap'd allocations) and trip the kernel's
                // vm.max_map_count limit on hosts with many per-file indices.
                let impg = self.load_sub_index_uncached(index_idx).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "degree pre-scan failed to load '{}': {e}",
                            self.index_paths[index_idx].display()
                        ),
                    )
                })?;
                let l2u = &self.local_to_unified[index_idx];

                // Iterate every local target in this file.
                let local_target_ids: Vec<u32> = impg.forest_map.entries.keys().copied().collect();

                // Per-thread reusable scratch buffers — drop only when the
                // closure returns, so a single allocation amortises across
                // every target processed by this rayon worker on this file.
                //
                // `local_seen` is a dense byte bitmap (Vec<u8>) so dedup is
                // O(1) per interval. We also remember which sample ids were
                // touched in `touched_samples` so we can clear `local_seen`
                // in O(unique-samples) instead of O(num_samples) between
                // targets — this matters when num_samples is much larger
                // than the actual fan-out of any one target.
                let mut local_seen: Vec<u8> = vec![0u8; num_samples];
                let mut touched_samples: Vec<u16> = Vec::new();

                for local_target_id in local_target_ids {
                    let unified_target_id = match l2u.get(local_target_id as usize) {
                        Some(&id) if id != u32::MAX => id,
                        _ => continue,
                    };
                    if !seq_included
                        .get(unified_target_id as usize)
                        .copied()
                        .unwrap_or(false)
                    {
                        continue;
                    }
                    let self_sample = seq_to_sample
                        .get(unified_target_id as usize)
                        .copied()
                        .unwrap_or(0);

                    let tree = match impg.get_or_load_tree(local_target_id) {
                        Some(t) => t,
                        None => continue,
                    };

                    // Reset only the bytes we touched on the previous target.
                    for &s in &touched_samples {
                        local_seen[s as usize] = 0;
                    }
                    touched_samples.clear();

                    for interval in tree.iter() {
                        let m = &interval.metadata;
                        let unified_query_id = match l2u.get(m.query_id() as usize) {
                            Some(&id) if id != u32::MAX => id,
                            _ => continue,
                        };
                        if !seq_included
                            .get(unified_query_id as usize)
                            .copied()
                            .unwrap_or(false)
                        {
                            continue;
                        }
                        let query_sample = seq_to_sample
                            .get(unified_query_id as usize)
                            .copied()
                            .unwrap_or(0);
                        if query_sample != self_sample {
                            // Safe: query_sample <= max_sample by construction.
                            let slot =
                                unsafe { local_seen.get_unchecked_mut(query_sample as usize) };
                            if *slot == 0 {
                                *slot = 1;
                                touched_samples.push(query_sample);
                            }
                        }
                    }

                    if !touched_samples.is_empty() {
                        let row_base = unified_target_id as usize * chunks_per_row;
                        for &sample in &touched_samples {
                            let chunk_idx = sample as usize / 64;
                            let bit_idx = sample as usize % 64;
                            // SAFETY: chunk_idx < chunks_per_row and
                            // row_base + chunks_per_row <= total_chunks.
                            unsafe {
                                bitset
                                    .get_unchecked(row_base + chunk_idx)
                                    .fetch_or(1u64 << bit_idx, Ordering::Relaxed);
                            }
                        }
                    }
                    // Tree Arc dropped here; with tree caching disabled on this
                    // transient Impg, the underlying COITree is freed immediately.
                }
                // Impg Arc dropped here — entire sub-index is freed.
                Ok(())
            },
        )?;

        // Final pass: popcount each row to recover the per-target degree.
        if chunks_per_row == 0 {
            return Ok(vec![0u16; num_unified]);
        }
        Ok((0..num_unified)
            .into_par_iter()
            .map(|t| {
                let row_base = t * chunks_per_row;
                let mut count: u32 = 0;
                for c in 0..chunks_per_row {
                    count += bitset[row_base + c].load(Ordering::Relaxed).count_ones();
                }
                count.min(u16::MAX as u32) as u16
            })
            .collect())
    }

    fn query_raw_overlapping(
        &self,
        unified_target_id: u32,
        start: i64,
        end: i64,
    ) -> Vec<RawAlignmentInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => return Vec::new(),
        };

        let mut results = Vec::new();
        for loc in locations {
            let impg = match self.get_sub_index(loc.index_idx()) {
                Ok(i) => i,
                Err(_) => continue,
            };
            if let Some(tree) = impg.get_or_load_tree(loc.local_target_id()) {
                let l2u = &self.local_to_unified[loc.index_idx()];
                tree.query(start, end, |interval| {
                    let m = &interval.metadata;
                    if let Some(&unified_query_id) = l2u.get(m.query_id() as usize) {
                        if unified_query_id != u32::MAX {
                            results.push(RawAlignmentInterval {
                                target_start: interval.first,
                                target_end: interval.last,
                                query_id: unified_query_id,
                                query_start: m.query_start(),
                                query_end: m.query_end(),
                                is_reverse: m.is_reverse_strand(),
                            });
                        }
                    }
                });
            }
        }
        results
    }

    /// Transient variant: loads each needed sub-index via `load_sub_index_transient`
    /// (which disables tree caching on the returned `Impg`) and drops it as soon as
    /// the relevant trees have been walked. The shared `self.sub_indices` vec is
    /// NEVER written to. Peak retained memory per call is bounded by one sub-index
    /// plus one COITree, regardless of how many alignment files contain the target.
    ///
    /// This is required by the depth command's non-transitive Phase 1/2 hot paths
    /// when running with ≫ 10⁴ per-file indices: the cached `query_raw_intervals`
    /// path would otherwise monotonically grow `sub_indices` across hub sequences
    /// until the process commit limit is hit.
    fn query_raw_intervals_transient(&self, unified_target_id: u32) -> Vec<RawAlignmentInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => return Vec::new(),
        };

        // Group by sub-index file so we load each file at most once per call even
        // if the unified forest map were to report multiple local target ids per
        // file (V2 bidirectional currently reports one, but we dedupe defensively).
        let mut by_index: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
        for loc in locations {
            by_index
                .entry(loc.index_idx())
                .or_default()
                .push(loc.local_target_id());
        }

        let mut results = Vec::new();
        for (index_idx, local_target_ids) in by_index {
            let impg = match self.load_sub_index_transient(index_idx) {
                Ok(i) => i,
                Err(e) => {
                    warn!(
                        "query_raw_intervals_transient: failed to load sub-index {:?}: {}",
                        self.index_paths[index_idx], e
                    );
                    continue;
                }
            };
            let l2u = &self.local_to_unified[index_idx];
            for local_target_id in local_target_ids {
                if let Some(tree) = impg.get_or_load_tree(local_target_id) {
                    for interval in tree.iter() {
                        let m = &interval.metadata;
                        let unified_query_id = match l2u.get(m.query_id() as usize) {
                            Some(&id) if id != u32::MAX => id,
                            _ => continue,
                        };
                        results.push(RawAlignmentInterval {
                            target_start: interval.first,
                            target_end: interval.last,
                            query_id: unified_query_id,
                            query_start: m.query_start(),
                            query_end: m.query_end(),
                            is_reverse: m.is_reverse_strand(),
                        });
                    }
                    // Tree Arc dropped here — tree caching is disabled on the
                    // transient Impg, so the underlying COITree is freed immediately.
                }
            }
            // Impg Arc dropped here — the entire sub-index is released.
        }
        results
    }

    /// Transient variant of `query_raw_overlapping`. Same memory-bounding
    /// rationale as `query_raw_intervals_transient`, using a coitree range query
    /// instead of iterating every interval.
    fn query_raw_overlapping_transient(
        &self,
        unified_target_id: u32,
        start: i64,
        end: i64,
    ) -> Vec<RawAlignmentInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => return Vec::new(),
        };

        let mut by_index: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
        for loc in locations {
            by_index
                .entry(loc.index_idx())
                .or_default()
                .push(loc.local_target_id());
        }

        let mut results = Vec::new();
        for (index_idx, local_target_ids) in by_index {
            let impg = match self.load_sub_index_transient(index_idx) {
                Ok(i) => i,
                Err(e) => {
                    warn!(
                        "query_raw_overlapping_transient: failed to load sub-index {:?}: {}",
                        self.index_paths[index_idx], e
                    );
                    continue;
                }
            };
            let l2u = &self.local_to_unified[index_idx];
            for local_target_id in local_target_ids {
                if let Some(tree) = impg.get_or_load_tree(local_target_id) {
                    tree.query(start, end, |interval| {
                        let m = &interval.metadata;
                        if let Some(&unified_query_id) = l2u.get(m.query_id() as usize) {
                            if unified_query_id != u32::MAX {
                                results.push(RawAlignmentInterval {
                                    target_start: interval.first,
                                    target_end: interval.last,
                                    query_id: unified_query_id,
                                    query_start: m.query_start(),
                                    query_end: m.query_end(),
                                    is_reverse: m.is_reverse_strand(),
                                });
                            }
                        }
                    });
                }
            }
        }
        results
    }

    /// Batch variant: groups all queries by sub-index file, then drives the
    /// per-file work in parallel via rayon. Each rayon worker loads one
    /// sub-index transiently, answers every query referencing it, and frees
    /// the sub-index before moving on. Peak retained memory is bounded to
    /// `rayon::current_num_threads()` sub-indices.
    ///
    /// History: previously processed files sequentially to keep peak at one
    /// sub-index. That serialised Phase 1 transitive depth to a single thread
    /// (`top -H` showed 1 worker running, 47 sleeping for 30+ minutes on a
    /// 200-PAF subset). Per-thread transient loads keep memory bounded while
    /// restoring the original goal of the batch design — "answer many queries
    /// per file load" — to actually run in parallel.
    fn batch_query_raw_overlapping(
        &self,
        queries: &[(u32, i64, i64)],
    ) -> io::Result<Vec<Vec<RawAlignmentInterval>>> {
        let n = queries.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        // Group: file_idx → Vec<(query_idx, local_target_id, start, end)>
        let mut by_file: FxHashMap<usize, Vec<(usize, u32, i64, i64)>> = FxHashMap::default();
        for (qi, &(unified_target_id, start, end)) in queries.iter().enumerate() {
            if let Some(locs) = self.forest_map.get(&unified_target_id) {
                for loc in locs {
                    by_file.entry(loc.index_idx()).or_default().push((
                        qi,
                        loc.local_target_id(),
                        start,
                        end,
                    ));
                }
            }
        }

        // Per-slot Mutex: when a unified target lives in multiple sub-index
        // files, two parallel workers may both push into the same `results[qi]`.
        // Lock contention is low because per-target file fan-out is typically 1
        // (per-PAF indices) and the critical section is a `Vec::extend` of a
        // small thread-local buffer.
        let results: Vec<PlMutex<Vec<RawAlignmentInterval>>> =
            (0..n).map(|_| PlMutex::new(Vec::new())).collect();

        // Sort file_order for deterministic scheduling (rayon may still steal
        // out of order). This keeps progress observable and makes reproduction
        // easier in case of regressions.
        let mut file_order: Vec<usize> = by_file.keys().copied().collect();
        file_order.sort_unstable_by_key(|&i| (std::cmp::Reverse(self.index_sizes[i]), i));
        let limiter = file_query_limiter();

        file_order
            .par_iter()
            .try_for_each(|&file_idx| -> io::Result<()> {
                let _permit = limiter.acquire(estimated_resident_bytes(self.index_sizes[file_idx]));
                let file_queries = &by_file[&file_idx];
                let impg = self.load_sub_index_transient(file_idx).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "batch raw depth query failed to load '{}': {e}",
                            self.index_paths[file_idx].display()
                        ),
                    )
                })?;
                let l2u = &self.local_to_unified[file_idx];

                // Per-query thread-local scratch buffer keeps the global
                // results[qi] mutex critical section to a single `extend`.
                let mut local: Vec<RawAlignmentInterval> = Vec::new();
                for &(qi, local_target_id, start, end) in file_queries {
                    if let Some(tree) = impg.get_or_load_tree(local_target_id) {
                        tree.query(start, end, |interval| {
                            let m = &interval.metadata;
                            if let Some(&unified_query_id) = l2u.get(m.query_id() as usize) {
                                if unified_query_id != u32::MAX {
                                    local.push(RawAlignmentInterval {
                                        target_start: interval.first,
                                        target_end: interval.last,
                                        query_id: unified_query_id,
                                        query_start: m.query_start(),
                                        query_end: m.query_end(),
                                        is_reverse: m.is_reverse_strand(),
                                    });
                                }
                            }
                        });
                    }
                    if !local.is_empty() {
                        results[qi].lock().append(&mut local);
                    }
                }
                // impg dropped here — sub-index freed immediately.
                Ok(())
            })?;

        Ok(results.into_iter().map(|m| m.into_inner()).collect())
    }

    /// CIGAR-carrying batch projected query. Same file-grouped, transient,
    /// load-once-per-file structure as `batch_query_raw_overlapping`, but runs
    /// the full single-hop projection (`Impg::query`, with CIGAR when
    /// `store_cigar`) on each loaded sub-index and translates results to unified
    /// IDs. Produces, per query, the same projected overlaps `Impg::query` would
    /// (minus the self-interval, which the BFS driver adds once), so a
    /// level-synchronised transitive driver can bound peak resident sub-indices
    /// to the file-parallel working set instead of the cached
    /// `T × working-set` of the per-chunk `query_all_indices` path.
    fn batch_query_overlapping_with_cigar(
        &self,
        queries: &[(u32, i64, i64)],
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> io::Result<Vec<Vec<AdjustedInterval>>> {
        let n = queries.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if n > u32::MAX as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CIGAR batch contains more than u32::MAX queries",
            ));
        }

        // Group: file_idx → Vec<(query_idx, local_target_id, start, end)>
        // query_idx is u32 to keep each tuple at 24 B instead of 32 B on
        // 64-bit targets; batches are explicitly bounded above.
        let mut by_file: FxHashMap<usize, Vec<(u32, u32, i64, i64)>> = FxHashMap::default();
        for (qi, &(unified_target_id, start, end)) in queries.iter().enumerate() {
            if let Some(locs) = self.forest_map.get(&unified_target_id) {
                for loc in locs {
                    by_file.entry(loc.index_idx()).or_default().push((
                        qi as u32,
                        loc.local_target_id(),
                        start,
                        end,
                    ));
                }
            }
        }

        // Retain file provenance while workers run so the final merge can
        // reproduce `query_all_indices` ordering. Visited-range acceptance is
        // order-sensitive; appending in rayon completion order made identical
        // graphs traverse different paths from run to run.
        let results: Vec<PlMutex<Vec<(u32, AdjustedInterval)>>> =
            (0..n).map(|_| PlMutex::new(Vec::new())).collect();

        let mut file_order: Vec<usize> = by_file.keys().copied().collect();
        file_order.sort_unstable_by_key(|&i| (std::cmp::Reverse(self.index_sizes[i]), i));
        let limiter = file_query_limiter();

        file_order
            .par_iter()
            .try_for_each(|&file_idx| -> io::Result<()> {
                let _permit = limiter.acquire(estimated_resident_bytes(self.index_sizes[file_idx]));
                let file_queries = &by_file[&file_idx];
                let impg = self.load_sub_index_transient(file_idx).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "batch CIGAR depth query failed to load '{}': {e}",
                            self.index_paths[file_idx].display()
                        ),
                    )
                })?;

                for &(qi, local_target_id, start, end) in file_queries {
                    let local_results = impg.query_without_self(
                        local_target_id,
                        start,
                        end,
                        store_cigar,
                        min_gap_compressed_identity,
                        sequence_index,
                        approximate_mode,
                    );
                    // Translate to unified IDs, dropping the self-interval: the
                    // driver owns self-interval handling (matches query_all_indices
                    // / batch_depth_bfs, where the self row is added once per hop).
                    let mut local: Vec<AdjustedInterval> = Vec::with_capacity(local_results.len());
                    for r in local_results {
                        if let Some(t) = self.translate_to_unified(r, file_idx) {
                            local.push(t);
                        }
                    }
                    if !local.is_empty() {
                        results[qi as usize]
                            .lock()
                            .extend(local.drain(..).map(|result| (file_idx as u32, result)));
                    }
                }
                // impg dropped here — sub-index freed immediately.
                Ok(())
            })?;

        Ok(results
            .into_iter()
            .map(|slot| {
                let mut tagged = slot.into_inner();
                // Stable sort preserves each sub-index tree's local query
                // order while restoring the original file-list order.
                tagged.sort_by_key(|(file_idx, _)| *file_idx);
                tagged.into_iter().map(|(_, result)| result).collect()
            })
            .collect())
    }
}

impl MultiImpg {
    /// Number of sub-indices currently resident in the shared `sub_indices`
    /// cache. Used by regression tests to assert that the transient query
    /// variants do not populate the cache, and by occasional diagnostics to
    /// report working-set size during long-running depth runs.
    ///
    /// **Do not call from hot paths.** This scans the entire `sub_indices`
    /// vec under a read lock and is O(num_alignment_files). At per-file scale
    /// (≥10⁴ files) the scan itself becomes non-trivial; one call per Phase 1
    /// iteration across 64 rayon workers would serialize everyone on the
    /// read lock and contend with the slow-path `get_sub_index` writers.
    pub fn loaded_sub_index_count(&self) -> usize {
        self.sub_indices
            .read()
            .unwrap()
            .iter()
            .filter(|s| s.is_some())
            .count()
    }

    /// Current soft cap on resident BFS/transitive sub-indices (`0` = unbounded).
    /// Exposed for tests of the adaptive [`MultiImpg::set_sub_index_cache_limit`].
    pub fn sub_index_cache_limit(&self) -> usize {
        self.sub_index_cache_limit.load(AtomicOrdering::Relaxed)
    }

    /// Current sub-index cache byte budget (estimated resident bytes; `0` =
    /// unbounded). Exposed for tests of [`MultiImpg::set_sub_index_cache_byte_budget`].
    pub fn sub_index_cache_byte_budget(&self) -> u64 {
        self.sub_index_cache_byte_budget
            .load(AtomicOrdering::Relaxed)
    }

    /// Sum of estimated resident bytes for currently-populated sub-index slots.
    /// Exposed for tests asserting the byte budget bounds residency.
    pub fn sub_index_cache_bytes(&self) -> u64 {
        self.sub_index_cache_bytes.load(AtomicOrdering::Relaxed)
    }

    /// Internal implementation of transitive queries.
    ///
    /// Matches the behavior of `Impg::query_transitive_dfs` and `Impg::query_transitive_bfs`,
    /// including min_distance_between_ranges checks and stack merging.
    fn transitive_query_impl(
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
        use_dfs: bool,
    ) -> io::Result<Vec<AdjustedInterval>> {
        // Initialize visited ranges.  Global depth --use-BFS calls this once
        // per 5 MB chunk; without a mask, prebuilding an entry for every
        // unified sequence is fixed O(num_sequences) work per chunk.  Keep
        // masked runs exact, but allocate the unmasked map lazily.
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, v.clone())).collect()
        } else {
            FxHashMap::default()
        };

        // Filter input range
        let target_len = self.seq_index.get_len_from_id(target_id).unwrap_or(0) as i64;
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_insert_with(|| SortedRanges::new(target_len, 0))
            .insert((range_start, range_end));

        let carry = store_cigar && !approximate_mode;
        let mut results = Vec::new();
        let mut initial_ranges = Vec::new();

        // Add filtered input ranges
        for (filtered_start, filtered_end) in filtered_input_range {
            results.push(self.make_self_interval(
                target_id,
                filtered_start,
                filtered_end,
                store_cigar,
            ));

            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                initial_ranges.push(TransitiveRange {
                    seq_id: target_id,
                    start: filtered_start,
                    end: filtered_end,
                    hub_to_anchor_cigar: if carry {
                        Arc::new(CigarOp::new_run((filtered_end - filtered_start).abs(), '='))
                    } else {
                        Arc::new(Vec::new())
                    },
                    anchor_strand: Strand::Forward,
                    anchor_id: target_id,
                    anchor_span: (
                        filtered_start.min(filtered_end),
                        filtered_start.max(filtered_end),
                    ),
                });
            }
        }

        if use_dfs {
            let mut stack: Vec<(TransitiveRange, u16)> =
                initial_ranges.into_iter().map(|range| (range, 0)).collect();
            while let Some((tr, current_depth)) = stack.pop() {
                if max_depth > 0 && current_depth >= max_depth {
                    continue;
                }
                let can_descend = max_depth == 0 || current_depth < max_depth.saturating_sub(1);
                let step_results = self.query_all_indices(
                    tr.seq_id,
                    tr.start,
                    tr.end,
                    store_cigar,
                    min_gap_compressed_identity,
                    sequence_index,
                    approximate_mode,
                )?;

                for result in step_results {
                    if !transitive_result_allowed(
                        &self.seq_index,
                        result.0.metadata,
                        target_id,
                        subset_filter,
                    ) {
                        continue;
                    }
                    let Some(hit) = transitive_hit_from_pairwise(&tr, result, carry, can_descend)
                    else {
                        continue;
                    };
                    if can_descend {
                        extend_multi_transitive_hit(
                            &self.seq_index,
                            &mut visited_ranges,
                            &hit,
                            min_distance_between_ranges,
                            min_transitive_len,
                            |range| stack.push((range, current_depth + 1)),
                        );
                    }
                    push_transitive_hit(&mut results, hit, min_output_length);
                }

                // Legacy non-CIGAR DFS merged pending ranges across parents.
                // A carried CIGAR belongs to one parent path and must never be
                // merged with another path.
                if !carry && !stack.is_empty() {
                    stack.sort_by_key(|(tr, _)| (tr.seq_id, tr.start));
                    let mut write = 0;
                    for read in 1..stack.len() {
                        if stack[write].0.seq_id == stack[read].0.seq_id
                            && stack[write].0.end >= stack[read].0.start
                        {
                            stack[write].0.end = stack[write].0.end.max(stack[read].0.end);
                        } else {
                            write += 1;
                            stack.swap(write, read);
                        }
                    }
                    stack.truncate(write + 1);
                }
            }
            return Ok(results);
        }

        // BFS is level-synchronised. All frontier regions at one depth are
        // grouped by alignment file in `batch_query_overlapping_with_cigar`,
        // so each file is loaded once for the whole level rather than once per
        // 5 Mb anchor chunk/frontier item.
        let mut current_depth = 0u16;
        let mut current_ranges = initial_ranges;
        while !current_ranges.is_empty() && (max_depth == 0 || current_depth < max_depth) {
            let can_descend = max_depth == 0 || current_depth < max_depth.saturating_sub(1);
            let queries: Vec<(u32, i64, i64)> = current_ranges
                .iter()
                .map(|tr| (tr.seq_id, tr.start, tr.end))
                .collect();
            let query_results = self.batch_query_overlapping_with_cigar(
                &queries,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            )?;
            debug_assert_eq!(query_results.len(), current_ranges.len());

            let mut next_depth_ranges = Vec::new();
            for (tr, step_results) in current_ranges.iter().zip(query_results) {
                for result in step_results {
                    if !transitive_result_allowed(
                        &self.seq_index,
                        result.0.metadata,
                        target_id,
                        subset_filter,
                    ) {
                        continue;
                    }
                    let Some(hit) = transitive_hit_from_pairwise(tr, result, carry, can_descend)
                    else {
                        continue;
                    };
                    if can_descend {
                        extend_multi_transitive_hit(
                            &self.seq_index,
                            &mut visited_ranges,
                            &hit,
                            min_distance_between_ranges,
                            min_transitive_len,
                            |range| next_depth_ranges.push(range),
                        );
                    }
                    push_transitive_hit(&mut results, hit, min_output_length);
                }
            }

            current_depth += 1;
            if !carry && !next_depth_ranges.is_empty() {
                next_depth_ranges.par_sort_by_key(|tr| (tr.seq_id, tr.start));
                let mut write = 0;
                for read in 1..next_depth_ranges.len() {
                    if next_depth_ranges[write].seq_id == next_depth_ranges[read].seq_id
                        && next_depth_ranges[write].end >= next_depth_ranges[read].start
                    {
                        next_depth_ranges[write].end = next_depth_ranges[write]
                            .end
                            .max(next_depth_ranges[read].end);
                    } else {
                        write += 1;
                        next_depth_ranges.swap(write, read);
                    }
                }
                next_depth_ranges.truncate(write + 1);
            }
            current_ranges = next_depth_ranges;
        }

        Ok(results)
    }
}

fn transitive_result_allowed(
    seq_index: &SequenceIndex,
    query_id: u32,
    anchor_id: u32,
    subset_filter: Option<&SubsetFilter>,
) -> bool {
    subset_filter.is_none_or(|filter| {
        query_id == anchor_id
            || seq_index
                .get_name(query_id)
                .is_some_and(|name| filter.matches(name))
    })
}

fn transitive_hit_from_pairwise(
    tr: &TransitiveRange,
    result: AdjustedInterval,
    carry: bool,
    can_descend: bool,
) -> Option<BfsHit> {
    let (query_interval, cigar_ops, target_interval) = result;
    let query_id = query_interval.metadata;
    if query_id == tr.seq_id {
        return None;
    }

    if carry {
        let strand_bc = if query_interval.first <= query_interval.last {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        let composed = compose_hop(
            tr,
            target_interval.first,
            target_interval.last,
            &cigar_ops,
            strand_bc,
        )?;
        let c_lo = query_interval.first.min(query_interval.last);
        let c_hi = query_interval.first.max(query_interval.last);
        let (query_first, query_last) = if composed.strand_ac == Strand::Forward {
            (c_lo, c_hi)
        } else {
            (c_hi, c_lo)
        };
        Some(BfsHit {
            query_id,
            query_first,
            query_last,
            result_cigar: composed.a_to_c,
            t_first: composed.anchor_lo,
            t_last: composed.anchor_hi,
            t_id: tr.anchor_id,
            parent_id: tr.seq_id,
            next_carry: can_descend.then_some(NextCarry {
                strand_ac: composed.strand_ac,
                anchor_id: tr.anchor_id,
                c_lo,
                c_hi,
                anchor_lo: composed.anchor_lo,
                anchor_hi: composed.anchor_hi,
            }),
        })
    } else {
        Some(BfsHit {
            query_id,
            query_first: query_interval.first,
            query_last: query_interval.last,
            result_cigar: cigar_ops,
            t_first: target_interval.first,
            t_last: target_interval.last,
            t_id: tr.seq_id,
            parent_id: tr.seq_id,
            next_carry: None,
        })
    }
}

fn extend_multi_transitive_hit(
    seq_index: &SequenceIndex,
    visited_ranges: &mut FxHashMap<u32, SortedRanges>,
    hit: &BfsHit,
    min_distance_between_ranges: i64,
    min_transitive_len: i64,
    push: impl FnMut(TransitiveRange),
) {
    if hit.query_id == hit.parent_id {
        return;
    }
    let ranges = visited_ranges.entry(hit.query_id).or_insert_with(|| {
        SortedRanges::new(
            seq_index.get_len_from_id(hit.query_id).unwrap_or(0) as i64,
            0,
        )
    });
    if min_distance_between_ranges > 0 {
        let new_min = hit.query_first.min(hit.query_last);
        let new_max = hit.query_first.max(hit.query_last);
        let idx = ranges
            .ranges
            .binary_search_by_key(&new_min, |&(start, _)| start)
            .unwrap_or_else(|i| i);
        if idx > 0 && (new_min - ranges.ranges[idx - 1].1).abs() < min_distance_between_ranges {
            return;
        }
        if idx < ranges.ranges.len()
            && (ranges.ranges[idx].0 - new_max).abs() < min_distance_between_ranges
        {
            return;
        }
    }
    let new_ranges = ranges.insert((hit.query_first, hit.query_last));
    extend_frontier_from_hit(hit, new_ranges, min_transitive_len, push);
}

fn push_transitive_hit(
    results: &mut Vec<AdjustedInterval>,
    hit: BfsHit,
    min_output_length: Option<i64>,
) {
    if min_output_length.is_some_and(|min_len| (hit.query_last - hit.query_first).abs() < min_len) {
        return;
    }
    results.push((
        Interval {
            first: hit.query_first,
            last: hit.query_last,
            metadata: hit.query_id,
        },
        hit.result_cigar,
        Interval {
            first: hit.t_first,
            last: hit.t_last,
            metadata: hit.t_id,
        },
    ));
}

impl MultiImpgCache {
    /// Load a cache from disk.
    pub fn load(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let cache: MultiImpgCache =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to decode cache: {e}"),
                    )
                })?;

        // Verify magic and version
        if &cache.magic != CACHE_MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid cache magic bytes",
            ));
        }
        if cache.version != CACHE_VERSION {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Unsupported cache version: {} (expected {})",
                    cache.version, CACHE_VERSION
                ),
            ));
        }

        Ok(cache)
    }

    /// Check if the cache is still valid (not stale).
    ///
    /// The cache is valid if:
    /// 1. All listed index files still exist
    /// 2. No index file has been modified (mtime + size check)
    /// 3. The list of files matches exactly
    pub fn is_valid(&self, index_paths: &[PathBuf], _list_file: &Path) -> std::io::Result<bool> {
        // Check file list length matches
        if self.manifest.len() != index_paths.len() {
            debug!(
                "Cache invalid: manifest has {} files, but {} index paths provided",
                self.manifest.len(),
                index_paths.len()
            );
            return Ok(false);
        }

        // Check each file's path, mtime, and size
        for (entry, path) in self.manifest.iter().zip(index_paths) {
            let path_str = path.to_string_lossy().to_string();
            if entry.path != path_str {
                debug!(
                    "Cache invalid: path mismatch '{}' vs '{}'",
                    entry.path, path_str
                );
                return Ok(false);
            }

            let metadata = match fs::metadata(path) {
                Ok(m) => m,
                Err(_) => {
                    debug!("Cache invalid: file not found '{}'", path_str);
                    return Ok(false);
                }
            };

            if metadata.len() != entry.size {
                debug!(
                    "Cache invalid: size mismatch for '{}' ({} vs {})",
                    path_str,
                    metadata.len(),
                    entry.size
                );
                return Ok(false);
            }

            let mtime = metadata
                .modified()?
                .duration_since(SystemTime::UNIX_EPOCH)
                .map_err(std::io::Error::other)?;

            if mtime.as_secs() != entry.mtime_secs {
                debug!(
                    "Cache invalid: mtime mismatch for '{}' ({} vs {})",
                    path_str,
                    mtime.as_secs(),
                    entry.mtime_secs
                );
                return Ok(false);
            }
        }

        Ok(true)
    }
}
