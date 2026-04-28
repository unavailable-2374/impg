// src/multi_impg.rs
//! Multi-file IMPG index implementation.
//!
//! `MultiImpg` coordinates queries across multiple per-file `.impg` indices,
//! presenting a unified view while internally managing ID translation.

use crate::forest_map::ForestMap;
use crate::impg::{AdjustedInterval, CigarOp, Impg, QueryMetadata, SortedRanges};
use crate::impg_index::{ImpgIndex, RawAlignmentInterval};
use crate::seqidx::SequenceIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use log::{debug, info, warn};
use parking_lot::Mutex as PlMutex;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, RwLock};
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

impl MultiImpg {
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

        Ok(Self {
            seq_index: unified_seq_index,
            forest_map: unified_forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            transient_header_cache: (0..num_indices).map(|_| PlMutex::new(None)).collect(),
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

        Ok(Self {
            seq_index: cache.unified_seq_index,
            forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified: cache.local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            transient_header_cache: (0..num_indices).map(|_| PlMutex::new(None)).collect(),
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
        let cache_enabled = self.tree_cache_enabled.load(std::sync::atomic::Ordering::Relaxed);
        impg.set_tree_cache_enabled(cache_enabled);

        // Store in sub-index cache
        {
            let mut indices = self.sub_indices.write().unwrap();
            indices[index_idx] = Some(Arc::clone(&impg));
        }

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

        // Cache miss. If the cache is already at its soft cap, evict every
        // slot before populating the new one. Doing this *before* the
        // expensive load/parse means we never temporarily hold N+1 entries
        // and never spike `vm.max_map_count` past the budget that produced
        // the cap in the first place.
        //
        // The eviction is "drop everything", not LRU: with chunked depth
        // workloads the per-thread access pattern over `index_idx` is
        // approximately uniform within a phase, so any per-slot priority
        // would have to be paid on every hit (expensive) for limited gain;
        // a periodic flush bounds memory at the cost of some extra header
        // re-parses, which empirically is a small fraction of the
        // chunk-processing time.
        if self.transient_cache_limit > 0
            && self.transient_cache_count.load(AtomicOrdering::Relaxed) >= self.transient_cache_limit
        {
            self.evict_transient_header_cache();
        }

        let arc = self.load_sub_index_uncached(index_idx)?;

        // Stash for next time. Race-tolerant: if another thread populated the
        // slot first, we drop our copy and use theirs (functionally identical).
        let mut slot = self.transient_header_cache[index_idx].lock();
        match *slot {
            Some(ref existing) => Ok(Arc::clone(existing)),
            None => {
                *slot = Some(Arc::clone(&arc));
                self.transient_cache_count
                    .fetch_add(1, AtomicOrdering::Relaxed);
                Ok(arc)
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
    ) -> Vec<AdjustedInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => {
                return vec![self.make_self_interval(
                    unified_target_id,
                    range_start,
                    range_end,
                    store_cigar,
                )]
            }
        };

        // Query all relevant sub-indices in parallel
        let results: Vec<Vec<AdjustedInterval>> = locations
            .par_iter()
            .filter_map(|loc| {
                let impg = match self.get_sub_index(loc.index_idx()) {
                    Ok(i) => i,
                    Err(e) => {
                        warn!("Failed to load sub-index {}: {}", loc.index_idx(), e);
                        return None;
                    }
                };

                // Query using local target ID
                let local_results = impg.query(
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

                Some(unified_results)
            })
            .collect();

        // Merge all results, but only include ONE self-interval
        let mut final_results = Vec::new();
        let mut seen_self = false;

        for result_set in results {
            for result in result_set {
                // Check if this is a self-interval (query == target == input)
                let is_self = result.0.metadata == unified_target_id
                    && result.2.metadata == unified_target_id
                    && result.0.first == range_start
                    && result.0.last == range_end;

                if is_self {
                    if !seen_self {
                        final_results.push(result);
                        seen_self = true;
                    }
                    // Skip duplicate self-intervals
                } else {
                    final_results.push(result);
                }
            }
        }

        // If we never saw a self-interval, add one
        if !seen_self {
            final_results.insert(
                0,
                self.make_self_interval(unified_target_id, range_start, range_end, store_cigar),
            );
        }

        // Sort results for deterministic ordering (excluding the self-interval which should stay first)
        // Sort by: query_id, query_start, query_end, target_start, target_end
        if final_results.len() > 1 {
            let self_interval = final_results.remove(0);
            final_results.sort_by(|a, b| {
                let a_key = (a.0.metadata, a.0.first, a.0.last, a.2.first, a.2.last);
                let b_key = (b.0.metadata, b.0.first, b.0.last, b.2.first, b.2.last);
                a_key.cmp(&b_key)
            });
            final_results.insert(0, self_interval);
        }

        final_results
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
    ) -> Vec<AdjustedInterval> {
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

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i64,
        range_end: i64,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        _cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
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
        _cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
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
    ) -> Vec<AdjustedInterval> {
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
    ) -> Vec<AdjustedInterval> {
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
                            results.push((query_start, query_end, target_start, target_end, target_id));
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
        self.tree_cache_enabled.store(enabled, std::sync::atomic::Ordering::Relaxed);
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
    ) -> Vec<u16> {
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

        (0..self.index_paths.len())
            .into_par_iter()
            .for_each(|index_idx| {
                // Pre-scan visits each file exactly once. Use the no-cache
                // loader so the transient header cache doesn't accumulate
                // hundreds of thousands of `Arc<Impg>` headers (each with its
                // own glibc-mmap'd allocations) and trip the kernel's
                // vm.max_map_count limit on hosts with many per-file indices.
                let impg = match self.load_sub_index_uncached(index_idx) {
                    Ok(i) => i,
                    Err(e) => {
                        warn!(
                            "Degree pre-scan: failed to load sub-index {:?}: {}",
                            self.index_paths[index_idx], e
                        );
                        return;
                    }
                };
                let l2u = &self.local_to_unified[index_idx];

                // Iterate every local target in this file.
                let local_target_ids: Vec<u32> =
                    impg.forest_map.entries.keys().copied().collect();

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
                            let slot = unsafe {
                                local_seen.get_unchecked_mut(query_sample as usize)
                            };
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
            });

        // Final pass: popcount each row to recover the per-target degree.
        if chunks_per_row == 0 {
            return vec![0u16; num_unified];
        }
        (0..num_unified)
            .into_par_iter()
            .map(|t| {
                let row_base = t * chunks_per_row;
                let mut count: u32 = 0;
                for c in 0..chunks_per_row {
                    count += bitset[row_base + c].load(Ordering::Relaxed).count_ones();
                }
                count.min(u16::MAX as u32) as u16
            })
            .collect()
    }

    fn query_raw_overlapping(&self, unified_target_id: u32, start: i64, end: i64) -> Vec<RawAlignmentInterval> {
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
    ) -> Vec<Vec<RawAlignmentInterval>> {
        let n = queries.len();
        if n == 0 {
            return Vec::new();
        }

        // Group: file_idx → Vec<(query_idx, local_target_id, start, end)>
        let mut by_file: FxHashMap<usize, Vec<(usize, u32, i64, i64)>> = FxHashMap::default();
        for (qi, &(unified_target_id, start, end)) in queries.iter().enumerate() {
            if let Some(locs) = self.forest_map.get(&unified_target_id) {
                for loc in locs {
                    by_file
                        .entry(loc.index_idx())
                        .or_default()
                        .push((qi, loc.local_target_id(), start, end));
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
        file_order.sort_unstable();

        file_order.par_iter().for_each(|&file_idx| {
            let file_queries = &by_file[&file_idx];
            let impg = match self.load_sub_index_transient(file_idx) {
                Ok(i) => i,
                Err(e) => {
                    warn!(
                        "batch_query_raw_overlapping: failed to load {:?}: {}",
                        self.index_paths[file_idx], e
                    );
                    return;
                }
            };
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
        });

        results.into_iter().map(|m| m.into_inner()).collect()
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
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, v.clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter()
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap_or(0);
                    (id, SortedRanges::new(len as i64, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();
        let mut stack: VecDeque<(u32, i64, i64, u16)> = VecDeque::new();

        // Add filtered input ranges
        for (filtered_start, filtered_end) in filtered_input_range {
            results.push(self.make_self_interval(
                target_id,
                filtered_start,
                filtered_end,
                store_cigar,
            ));

            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                stack.push_back((target_id, filtered_start, filtered_end, 0));
            }
        }

        while let Some((current_target_id, current_start, current_end, current_depth)) = if use_dfs
        {
            stack.pop_back() // DFS: pop from back (LIFO) - O(1)
        } else {
            stack.pop_front() // BFS: pop from front (FIFO) - O(1) with VecDeque
        } {
            if max_depth > 0 && current_depth >= max_depth {
                continue;
            }

            debug!(
                "Transitive query: {}:{}-{}, depth={}",
                self.seq_index.get_name(current_target_id).unwrap_or("?"),
                current_start,
                current_end,
                current_depth
            );

            // Query all sub-indices for this region
            let step_results = self.query_all_indices(
                current_target_id,
                current_start,
                current_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            );

            // Process results - matches original Impg behavior
            for result in step_results {
                let query_id = result.0.metadata;

                // Skip self-referential results
                if query_id == current_target_id {
                    continue;
                }

                // Apply subset filter if provided (always keep original target)
                if let Some(filter) = subset_filter {
                    if query_id != target_id {
                        if let Some(name) = self.seq_index.get_name(query_id) {
                            if !filter.matches(name) {
                                continue;
                            }
                        }
                    }
                }

                // Normalize coordinates
                let (adjusted_query_start, adjusted_query_end) = if result.0.first <= result.0.last
                {
                    (result.0.first, result.0.last)
                } else {
                    (result.0.last, result.0.first)
                };

                let length = (result.0.last - result.0.first).abs();

                // Add to results only if it passes min_output_length filter
                let should_add_to_output = if let Some(min_len) = min_output_length {
                    length >= min_len
                } else {
                    true
                };
                if should_add_to_output {
                    results.push(result);
                }

                // Only add non-overlapping portions to the stack for further exploration
                let ranges = visited_ranges.entry(query_id).or_insert_with(|| {
                    let len = self.seq_index.get_len_from_id(query_id).unwrap_or(0);
                    SortedRanges::new(len as i64, 0)
                });

                let mut should_add = true;

                // Check if the range is too close to any existing ranges
                if min_distance_between_ranges > 0 {
                    let (new_min, new_max) = (adjusted_query_start, adjusted_query_end);

                    // Find insertion point in sorted ranges
                    let idx = match ranges
                        .ranges
                        .binary_search_by_key(&new_min, |&(start, _)| start)
                    {
                        Ok(i) => i,
                        Err(i) => i,
                    };

                    // Only need to check adjacent ranges due to sorting
                    if idx > 0 {
                        // Check previous range
                        let (_, prev_end) = ranges.ranges[idx - 1];
                        if (new_min - prev_end).abs() < min_distance_between_ranges {
                            should_add = false;
                        }
                    }
                    if idx < ranges.ranges.len() {
                        // Check next range
                        let (next_start, _) = ranges.ranges[idx];
                        if (next_start - new_max).abs() < min_distance_between_ranges {
                            should_add = false;
                        }
                    }
                }

                if should_add {
                    let new_ranges = ranges.insert((adjusted_query_start, adjusted_query_end));

                    // Add non-overlapping portions to stack
                    for (new_start, new_end) in new_ranges {
                        if (new_end - new_start).abs() >= min_transitive_len {
                            stack.push_back((query_id, new_start, new_end, current_depth + 1));
                        }
                    }
                }
            }

            // Merge contiguous/overlapping ranges with same sequence_id (matches original)
            let slice = stack.make_contiguous();
            slice.sort_by_key(|(id, start, _, _)| (*id, *start));

            let mut write = 0;
            for read in 1..slice.len() {
                if slice[write].0 == slice[read].0 &&   // Same sequence_id
                    slice[write].2 >= slice[read].1
                // Overlapping or contiguous
                {
                    // Merge by extending end
                    slice[write].2 = slice[write].2.max(slice[read].2);
                } else {
                    write += 1;
                    slice.swap(write, read);
                }
            }
            if !stack.is_empty() {
                stack.truncate(write + 1);
            }
        }

        results
    }
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
