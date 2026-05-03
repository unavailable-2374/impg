use crate::impg::{CigarOp, SortedRanges};
use crate::impg_index::{ImpgIndex, RawAlignmentInterval};
use crate::sequence_index::UnifiedSequenceIndex;
use bitvec::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info};
use parking_lot::Mutex;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeMap;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread::JoinHandle;

use crate::commands::depth_checkpoint::{
    compute_invalidation_hash, encode_chunk_id_phase1, encode_chunk_id_phase2,
    encode_work_record, encode_worklog_header, replay_work_log, truncate_to,
    DepthCheckpoint, HashInputs, ResumeState, CHUNK_ID_BOUNDARY, CHUNK_ID_FAI_DONE, TSV_SUFFIX,
    WORKLOG_SUFFIX,
};

// ============================================================================
// `depth-trace` cargo feature — diagnostic eprintln! probes for triaging
// `impg depth` ↔ `impg query -x` divergence (PLAN_depth_100pct.md §3.3).
//
// Off by default: the macro expands to a no-op so there is zero runtime cost.
// Enable with `cargo build --features depth-trace` to dump per-stage trace
// records to stderr. Output is emitted as `[depth-trace] ...` lines.
//
// Stages emitted (greppable tags):
//   FALLBACK  project_hop0_coords fell through to (region_start, region_end)
//   PASS1     pass-1 finished building seq_anchor_coverage
//   HOP2      pass-2 hop ≥ 1 hit projected through `project_hop0_coords`
//   SWEEP     alignments handed to `sweep_line_depth` (count + unique samples)
// ============================================================================
#[cfg(feature = "depth-trace")]
macro_rules! depth_trace {
    ($($arg:tt)*) => {{
        eprintln!("[depth-trace] {}", format_args!($($arg)*));
    }};
}
#[cfg(not(feature = "depth-trace"))]
macro_rules! depth_trace {
    ($($arg:tt)*) => {};
}

/// Configuration for the depth command
pub struct DepthConfig {
    pub transitive: bool,
    pub transitive_dfs: bool,
    pub max_depth: u16,
    pub min_transitive_len: i64,
    pub min_distance_between_ranges: i64,
    pub merge_adjacent: bool,
    /// When true, use CIGAR-precise BFS for transitive depth (--use-BFS).
    /// When false (default), use raw-interval BFS with linear interpolation.
    pub use_cigar_bfs: bool,
}

// ============================================================================
// Sample filtering for depth computation
// ============================================================================

/// Filter for samples to include in depth calculation
#[derive(Debug, Clone)]
pub struct SampleFilter {
    /// Set of sample names to include (empty = include all)
    samples: FxHashSet<String>,
}

impl SampleFilter {
    /// Create a new filter with no restrictions (all samples included)
    pub fn new() -> Self {
        Self {
            samples: FxHashSet::default(),
        }
    }

    /// Create from a list of sample names
    pub fn from_samples(samples: Vec<String>) -> Self {
        Self {
            samples: samples.into_iter().collect(),
        }
    }

    /// Load from a file (one sample name per line)
    pub fn from_file(path: &str) -> io::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let samples: FxHashSet<String> = content
            .lines()
            .map(|line| line.trim())
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|s| s.to_string())
            .collect();
        Ok(Self { samples })
    }

    /// Check if a sample should be included
    pub fn includes(&self, sample: &str) -> bool {
        self.samples.is_empty() || self.samples.contains(sample)
    }

    /// Check if the filter is active (has restrictions)
    pub fn is_active(&self) -> bool {
        !self.samples.is_empty()
    }

    /// Get the number of samples in the filter
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Check if the filter is empty (no restrictions)
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Get the set of sample names
    pub fn get_samples(&self) -> &FxHashSet<String> {
        &self.samples
    }
}

impl Default for SampleFilter {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Depth statistics for stats mode
// ============================================================================

/// Statistics for depth distribution across all sequences
#[derive(Debug, Clone, Default)]
pub struct DepthStats {
    /// Total bases analyzed (anchor perspective: sum of anchor interval lengths)
    pub total_bases: i64,
    /// Distribution: depth -> total bases at that depth (anchor perspective)
    pub depth_distribution: FxHashMap<usize, i64>,
    /// Per-depth intervals for output files: depth -> Vec<(seq_name, start, end)>
    pub depth_intervals: FxHashMap<usize, Vec<(String, i64, i64)>>,
    /// Total pangenome bases: sum of all samples' actual query bases (union per sample)
    pub pangenome_total_bases: i64,
    /// Pangenome distribution: depth -> total pangenome bases at that depth
    pub pangenome_depth_distribution: FxHashMap<usize, i64>,
}

impl DepthStats {
    /// Create a new empty stats container
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an interval with a specific depth and pangenome bases
    pub fn add_interval(
        &mut self,
        seq_name: &str,
        start: i64,
        end: i64,
        depth: usize,
        pangenome_bases: i64,
    ) {
        let length = end - start;
        if length <= 0 {
            return;
        }
        self.total_bases += length;
        *self.depth_distribution.entry(depth).or_insert(0) += length;
        self.depth_intervals
            .entry(depth)
            .or_default()
            .push((seq_name.to_string(), start, end));
        self.pangenome_total_bases += pangenome_bases;
        *self.pangenome_depth_distribution.entry(depth).or_insert(0) += pangenome_bases;
    }

    /// Merge another DepthStats into this one
    pub fn merge(&mut self, other: &DepthStats) {
        self.total_bases += other.total_bases;
        for (&depth, &bases) in &other.depth_distribution {
            *self.depth_distribution.entry(depth).or_insert(0) += bases;
        }
        for (depth, intervals) in &other.depth_intervals {
            self.depth_intervals
                .entry(*depth)
                .or_default()
                .extend(intervals.iter().cloned());
        }
        self.pangenome_total_bases += other.pangenome_total_bases;
        for (&depth, &bases) in &other.pangenome_depth_distribution {
            *self.pangenome_depth_distribution.entry(depth).or_insert(0) += bases;
        }
    }

    /// Get sorted list of depths with their base counts
    pub fn get_sorted_distribution(&self) -> Vec<(usize, i64)> {
        let mut dist: Vec<_> = self
            .depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        dist.sort_by_key(|&(d, _)| d);
        dist
    }

    /// Get maximum depth
    pub fn max_depth(&self) -> usize {
        self.depth_distribution.keys().copied().max().unwrap_or(0)
    }

    /// Write summary to a writer
    pub fn write_summary<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "Anchor-perspective depth distribution:")?;
        writeln!(writer, "Total anchor bases: {}", self.total_bases)?;
        for (depth, bases) in self.get_sorted_distribution() {
            let pct = if self.total_bases > 0 {
                100.0 * bases as f64 / self.total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        writeln!(writer)?;
        writeln!(writer, "Pangenome-wide depth distribution:")?;
        writeln!(
            writer,
            "Total pangenome bases: {}",
            self.pangenome_total_bases
        )?;
        let mut pg_dist: Vec<_> = self
            .pangenome_depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        pg_dist.sort_by_key(|&(d, _)| d);
        for (depth, bases) in pg_dist {
            let pct = if self.pangenome_total_bases > 0 {
                100.0 * bases as f64 / self.pangenome_total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        Ok(())
    }

    /// Write per-depth BED files
    pub fn write_depth_bed_files(&self, prefix: &str) -> io::Result<()> {
        for (&depth, intervals) in &self.depth_intervals {
            let path = format!("{}.depth{}.bed", prefix, depth);
            let file = std::fs::File::create(&path)?;
            let mut writer = BufWriter::new(file);
            for (seq_name, start, end) in intervals {
                writeln!(writer, "{}\t{}\t{}", seq_name, start, end)?;
            }
            writer.flush()?;
            info!("Wrote {} intervals to {}", intervals.len(), path);
        }
        Ok(())
    }
}

// ============================================================================
// Combined depth output with sample lists
// ============================================================================

/// Interval with depth and list of covering samples
#[derive(Debug, Clone)]
pub struct DepthIntervalWithSamples {
    /// Sequence name
    pub seq_name: String,
    /// Start position (0-based)
    pub start: i64,
    /// End position (exclusive)
    pub end: i64,
    /// Depth (number of unique samples)
    pub depth: usize,
    /// Sorted list of sample names covering this interval
    pub samples: Vec<String>,
}

impl DepthIntervalWithSamples {
    pub fn new(seq_name: String, start: i64, end: i64, depth: usize, samples: Vec<String>) -> Self {
        Self {
            seq_name,
            start,
            end,
            depth,
            samples,
        }
    }
}

/// Statistics with sample tracking for combined output
#[derive(Debug, Clone, Default)]
pub struct DepthStatsWithSamples {
    /// Total bases analyzed (anchor perspective)
    pub total_bases: i64,
    /// Distribution: depth -> total bases at that depth (anchor perspective)
    pub depth_distribution: FxHashMap<usize, i64>,
    /// All intervals with samples (unsorted, to be sorted at output time)
    pub intervals: Vec<DepthIntervalWithSamples>,
    /// Total pangenome bases: sum of all samples' actual query bases
    pub pangenome_total_bases: i64,
    /// Pangenome distribution: depth -> total pangenome bases at that depth
    pub pangenome_depth_distribution: FxHashMap<usize, i64>,
}

impl DepthStatsWithSamples {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an interval with depth, sample list, and pangenome bases
    pub fn add_interval(
        &mut self,
        seq_name: &str,
        start: i64,
        end: i64,
        depth: usize,
        samples: Vec<String>,
        pangenome_bases: i64,
    ) {
        let length = end - start;
        if length <= 0 {
            return;
        }
        self.total_bases += length;
        *self.depth_distribution.entry(depth).or_insert(0) += length;
        self.intervals.push(DepthIntervalWithSamples::new(
            seq_name.to_string(),
            start,
            end,
            depth,
            samples,
        ));
        self.pangenome_total_bases += pangenome_bases;
        *self.pangenome_depth_distribution.entry(depth).or_insert(0) += pangenome_bases;
    }

    /// Merge another DepthStatsWithSamples into this one
    pub fn merge(&mut self, other: DepthStatsWithSamples) {
        self.total_bases += other.total_bases;
        for (&depth, &bases) in &other.depth_distribution {
            *self.depth_distribution.entry(depth).or_insert(0) += bases;
        }
        self.intervals.extend(other.intervals);
        self.pangenome_total_bases += other.pangenome_total_bases;
        for (&depth, &bases) in &other.pangenome_depth_distribution {
            *self.pangenome_depth_distribution.entry(depth).or_insert(0) += bases;
        }
    }

    /// Get sorted list of depths with their base counts
    pub fn get_sorted_distribution(&self) -> Vec<(usize, i64)> {
        let mut dist: Vec<_> = self
            .depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        dist.sort_by_key(|&(d, _)| d);
        dist
    }

    /// Get maximum depth
    pub fn max_depth(&self) -> usize {
        self.depth_distribution.keys().copied().max().unwrap_or(0)
    }

    /// Write summary to a writer
    pub fn write_summary<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "Anchor-perspective depth distribution:")?;
        writeln!(writer, "Total anchor bases: {}", self.total_bases)?;
        for (depth, bases) in self.get_sorted_distribution() {
            let pct = if self.total_bases > 0 {
                100.0 * bases as f64 / self.total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        writeln!(writer)?;
        writeln!(writer, "Pangenome-wide depth distribution:")?;
        writeln!(
            writer,
            "Total pangenome bases: {}",
            self.pangenome_total_bases
        )?;
        let mut pg_dist: Vec<_> = self
            .pangenome_depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        pg_dist.sort_by_key(|&(d, _)| d);
        for (depth, bases) in pg_dist {
            let pct = if self.pangenome_total_bases > 0 {
                100.0 * bases as f64 / self.pangenome_total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        Ok(())
    }

    /// Write combined output file sorted by chromosome and position
    /// Write combined output file sorted by chromosome and position.
    /// Intervals shorter than `min_interval_len` bp are absorbed into adjacent neighbors.
    pub fn write_combined_output(&mut self, prefix: &str, min_interval_len: i64) -> io::Result<()> {
        let path = format!("{}.combined.bed", prefix);
        let file = std::fs::File::create(&path)?;
        let mut writer = BufWriter::new(file);

        // Sort by sequence name then start position
        self.intervals.sort_by(|a, b| {
            a.seq_name
                .cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
                .then_with(|| a.end.cmp(&b.end))
        });

        writeln!(writer, "#seq_name\tlength\tdepth\tsamples")?;

        let intervals = if min_interval_len > 0 {
            merge_short_depth_intervals(std::mem::take(&mut self.intervals), min_interval_len)
        } else {
            std::mem::take(&mut self.intervals)
        };

        for interval in &intervals {
            let samples_str = interval.samples.join(",");
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                interval.seq_name,
                interval.end - interval.start,
                interval.depth,
                samples_str
            )?;
        }

        writer.flush()?;
        info!("Wrote {} intervals to {} (combined output)", intervals.len(), path);
        Ok(())
    }
}

// ============================================================================
// Region query result for region query mode
// ============================================================================

/// Result of a region query with per-sample position tracking
#[derive(Debug, Clone)]
pub struct RegionDepthResult {
    /// Reference sequence name
    pub ref_seq: String,
    /// Reference start position
    pub ref_start: i64,
    /// Reference end position
    pub ref_end: i64,
    /// Depth (number of samples covering this region)
    pub depth: usize,
    /// Per-sample positions: sample -> Vec<(seq_name, start, end)>
    /// Multiple alignments per sample are possible
    pub sample_positions: FxHashMap<String, Vec<(String, i64, i64)>>,
}

impl RegionDepthResult {
    /// Create a new result
    pub fn new(ref_seq: String, ref_start: i64, ref_end: i64) -> Self {
        Self {
            ref_seq,
            ref_start,
            ref_end,
            depth: 0,
            sample_positions: FxHashMap::default(),
        }
    }

    /// Add a sample position. No-op if the exact `(seq_name, start, end)` tuple
    /// is already present for this sample — this keeps the per-sample column
    /// from emitting `pos;pos` when several distinct raw alignments project
    /// (via linear interpolation) to coincident coordinates within the sweep
    /// window. Distinct projected ranges are still preserved.
    pub fn add_sample_position(&mut self, sample: &str, seq_name: &str, start: i64, end: i64) {
        let positions = self.sample_positions.entry(sample.to_string()).or_default();
        if positions
            .iter()
            .any(|(s, ds, de)| s == seq_name && *ds == start && *de == end)
        {
            return;
        }
        positions.push((seq_name.to_string(), start, end));
    }

    /// Update depth based on sample_positions
    pub fn update_depth(&mut self) {
        self.depth = self.sample_positions.len();
    }

    /// Format sample positions as semicolon-separated strings.
    ///
    /// Positions within a sample's cell are sorted by `(seq_name, start, end)`
    /// to give a deterministic, sweep-emission-order-independent output.
    /// Stage 4 unified the region sweep onto compact `u16` sample IDs / `u32`
    /// seq IDs, which means the in-cell order produced directly by the sweep
    /// now depends on `query_id` rather than the lexicographic `query_name`
    /// order the previous String-based path happened to produce. Sorting at
    /// the formatter restores byte-identical TSV output for downstream tools
    /// that pattern-match on column text — at the cost of one stable sort
    /// per cell, which is negligible compared to the sweep itself.
    pub fn format_sample_positions(&self, sample: &str) -> String {
        match self.sample_positions.get(sample) {
            Some(positions) if !positions.is_empty() => {
                let mut sorted: Vec<&(String, i64, i64)> = positions.iter().collect();
                sorted.sort_unstable_by(|a, b| {
                    a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2))
                });
                sorted
                    .iter()
                    .map(|(seq, start, end)| format!("{}:{}-{}", seq, start, end))
                    .collect::<Vec<_>>()
                    .join(";")
            }
            _ => "NA".to_string(),
        }
    }
}

/// Write region depth results in tabular format
pub fn write_region_depth_output<W: Write>(
    writer: &mut W,
    results: &[RegionDepthResult],
    sample_order: &[String],
) -> io::Result<()> {
    // Write header
    write!(writer, "#ref_seq\tref_start\tref_end\tdepth")?;
    for sample in sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Write data rows
    for result in results {
        write!(
            writer,
            "{}\t{}\t{}\t{}",
            result.ref_seq, result.ref_start, result.ref_end, result.depth
        )?;
        for sample in sample_order {
            write!(writer, "\t{}", result.format_sample_positions(sample))?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

// ============================================================================
// New data structures for sample-based depth computation
// ============================================================================

/// Stores sequence lengths and sample-to-sequence mappings
#[derive(Debug, Clone)]
pub(crate) struct SequenceLengths {
    /// seq_name -> length
    pub lengths: FxHashMap<String, i64>,
    /// sample_name -> list of seq_names belonging to this sample
    pub sample_to_seqs: FxHashMap<String, Vec<String>>,
}

impl SequenceLengths {
    /// Build from FAI file list
    pub fn from_fai_list(fai_list_path: &str, separator: &str) -> io::Result<Self> {
        let mut lengths = FxHashMap::default();
        let mut sample_to_seqs: FxHashMap<String, Vec<String>> = FxHashMap::default();

        let content = std::fs::read_to_string(fai_list_path).map_err(|e| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Failed to read FAI list file '{}': {}", fai_list_path, e),
            )
        })?;

        for line in content.lines() {
            let fai_path = line.trim();
            if fai_path.is_empty() || fai_path.starts_with('#') {
                continue;
            }

            let fai_content = std::fs::read_to_string(fai_path).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Failed to read FAI file '{}': {}", fai_path, e),
                )
            })?;

            for fai_line in fai_content.lines() {
                let parts: Vec<&str> = fai_line.split('\t').collect();
                if parts.len() >= 2 {
                    let seq_name = parts[0].to_string();
                    let seq_len: i64 = parts[1].parse().unwrap_or(0);

                    let sample = extract_sample(&seq_name, separator);
                    lengths.insert(seq_name.clone(), seq_len);
                    sample_to_seqs
                        .entry(sample)
                        .or_default()
                        .push(seq_name);
                }
            }
        }

        // Sort sequences within each sample for deterministic ordering
        for seqs in sample_to_seqs.values_mut() {
            seqs.sort();
        }

        Ok(SequenceLengths {
            lengths,
            sample_to_seqs,
        })
    }
}

// ============================================================================
// ID-based compact data structures for memory-efficient depth computation
// ============================================================================

/// Sample index: maps sample names to compact IDs (u16) and vice versa
/// Max 65535 samples, which is sufficient for most pangenome analyses
#[derive(Debug, Clone)]
pub(crate) struct SampleIndex {
    name_to_id: FxHashMap<String, u16>,
    id_to_name: Vec<String>, // Use Vec for O(1) lookup
}

impl SampleIndex {
    /// Create a new empty sample index
    pub fn new() -> Self {
        Self {
            name_to_id: FxHashMap::default(),
            id_to_name: Vec::new(),
        }
    }

    /// Get sample ID from name
    #[inline]
    pub fn get_id(&self, name: &str) -> Option<u16> {
        self.name_to_id.get(name).copied()
    }

    /// Get sample name from ID
    #[inline]
    pub fn get_name(&self, id: u16) -> Option<&str> {
        self.id_to_name.get(id as usize).map(|s| s.as_str())
    }

    /// Get number of samples
    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }
}

impl Default for SampleIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Compact sequence lengths using numeric IDs
/// Memory usage: O(num_sequences) instead of O(num_sequences × avg_name_length)
#[derive(Debug, Clone)]
pub(crate) struct CompactSequenceLengths {
    /// Sequence lengths indexed by seq_id (from impg's SequenceIndex)
    lengths: Vec<i64>,
    /// seq_id -> sample_id (reverse lookup)
    seq_to_sample: Vec<u16>,
    /// Sample index for name<->ID conversion
    sample_index: SampleIndex,
}

impl CompactSequenceLengths {
    /// Build from impg sequence index and separator
    pub fn from_impg(impg: &impl ImpgIndex, separator: &str) -> Self {
        // First pass: collect all samples
        let mut sample_names: FxHashSet<String> = FxHashSet::default();
        for id in 0..impg.seq_index().len() as u32 {
            if let Some(name) = impg.seq_index().get_name(id) {
                sample_names.insert(extract_sample(name, separator));
            }
        }

        // Build sample index
        let mut sorted_samples: Vec<_> = sample_names.into_iter().collect();
        sorted_samples.sort();
        if sorted_samples.len() > u16::MAX as usize {
            panic!(
                "Too many samples ({}) — the depth command supports at most {} distinct samples",
                sorted_samples.len(),
                u16::MAX
            );
        }
        let mut name_to_id = FxHashMap::default();
        let mut id_to_name = Vec::with_capacity(sorted_samples.len());
        for (i, sample) in sorted_samples.into_iter().enumerate() {
            name_to_id.insert(sample.clone(), i as u16);
            id_to_name.push(sample);
        }
        let sample_index = SampleIndex {
            name_to_id,
            id_to_name,
        };

        // Second pass: build mappings
        let num_seqs = impg.seq_index().len();
        let mut lengths = vec![0i64; num_seqs];
        let mut seq_to_sample = vec![0u16; num_seqs];

        for seq_id in 0..num_seqs as u32 {
            if let (Some(name), Some(len)) = (
                impg.seq_index().get_name(seq_id),
                impg.seq_index().get_len_from_id(seq_id),
            ) {
                let sample_name = extract_sample(name, separator);
                if let Some(sample_id) = sample_index.get_id(&sample_name) {
                    lengths[seq_id as usize] = len as i64;
                    seq_to_sample[seq_id as usize] = sample_id;
                }
            }
        }

        CompactSequenceLengths {
            lengths,
            seq_to_sample,
            sample_index,
        }
    }

    /// Get sequence length by ID
    #[inline]
    pub fn get_length(&self, seq_id: u32) -> i64 {
        self.lengths.get(seq_id as usize).copied().unwrap_or(0)
    }

    /// Get sample ID for a sequence
    #[inline]
    pub fn get_sample_id(&self, seq_id: u32) -> u16 {
        self.seq_to_sample
            .get(seq_id as usize)
            .copied()
            .unwrap_or(0)
    }

    /// Get sample index
    #[inline]
    pub fn sample_index(&self) -> &SampleIndex {
        &self.sample_index
    }

    /// Borrow the full seq_id -> sample_id slice.
    /// Used by the pre-scan file-parallel degree computation in `MultiImpg`.
    #[inline]
    pub fn seq_to_sample_slice(&self) -> &[u16] {
        &self.seq_to_sample
    }
}

/// Compact alignment info using numeric IDs
/// Memory usage: ~60% less than SampleAlignmentInfo for large datasets
#[derive(Debug, Clone)]
pub(crate) struct CompactAlignmentInfo {
    /// Sample ID (from SampleIndex)
    pub sample_id: u16,
    /// Query sequence ID (from impg's SequenceIndex)
    pub query_id: u32,
    /// Query coordinates
    pub query_start: i64,
    pub query_end: i64,
    /// Target coordinates (on anchor sequence)
    pub target_start: i64,
    pub target_end: i64,
    /// True if alignment is on reverse strand
    pub is_reverse: bool,
}

impl CompactAlignmentInfo {
    /// Create a new compact alignment info
    pub fn new(
        sample_id: u16,
        query_id: u32,
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        is_reverse: bool,
    ) -> Self {
        Self {
            sample_id,
            query_id,
            query_start,
            query_end,
            target_start,
            target_end,
            is_reverse,
        }
    }
}

/// Bitmap for tracking active samples in sweep-line algorithm
/// Memory efficient: uses 1 bit per sample instead of hash map entries
/// Supports counting multiple overlaps per sample
#[derive(Debug, Clone)]
pub(crate) struct SampleBitmap {
    /// Bit set for sample presence (sample_id -> present)
    bits: BitVec,
    /// Count of active alignments per sample (for overlapping alignments from same sample)
    counts: Vec<u32>,
    /// Number of unique active samples (cached for O(1) depth query)
    active_count: usize,
}

impl SampleBitmap {
    /// Create a new bitmap for the given number of samples
    pub fn new(num_samples: usize) -> Self {
        Self {
            bits: bitvec![0; num_samples],
            counts: vec![0; num_samples],
            active_count: 0,
        }
    }

    /// Add an alignment for a sample (increment count)
    #[inline]
    pub fn add(&mut self, sample_id: u16) {
        let idx = sample_id as usize;
        if idx < self.counts.len() {
            if self.counts[idx] == 0 {
                self.bits.set(idx, true);
                self.active_count += 1;
            }
            self.counts[idx] += 1;
        }
    }

    /// Remove an alignment for a sample (decrement count)
    #[inline]
    pub fn remove(&mut self, sample_id: u16) {
        let idx = sample_id as usize;
        if idx < self.counts.len() && self.counts[idx] > 0 {
            self.counts[idx] -= 1;
            if self.counts[idx] == 0 {
                self.bits.set(idx, false);
                self.active_count -= 1;
            }
        }
    }

    /// Get the number of unique active samples (depth)
    #[inline]
    pub fn depth(&self) -> usize {
        self.active_count
    }

    /// Iterate over active sample IDs
    pub fn active_samples(&self) -> impl Iterator<Item = u16> + '_ {
        self.bits.iter_ones().map(|idx| idx as u16)
    }
}

impl Default for SampleBitmap {
    fn default() -> Self {
        Self::new(0)
    }
}

/// High-performance interval set with efficient intersection and subtraction.
///
/// Backed by `BTreeMap<start, end>` (half-open [start, end)), invariant:
/// intervals are sorted, non-overlapping, and non-adjacent (touching intervals
/// are merged on insertion).
///
/// All mutating ops are O(log N + k) where k = number of intervals affected
/// (typically 0-2 in `ConcurrentProcessedTracker` usage). The previous sorted
/// `Vec` backend was O(N) per insert due to `Vec::insert`/`splice` memmove;
/// under VGP 581-sample workloads that pushed claim_unprocessed to >90% of
/// Phase 2 CPU as the per-sequence IntervalSet grew into the millions.
#[derive(Debug, Clone, Default)]
pub(crate) struct IntervalSet {
    /// Sorted, non-overlapping, non-adjacent intervals keyed by start
    intervals: BTreeMap<i64, i64>,
    /// Total covered length (for quick empty check)
    total_length: i64,
}

impl IntervalSet {
    /// Create a new empty interval set
    pub fn new() -> Self {
        Self {
            intervals: BTreeMap::new(),
            total_length: 0,
        }
    }

    /// Create interval set with a single interval
    pub fn new_single(start: i64, end: i64) -> Self {
        let mut s = Self::new();
        if start < end {
            s.intervals.insert(start, end);
            s.total_length = end - start;
        }
        s
    }

    /// Add an interval to this set (merging with existing if overlapping or touching).
    /// O(log N + k), where k = number of intervals merged (typically 0-2).
    pub fn add(&mut self, start: i64, end: i64) {
        if start >= end {
            return;
        }

        let mut merge_start = start;
        let mut merge_end = end;
        let mut removed_length: i64 = 0;

        // Check the interval immediately before `start`: it may reach into us
        // if prev_end >= start (touching or overlapping).
        let prev = self
            .intervals
            .range(..start)
            .next_back()
            .map(|(&k, &v)| (k, v));
        if let Some((prev_start, prev_end)) = prev {
            if prev_end >= start {
                merge_start = prev_start;
                merge_end = merge_end.max(prev_end);
                removed_length += prev_end - prev_start;
                self.intervals.remove(&prev_start);
            }
        }

        // Any interval whose start is in [start, end] is touching or overlapping.
        // Collect keys first to avoid iterator-invalidation during removal.
        let keys_to_merge: Vec<i64> = self
            .intervals
            .range(start..=end)
            .map(|(&k, _)| k)
            .collect();
        for k in keys_to_merge {
            let v = self.intervals.remove(&k).unwrap();
            merge_end = merge_end.max(v);
            removed_length += v - k;
        }

        self.intervals.insert(merge_start, merge_end);
        let added_length = merge_end - merge_start;
        self.total_length = self.total_length - removed_length + added_length;
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.total_length == 0
    }

    /// Get total covered length
    pub fn total_length(&self) -> i64 {
        self.total_length
    }

    /// Iterator over all intervals, sorted by start.
    pub fn intervals(&self) -> impl Iterator<Item = (i64, i64)> + '_ {
        self.intervals.iter().map(|(&s, &e)| (s, e))
    }

    /// Iterator over intervals that overlap `[range_start, range_end)`, in order.
    /// O(log N + k) where k is the number of overlapping intervals — much
    /// cheaper than iterating all intervals when the set is dense.
    pub fn iter_overlapping(
        &self,
        range_start: i64,
        range_end: i64,
    ) -> impl Iterator<Item = (i64, i64)> + '_ {
        // An interval (s, e) overlaps [range_start, range_end) iff
        //   s < range_end AND e > range_start.
        // Start scanning from the predecessor of range_start (its key < range_start
        // but it may reach forward), then through all keys < range_end.
        let scan_from = self
            .intervals
            .range(..range_start)
            .next_back()
            .map(|(&k, _)| k)
            .unwrap_or(i64::MIN);
        self.intervals
            .range(scan_from..range_end)
            .filter_map(move |(&s, &e)| {
                if e > range_start {
                    Some((s, e))
                } else {
                    None
                }
            })
    }

    /// Subtract a single interval from this set.
    /// O(log N + k), where k = number of intervals overlapping [sub_start, sub_end).
    pub fn subtract(&mut self, sub_start: i64, sub_end: i64) {
        if sub_start >= sub_end || self.intervals.is_empty() {
            return;
        }

        let mut removed_length: i64 = 0;
        let mut to_insert: Vec<(i64, i64)> = Vec::new();

        // Predecessor (key < sub_start): may partially overlap if prev_end > sub_start.
        let prev = self
            .intervals
            .range(..sub_start)
            .next_back()
            .map(|(&k, &v)| (k, v));
        if let Some((prev_start, prev_end)) = prev {
            if prev_end > sub_start {
                // Split / trim the predecessor.
                self.intervals.remove(&prev_start);
                let overlap_end = prev_end.min(sub_end);
                removed_length += overlap_end - sub_start;
                if prev_start < sub_start {
                    to_insert.push((prev_start, sub_start));
                }
                if prev_end > sub_end {
                    to_insert.push((sub_end, prev_end));
                }
            }
        }

        // All intervals whose start is in [sub_start, sub_end) are affected.
        // Invariant guarantees they're disjoint from the predecessor case above
        // (the predecessor can only span past sub_end if there are NO intervals
        // with start in [sub_start, sub_end), by non-overlap).
        let keys_in_range: Vec<i64> = self
            .intervals
            .range(sub_start..sub_end)
            .map(|(&k, _)| k)
            .collect();
        for k in keys_in_range {
            let v = self.intervals.remove(&k).unwrap();
            let overlap_end = v.min(sub_end);
            removed_length += overlap_end - k;
            if v > sub_end {
                to_insert.push((sub_end, v));
            }
        }

        for (s, e) in to_insert {
            self.intervals.insert(s, e);
        }
        self.total_length -= removed_length;
    }
}

/// Map a target-side offset to query-side offset using integer mul-div with
/// round-half-up semantics. Replaces `(off as f64 * (q_len/t_len)).round() as i64`.
///
/// Safety: requires `t_len > 0`. Caller's responsibility to pre-check; we do NOT
/// branch here for hot-path reasons. With `t_len <= 3e8` and `q_len <= 3e8`, the
/// product fits i64 (max ~9e16 << i64::MAX ~9.2e18). `saturating_mul` guards
/// against pathological overflow without changing results in normal ranges.
#[inline(always)]
fn proj_offset(off: i64, q_len: i64, t_len: i64) -> i64 {
    debug_assert!(t_len > 0);
    debug_assert!(off >= 0);
    debug_assert!(q_len >= 0);
    let prod = off.saturating_mul(q_len);
    (prod + (t_len >> 1)) / t_len
}

/// Map a target coordinate range to query coordinates using linear interpolation.
/// This provides smooth coordinate transitions between adjacent windows.
///
/// Note: CIGAR-precise mapping was considered but produces more discontinuities
/// because it reveals actual alignment structure differences. Linear interpolation
/// smooths these over at the cost of some coordinate accuracy.
#[allow(unused_variables)]
#[inline]
fn map_target_to_query_linear(
    cigar: &[CigarOp],
    aln_target_start: i64,
    aln_target_end: i64,
    aln_query_start: i64,
    aln_query_end: i64,
    window_target_start: i64,
    window_target_end: i64,
    is_reverse: bool,
) -> (i64, i64) {
    let target_len = aln_target_end - aln_target_start;
    let query_len = aln_query_end - aln_query_start;

    if target_len == 0 {
        return (aln_query_start, aln_query_end);
    }

    let window_target_start_clamped = window_target_start.max(aln_target_start);
    let window_target_end_clamped = window_target_end.min(aln_target_end);

    let start_offset = window_target_start_clamped - aln_target_start;
    let end_offset = window_target_end_clamped - aln_target_start;

    if is_reverse {
        let q_end = aln_query_end - proj_offset(start_offset, query_len, target_len);
        let q_start = aln_query_end - proj_offset(end_offset, query_len, target_len);
        (q_start.min(q_end), q_start.max(q_end))
    } else {
        let q_start = aln_query_start + proj_offset(start_offset, query_len, target_len);
        let q_end = aln_query_start + proj_offset(end_offset, query_len, target_len);
        (q_start.min(q_end), q_start.max(q_end))
    }
}

/// Map query sub-range [q_start, q_end] back to target coordinates using linear interpolation.
/// Inverse of the target→query mapping used in sweep construction.
#[inline]
fn inverse_map_query_to_target(
    aln_target_start: i64,
    aln_target_end: i64,
    aln_query_start: i64,
    aln_query_end: i64,
    q_start: i64,
    q_end: i64,
    is_reverse: bool,
) -> (i64, i64) {
    let target_len = aln_target_end - aln_target_start;
    let query_len = aln_query_end - aln_query_start;
    if query_len == 0 {
        return (aln_target_start, aln_target_end);
    }
    if is_reverse {
        let end_off = aln_query_end - q_start;
        let start_off = aln_query_end - q_end;
        let t_start = aln_target_start + proj_offset(start_off, target_len, query_len);
        let t_end = aln_target_start + proj_offset(end_off, target_len, query_len);
        (t_start.min(t_end).max(aln_target_start), t_start.max(t_end).min(aln_target_end))
    } else {
        let start_off = q_start - aln_query_start;
        let end_off = q_end - aln_query_start;
        let t_start = aln_target_start + proj_offset(start_off, target_len, query_len);
        let t_end = aln_target_start + proj_offset(end_off, target_len, query_len);
        (t_start.min(t_end).max(aln_target_start), t_start.max(t_end).min(aln_target_end))
    }
}

/// Extract sample from PanSN format sequence name (sample#haplotype#chr -> sample)
pub fn extract_sample(seq_name: &str, separator: &str) -> String {
    let parts: Vec<&str> = seq_name.split(separator).collect();
    if !parts.is_empty() {
        parts[0].to_string()
    } else {
        seq_name.to_string()
    }
}

/// (inter_start, inter_end, anchor_start, anchor_end) for one hop-0 segment.
type HopZeroSeg = (i64, i64, i64, i64);

/// Project `[t_start, t_end)` on an intermediate sequence back onto the
/// anchor using the containing hop-0 segment. Falls back to `[region_start,
/// region_end)` when no segment contains the range (gap between hop-0
/// hits, or strand-specific coverage that does not match).
#[inline]
/// Project a hop≥1 BFS hit's target-side range `(t_start, t_end)` (lying on
/// some intermediate hub sequence) back onto the anchor's coordinate system,
/// using the pass-1 `seq_anchor_coverage` map of `(q_start, q_end, anc_start,
/// anc_end)` segments registered for that hub.
///
/// Clip-then-project semantics:
///
/// - For every segment whose `(seq_start, seq_end)` overlaps the hit's
///   `(t_start, t_end)`, clip the hit to the segment, linearly interpolate the
///   clipped range to anchor coordinates, and union the projections.
/// - Multiple overlapping segments contribute via min/max — equivalent to
///   "take the convex hull of the projected ranges on anchor".
/// - If no segment overlaps the hit, return an empty range
///   `(region_start, region_start)` so the caller's
///   `if a_start >= a_end { continue }` drops the hit instead of crediting the
///   entire anchor span. This is the critical correctness fix
///   (PLAN_depth_100pct.md §3): the previous strict-containment + whole-region
///   fallback was the source of raw-BFS over-reports at high transitive depth.
fn project_hop0_coords(
    segments: Option<&Vec<HopZeroSeg>>,
    t_start: i64,
    t_end: i64,
    region_start: i64,
    region_end: i64,
) -> (i64, i64) {
    if let Some(segs) = segments {
        let mut proj_min: i64 = i64::MAX;
        let mut proj_max: i64 = i64::MIN;
        for &(seq_start, seq_end, anc_start, anc_end) in segs {
            // Clip the hit's t-range into this segment.
            let clip_s = t_start.max(seq_start);
            let clip_e = t_end.min(seq_end);
            if clip_s >= clip_e {
                continue; // no overlap with this segment
            }
            let seq_len = seq_end - seq_start;
            let anc_len = anc_end - anc_start;
            if seq_len <= 0 || anc_len <= 0 {
                continue;
            }
            let proj_s = anc_start + proj_offset(clip_s - seq_start, anc_len, seq_len);
            let proj_e = anc_start + proj_offset(clip_e - seq_start, anc_len, seq_len);
            let lo = proj_s.min(proj_e);
            let hi = proj_s.max(proj_e);
            if lo < proj_min {
                proj_min = lo;
            }
            if hi > proj_max {
                proj_max = hi;
            }
        }
        if proj_min < proj_max {
            return (
                proj_min.max(region_start),
                proj_max.min(region_end),
            );
        }
    }
    // No segment overlaps the hit: drop it. Empty range routes to the caller's
    // `if a_start >= a_end { continue }` check.
    depth_trace!(
        "FALLBACK project_hop0_coords t={}-{} segs={} → empty [no overlap; hit dropped]",
        t_start,
        t_end,
        segments.map_or(0, |v| v.len())
    );
    (region_start, region_start)
}
// ============================================================================
// Compact depth computation using ID-based data structures
// ============================================================================

/// Compact depth event for sweep-line algorithm using numeric IDs
#[derive(Debug, Clone, PartialEq, Eq)]
struct CompactDepthEvent {
    position: i64,
    is_start: bool,
    sample_id: u16,
    /// Index into the alignment info array
    alignment_idx: usize,
}

impl CompactDepthEvent {
    /// Pack (position, !is_start, sample_id) into a single u64 sort key matching
    /// the Ord impl: position ascending, then starts before ends, then sample_id.
    /// Assumes 0 <= position < 2^47 (chr length << 2^47).
    #[inline]
    fn packed_sort_key(&self) -> u64 {
        debug_assert!(self.position >= 0 && self.position < (1i64 << 47));
        ((self.position as u64) << 17)
            | ((!self.is_start as u64) << 16)
            | (self.sample_id as u64)
    }
}

impl Ord for CompactDepthEvent {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| other.is_start.cmp(&self.is_start)) // Starts before ends at same position
            .then_with(|| self.sample_id.cmp(&other.sample_id))
    }
}

impl PartialOrd for CompactDepthEvent {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Check if an alignment is a same-sample alignment (intra-genome duplication).
/// Returns true if the alignment should be filtered out.
#[inline]
fn is_self_alignment(
    query_sample_id: u16,
    anchor_sample_id: u16,
    _query_id: u32,
    _target_seq_id: u32,
) -> bool {
    query_sample_id == anchor_sample_id
}

// ============================================================================
// Windowed depth computation - fixed window size, sample by sample
// ============================================================================

/// Default window size for windowed depth computation (50kb)
const DEFAULT_WINDOW_SIZE: i64 = 50_000;

// ============================================================================
// Windowed depth structures - sparse storage, streaming output
// ============================================================================

/// Sample position entry (sparse representation)
/// (sample_id, query_seq_id, query_start, query_end)
type SamplePosition = (u16, u32, i64, i64);

/// Depth interval with sparse sample storage
/// Memory: 24 bytes base + 16 bytes per actual sample (vs 16 bytes per ALL samples)
#[derive(Debug, Clone)]
struct SparseDepthInterval {
    start: i64,
    end: i64,
    /// Sparse sample positions - only stores samples that are present
    /// Sorted by sample_id for binary search lookup
    samples: Vec<SamplePosition>,
    /// Total pangenome bases: sum of all samples' union query lengths for this interval.
    /// Accounts for multiple overlapping alignments from the same sample by computing
    /// the union of their query-side projections (grouped by contig).
    pangenome_bases: i64,
}

impl SparseDepthInterval {
    #[inline]
    fn depth(&self) -> usize {
        self.samples.len()
    }

    /// Get position for a specific sample (binary search)
    #[inline]
    fn get_sample(&self, sample_id: u16) -> Option<(u32, i64, i64)> {
        self.samples
            .binary_search_by_key(&sample_id, |s| s.0)
            .ok()
            .map(|idx| {
                (
                    self.samples[idx].1,
                    self.samples[idx].2,
                    self.samples[idx].3,
                )
            })
    }
}

// ============================================================================
// Streaming depth output pipeline
// ============================================================================

/// Streaming depth output pipeline.
///
/// Processes intervals one at a time through:
///   sweep_line → window_split → merge → format/stats → buffered_write
///
/// This eliminates the O(K × N_samples) peak memory from collecting all
/// SparseDepthIntervals into a Vec, reducing it to O(N_samples) per interval.
struct StreamingDepthEmitter<'a> {
    // ---- Output buffer with periodic flushing ----
    buf: Vec<u8>,
    writer: &'a Option<DepthWriter>,

    // ---- Formatting context ----
    seq_name: &'a str,
    anchor_sample_id: u16,
    num_samples: usize,
    seq_index: &'a crate::seqidx::SequenceIndex,
    row_counter: &'a AtomicUsize,
    intervals_counter: &'a AtomicUsize,
    should_output: bool,

    // ---- Stats accumulators (thread-local, merged into global at end) ----
    stats_mode: bool,
    local_stats: Option<DepthStats>,
    local_combined: Option<DepthStatsWithSamples>,
    sample_index: &'a SampleIndex,

    // ---- Short-interval merging ----
    /// Intervals shorter than this (bp) are absorbed into adjacent neighbors. 0 = disabled.
    min_interval_len: i64,
    /// Buffer for windowed output (only used when `window_size.is_some()`).
    /// The windowed path retains buffered semantics so cross-window merging stays
    /// equivalent to the pre-streaming implementation. Non-windowed paths bypass
    /// this entirely and stream through `pending` (or directly when no merge is
    /// requested), so this stays empty in the hot Phase 1/2 case.
    seq_intervals: Vec<SparseDepthInterval>,
    /// Two-element sliding window for `min_interval_len > 0` streaming.
    ///
    /// State machine (see `merge_into_pending`): equivalent to the two-pass
    /// `merge_short_intervals` algorithm, but holds only the currently-growing
    /// interval rather than the full per-region Vec. Peak memory drops from
    /// `O(K × N_samples)` to `O(N_samples)`.
    pending: Option<SparseDepthInterval>,

    // ---- Windowing ----
    window_size: Option<i64>,
}

impl<'a> StreamingDepthEmitter<'a> {
    /// Process one interval from the sweep line through the full pipeline.
    ///
    /// Three paths:
    ///   1. Windowed (`window_size.is_some()`): split at window boundaries, buffer
    ///      in `seq_intervals` for cross-window `merge_short_intervals` at flush.
    ///      This preserves the pre-streaming windowed semantics.
    ///   2. Non-windowed, no short-interval merging (`min_interval_len <= 0`):
    ///      true streaming — emit immediately, drop the interval and its samples
    ///      Vec right away. This is the hot Phase 1 path; previous behavior was
    ///      to accumulate the full chromosome's intervals before flushing.
    ///   3. Non-windowed with merging (`min_interval_len > 0`): two-element
    ///      sliding window via `pending`. Equivalent to the two-pass
    ///      `merge_short_intervals` algorithm but with O(1) buffered intervals.
    fn emit(&mut self, interval: SparseDepthInterval) {
        if let Some(ws) = self.window_size {
            self.split_and_forward(interval, ws);
        } else if self.min_interval_len <= 0 {
            self.emit_final(interval);
        } else {
            self.merge_into_pending(interval);
        }
    }

    /// Sliding-window equivalent of `merge_short_intervals` for the streaming path.
    ///
    /// State machine derivation:
    ///
    ///   - Original Pass 1 (left-to-right): a short `new` is absorbed into its left
    ///     neighbor if one exists; otherwise it is pushed (becomes a "leading short").
    ///     A long `new` is always pushed.
    ///
    ///   - Original Pass 2: any leading items shorter than the first item that is
    ///     itself ≥ `min_len` are absorbed into that first long item.
    ///
    /// Streaming reformulation (only `pending` is held):
    ///
    ///   * `new` is short → right-extend `pending` (Pass 1 left-absorb).
    ///   * `new` is long and `pending` is already long (length ≥ min_len) → emit
    ///     `pending`, then `pending = new` (no Pass-2 absorption needed).
    ///   * `new` is long and `pending` is still short → left-absorb `pending` into
    ///     `new` (Pass-2 leading-shorts-into-first-long), then `pending = new`.
    ///
    /// At gap/seq end, `flush_seq_intervals` emits whatever `pending` holds.
    fn merge_into_pending(&mut self, new: SparseDepthInterval) {
        let new_is_short = (new.end - new.start) < self.min_interval_len;
        match self.pending.take() {
            None => {
                self.pending = Some(new);
            }
            Some(mut pending) => {
                if new_is_short {
                    Self::absorb_right(&mut pending, new);
                    self.pending = Some(pending);
                } else if (pending.end - pending.start) >= self.min_interval_len {
                    self.emit_final(pending);
                    self.pending = Some(new);
                } else {
                    let mut absorbed = new;
                    Self::absorb_left(pending, &mut absorbed);
                    self.pending = Some(absorbed);
                }
            }
        }
    }

    /// Right-extend `left` by absorbing `right`: end advances to `right.end`,
    /// pangenome bases sum, and `right`'s samples are unioned into `left`'s
    /// (per-sample query coordinates expand to the union of overlapping ranges).
    fn absorb_right(left: &mut SparseDepthInterval, right: SparseDepthInterval) {
        left.end = right.end;
        left.pangenome_bases += right.pangenome_bases;
        let mut sample_idx: FxHashMap<u16, usize> = left
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.0, i))
            .collect();
        for sp in right.samples {
            if let Some(&i) = sample_idx.get(&sp.0) {
                left.samples[i].2 = left.samples[i].2.min(sp.2);
                left.samples[i].3 = left.samples[i].3.max(sp.3);
            } else {
                sample_idx.insert(sp.0, left.samples.len());
                left.samples.push(sp);
            }
        }
        left.samples.sort_by_key(|s| s.0);
    }

    /// Left-extend `right` by absorbing `left`: start retreats to `left.start`,
    /// pangenome bases sum, and `left`'s samples are unioned into `right`'s.
    /// Mirrors the Pass-2 leading-shorts-into-first-long path of
    /// `merge_short_intervals`.
    fn absorb_left(left: SparseDepthInterval, right: &mut SparseDepthInterval) {
        right.start = left.start;
        right.pangenome_bases += left.pangenome_bases;
        let mut sample_idx: FxHashMap<u16, usize> = right
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.0, i))
            .collect();
        for sp in left.samples {
            if let Some(&i) = sample_idx.get(&sp.0) {
                right.samples[i].2 = right.samples[i].2.min(sp.2);
                right.samples[i].3 = right.samples[i].3.max(sp.3);
            } else {
                sample_idx.insert(sp.0, right.samples.len());
                right.samples.push(sp);
            }
        }
        right.samples.sort_by_key(|s| s.0);
    }

    /// Split an interval at window boundaries and forward each piece.
    fn split_and_forward(&mut self, interval: SparseDepthInterval, window_size: i64) {
        let mut pos = interval.start;
        let interval_len = interval.end - interval.start;
        while pos < interval.end {
            let next_boundary = ((pos / window_size) + 1) * window_size;
            let chunk_end = interval.end.min(next_boundary);
            if pos >= chunk_end {
                break;
            }

            let samples: Vec<SamplePosition> = interval
                .samples
                .iter()
                .map(|&(sid, qid, qs, qe)| {
                    if interval_len <= 0 {
                        return (sid, qid, qs, qe);
                    }
                    let q_len = qe - qs;
                    let off_s = pos - interval.start;
                    let off_e = chunk_end - interval.start;
                    let new_qs = qs + proj_offset(off_s, q_len, interval_len);
                    let new_qe = qs + proj_offset(off_e, q_len, interval_len);
                    (sid, qid, new_qs, new_qe)
                })
                .collect();

            let chunk_pangenome_bases = if interval_len > 0 {
                proj_offset(chunk_end - pos, interval.pangenome_bases, interval_len)
            } else {
                interval.pangenome_bases
            };

            self.seq_intervals.push(SparseDepthInterval {
                start: pos,
                end: chunk_end,
                samples,
                pangenome_bases: chunk_pangenome_bases,
            });
            pos = chunk_end;
        }
    }



    /// Terminal stage: format the interval (TSV or stats) and buffer the output.
    fn emit_final(&mut self, interval: SparseDepthInterval) {
        if !self.should_output {
            return;
        }

        if self.stats_mode {
            let depth = interval.depth();
            if let Some(ref mut ls) = self.local_stats {
                ls.add_interval(
                    self.seq_name,
                    interval.start,
                    interval.end,
                    depth,
                    interval.pangenome_bases,
                );
            }
            if let Some(ref mut lc) = self.local_combined {
                let sample_names: Vec<String> = interval
                    .samples
                    .iter()
                    .filter_map(|&(sid, _, _, _)| {
                        self.sample_index.get_name(sid).map(|s| s.to_string())
                    })
                    .collect();
                lc.add_interval(
                    self.seq_name,
                    interval.start,
                    interval.end,
                    depth,
                    sample_names,
                    interval.pangenome_bases,
                );
            }
        } else {
            let rid = self.row_counter.fetch_add(1, Ordering::Relaxed) + 1;
            let _ = write!(
                self.buf,
                "{}\t{}\t{}",
                rid,
                interval.end - interval.start,
                interval.depth()
            );
            // Anchor position first
            let _ = write!(self.buf, "\t{}:{}-{}", self.seq_name, interval.start, interval.end);
            // Then other samples (non-anchor) sorted by sample_id
            for &(sid, query_id, q_start, q_end) in &interval.samples {
                if sid != self.anchor_sample_id {
                    if let Some(name) = self.seq_index.get_name(query_id) {
                        let _ = write!(self.buf, ";{}:{}-{}", name, q_start, q_end);
                    }
                }
            }
            let _ = writeln!(self.buf);
            self.intervals_counter.fetch_add(1, Ordering::Relaxed);

            // Periodic flush to bound per-thread buffer memory.
            // Hand ownership of the buffer to the writer thread (move, not
            // copy) and start a fresh Vec so the worker keeps zero-copy
            // ownership semantics on the next flush.
            if self.buf.len() >= STREAMING_FLUSH_THRESHOLD {
                if let Some(ref w) = self.writer {
                    let buf = std::mem::take(&mut self.buf);
                    let _ = w.send(buf);
                } else {
                    self.buf.clear();
                }
            }
        }
        // interval is DROPPED here — Vec<SamplePosition> freed immediately
    }

    /// Flush buffered intervals for the current contiguous depth region.
    ///
    /// Two independent buffers may carry state across this call:
    ///   - `seq_intervals` (windowed path only): drained through
    ///     `merge_short_intervals` to preserve cross-window merging semantics.
    ///   - `pending` (non-windowed `min_interval_len > 0` streaming path):
    ///     emitted as-is. If `pending` is still shorter than `min_interval_len`
    ///     at flush time, no long interval ever followed in this gap, so it is
    ///     emitted unchanged — matching `merge_short_intervals`'s "boundary is
    ///     None → keep all leading shorts as-is" branch.
    ///
    /// Callers invoke this between Phase 2 gaps so leading shorts of one gap
    /// never absorb into the trailing tail of the previous (non-contiguous) gap.
    fn flush_seq_intervals(&mut self) {
        if !self.seq_intervals.is_empty() {
            let intervals = std::mem::take(&mut self.seq_intervals);
            let merged = merge_short_intervals(intervals, self.min_interval_len);
            for interval in merged {
                self.emit_final(interval);
            }
        }
        if let Some(p) = self.pending.take() {
            self.emit_final(p);
        }
    }

    /// Flush everything: pending intervals, remaining buffer, return stats for merging.
    fn flush(&mut self) -> io::Result<()> {
        self.flush_seq_intervals();
        if !self.buf.is_empty() {
            if let Some(ref w) = self.writer {
                let buf = std::mem::take(&mut self.buf);
                w.send(buf)?;
            } else {
                self.buf.clear();
            }
        }
        Ok(())
    }

    /// Merge thread-local stats into global accumulators.
    fn merge_stats_into(
        self,
        stats_accumulator: &Option<Mutex<DepthStats>>,
        stats_combined_acc: &Option<Mutex<DepthStatsWithSamples>>,
    ) {
        if let Some(ls) = self.local_stats {
            if let Some(ref acc) = stats_accumulator {
                acc.lock().merge(&ls);
            }
        }
        if let Some(lc) = self.local_combined {
            if let Some(ref acc) = stats_combined_acc {
                acc.lock().merge(lc);
            }
        }
    }
}

/// Split intervals at fixed window boundaries.
/// Each interval crossing a window edge is cut into pieces,
/// with sample query coordinates proportionally adjusted.
fn split_intervals_by_window(
    intervals: Vec<SparseDepthInterval>,
    window_size: i64,
) -> Vec<SparseDepthInterval> {
    let mut result: Vec<SparseDepthInterval> = Vec::new();
    for interval in intervals {
        let mut pos = interval.start;
        while pos < interval.end {
            // Next window boundary after pos
            let next_boundary = ((pos / window_size) + 1) * window_size;
            let chunk_end = interval.end.min(next_boundary);

            if pos >= chunk_end {
                break;
            }

            let interval_len = interval.end - interval.start;
            // Proportionally adjust sample query coordinates for this chunk
            let samples: Vec<SamplePosition> = interval
                .samples
                .iter()
                .map(|&(sid, qid, qs, qe)| {
                    if interval_len <= 0 {
                        return (sid, qid, qs, qe);
                    }
                    let q_len = qe - qs;
                    let off_s = pos - interval.start;
                    let off_e = chunk_end - interval.start;
                    let new_qs = qs + proj_offset(off_s, q_len, interval_len);
                    let new_qe = qs + proj_offset(off_e, q_len, interval_len);
                    (sid, qid, new_qs, new_qe)
                })
                .collect();

            // Proportionally split pangenome_bases
            let chunk_pangenome_bases = if interval_len > 0 {
                proj_offset(chunk_end - pos, interval.pangenome_bases, interval_len)
            } else {
                interval.pangenome_bases
            };

            result.push(SparseDepthInterval {
                start: pos,
                end: chunk_end,
                samples,
                pangenome_bases: chunk_pangenome_bases,
            });
            pos = chunk_end;
        }
    }
    result
}

/// Absorb intervals shorter than `min_len` bp into adjacent neighbors.
///
/// Short intervals caused by transient alignment dropouts (e.g., one sample
/// briefly loses coverage) are absorbed into their left neighbor; if no left
/// neighbor exists, they are absorbed into the right neighbor.  The absorbing
/// neighbor's depth and sample set are preserved — only its coordinate span
/// grows and its `pangenome_bases` accumulates the absorbed interval's bases.
///
/// Two passes are sufficient: the left-to-right pass handles all interior and
/// trailing short intervals; the second pass absorbs any leading short
/// intervals into the first long interval.
fn merge_short_intervals(
    intervals: Vec<SparseDepthInterval>,
    min_len: i64,
) -> Vec<SparseDepthInterval> {
    if min_len <= 0 || intervals.len() <= 1 {
        return intervals;
    }

    // Pass 1: left-to-right — absorb each short interval into its left neighbor.
    let mut result: Vec<SparseDepthInterval> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if interval.end - interval.start < min_len {
            if let Some(left) = result.last_mut() {
                left.end = interval.end;
                left.pangenome_bases += interval.pangenome_bases;
                // Build a position index for O(1) lookup per absorbed sample.
                let mut sample_idx: FxHashMap<u16, usize> = left
                    .samples
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.0, i))
                    .collect();
                for sp in interval.samples {
                    if let Some(&i) = sample_idx.get(&sp.0) {
                        left.samples[i].2 = left.samples[i].2.min(sp.2);
                        left.samples[i].3 = left.samples[i].3.max(sp.3);
                    } else {
                        sample_idx.insert(sp.0, left.samples.len());
                        left.samples.push(sp);
                    }
                }
                // samples remains sorted because active_samples() yielded sorted IDs;
                // newly pushed entries from the absorbed interval are already sorted
                // relative to existing ones only if the absorbed set is a superset disjoint
                // from left — in the general case we must re-sort.
                left.samples.sort_by_key(|s| s.0);
                continue;
            }
            // No left neighbor yet — fall through; pass 2 will absorb into right.
        }
        result.push(interval);
    }

    // Pass 2: absorb any remaining leading short intervals into the first long one.
    let boundary = result.iter().position(|iv| iv.end - iv.start >= min_len);
    if let Some(b) = boundary {
        if b > 0 {
            let extra_pb: i64 = result[..b].iter().map(|iv| iv.pangenome_bases).sum();
            let new_start = result[0].start;
            // Union samples from all absorbed leading intervals into result[b]
            let absorbed_samples: Vec<SamplePosition> = result[..b]
                .iter()
                .flat_map(|iv| iv.samples.iter().copied())
                .collect();
            result[b].start = new_start;
            result[b].pangenome_bases += extra_pb;
            let mut sample_idx: FxHashMap<u16, usize> = result[b]
                .samples
                .iter()
                .enumerate()
                .map(|(i, s)| (s.0, i))
                .collect();
            for sp in absorbed_samples {
                if let Some(&i) = sample_idx.get(&sp.0) {
                    result[b].samples[i].2 = result[b].samples[i].2.min(sp.2);
                    result[b].samples[i].3 = result[b].samples[i].3.max(sp.3);
                } else {
                    sample_idx.insert(sp.0, result[b].samples.len());
                    result[b].samples.push(sp);
                }
            }
            result[b].samples.sort_by_key(|s| s.0);
            result.drain(..b);
        }
    }
    // If boundary is None every interval is shorter than min_len — keep as-is.

    result
}

/// Same absorption logic as `merge_short_intervals` but for the stats `DepthIntervalWithSamples`
/// type. Short intervals are absorbed into their left neighbor (or right if no left exists);
/// the absorbing neighbor's depth is kept and the sample sets are unioned.
fn merge_short_depth_intervals(
    intervals: Vec<DepthIntervalWithSamples>,
    min_len: i64,
) -> Vec<DepthIntervalWithSamples> {
    if min_len <= 0 || intervals.len() <= 1 {
        return intervals;
    }

    // Pass 1: left-to-right absorption
    let mut result: Vec<DepthIntervalWithSamples> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if interval.end - interval.start < min_len {
            if let Some(left) = result.last_mut() {
                left.end = interval.end;
                for s in interval.samples {
                    if !left.samples.contains(&s) {
                        left.samples.push(s);
                    }
                }
                continue;
            }
        }
        result.push(interval);
    }

    // Pass 2: absorb leading short intervals into the first long one
    let boundary = result.iter().position(|iv| iv.end - iv.start >= min_len);
    if let Some(b) = boundary {
        if b > 0 {
            let new_start = result[0].start;
            let extra_samples: Vec<String> = result[..b]
                .iter()
                .flat_map(|iv| iv.samples.iter().cloned())
                .collect();
            result[b].start = new_start;
            for s in extra_samples {
                if !result[b].samples.contains(&s) {
                    result[b].samples.push(s);
                }
            }
            result.drain(..b);
        }
    }

    result
}

/// Thread-safe processed region tracker using per-sequence locks.
/// Each sequence has its own Mutex, so concurrent access to different
/// sequences has zero contention. This enables fully parallel depth
/// computation across all sequences without batch synchronization barriers.
#[derive(Debug)]
struct ConcurrentProcessedTracker {
    /// Per-sequence IntervalSet behind a Mutex for thread-safe access
    processed: Vec<Mutex<IntervalSet>>,
}

impl ConcurrentProcessedTracker {
    /// Create a new tracker for the given number of sequences
    fn new(num_sequences: usize) -> Self {
        Self {
            processed: (0..num_sequences)
                .map(|_| Mutex::new(IntervalSet::new()))
                .collect(),
        }
    }

    /// Get unprocessed regions for a sequence within [start, end)
    fn get_unprocessed(&self, seq_id: u32, start: i64, end: i64) -> Vec<(i64, i64)> {
        let lock = self.processed[seq_id as usize].lock();
        if lock.is_empty() {
            return vec![(start, end)];
        }
        // Subtract only those processed intervals that actually overlap [start, end).
        // `iter_overlapping` is O(log N + k) instead of O(N).
        let mut unprocessed = IntervalSet::new_single(start, end);
        for (s, e) in lock.iter_overlapping(start, end) {
            unprocessed.subtract(s, e);
        }
        unprocessed.intervals().collect()
    }

    /// Mark a region as processed
    fn mark_processed(&self, seq_id: u32, start: i64, end: i64) {
        let mut lock = self.processed[seq_id as usize].lock();
        lock.add(start, end);
    }

    /// Mark a batch of regions as processed: Vec<(seq_id, start, end)>
    fn mark_processed_batch(&self, regions: &[(u32, i64, i64)]) {
        for &(seq_id, start, end) in regions {
            self.mark_processed(seq_id, start, end);
        }
    }

    /// Atomically claim unprocessed ranges within [start, end): returns what was claimed and marks it processed.
    fn claim_unprocessed(&self, seq_id: u32, start: i64, end: i64) -> Vec<(i64, i64)> {
        let mut lock = self.processed[seq_id as usize].lock();
        let unprocessed = if lock.is_empty() {
            vec![(start, end)]
        } else {
            // Only the processed intervals that actually overlap [start, end)
            // can affect the result. Bounded range iteration keeps each claim
            // O(log N + k) instead of O(N).
            let mut result = IntervalSet::new_single(start, end);
            for (s, e) in lock.iter_overlapping(start, end) {
                result.subtract(s, e);
            }
            result.intervals().collect::<Vec<_>>()
        };
        for &(s, e) in &unprocessed {
            lock.add(s, e);
        }
        unprocessed
    }

    /// Bool-only variant of [`claim_unprocessed`].
    ///
    /// Returns `true` iff any portion of `[start, end)` was previously
    /// unprocessed; in either case, the entire range is marked processed
    /// afterwards. The post-state is identical to `claim_unprocessed`'s
    /// (`IntervalSet::add` is idempotent), but no `Vec<(i64, i64)>` is ever
    /// allocated.
    ///
    /// Use at hot-loop call sites that only need the empty/non-empty bit
    /// (e.g., `process_anchor_region_raw_streaming` calls this once per query
    /// alignment — at CHM13 / 580-sample scale that is tens of millions of
    /// invocations and the original Vec allocation showed up as significant
    /// allocator pressure under jemalloc profiling).
    #[allow(dead_code)] // retained for parity tests vs claim_any_unprocessed_batch
    fn claim_any_unprocessed(&self, seq_id: u32, start: i64, end: i64) -> bool {
        if start >= end {
            return false;
        }
        let mut lock = self.processed[seq_id as usize].lock();
        let span = end - start;
        let any_unprocessed = if lock.is_empty() {
            true
        } else {
            // IntervalSet maintains disjoint, sorted, merged intervals, so
            // summing clipped overlap lengths gives the exact covered area.
            let mut covered: i64 = 0;
            for (s, e) in lock.iter_overlapping(start, end) {
                let cs = s.max(start);
                let ce = e.min(end);
                if ce > cs {
                    covered += ce - cs;
                    if covered >= span {
                        break;
                    }
                }
            }
            covered < span
        };
        lock.add(start, end);
        any_unprocessed
    }

    /// Batched [`claim_any_unprocessed`] for one `seq_id`: takes the per-seq
    /// mutex once and processes `ranges` in input order, with the same
    /// per-range semantics as the single-shot version. The output is appended
    /// to `out` in input order so callers can map results back to range indices.
    /// End state and per-range bool sequence are identical to N serial calls.
    fn claim_any_unprocessed_batch(
        &self,
        seq_id: u32,
        ranges: &[(i64, i64)],
        out: &mut Vec<bool>,
    ) {
        out.reserve(ranges.len());
        if ranges.is_empty() {
            return;
        }
        let mut lock = self.processed[seq_id as usize].lock();
        for &(start, end) in ranges {
            if start >= end {
                out.push(false);
                continue;
            }
            let span = end - start;
            let any_unprocessed = if lock.is_empty() {
                true
            } else {
                let mut covered: i64 = 0;
                for (s, e) in lock.iter_overlapping(start, end) {
                    let cs = s.max(start);
                    let ce = e.min(end);
                    if ce > cs {
                        covered += ce - cs;
                        if covered >= span {
                            break;
                        }
                    }
                }
                covered < span
            };
            lock.add(start, end);
            out.push(any_unprocessed);
        }
    }
}

// ============================================================================
// Global depth computation - connected-component traversal
// ============================================================================

/// Sidecar cache for the alignment-degree pre-scan.
///
/// `compute_sample_degrees` walks every per-file index header on every
/// invocation; for ≥10⁵-file (per-file index) workloads at 580+ samples this
/// is the dominant Phase-0 cost on repeat runs. Persisting the result to a
/// bin file next to the input alignment list lets repeat runs against the
/// same input load degrees from disk in well under a second instead of
/// re-walking every sub-index.
///
/// **Invalidation key** combines:
/// - hash of every alignment file path + its mtime (seconds) + size, so any
///   on-disk change forces a rebuild;
/// - hash of `seq_included` (the inclusion bitmap), since a different
///   `--min-seq-length` produces a different result;
/// - `min_seq_length` itself (defensive — equivalent to `seq_included` hash
///   for the canonical caller, but keeps the key meaningful if the caller
///   ever derives `seq_included` from something other than length).
///
/// The cache is intentionally per-input rather than per-output: multiple
/// `impg depth` invocations with different `--output-prefix` against the
/// same alignment list should all share the same cached degrees.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
struct DegreesCache {
    /// Schema version — bump when the on-disk layout changes incompatibly.
    schema_version: u32,
    /// Combined invalidation hash (see DegreesCache docs).
    invalidation_hash: u64,
    /// Number of unified sequences at cache time. Validated post-load.
    num_sequences: u32,
    /// Per-sequence degree (indexed by unified seq_id), length =
    /// `num_sequences`.
    degrees: Vec<u16>,
}

const DEGREES_CACHE_SCHEMA_VERSION: u32 = 1;

/// Compute the invalidation hash from alignment file mtime/size tuples,
/// `seq_included`, and `min_seq_length`. Returns None if any alignment file
/// metadata cannot be read (cache cannot be reliably keyed → bypass).
fn compute_degrees_invalidation_hash(
    alignment_files: &[String],
    seq_included: &[bool],
    min_seq_length: i64,
) -> Option<u64> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    if alignment_files.is_empty() {
        return None;
    }

    let mut hasher = DefaultHasher::new();
    DEGREES_CACHE_SCHEMA_VERSION.hash(&mut hasher);
    alignment_files.len().hash(&mut hasher);
    for path in alignment_files {
        path.hash(&mut hasher);
        let meta = match std::fs::metadata(path) {
            Ok(m) => m,
            Err(_) => return None,
        };
        meta.len().hash(&mut hasher);
        // mtime in seconds — we accept second-level granularity; sub-second
        // edits within the same second are uncommon for alignment input.
        if let Ok(modified) = meta.modified() {
            if let Ok(dur) = modified.duration_since(std::time::UNIX_EPOCH) {
                dur.as_secs().hash(&mut hasher);
            }
        }
    }
    seq_included.len().hash(&mut hasher);
    // Hash the seq_included bitmap as packed u64 chunks for compactness.
    let mut chunk: u64 = 0;
    let mut bit: u32 = 0;
    for &v in seq_included {
        if v {
            chunk |= 1u64 << bit;
        }
        bit += 1;
        if bit == 64 {
            chunk.hash(&mut hasher);
            chunk = 0;
            bit = 0;
        }
    }
    if bit > 0 {
        chunk.hash(&mut hasher);
    }
    min_seq_length.hash(&mut hasher);
    Some(hasher.finish())
}

/// Choose where to store the degrees sidecar.
///
/// Strategy: place it next to the first alignment file (a stable, input-
/// derived location that survives across runs without polluting `~/.cache`).
/// The filename incorporates the invalidation hash so concurrent runs against
/// different `seq_included` / `min_seq_length` configurations don't trample
/// each other.
fn degrees_cache_path(alignment_files: &[String], hash: u64) -> Option<std::path::PathBuf> {
    use std::path::Path;
    let first = alignment_files.first()?;
    let parent = Path::new(first).parent().unwrap_or(Path::new("."));
    Some(parent.join(format!("impg_depth_degrees_{:016x}.bin", hash)))
}

/// Try to load a previously-saved degrees cache.
fn try_load_degrees_cache(
    cache_path: &std::path::Path,
    expected_hash: u64,
    expected_num_sequences: u32,
) -> Option<Vec<u16>> {
    let file = std::fs::File::open(cache_path).ok()?;
    let mut reader = std::io::BufReader::new(file);
    let cache: DegreesCache =
        bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard()).ok()?;
    if cache.schema_version != DEGREES_CACHE_SCHEMA_VERSION {
        return None;
    }
    if cache.invalidation_hash != expected_hash {
        return None;
    }
    if cache.num_sequences != expected_num_sequences {
        return None;
    }
    if cache.degrees.len() != expected_num_sequences as usize {
        return None;
    }
    Some(cache.degrees)
}

/// Save degrees to the sidecar file atomically (write to temp file then
/// rename). Returns Ok even if the rename succeeds; logs but does not
/// surface write errors to the caller (cache is best-effort).
fn save_degrees_cache(
    cache_path: &std::path::Path,
    invalidation_hash: u64,
    degrees: &[u16],
) -> io::Result<()> {
    let cache = DegreesCache {
        schema_version: DEGREES_CACHE_SCHEMA_VERSION,
        invalidation_hash,
        num_sequences: degrees.len() as u32,
        degrees: degrees.to_vec(),
    };
    let tmp_path = cache_path.with_extension("bin.tmp");
    {
        let file = std::fs::File::create(&tmp_path)?;
        let mut writer = std::io::BufWriter::new(file);
        bincode::serde::encode_into_std_write(&cache, &mut writer, bincode::config::standard())
            .map_err(|e| io::Error::other(format!("bincode encode: {}", e)))?;
        writer.flush()?;
    }
    std::fs::rename(&tmp_path, cache_path)?;
    Ok(())
}

/// Pre-scan: compute alignment degree for each included sequence.
/// Degree = number of unique OTHER samples with direct alignments to this sequence.
/// Used to automatically identify hub sequences for Phase 1 when --ref is not specified.
///
/// On repeat runs against the same alignment input + `seq_included`
/// configuration, loads the result from a sidecar bin file (next to the
/// first alignment file). On any cache miss / hash mismatch / IO error, falls
/// back to a full pre-scan and writes a fresh sidecar.
fn compute_alignment_degrees(
    impg: &(impl ImpgIndex + Sync),
    compact_lengths: &CompactSequenceLengths,
    seq_included: &[bool],
    min_seq_length: i64,
) -> Vec<u16> {
    let alignment_files = impg.alignment_files();
    let cache_key = compute_degrees_invalidation_hash(alignment_files, seq_included, min_seq_length);
    let cache_path = cache_key.and_then(|h| degrees_cache_path(alignment_files, h));

    if let (Some(hash), Some(path)) = (cache_key, cache_path.as_ref()) {
        if let Some(cached) =
            try_load_degrees_cache(path, hash, seq_included.len() as u32)
        {
            info!(
                "Loaded alignment degrees from sidecar cache: {} ({} entries)",
                path.display(),
                cached.len()
            );
            return cached;
        }
    }

    // Cache miss / disabled → full pre-scan via the trait method.
    // - Single `Impg`: default parallel-by-target impl (fine: one file, one cache).
    // - `MultiImpg`: file-parallel override loading each sub-index transiently
    //   and dropping it immediately, bounding peak retained sub-indices to the
    //   rayon worker count.
    let degrees =
        impg.compute_sample_degrees(seq_included, compact_lengths.seq_to_sample_slice());

    // Best-effort sidecar save. Failure is logged but not propagated — the
    // cache is purely an optimisation.
    if let (Some(hash), Some(path)) = (cache_key, cache_path) {
        if let Err(e) = save_degrees_cache(&path, hash, &degrees) {
            log::warn!(
                "Failed to write degrees sidecar cache to {}: {}",
                path.display(),
                e
            );
        } else {
            debug!("Saved alignment degrees sidecar cache: {}", path.display());
        }
    }

    degrees
}

/// Build sequence processing order: sorted by degree descending, then length descending.
/// If ref_sample is specified, that sample's sequences are moved to the front.
fn build_sequence_order(
    impg: &impl ImpgIndex,
    compact_lengths: &CompactSequenceLengths,
    ref_sample_id: Option<u16>,
    seq_included: &[bool],
    degrees: &[u16],
) -> Vec<u32> {
    let num_seqs = impg.seq_index().len();
    let mut seq_order: Vec<(u32, i64, bool, u16)> = (0..num_seqs as u32)
        .filter_map(|seq_id| {
            if !seq_included.get(seq_id as usize).copied().unwrap_or(false) {
                return None;
            }
            let len = compact_lengths.get_length(seq_id);
            if len <= 0 {
                return None;
            }
            let is_ref = ref_sample_id
                .map(|ref_id| compact_lengths.get_sample_id(seq_id) == ref_id)
                .unwrap_or(false);
            let degree = degrees.get(seq_id as usize).copied().unwrap_or(0);
            Some((seq_id, len, is_ref, degree))
        })
        .collect();

    // Sort: ref sample first, then by degree descending, then by length descending
    seq_order.sort_by(|a, b| {
        b.2.cmp(&a.2) // ref sequences first
            .then_with(|| b.3.cmp(&a.3)) // then by degree descending
            .then_with(|| b.1.cmp(&a.1)) // then by length descending
    });

    seq_order
        .into_iter()
        .map(|(seq_id, _, _, _)| seq_id)
        .collect()
}

/// Result of processing a single anchor region
struct AnchorRegionResult {
    /// Depth intervals for output
    intervals: Vec<SparseDepthInterval>,
    /// All discovered regions to mark as processed: (seq_id, start, end)
    discovered_regions: Vec<(u32, i64, i64)>,
    /// The anchor seq_id
    anchor_seq_id: u32,
    /// The anchor sample_id
    anchor_sample_id: u16,
}

/// A single hit from the depth-specific raw-interval BFS.
/// Stores query and target coordinates plus orientation.
struct DepthBfsHit {
    query_id: u32,
    query_start: i64,
    query_end: i64,
    target_id: u32,
    target_start: i64,
    target_end: i64,
    is_reverse: bool,
}

/// Per-chunk BFS state used by `batch_depth_bfs`.
///
/// Each chunk maintains its own frontier queue and visited-ranges dedup table
/// independently. `batch_depth_bfs` drives all chunk states together, loading
/// each sub-index file only once per round regardless of how many chunks need
/// it, so peak memory is bounded to one sub-index at a time.
struct BfsChunkState {
    queue: std::collections::VecDeque<(u32, i64, i64, u16)>, // (target_id, start, end, depth)
    visited_ranges: FxHashMap<u32, SortedRanges>,
    results: Vec<DepthBfsHit>,
}

impl BfsChunkState {
    fn new(
        anchor_seq_id: u32,
        region_start: i64,
        region_end: i64,
        anchor_len: i64,
        min_transitive_len: i64,
    ) -> Self {
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = FxHashMap::default();
        let anchor_sorted = visited_ranges
            .entry(anchor_seq_id)
            .or_insert_with(|| SortedRanges::new(anchor_len, 0));
        let filtered = anchor_sorted.insert((region_start, region_end));

        let mut queue = std::collections::VecDeque::new();
        for (s, e) in filtered {
            if (e - s).abs() >= min_transitive_len {
                queue.push_back((anchor_seq_id, s, e, 0u16));
            }
        }

        Self {
            queue,
            visited_ranges,
            results: Vec::with_capacity(64),
        }
    }

    /// Process raw alignment results for one BFS hop, record hits, and enqueue
    /// next-level frontier items. Mirrors the inner loop of `depth_transitive_bfs`.
    fn process_hop(
        &mut self,
        current_target_id: u32,
        current_start: i64,
        current_end: i64,
        current_depth: u16,
        raw_alns: &[RawAlignmentInterval],
        seq_len_fn: &impl Fn(u32) -> i64,
        min_transitive_len: i64,
        min_distance_between_ranges: i64,
    ) {
        let mut next_ranges: Vec<(u32, i64, i64)> = Vec::new();

        for aln in raw_alns {
            if aln.query_id == current_target_id {
                continue;
            }

            let aln_target_start = aln.target_start;
            let aln_target_end = aln.target_end;

            let clipped_target_start = aln_target_start.max(current_start);
            let clipped_target_end = aln_target_end.min(current_end);
            if clipped_target_start >= clipped_target_end {
                continue;
            }

            let target_len_aln = aln_target_end - aln_target_start;
            let query_len_aln = aln.query_end - aln.query_start;

            let (clipped_query_start, clipped_query_end) = if target_len_aln > 0 {
                let off_s = clipped_target_start - aln_target_start;
                let off_e = clipped_target_end - aln_target_start;
                if aln.is_reverse {
                    let cqe = aln.query_end - proj_offset(off_s, query_len_aln, target_len_aln);
                    let cqs = aln.query_end - proj_offset(off_e, query_len_aln, target_len_aln);
                    (cqs.min(cqe), cqs.max(cqe))
                } else {
                    let cqs = aln.query_start + proj_offset(off_s, query_len_aln, target_len_aln);
                    let cqe = aln.query_start + proj_offset(off_e, query_len_aln, target_len_aln);
                    (cqs.min(cqe), cqs.max(cqe))
                }
            } else {
                (aln.query_start, aln.query_end)
            };

            self.results.push(DepthBfsHit {
                query_id: aln.query_id,
                query_start: clipped_query_start,
                query_end: clipped_query_end,
                target_id: current_target_id,
                target_start: clipped_target_start,
                target_end: clipped_target_end,
                is_reverse: aln.is_reverse,
            });

            let explore_start = aln.query_start.min(aln.query_end);
            let explore_end = aln.query_start.max(aln.query_end);
            let query_seq_len = seq_len_fn(aln.query_id);
            let ranges = self
                .visited_ranges
                .entry(aln.query_id)
                .or_insert_with(|| SortedRanges::new(query_seq_len, 0));

            let mut should_add = true;
            if min_distance_between_ranges > 0 {
                let (new_min, new_max) = (explore_start, explore_end);
                let idx = match ranges
                    .ranges
                    .binary_search_by_key(&new_min, |&(start, _)| start)
                {
                    Ok(i) => i,
                    Err(i) => i,
                };
                if idx > 0 {
                    let (_, prev_end) = ranges.ranges[idx - 1];
                    if (new_min - prev_end).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
                if should_add && idx < ranges.ranges.len() {
                    let (next_start, _) = ranges.ranges[idx];
                    if (next_start - new_max).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
            }

            if should_add {
                for (ns, ne) in ranges.insert((explore_start, explore_end)) {
                    if (ne - ns).abs() >= min_transitive_len {
                        next_ranges.push((aln.query_id, ns, ne));
                    }
                }
            }
        }

        // Sort and merge contiguous ranges before enqueueing
        if !next_ranges.is_empty() {
            next_ranges.sort_by_key(|(id, start, _)| (*id, *start));
            let mut write = 0;
            for read in 1..next_ranges.len() {
                if next_ranges[write].0 == next_ranges[read].0
                    && next_ranges[write].2 >= next_ranges[read].1
                {
                    next_ranges[write].2 = next_ranges[write].2.max(next_ranges[read].2);
                } else {
                    write += 1;
                    next_ranges.swap(write, read);
                }
            }
            next_ranges.truncate(write + 1);
            for (seq_id, start, end) in next_ranges {
                self.queue.push_back((seq_id, start, end, current_depth + 1));
            }
        }
    }
}

/// Batch BFS across all hub chunks simultaneously.
///
/// Instead of each chunk independently loading sub-index files (causing
/// T × sub_index_size peak memory for T concurrent rayon threads), this
/// function drives all chunk frontiers together: in each round, every pending
/// query is collected, grouped by alignment file via
/// `ImpgIndex::batch_query_raw_overlapping`, and each file is loaded exactly
/// once to answer all queries referencing it before being freed. Sequential
/// file processing keeps peak memory at one sub-index + all chunk states.
///
/// For hub-heavy workloads where many chunks share the same alignment files
/// (e.g., 50 × 5 MB chunks of the same hub chromosome all need sample F's
/// alignment file), the batch approach also reduces total I/O: F is read once
/// rather than 50 times.
///
/// Returns one `Vec<DepthBfsHit>` per input chunk, in the same order.
fn batch_depth_bfs(
    impg: &impl ImpgIndex,
    chunks: &[(u32, i64, i64)], // (anchor_seq_id, chunk_start, chunk_end)
    max_depth: u16,
    min_transitive_len: i64,
    min_distance_between_ranges: i64,
) -> Vec<Vec<DepthBfsHit>> {
    if chunks.is_empty() {
        return Vec::new();
    }

    let seq_len_fn = |id: u32| impg.seq_index().get_len_from_id(id).unwrap_or(0) as i64;

    let mut states: Vec<BfsChunkState> = chunks
        .iter()
        .map(|&(anchor_id, start, end)| {
            BfsChunkState::new(anchor_id, start, end, seq_len_fn(anchor_id), min_transitive_len)
        })
        .collect();

    let mut pending: Vec<(usize, u32, i64, i64, u16)> = Vec::new();
    let mut queries: Vec<(u32, i64, i64)> = Vec::new();
    loop {
        // Drain all active queue items across all chunk states.
        // This is safe because process_hop only enqueues new items at depth+1,
        // so everything in the queue right now belongs to the same BFS level
        // for each chunk. Draining the full queue then calling process_hop
        // preserves level-by-level BFS ordering per chunk.
        // Items that exceed max_depth are discarded (same semantics as the
        // sequential depth_transitive_bfs which calls `continue` on depth check).
        pending.clear();
        for (ci, state) in states.iter_mut().enumerate() {
            while let Some((target_id, start, end, depth)) = state.queue.pop_front() {
                if max_depth > 0 && depth >= max_depth {
                    continue;
                }
                pending.push((ci, target_id, start, end, depth));
            }
        }

        if pending.is_empty() {
            break;
        }

        // Build query list (one entry per pending BFS item)
        queries.clear();
        queries.extend(pending.iter().map(|&(_, target_id, start, end, _)| (target_id, start, end)));

        // Batch query: each alignment file is loaded at most once to serve
        // all queries that reference it, then immediately freed.
        let batch_results = impg.batch_query_raw_overlapping(&queries);

        // Distribute results back and advance each chunk's frontier
        for ((ci, target_id, start, end, depth), raw_alns) in
            pending.iter().zip(batch_results.iter())
        {
            states[*ci].process_hop(
                *target_id,
                *start,
                *end,
                *depth,
                raw_alns,
                &seq_len_fn,
                min_transitive_len,
                min_distance_between_ranges,
            );
        }
    }

    states.into_iter().map(|s| s.results).collect()
}

/// Sweep-line half of `process_anchor_region_transitive_raw`, operating on
/// pre-computed BFS hits instead of running BFS internally.
///
/// Used by the Phase 1 batch path: `batch_depth_bfs` produces hits for all
/// chunks in one coordinated pass (sequential file loading, O(sub_index_size)
/// peak memory), then this function runs the sweep-line in parallel across
/// all chunks (CPU-only, no sub-index loading).
fn process_anchor_region_transitive_raw_with_hits(
    hits: Vec<DepthBfsHit>,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
    global_used: &ConcurrentProcessedTracker,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::with_capacity(hits.len() + 1);
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::with_capacity(hits.len() + 1);
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    // Pass 1: build anchor coverage map from hop-0 hits (target on anchor)
    let mut seq_anchor_coverage: FxHashMap<u32, Vec<HopZeroSeg>> = FxHashMap::default();
    for hit in &hits {
        if hit.target_id == anchor_seq_id {
            let q_start = hit.query_start.min(hit.query_end);
            let q_end = hit.query_start.max(hit.query_end);
            let t_start = hit.target_start.min(hit.target_end);
            let t_end = hit.target_start.max(hit.target_end);
            seq_anchor_coverage
                .entry(hit.query_id)
                .or_default()
                .push((q_start, q_end, t_start, t_end));
        }
    }
    depth_trace!(
        "PASS1 stage=raw_with_hits anchor={} region={}-{} hits={} keys={}",
        anchor_seq_id,
        region_start,
        region_end,
        hits.len(),
        seq_anchor_coverage.len()
    );

    // Pass 2: process all hits for depth (identical to process_anchor_region_transitive_raw)
    for hit in &hits {
        if min_seq_length > 0
            && !seq_included
                .get(hit.query_id as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(hit.query_id);
        if query_sample_id == anchor_sample_id {
            continue;
        }

        let query_start = hit.query_start.min(hit.query_end);
        let query_end = hit.query_start.max(hit.query_end);

        let hit_t_start = hit.target_start.min(hit.target_end);
        let hit_t_end = hit.target_start.max(hit.target_end);
        let (full_a_start, full_a_end) = if hit.target_id == anchor_seq_id {
            (hit_t_start.max(region_start), hit_t_end.min(region_end))
        } else {
            let p = project_hop0_coords(
                seq_anchor_coverage.get(&hit.target_id),
                hit_t_start,
                hit_t_end,
                region_start,
                region_end,
            );
            depth_trace!(
                "HOP2 stage=raw_with_hits sample={} q_id={} t_id={} t={}-{} a={}-{}",
                query_sample_id,
                hit.query_id,
                hit.target_id,
                hit_t_start,
                hit_t_end,
                p.0,
                p.1
            );
            p
        };
        if full_a_start >= full_a_end {
            continue;
        }

        // Always add the hit to the sweep-line for correct depth counting. See
        // process_anchor_region_raw_streaming for the full rationale — gating
        // the sweep on claim_unprocessed undercounts samples when another anchor
        // has already claimed the query range.
        let (ua_start, ua_end) = inverse_map_query_to_target(
            full_a_start,
            full_a_end,
            query_start,
            query_end,
            query_start,
            query_end,
            hit.is_reverse,
        );
        if ua_start < ua_end {
            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                hit.query_id,
                query_start,
                query_end,
                ua_start,
                ua_end,
                hit.is_reverse,
            ));
        }

        let claimed = global_used.claim_unprocessed(hit.query_id, query_start, query_end);
        for (uq_start, uq_end) in claimed {
            discovered_regions.push((hit.query_id, uq_start, uq_end));
        }
    }

    #[cfg(feature = "depth-trace")]
    {
        let unique: FxHashSet<u16> = alignments.iter().map(|a| a.sample_id).collect();
        depth_trace!(
            "SWEEP stage=raw_with_hits anchor={} region={}-{} alignments={} unique_samples={}",
            anchor_seq_id,
            region_start,
            region_end,
            alignments.len(),
            unique.len()
        );
    }
    let seq_intervals = sweep_line_depth(&alignments, num_samples, region_start, region_end);

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Depth-specific raw-interval BFS using `query_raw_overlapping()` at each hop.
///
/// Unlike the CIGAR-based BFS in `impg.rs`, this uses raw alignment extents with
/// linear interpolation to clip query coordinates. This is faster (no CIGAR walk)
/// and avoids chunk boundary gaps from CIGAR-precise projection, making it better
/// suited for depth computation where base-level precision is not needed.
///
/// Returns all discovered hits (excluding self-referential alignments).
fn depth_transitive_bfs(
    impg: &impl ImpgIndex,
    target_id: u32,
    range_start: i64,
    range_end: i64,
    max_depth: u16,
    min_transitive_len: i64,
    min_distance_between_ranges: i64,
    _use_dfs: bool,
) -> Vec<DepthBfsHit> {
    use std::collections::VecDeque;

    let mut results: Vec<DepthBfsHit> = Vec::new();

    // Lazily allocated visited ranges per sequence (not pre-allocated for all sequences)
    let mut visited_ranges: FxHashMap<u32, SortedRanges> = FxHashMap::default();

    // Initialize visited for the starting target
    let target_len = impg.seq_index().get_len_from_id(target_id).unwrap_or(0) as i64;
    let target_sorted = visited_ranges
        .entry(target_id)
        .or_insert_with(|| SortedRanges::new(target_len, 0));
    let filtered_input_range = target_sorted.insert((range_start, range_end));

    // BFS queue: (seq_id, start, end, current_depth)
    let mut queue: VecDeque<(u32, i64, i64, u16)> = VecDeque::new();
    for (s, e) in filtered_input_range {
        if (e - s).abs() >= min_transitive_len {
            queue.push_back((target_id, s, e, 0));
        }
    }

    while let Some((current_target_id, current_start, current_end, current_depth)) =
        queue.pop_front()
    {
        if max_depth > 0 && current_depth >= max_depth {
            continue;
        }

        let raw_alns = impg.query_raw_overlapping_transient(current_target_id, current_start, current_end);

        // Collect next-depth ranges to sort and merge before adding to queue
        let mut next_ranges: Vec<(u32, i64, i64)> = Vec::new();

        for aln in &raw_alns {
            // Skip self-referential (same sequence)
            if aln.query_id == current_target_id {
                continue;
            }

            let aln_target_start = aln.target_start;
            let aln_target_end = aln.target_end;

            // Clip target to current range
            let clipped_target_start = aln_target_start.max(current_start);
            let clipped_target_end = aln_target_end.min(current_end);
            if clipped_target_start >= clipped_target_end {
                continue;
            }

            // Linear interpolation to compute query coordinates
            let target_len_aln = aln_target_end - aln_target_start;
            let query_len_aln = aln.query_end - aln.query_start;

            let (clipped_query_start, clipped_query_end) = if target_len_aln > 0 {
                let off_s = clipped_target_start - aln_target_start;
                let off_e = clipped_target_end - aln_target_start;
                if aln.is_reverse {
                    let cqe = aln.query_end - proj_offset(off_s, query_len_aln, target_len_aln);
                    let cqs = aln.query_end - proj_offset(off_e, query_len_aln, target_len_aln);
                    (cqs.min(cqe), cqs.max(cqe))
                } else {
                    let cqs = aln.query_start + proj_offset(off_s, query_len_aln, target_len_aln);
                    let cqe = aln.query_start + proj_offset(off_e, query_len_aln, target_len_aln);
                    (cqs.min(cqe), cqs.max(cqe))
                }
            } else {
                (aln.query_start, aln.query_end)
            };

            // Record this hit (clipped coordinates for depth calculation)
            results.push(DepthBfsHit {
                query_id: aln.query_id,
                query_start: clipped_query_start,
                query_end: clipped_query_end,
                target_id: current_target_id,
                target_start: clipped_target_start,
                target_end: clipped_target_end,
                is_reverse: aln.is_reverse,
            });

            // For BFS exploration: use FULL alignment query extent (not clipped).
            // Linear interpolation can underestimate query range when indels are present,
            // causing hop 2+ to miss alignments at range boundaries. Using the full
            // extent guarantees we explore at least as widely as CIGAR-precise BFS.
            // visited_ranges prevents re-exploration, bounding the cost.
            let explore_start = aln.query_start.min(aln.query_end);
            let explore_end = aln.query_start.max(aln.query_end);

            let query_seq_len = impg.seq_index().get_len_from_id(aln.query_id).unwrap_or(0) as i64;
            let ranges = visited_ranges
                .entry(aln.query_id)
                .or_insert_with(|| SortedRanges::new(query_seq_len, 0));

            // Check proximity to existing ranges
            let mut should_add = true;
            if min_distance_between_ranges > 0 {
                let (new_min, new_max) = (explore_start, explore_end);
                let idx = match ranges
                    .ranges
                    .binary_search_by_key(&new_min, |&(start, _)| start)
                {
                    Ok(i) => i,
                    Err(i) => i,
                };

                if idx > 0 {
                    let (_, prev_end) = ranges.ranges[idx - 1];
                    if (new_min - prev_end).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
                if should_add && idx < ranges.ranges.len() {
                    let (next_start, _) = ranges.ranges[idx];
                    if (next_start - new_max).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
            }

            if should_add {
                let new_ranges = ranges.insert((explore_start, explore_end));
                for (new_start, new_end) in new_ranges {
                    if (new_end - new_start).abs() >= min_transitive_len {
                        next_ranges.push((aln.query_id, new_start, new_end));
                    }
                }
            }
        }

        // Sort and merge contiguous ranges before adding to queue
        if !next_ranges.is_empty() {
            next_ranges.sort_by_key(|(id, start, _)| (*id, *start));

            let mut write = 0;
            for read in 1..next_ranges.len() {
                if next_ranges[write].0 == next_ranges[read].0
                    && next_ranges[write].2 >= next_ranges[read].1
                {
                    next_ranges[write].2 = next_ranges[write].2.max(next_ranges[read].2);
                } else {
                    write += 1;
                    next_ranges.swap(write, read);
                }
            }
            next_ranges.truncate(write + 1);

            for (seq_id, start, end) in next_ranges {
                queue.push_back((seq_id, start, end, current_depth + 1));
            }
        }
    }

    results
}

/// Process a single anchor region using raw-interval BFS for transitive depth.
///
/// This is the default transitive path. Uses `depth_transitive_bfs()` which
/// operates on raw alignment extents with linear interpolation, avoiding
/// CIGAR walk overhead and chunk boundary gaps.
fn process_anchor_region_transitive_raw(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
    global_used: &ConcurrentProcessedTracker,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let hits = depth_transitive_bfs(
        impg,
        anchor_seq_id,
        region_start,
        region_end,
        config.max_depth,
        config.min_transitive_len,
        config.min_distance_between_ranges,
        config.transitive_dfs,
    );

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    // Pass 1: Build anchor coverage map from hop 0 results (target on anchor)
    let mut seq_anchor_coverage: FxHashMap<u32, Vec<HopZeroSeg>> = FxHashMap::default();
    for hit in &hits {
        if hit.target_id == anchor_seq_id {
            let q_start = hit.query_start.min(hit.query_end) as i64;
            let q_end = hit.query_start.max(hit.query_end) as i64;
            let t_start = hit.target_start.min(hit.target_end) as i64;
            let t_end = hit.target_start.max(hit.target_end) as i64;
            seq_anchor_coverage
                .entry(hit.query_id)
                .or_default()
                .push((q_start, q_end, t_start, t_end));
        }
    }
    depth_trace!(
        "PASS1 stage=transitive_raw anchor={} region={}-{} hits={} keys={}",
        anchor_seq_id,
        region_start,
        region_end,
        hits.len(),
        seq_anchor_coverage.len()
    );

    // Pass 2: Process all hits for depth
    for hit in &hits {
        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0
            && !seq_included
                .get(hit.query_id as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(hit.query_id);

        // Skip self-alignment (same sample)
        if query_sample_id == anchor_sample_id {
            continue;
        }

        let query_start = hit.query_start.min(hit.query_end) as i64;
        let query_end = hit.query_start.max(hit.query_end) as i64;

        // Compute full anchor coordinates for this hit (existing logic)
        let hit_t_start = hit.target_start.min(hit.target_end) as i64;
        let hit_t_end = hit.target_start.max(hit.target_end) as i64;
        let (full_a_start, full_a_end) = if hit.target_id == anchor_seq_id {
            (hit_t_start.max(region_start), hit_t_end.min(region_end))
        } else {
            let p = project_hop0_coords(
                seq_anchor_coverage.get(&hit.target_id),
                hit_t_start, hit_t_end,
                region_start, region_end,
            );
            depth_trace!(
                "HOP2 stage=transitive_raw sample={} q_id={} t_id={} t={}-{} a={}-{}",
                query_sample_id, hit.query_id, hit.target_id,
                hit_t_start, hit_t_end, p.0, p.1
            );
            p
        };
        if full_a_start >= full_a_end {
            continue;
        }

        // Always add the hit to the sweep-line for correct depth counting. See
        // process_anchor_region_raw_streaming for the full rationale.
        let (ua_start, ua_end) = inverse_map_query_to_target(
            full_a_start, full_a_end,
            query_start, query_end,
            query_start, query_end,
            hit.is_reverse,
        );
        if ua_start < ua_end {
            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                hit.query_id,
                query_start,
                query_end,
                ua_start,
                ua_end,
                hit.is_reverse,
            ));
        }

        let claimed = global_used.claim_unprocessed(hit.query_id, query_start, query_end);
        for (uq_start, uq_end) in claimed {
            discovered_regions.push((hit.query_id, uq_start, uq_end));
        }
    }

    #[cfg(feature = "depth-trace")]
    {
        let unique: FxHashSet<u16> = alignments.iter().map(|a| a.sample_id).collect();
        depth_trace!(
            "SWEEP stage=transitive_raw anchor={} region={}-{} alignments={} unique_samples={}",
            anchor_seq_id, region_start, region_end,
            alignments.len(), unique.len()
        );
    }
    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(&alignments, num_samples, region_start, region_end);

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region: query alignments, compute depth via sweep-line.
/// Returns depth intervals and all discovered regions (for marking as processed).
///
/// Query strategy:
/// - Non-transitive (default): direct query — O(log n), bounded results.
///   impg's bidirectional index means 1-hop already captures both alignment directions.
/// - Transitive (-x): routes to raw BFS (default) or CIGAR BFS (--use-BFS).
///   Caller should pass bounded-size regions (see TRANSITIVE_CHUNK_SIZE).
fn process_anchor_region(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
    global_used: &ConcurrentProcessedTracker,
) -> AnchorRegionResult {
    let is_transitive = config.transitive || config.transitive_dfs;

    // Transitive mode: route to raw BFS (default) or CIGAR BFS (--use-BFS)
    if is_transitive {
        if config.use_cigar_bfs {
            return process_anchor_region_transitive_cigar(
                impg,
                config,
                compact_lengths,
                num_samples,
                anchor_seq_id,
                anchor_sample_id,
                region_start,
                region_end,
                seq_included,
                min_seq_length,
                global_used,
            );
        } else {
            return process_anchor_region_transitive_raw(
                impg,
                config,
                compact_lengths,
                num_samples,
                anchor_seq_id,
                anchor_sample_id,
                region_start,
                region_end,
                seq_included,
                min_seq_length,
                global_used,
            );
        }
    }

    // Non-transitive mode: direct 1-hop query + sweep-line
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let overlaps = impg.query(
        anchor_seq_id,
        region_start,
        region_end,
        false,
        None,
        None,
        false,
    );

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Add self (anchor sample covers the range)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_id = query_interval.metadata;

        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0
            && !seq_included
                .get(query_id as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(query_id);

        if is_self_alignment(query_sample_id, anchor_sample_id, query_id, anchor_seq_id) {
            continue;
        }

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;
        let target_start = target_interval.first.min(target_interval.last) as i64;
        let target_end = target_interval.first.max(target_interval.last) as i64;

        // Clip target to range boundaries
        let clipped_target_start = target_start.max(region_start);
        let clipped_target_end = target_end.min(region_end);

        if clipped_target_start >= clipped_target_end {
            continue;
        }

        // Proportionally adjust query coordinates for clipped target
        let target_len = target_end - target_start;
        let query_len = query_end - query_start;
        let (cq_start, cq_end) = if target_len > 0 {
            let off_s = clipped_target_start - target_start;
            let off_e = clipped_target_end - target_start;
            if is_reverse {
                let cqe = query_end - proj_offset(off_s, query_len, target_len);
                let cqs = query_end - proj_offset(off_e, query_len, target_len);
                (cqs.min(cqe), cqs.max(cqe))
            } else {
                let cqs = query_start + proj_offset(off_s, query_len, target_len);
                let cqe = query_start + proj_offset(off_e, query_len, target_len);
                (cqs.min(cqe), cqs.max(cqe))
            }
        } else {
            (query_start, query_end)
        };

        alignments.push(CompactAlignmentInfo::new(
            query_sample_id,
            query_id,
            cq_start,
            cq_end,
            clipped_target_start,
            clipped_target_end,
            is_reverse,
        ));

        discovered_regions.push((query_id, query_start, query_end));
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(&alignments, num_samples, region_start, region_end);

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region using pre-scanned raw alignment data.
/// No CIGAR projection — uses raw coordinates from interval tree metadata.
/// This is the fast path for non-transitive depth computation.
fn process_anchor_region_raw(
    raw_intervals: &[RawAlignmentInterval],
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    global_used: &ConcurrentProcessedTracker,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::with_capacity(raw_intervals.len() + 1);
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let mut alignments: Vec<CompactAlignmentInfo> =
        Vec::with_capacity(raw_intervals.len() + 1);

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    // Per-task buffer: collects (orig_idx, query_id, cq_start, cq_end). After
    // the loop we group runs of equal query_id and call claim_*_batch once per
    // group, holding the per-seq mutex for the whole group instead of per-aln.
    let mut claim_buf: Vec<(usize, u32, i64, i64)> = Vec::with_capacity(raw_intervals.len());

    // Binary search: only scan intervals where target_start < region_end
    // (alignment_table is sorted by target_start after pre-scan)
    let end_idx = raw_intervals.partition_point(|aln| (aln.target_start as i64) < region_end);
    for aln in &raw_intervals[..end_idx] {
        let target_start = aln.target_start as i64;
        let target_end = aln.target_end as i64;

        // Skip intervals that end before region starts
        if target_end <= region_start {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(aln.query_id);

        // Filter self-alignments (same sample, different sequence)
        if is_self_alignment(
            query_sample_id,
            anchor_sample_id,
            aln.query_id,
            anchor_seq_id,
        ) {
            continue;
        }

        let query_start = aln.query_start as i64;
        let query_end = aln.query_end as i64;

        // Stage 2 J: do NOT pre-clip target/query for storage. The sweep-line
        // clips per-emit-interval, performing a single projection + single
        // rounding (instead of two compounded roundings). We compute the
        // clipped query range here ONLY for the ownership claim below.
        let clipped_t_start = target_start.max(region_start);
        let clipped_t_end = target_end.min(region_end);
        let target_len = target_end - target_start;
        let query_len = query_end - query_start;
        let (cq_start, cq_end) = if target_len > 0 {
            let off_s = clipped_t_start - target_start;
            let off_e = clipped_t_end - target_start;
            if aln.is_reverse {
                let qe = query_end - proj_offset(off_s, query_len, target_len);
                let qs = query_end - proj_offset(off_e, query_len, target_len);
                (qs.min(qe), qs.max(qe))
            } else {
                let qs = query_start + proj_offset(off_s, query_len, target_len);
                let qe = query_start + proj_offset(off_e, query_len, target_len);
                (qs.min(qe), qs.max(qe))
            }
        } else {
            (query_start, query_end)
        };

        // Always add the alignment to the sweep-line for correct depth counting.
        // See the streaming variant for the full rationale — gating the sweep on
        // claim_unprocessed causes samples that genuinely cover the anchor region
        // to be omitted when their query range was claimed by a prior anchor.
        //
        // Stage 2 J: store ORIGINAL alignment coords (not pre-clipped). Sweep-line
        // events at aln.target_start / aln.target_end may fall outside
        // [region_start, region_end), but the sweep already clips per-emit interval,
        // so out-of-region events trigger no emit while still toggling sample-active
        // state correctly across the region boundary.
        let aln_idx = alignments.len();
        alignments.push(CompactAlignmentInfo::new(
            query_sample_id,
            aln.query_id,
            query_start,
            query_end,
            target_start,
            target_end,
            aln.is_reverse,
        ));

        claim_buf.push((aln_idx, aln.query_id, cq_start, cq_end));
    }

    // Stage 1 F: batch claims by query_id. Stable sort_by_key keeps secondary
    // idx ascending so within-group call order matches iteration order; the
    // post-pass sort by idx restores discovered_regions to original order.
    if !claim_buf.is_empty() {
        claim_buf.sort_by_key(|&(idx, qid, _, _)| (qid, idx));
        let mut group_buf: Vec<(i64, i64)> = Vec::new();
        let mut bool_buf: Vec<bool> = Vec::new();
        let mut newly_owned: Vec<(usize, u32, i64, i64)> = Vec::new();
        let mut i = 0;
        while i < claim_buf.len() {
            let qid = claim_buf[i].1;
            let mut j = i;
            while j < claim_buf.len() && claim_buf[j].1 == qid {
                j += 1;
            }
            group_buf.clear();
            group_buf.extend(claim_buf[i..j].iter().map(|&(_, _, s, e)| (s, e)));
            bool_buf.clear();
            global_used.claim_any_unprocessed_batch(qid, &group_buf, &mut bool_buf);
            for (k, &b) in bool_buf.iter().enumerate() {
                if b {
                    let (idx, _, s, e) = claim_buf[i + k];
                    newly_owned.push((idx, qid, s, e));
                }
            }
            i = j;
        }
        newly_owned.sort_by_key(|&(idx, _, _, _)| idx);
        for (_, qid, s, e) in newly_owned {
            discovered_regions.push((qid, s, e));
        }
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(&alignments, num_samples, region_start, region_end);

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region using query_transitive_bfs/dfs for discovery,
/// then map results to anchor coordinates for the depth sweep-line.
///
/// Uses the same BFS as the query command (guaranteeing identical sample discovery).
/// Anchor coordinate mapping:
/// - Hop 0 results (target on anchor): precise CIGAR-projected coordinates
/// - Hop 1+ results (target on intermediate): projected back via the intermediate's
///   known anchor coverage from hop 0 results (linear interpolation)
/// - Deeper hops with unknown intermediate: full anchor region as fallback
fn process_anchor_region_transitive_cigar(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
    global_used: &ConcurrentProcessedTracker,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    // Use the same CIGAR-precise BFS/DFS as the query command (--use-BFS path)
    let overlaps = if config.transitive_dfs {
        impg.query_transitive_dfs(
            anchor_seq_id,
            region_start,
            region_end,
            None,
            config.max_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            None,
            false,
            None,
            None,
            false,
            None,
        )
    } else {
        impg.query_transitive_bfs(
            anchor_seq_id,
            region_start,
            region_end,
            None,
            config.max_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            None,
            false,
            None,
            None,
            false,
            None,
        )
    };

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    // Pass 1: Build anchor coverage map from hop 0 results (target on anchor)
    let mut seq_anchor_coverage: FxHashMap<u32, Vec<HopZeroSeg>> = FxHashMap::default();
    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;
        if target_interval.metadata == anchor_seq_id {
            let query_id = query_interval.metadata;
            let q_start = query_interval.first.min(query_interval.last) as i64;
            let q_end = query_interval.first.max(query_interval.last) as i64;
            let t_start = target_interval.first.min(target_interval.last) as i64;
            let t_end = target_interval.first.max(target_interval.last) as i64;
            seq_anchor_coverage
                .entry(query_id)
                .or_default()
                .push((q_start, q_end, t_start, t_end));
        }
    }
    depth_trace!(
        "PASS1 stage=transitive_cigar anchor={} region={}-{} overlaps={} keys={}",
        anchor_seq_id,
        region_start,
        region_end,
        overlaps.len(),
        seq_anchor_coverage.len()
    );

    // Pass 2: Process all results for depth
    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_id = query_interval.metadata;

        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0
            && !seq_included
                .get(query_id as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(query_id);

        // Skip self-alignment (same sample)
        if query_sample_id == anchor_sample_id {
            continue;
        }

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;

        // Compute full anchor coordinates for this overlap (existing logic)
        let hit_t_start = target_interval.first.min(target_interval.last) as i64;
        let hit_t_end = target_interval.first.max(target_interval.last) as i64;
        let (full_a_start, full_a_end) = if target_interval.metadata == anchor_seq_id {
            (hit_t_start.max(region_start), hit_t_end.min(region_end))
        } else {
            let p = project_hop0_coords(
                seq_anchor_coverage.get(&target_interval.metadata),
                hit_t_start, hit_t_end,
                region_start, region_end,
            );
            depth_trace!(
                "HOP2 stage=transitive_cigar sample={} q_id={} t_id={} t={}-{} a={}-{}",
                query_sample_id, query_id, target_interval.metadata,
                hit_t_start, hit_t_end, p.0, p.1
            );
            p
        };
        if full_a_start >= full_a_end {
            continue;
        }

        // Always add the hit to the sweep-line for correct depth counting. See
        // process_anchor_region_raw_streaming for the full rationale.
        let (ua_start, ua_end) = inverse_map_query_to_target(
            full_a_start, full_a_end,
            query_start, query_end,
            query_start, query_end,
            is_reverse,
        );
        if ua_start < ua_end {
            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                query_id,
                query_start,
                query_end,
                ua_start,
                ua_end,
                is_reverse,
            ));
        }

        let claimed = global_used.claim_unprocessed(query_id, query_start, query_end);
        for (uq_start, uq_end) in claimed {
            discovered_regions.push((query_id, uq_start, uq_end));
        }
    }

    // Augment discovered_regions with raw alignment extents for direct (hop 0) overlaps.
    // The BFS loop above records CIGAR-projected sub-ranges (line 1899) which may leave
    // gaps at chunk boundaries due to indels. Raw extents ensure the full alignment
    // coverage is marked as processed, preventing Phase 2 from re-processing these
    // regions and producing duplicate output (e.g., CHM13 appearing in Phase 2 rows).
    let raw_hop0 =
        impg.query_raw_overlapping_transient(anchor_seq_id, region_start, region_end);
    for aln in &raw_hop0 {
        if min_seq_length > 0
            && !seq_included
                .get(aln.query_id as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }
        let query_sample_id = compact_lengths.get_sample_id(aln.query_id);
        if is_self_alignment(
            query_sample_id,
            anchor_sample_id,
            aln.query_id,
            anchor_seq_id,
        ) {
            continue;
        }
        discovered_regions.push((aln.query_id, aln.query_start as i64, aln.query_end as i64));
    }

    #[cfg(feature = "depth-trace")]
    {
        let unique: FxHashSet<u16> = alignments.iter().map(|a| a.sample_id).collect();
        depth_trace!(
            "SWEEP stage=transitive_cigar anchor={} region={}-{} alignments={} unique_samples={}",
            anchor_seq_id, region_start, region_end,
            alignments.len(), unique.len()
        );
    }
    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(&alignments, num_samples, region_start, region_end);

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Compute the total length of the union of a set of intervals.
/// Intervals may overlap; overlapping regions are counted only once.
fn interval_union_length(intervals: &mut Vec<(i64, i64)>) -> i64 {
    if intervals.is_empty() {
        return 0;
    }
    if intervals.len() == 1 {
        return (intervals[0].1 - intervals[0].0).max(0);
    }
    intervals.sort_unstable_by_key(|&(s, _)| s);
    let mut total = 0i64;
    let mut merged_start = intervals[0].0;
    let mut merged_end = intervals[0].1;
    for &(s, e) in &intervals[1..] {
        if s <= merged_end {
            merged_end = merged_end.max(e);
        } else {
            total += merged_end - merged_start;
            merged_start = s;
            merged_end = e;
        }
    }
    total += merged_end - merged_start;
    total.max(0)
}

/// Sweep-line algorithm: given alignments, produce depth intervals with sample tracking.
/// Callers must include the anchor sample as an alignment covering [region_start, region_end]
/// so that depth naturally counts all samples (including the anchor itself).
fn sweep_line_depth(
    alignments: &[CompactAlignmentInfo],
    num_samples: usize,
    region_start: i64,
    region_end: i64,
) -> Vec<SparseDepthInterval> {
    // Build sweep-line events
    let mut events: Vec<CompactDepthEvent> = Vec::with_capacity(alignments.len() * 2);
    for (idx, aln) in alignments.iter().enumerate() {
        events.push(CompactDepthEvent {
            position: aln.target_start,
            is_start: true,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
        events.push(CompactDepthEvent {
            position: aln.target_end,
            is_start: false,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
    }
    events.sort_by_key(|e| e.packed_sort_key());

    // Sweep-line to compute depth intervals with SPARSE storage
    let mut seq_intervals: Vec<SparseDepthInterval> = Vec::new();
    let mut active_bitmap = SampleBitmap::new(num_samples);
    let mut active_alns: Vec<Vec<usize>> = vec![Vec::new(); num_samples];
    let mut prev_pos: Option<i64> = None;
    // Hoisted to avoid repeated HashMap allocation; cleared before each per-sample use.
    let mut query_intervals_by_contig: FxHashMap<u32, Vec<(i64, i64)>> = FxHashMap::default();

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev && active_bitmap.depth() > 0 {
                let interval_start = prev;
                let interval_end = event.position;

                // Clip to anchor region
                let clipped_start = interval_start.max(region_start);
                let clipped_end = interval_end.min(region_end);

                if clipped_start < clipped_end {
                    // Build SPARSE sample positions + compute pangenome bases
                    let mut samples: Vec<SamplePosition> =
                        Vec::with_capacity(active_bitmap.depth());
                    let mut pangenome_bases: i64 = 0;

                    for sample_id in active_bitmap.active_samples() {
                        let alns = &active_alns[sample_id as usize];

                        // Compute query projections for ALL active alignments of this sample,
                        // grouped by query_id (contig) to correctly compute union lengths.
                        query_intervals_by_contig.clear();
                        let mut best_idx: Option<usize> = None;
                        let mut best_overlap: i64 = -1;
                        let mut best_qs: i64 = 0;
                        let mut best_qe: i64 = 0;

                        for &idx in alns {
                            let aln = &alignments[idx];
                            let overlap_start = clipped_start.max(aln.target_start);
                            let overlap_end = clipped_end.min(aln.target_end);
                            let overlap = overlap_end - overlap_start;

                            if overlap <= 0 {
                                continue;
                            }

                            let (q_start, q_end) = map_target_to_query_linear(
                                &[],
                                aln.target_start,
                                aln.target_end,
                                aln.query_start,
                                aln.query_end,
                                clipped_start,
                                clipped_end,
                                aln.is_reverse,
                            );
                            query_intervals_by_contig
                                .entry(aln.query_id)
                                .or_default()
                                .push((q_start, q_end));

                            if overlap > best_overlap {
                                best_overlap = overlap;
                                best_idx = Some(idx);
                                best_qs = q_start;
                                best_qe = q_end;
                            }
                        }

                        // Best alignment for TSV output: reuse the (q_start, q_end)
                        // already computed above instead of mapping a second time.
                        if let Some(idx) = best_idx {
                            let aln = &alignments[idx];
                            samples.push((sample_id, aln.query_id, best_qs, best_qe));
                        }

                        // Pangenome bases: union of query projections per contig
                        for (_, intervals) in query_intervals_by_contig.iter_mut() {
                            pangenome_bases += interval_union_length(intervals);
                        }
                    }

                    if !samples.is_empty() {
                        // active_samples() yields IDs via iter_ones() in ascending order,
                        // so samples is already sorted — no sort needed.
                        seq_intervals.push(SparseDepthInterval {
                            start: clipped_start,
                            end: clipped_end,
                            samples,
                            pangenome_bases,
                        });
                    }
                }
            }
        }

        // Update active samples
        if event.is_start {
            active_bitmap.add(event.sample_id);
            active_alns[event.sample_id as usize].push(event.alignment_idx);
        } else {
            active_bitmap.remove(event.sample_id);
            let v = &mut active_alns[event.sample_id as usize];
            if let Some(pos) = v.iter().position(|&idx| idx == event.alignment_idx) {
                v.swap_remove(pos);
            }
        }

        prev_pos = Some(event.position);
    }

    seq_intervals
}

/// Streaming variant of [`sweep_line_depth`]. Calls `emit` for each depth interval
/// instead of collecting into a Vec. Peak memory is O(N_samples) per interval
/// rather than O(K × N_samples) for the full result set.
fn sweep_line_depth_streaming(
    alignments: &[CompactAlignmentInfo],
    num_samples: usize,
    region_start: i64,
    region_end: i64,
    emit: &mut impl FnMut(SparseDepthInterval),
) {
    // Build sweep-line events
    let mut events: Vec<CompactDepthEvent> = Vec::with_capacity(alignments.len() * 2);
    for (idx, aln) in alignments.iter().enumerate() {
        events.push(CompactDepthEvent {
            position: aln.target_start,
            is_start: true,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
        events.push(CompactDepthEvent {
            position: aln.target_end,
            is_start: false,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
    }
    events.sort_by_key(|e| e.packed_sort_key());

    let mut active_bitmap = SampleBitmap::new(num_samples);
    let mut active_alns: Vec<Vec<usize>> = vec![Vec::new(); num_samples];
    let mut prev_pos: Option<i64> = None;
    // Hoisted to avoid repeated HashMap allocation; cleared before each per-sample use.
    let mut query_intervals_by_contig: FxHashMap<u32, Vec<(i64, i64)>> = FxHashMap::default();

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev && active_bitmap.depth() > 0 {
                let clipped_start = prev.max(region_start);
                let clipped_end = event.position.min(region_end);

                if clipped_start < clipped_end {
                    let mut samples: Vec<SamplePosition> =
                        Vec::with_capacity(active_bitmap.depth());
                    let mut pangenome_bases: i64 = 0;

                    for sample_id in active_bitmap.active_samples() {
                        let alns = &active_alns[sample_id as usize];
                        query_intervals_by_contig.clear();
                        let mut best_idx: Option<usize> = None;
                        let mut best_overlap: i64 = -1;
                        let mut best_qs: i64 = 0;
                        let mut best_qe: i64 = 0;

                        for &idx in alns {
                            let aln = &alignments[idx];
                            let overlap_start = clipped_start.max(aln.target_start);
                            let overlap_end = clipped_end.min(aln.target_end);
                            let overlap = overlap_end - overlap_start;
                            if overlap <= 0 {
                                continue;
                            }

                            let (q_start, q_end) = map_target_to_query_linear(
                                &[],
                                aln.target_start,
                                aln.target_end,
                                aln.query_start,
                                aln.query_end,
                                clipped_start,
                                clipped_end,
                                aln.is_reverse,
                            );
                            query_intervals_by_contig
                                .entry(aln.query_id)
                                .or_default()
                                .push((q_start, q_end));
                            if overlap > best_overlap {
                                best_overlap = overlap;
                                best_idx = Some(idx);
                                best_qs = q_start;
                                best_qe = q_end;
                            }
                        }

                        if let Some(idx) = best_idx {
                            let aln = &alignments[idx];
                            samples.push((sample_id, aln.query_id, best_qs, best_qe));
                        }

                        for (_, intervals) in query_intervals_by_contig.iter_mut() {
                            pangenome_bases += interval_union_length(intervals);
                        }
                    }

                    if !samples.is_empty() {
                        // active_samples() yields IDs via iter_ones() in ascending order,
                        // so samples is already sorted — no sort needed.
                        emit(SparseDepthInterval {
                            start: clipped_start,
                            end: clipped_end,
                            samples,
                            pangenome_bases,
                        });
                    }
                }
            }
        }

        if event.is_start {
            active_bitmap.add(event.sample_id);
            active_alns[event.sample_id as usize].push(event.alignment_idx);
        } else {
            active_bitmap.remove(event.sample_id);
            let v = &mut active_alns[event.sample_id as usize];
            if let Some(pos) = v.iter().position(|&idx| idx == event.alignment_idx) {
                v.swap_remove(pos);
            }
        }
        prev_pos = Some(event.position);
    }
}

/// Streaming variant of [`process_anchor_region_raw`]. Returns only discovered_regions,
/// streaming depth intervals through `emit` for immediate processing.
fn process_anchor_region_raw_streaming(
    raw_intervals: &[RawAlignmentInterval],
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    global_used: &ConcurrentProcessedTracker,
    emit: &mut impl FnMut(SparseDepthInterval),
) -> Vec<(u32, i64, i64)> {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::with_capacity(raw_intervals.len() + 1);
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let mut alignments: Vec<CompactAlignmentInfo> =
        Vec::with_capacity(raw_intervals.len() + 1);
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id,
        anchor_seq_id,
        region_start,
        region_end,
        region_start,
        region_end,
        false,
    ));

    // Per-task buffer: see process_anchor_region_raw for rationale (Stage 1 F).
    let mut claim_buf: Vec<(usize, u32, i64, i64)> = Vec::with_capacity(raw_intervals.len());

    let end_idx = raw_intervals.partition_point(|aln| (aln.target_start as i64) < region_end);
    for aln in &raw_intervals[..end_idx] {
        let target_start = aln.target_start as i64;
        let target_end = aln.target_end as i64;
        if target_end <= region_start {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(aln.query_id);
        if is_self_alignment(
            query_sample_id,
            anchor_sample_id,
            aln.query_id,
            anchor_seq_id,
        ) {
            continue;
        }

        let query_start = aln.query_start as i64;
        let query_end = aln.query_end as i64;

        // Stage 2 J: clipped query coords used ONLY for the ownership claim;
        // alignment is stored with original coords so sweep-line performs a
        // single projection + single rounding per emit.
        let clipped_t_start = target_start.max(region_start);
        let clipped_t_end = target_end.min(region_end);
        let target_len = target_end - target_start;
        let query_len = query_end - query_start;
        let (cq_start, cq_end) = if target_len > 0 {
            let off_s = clipped_t_start - target_start;
            let off_e = clipped_t_end - target_start;
            if aln.is_reverse {
                let qe = query_end - proj_offset(off_s, query_len, target_len);
                let qs = query_end - proj_offset(off_e, query_len, target_len);
                (qs.min(qe), qs.max(qe))
            } else {
                let qs = query_start + proj_offset(off_s, query_len, target_len);
                let qe = query_start + proj_offset(off_e, query_len, target_len);
                (qs.min(qe), qs.max(qe))
            }
        } else {
            (query_start, query_end)
        };

        // Always add the alignment to the sweep-line for correct depth counting.
        // Depth at the anchor's target position is a function of sample coverage
        // and is independent of which anchor "owns" the query sub-range for output.
        // Using claim_unprocessed to gate the sweep causes samples that genuinely
        // cover the anchor region to be omitted when their query range was claimed
        // by a prior anchor, producing a silent undercount.
        //
        // Stage 2 J: store ORIGINAL alignment coords (not pre-clipped). Out-of-region
        // events toggle sample-active state but produce no emit (sweep clips per
        // emit interval), so the correct samples are active across [region_start,
        // region_end) without compounded clipping rounding error.
        let aln_idx = alignments.len();
        alignments.push(CompactAlignmentInfo::new(
            query_sample_id,
            aln.query_id,
            query_start,
            query_end,
            target_start,
            target_end,
            aln.is_reverse,
        ));

        claim_buf.push((aln_idx, aln.query_id, cq_start, cq_end));
    }

    // Stage 1 F: batched claim by query_id. See process_anchor_region_raw.
    if !claim_buf.is_empty() {
        claim_buf.sort_by_key(|&(idx, qid, _, _)| (qid, idx));
        let mut group_buf: Vec<(i64, i64)> = Vec::new();
        let mut bool_buf: Vec<bool> = Vec::new();
        let mut newly_owned: Vec<(usize, u32, i64, i64)> = Vec::new();
        let mut i = 0;
        while i < claim_buf.len() {
            let qid = claim_buf[i].1;
            let mut j = i;
            while j < claim_buf.len() && claim_buf[j].1 == qid {
                j += 1;
            }
            group_buf.clear();
            group_buf.extend(claim_buf[i..j].iter().map(|&(_, _, s, e)| (s, e)));
            bool_buf.clear();
            global_used.claim_any_unprocessed_batch(qid, &group_buf, &mut bool_buf);
            for (k, &b) in bool_buf.iter().enumerate() {
                if b {
                    let (idx, _, s, e) = claim_buf[i + k];
                    newly_owned.push((idx, qid, s, e));
                }
            }
            i = j;
        }
        newly_owned.sort_by_key(|&(idx, _, _, _)| idx);
        for (_, qid, s, e) in newly_owned {
            discovered_regions.push((qid, s, e));
        }
    }

    sweep_line_depth_streaming(&alignments, num_samples, region_start, region_end, emit);

    discovered_regions
}

/// Maximum region size for a single transitive BFS query (5 MB).
/// Prevents BFS from visiting the entire pangenome when starting from a large chromosome.
/// Non-transitive queries don't need this limit (they're O(log n) per lookup).
const TRANSITIVE_CHUNK_SIZE: i64 = 5_000_000;

/// Flush threshold for streaming depth output buffers (4 MB).
/// Each thread accumulates TSV text up to this limit before sending the
/// buffer over the writer channel, balancing memory usage against per-send
/// channel overhead.
const STREAMING_FLUSH_THRESHOLD: usize = 4 * 1024 * 1024;

/// Bounded MPSC channel capacity for the dedicated TSV writer thread.
///
/// 256 buffers × ~4 MB worker buffers ≈ 1 GB worst-case backpressure ceiling.
/// Workers block on `send` when the channel is full, providing natural
/// throttling without losing any rows.
const WRITER_CHANNEL_CAPACITY: usize = 256;

/// Dedicated writer thread + bounded MPSC channel for TSV (and, when
/// resumable, work-log) output.
///
/// Replaces the previous `Mutex<BufWriter>` design: at high thread counts the
/// shared mutex serialised every `write_all` call across all workers, even
/// though each worker had already done the formatting on its own 4 MB buffer.
/// Now workers `send_*(...)` (blocks if the bounded channel is full) and the
/// dedicated thread drains the channel into the underlying writer in order.
///
/// **Invariant**: the writer thread MUST be joined via `finish()` at end of
/// run. Dropping `tx` flushes the channel; the writer thread then exits, and
/// `finish()` propagates any IO error from the writer thread back to the
/// caller. Without this join the program may exit before the underlying file
/// is fully drained, producing a truncated TSV.
///
/// When checkpoint/resume is enabled the writer additionally drives a
/// `<prefix>.depth.work.bin` append-only log. The writer is the only thread
/// touching either file, so a per-chunk `WriterMsg::Chunk { tsv, work }` is
/// applied atomically: TSV bytes and the matching work-log record advance
/// together, and a subsequent `Barrier` returns the post-flush stream
/// positions of *both* files. That coupling is what lets the ckpt offset
/// pair stay consistent — anything past the captured offsets gets truncated
/// on resume regardless of whether it was a half-flushed TSV row or a
/// half-flushed work-log record.
enum WriterMsg {
    /// Plain TSV bytes with no matching work-log record. Used for the
    /// non-resume path and for boundary writes that don't correspond to a
    /// chunk.
    TsvOnly(Vec<u8>),
    /// Atomic per-chunk write. Both halves are pushed in lockstep.
    Chunk { tsv: Vec<u8>, work: Vec<u8> },
    /// Drain queued messages, flush + `sync_data` both files, return the
    /// (tsv_offset, work_offset) pair to the requester.
    Barrier(SyncSender<io::Result<(u64, u64)>>),
}

/// TSV destination: a file we can `sync_data` + query stream position on, or
/// stdout (no checkpointing — barrier returns dummies).
enum TsvSink {
    File(BufWriter<std::fs::File>),
    Stdout(BufWriter<std::io::Stdout>),
}

impl TsvSink {
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::File(w) => w.write_all(buf),
            Self::Stdout(w) => w.write_all(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::File(w) => w.flush(),
            Self::Stdout(w) => w.flush(),
        }
    }
    /// Returns the post-flush byte offset on a real file; 0 for stdout.
    fn flush_and_sync(&mut self) -> io::Result<u64> {
        match self {
            Self::File(w) => {
                w.flush()?;
                w.get_ref().sync_data()?;
                // BufWriter<File> implements Seek (via inner), but only after
                // a flush; we just flushed.
                w.stream_position()
            }
            Self::Stdout(w) => {
                w.flush()?;
                Ok(0)
            }
        }
    }
}

struct DepthWriter {
    tx: SyncSender<WriterMsg>,
    handle: JoinHandle<io::Result<()>>,
    /// True iff a work-log file backs this writer; gates `send_chunk_bundle`
    /// and meaningful barriers.
    has_work_log: bool,
}

/// Options for opening the writer; controls header-write and work-log paths.
struct WriterOpts<'a> {
    /// Output prefix; when None, TSV goes to stdout and resume is disabled.
    output_prefix: Option<&'a str>,
    /// Append to existing TSV (resume path) instead of truncating.
    append_tsv: bool,
    /// Append to existing work-log (resume path) instead of truncating.
    append_worklog: bool,
    /// If true, allocate and drive a work-log alongside the TSV. Implies
    /// `output_prefix.is_some()`.
    with_work_log: bool,
    /// Skip writing the TSV header. Set on resume so we don't duplicate it
    /// after truncating the file to `tsv_byte_offset`.
    skip_header: bool,
}

impl DepthWriter {
    fn open(opts: WriterOpts<'_>) -> io::Result<Self> {
        let WriterOpts {
            output_prefix,
            append_tsv,
            append_worklog,
            with_work_log,
            skip_header,
        } = opts;

        if with_work_log && output_prefix.is_none() {
            return Err(io::Error::other(
                "work-log requires --output-prefix (cannot resume to stdout)",
            ));
        }

        // Open TSV sink.
        //
        // Resume mode opens the existing file in `write` mode (no `truncate`,
        // no `O_APPEND`) and seeks to EOF so subsequent writes extend it.
        // We deliberately do *not* set the OS-level append flag: combined
        // with explicit `seek` it adds ambiguity on networked filesystems
        // and the single-writer-thread invariant already gives us linear
        // append semantics.
        let mut tsv = if let Some(prefix) = output_prefix {
            let path = format!("{}{}", prefix, TSV_SUFFIX);
            let mut f = std::fs::OpenOptions::new();
            f.write(true).create(true);
            if !append_tsv {
                f.truncate(true);
            }
            let file = f.open(&path)?;
            let mut bw = BufWriter::with_capacity(8 * 1024 * 1024, file);
            if append_tsv {
                bw.seek(SeekFrom::End(0))?;
            }
            TsvSink::File(bw)
        } else {
            TsvSink::Stdout(BufWriter::with_capacity(1024 * 1024, std::io::stdout()))
        };

        // Optionally open work-log. Same rationale for skipping `O_APPEND`.
        let mut work: Option<BufWriter<std::fs::File>> = if with_work_log {
            let prefix = output_prefix.expect("with_work_log implies output_prefix");
            let path = format!("{}{}", prefix, WORKLOG_SUFFIX);
            let mut f = std::fs::OpenOptions::new();
            f.write(true).create(true).read(true);
            if !append_worklog {
                f.truncate(true);
            }
            let file = f.open(&path)?;
            let mut bw = BufWriter::with_capacity(2 * 1024 * 1024, file);
            if append_worklog {
                bw.seek(SeekFrom::End(0))?;
            }
            Some(bw)
        } else {
            None
        };

        let has_work_log = work.is_some();

        // Header logic: TSV header on fresh files only. Work-log header on
        // fresh files; on append we trust the existing one.
        if !skip_header {
            tsv.write_all(b"#id\tlength\tdepth\tpositions\n")?;
        }
        if let Some(w) = work.as_mut() {
            if !append_worklog {
                w.write_all(&encode_worklog_header())?;
            }
        }

        let (tx, rx) = sync_channel::<WriterMsg>(WRITER_CHANNEL_CAPACITY);
        let handle = std::thread::Builder::new()
            .name("depth-writer".to_string())
            .spawn(move || -> io::Result<()> {
                let mut tsv = tsv;
                let mut work = work;
                while let Ok(msg) = rx.recv() {
                    match msg {
                        WriterMsg::TsvOnly(buf) => {
                            tsv.write_all(&buf)?;
                        }
                        WriterMsg::Chunk { tsv: tbuf, work: wbuf } => {
                            // TSV first, then work-log. Both are append-only;
                            // ordering between the two files doesn't matter for
                            // correctness, but doing them back-to-back keeps
                            // the per-chunk advance together when a barrier
                            // captures offsets.
                            if !tbuf.is_empty() {
                                tsv.write_all(&tbuf)?;
                            }
                            if let Some(w) = work.as_mut() {
                                if !wbuf.is_empty() {
                                    w.write_all(&wbuf)?;
                                }
                            }
                        }
                        WriterMsg::Barrier(reply) => {
                            let result: io::Result<(u64, u64)> = (|| {
                                let tsv_off = tsv.flush_and_sync()?;
                                let work_off = if let Some(w) = work.as_mut() {
                                    w.flush()?;
                                    w.get_ref().sync_data()?;
                                    w.stream_position()?
                                } else {
                                    0u64
                                };
                                Ok((tsv_off, work_off))
                            })();
                            // Best-effort send; if the receiver is gone the
                            // requester crashed — let the writer keep going
                            // until the channel closes for real.
                            let _ = reply.send(result);
                        }
                    }
                }
                tsv.flush()?;
                if let Some(w) = work.as_mut() {
                    w.flush()?;
                }
                Ok(())
            })?;

        Ok(DepthWriter {
            tx,
            handle,
            has_work_log,
        })
    }

    /// Send TSV bytes only. Backward-compatible with the pre-resume API.
    fn send(&self, buf: Vec<u8>) -> io::Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        self.tx.send(WriterMsg::TsvOnly(buf)).map_err(|e| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("depth writer channel closed: {}", e),
            )
        })
    }

    /// Atomic per-chunk write: TSV bytes + matching work-log record. Use this
    /// in CIGAR Phase 1/2 hot loops when checkpoint/resume is enabled. When
    /// the writer was opened without a work-log this falls back to a TSV-only
    /// send and the `work` buffer is discarded (the caller has decided not to
    /// checkpoint).
    fn send_chunk_bundle(&self, tsv: Vec<u8>, work: Vec<u8>) -> io::Result<()> {
        if !self.has_work_log {
            return self.send(tsv);
        }
        if tsv.is_empty() && work.is_empty() {
            return Ok(());
        }
        self.tx
            .send(WriterMsg::Chunk { tsv, work })
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("depth writer channel closed: {}", e),
                )
            })
    }

    /// Issue a barrier and wait for the writer thread to flush + `sync_data`
    /// both backing files. Returns `(tsv_byte_offset, work_byte_offset)` post
    /// flush. Stdout-mode TSV always returns 0.
    fn barrier_and_offsets(&self) -> io::Result<(u64, u64)> {
        let (reply_tx, reply_rx) = sync_channel::<io::Result<(u64, u64)>>(1);
        self.tx.send(WriterMsg::Barrier(reply_tx)).map_err(|e| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("depth writer channel closed: {}", e),
            )
        })?;
        match reply_rx.recv() {
            Ok(res) => res,
            Err(_) => Err(io::Error::other(
                "depth writer thread dropped barrier reply channel",
            )),
        }
    }

    /// Drop the sender (signals the writer to flush + exit) and join the
    /// writer thread. Propagates any IO error from the thread.
    fn finish(self) -> io::Result<()> {
        let DepthWriter { tx, handle, .. } = self;
        drop(tx);
        match handle.join() {
            Ok(res) => res,
            Err(_) => Err(io::Error::other("depth writer thread panicked")),
        }
    }
}

/// Number of completed chunks between two automatic ckpt commits.
///
/// Phase 1/2 chunks are 5 MB of sequence each; CIGAR BFS through the
/// alignment network typically runs each chunk in the seconds–minutes range.
/// At 64 chunks per commit we're committing on the order of every few
/// minutes, which strikes the balance between (a) ckpt-write IO amortising
/// to noise and (b) bounded re-work after a crash.
const CHECKPOINT_CHUNK_INTERVAL: usize = 64;

/// Drives periodic ckpt commits without requiring a dedicated thread.
///
/// One worker out of the rayon par_iter wins the `try_lock` whenever the
/// shared chunk counter rolls past a multiple of `CHECKPOINT_CHUNK_INTERVAL`,
/// issues a writer barrier, captures the post-flush byte offsets and current
/// counter snapshot, and persists a fresh ckpt via `save_atomic`. Other
/// workers continue pumping chunk bundles into the writer channel.
///
/// When `state` is `None` (resume disabled or a stats-mode run) every method
/// is a no-op so call sites can stay branchless.
struct CheckpointController<'a> {
    state: Option<CheckpointState<'a>>,
}

struct CheckpointState<'a> {
    prefix: &'a str,
    invalidation_hash: u64,
    interval: usize,
    counter: AtomicUsize,
    commit_lock: Mutex<()>,
    row_counter: &'a AtomicUsize,
    intervals_counter: &'a AtomicUsize,
    writer: &'a DepthWriter,
}

impl<'a> CheckpointController<'a> {
    fn disabled() -> Self {
        Self { state: None }
    }

    fn new(
        prefix: &'a str,
        invalidation_hash: u64,
        writer: &'a DepthWriter,
        row_counter: &'a AtomicUsize,
        intervals_counter: &'a AtomicUsize,
    ) -> Self {
        // Tests can set `IMPG_CHECKPOINT_INTERVAL_OVERRIDE` to a small number
        // (e.g. 1) so a tiny synthetic dataset still exercises the periodic
        // commit path. Production runs use the constant default.
        let interval = std::env::var("IMPG_CHECKPOINT_INTERVAL_OVERRIDE")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|n| *n > 0)
            .unwrap_or(CHECKPOINT_CHUNK_INTERVAL);
        Self {
            state: Some(CheckpointState {
                prefix,
                invalidation_hash,
                interval,
                counter: AtomicUsize::new(0),
                commit_lock: Mutex::new(()),
                row_counter,
                intervals_counter,
                writer,
            }),
        }
    }

    fn enabled(&self) -> bool {
        self.state.is_some()
    }

    /// Called by every worker after each completed chunk. Increments the
    /// shared counter and triggers an opportunistic commit when the
    /// configured interval is reached. Only one worker per crossing actually
    /// runs the commit (the rest fall through on `try_lock`).
    fn note_chunk_done(&self) -> io::Result<()> {
        let Some(s) = self.state.as_ref() else {
            return Ok(());
        };
        let n = s.counter.fetch_add(1, Ordering::Relaxed) + 1;
        if n % s.interval == 0 {
            if let Some(_g) = s.commit_lock.try_lock() {
                self.commit_inner(s)?;
            }
        }
        Ok(())
    }

    /// Force a commit (waits for the lock if another worker is currently
    /// committing). Used at phase boundaries and at the end of the run.
    fn force_commit(&self) -> io::Result<()> {
        let Some(s) = self.state.as_ref() else {
            return Ok(());
        };
        let _g = s.commit_lock.lock();
        self.commit_inner(s)
    }

    fn commit_inner(&self, s: &CheckpointState<'_>) -> io::Result<()> {
        let (tsv_off, work_off) = s.writer.barrier_and_offsets()?;
        let ck = DepthCheckpoint {
            schema_version: crate::commands::depth_checkpoint::CKPT_SCHEMA_VERSION,
            invalidation_hash: s.invalidation_hash,
            tsv_byte_offset: tsv_off,
            work_byte_offset: work_off,
            row_counter: s.row_counter.load(Ordering::Relaxed) as u64,
            intervals_counter: s.intervals_counter.load(Ordering::Relaxed) as u64,
        };
        ck.save_atomic(s.prefix)?;
        debug!(
            "ckpt committed: tsv={} work={} row={} intervals={}",
            ck.tsv_byte_offset,
            ck.work_byte_offset,
            ck.row_counter,
            ck.intervals_counter,
        );
        // Test hook: when set, abort the process immediately after a successful
        // commit. The next run with `--resume` should pick up exactly here.
        // Production runs never set this env var.
        if std::env::var_os("IMPG_TEST_EXIT_AFTER_COMMIT").is_some() {
            // The writer thread, channel buffers and any pending Phase 1/2
            // messages are intentionally abandoned — that's the whole point
            // of the simulated crash.
            eprintln!(
                "[impg test hook] IMPG_TEST_EXIT_AFTER_COMMIT set; aborting after \
                 ckpt commit (tsv={}, work={})",
                ck.tsv_byte_offset, ck.work_byte_offset,
            );
            std::process::exit(99);
        }
        Ok(())
    }
}

/// Reference-free global depth computation.
///
/// Algorithm:
/// 1. Sort all sequences by length descending (--ref sample gets priority)
/// 2. For each sequence, find unprocessed regions
/// 3. Query alignments for each unprocessed region (direct or transitive)
/// 4. Compute depth via sweep-line, mark all discovered regions as processed
/// 5. Every base in the pangenome is assigned to exactly one output row
///
/// Query modes:
/// - Default (non-transitive): direct 1-hop queries. Fast, bounded memory.
///   impg's bidirectional index captures both alignment directions in one lookup.
/// - Transitive (-x): BFS with 5MB region chunking to bound memory.
///
/// Key properties:
/// - Every base processed exactly once (no double-counting)
/// - Hub-first ordering: high-connectivity sequences anchor first (auto-detected or via --ref)
/// - --ref: ref-anchored mode, guarantees ref sample's coordinate system for covered regions
/// - --ref-only: ref-only mode, output filtered to ref sample's anchored regions only
#[allow(clippy::too_many_arguments)]
pub fn compute_depth_global(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    separator: &str,
    output_prefix: Option<&str>,
    fai_list: Option<&str>,
    window_size: Option<i64>,
    min_interval_len: i64,
    ref_sample: Option<&str>,
    ref_only: bool,
    min_seq_length: i64,
    stats_mode: bool,
    stats_combined: bool,
    resume: bool,
    alignment_files: &[String],
) -> io::Result<()> {
    // ------------------------------------------------------------------
    // Checkpoint / resume bootstrap.
    //
    // `resume` is the sole user-facing toggle:
    //   - false: legacy behavior. Writer is TSV-only; no ckpt or work-log.
    //   - true:  CIGAR-precise BFS path only (validated upstream). On a fresh
    //            run we still write `<prefix>.depth.work.bin` + `<prefix>
    //            .depth.ckpt` so the *next* invocation can resume. On a
    //            subsequent run the work-log is replayed to rebuild
    //            `tracker` / `global_used` and any chunk_id present in the
    //            replayed records is short-circuited.
    //
    // The caller (main.rs) has already rejected combinations we don't
    // support (--stats, raw-BFS, stdout output, region query mode). The
    // assertions below are a defence-in-depth check.
    // ------------------------------------------------------------------
    if resume {
        if !config.use_cigar_bfs {
            return Err(io::Error::other(
                "internal: --resume reached compute_depth_global without --use-BFS",
            ));
        }
        if !(config.transitive || config.transitive_dfs) {
            return Err(io::Error::other(
                "internal: --resume requires transitive depth (-x or --transitive-dfs)",
            ));
        }
        if stats_mode || stats_combined {
            return Err(io::Error::other(
                "internal: --resume is not yet supported with --stats",
            ));
        }
        if output_prefix.is_none() {
            return Err(io::Error::other(
                "internal: --resume requires an --output-prefix",
            ));
        }
    }

    // Compute invalidation hash early so we can compare against any existing
    // ckpt before consuming wall time on a stale state.
    let invalidation_hash = compute_invalidation_hash(
        alignment_files,
        &HashInputs {
            ref_sample,
            ref_only,
            min_seq_length,
            min_interval_len,
            window_size,
            merge_adjacent: config.merge_adjacent,
            stats_mode,
            stats_combined,
            separator,
            fai_list,
            use_cigar_bfs: config.use_cigar_bfs,
            transitive: config.transitive,
            transitive_dfs: config.transitive_dfs,
            max_depth: config.max_depth,
            min_transitive_len: config.min_transitive_len,
            min_distance_between_ranges: config.min_distance_between_ranges,
        },
    );

    let resume_state: Option<ResumeState> = if resume {
        // Safe to unwrap: validated above.
        let prefix = output_prefix.unwrap();
        match DepthCheckpoint::try_load(prefix)? {
            Some(ck) => {
                if ck.invalidation_hash != invalidation_hash {
                    return Err(io::Error::other(format!(
                        "checkpoint at {prefix}.depth.ckpt was produced with a different \
                         configuration or alignment fingerprint (hash {} != {}); refusing to \
                         resume. Delete the .ckpt + .work.bin files to start fresh.",
                        ck.invalidation_hash, invalidation_hash,
                    )));
                }
                info!(
                    "Resuming depth: ckpt @ tsv={} bytes, work={} bytes, row_counter={}",
                    ck.tsv_byte_offset, ck.work_byte_offset, ck.row_counter
                );

                // Truncate live files to ckpt offsets. After this point any
                // post-commit dirty bytes are physically gone; the writer
                // thread (when we start it) will append cleanly.
                let tsv_path: std::path::PathBuf =
                    format!("{}{}", prefix, TSV_SUFFIX).into();
                let work_path: std::path::PathBuf =
                    format!("{}{}", prefix, WORKLOG_SUFFIX).into();
                truncate_to(&tsv_path, ck.tsv_byte_offset)?;
                truncate_to(&work_path, ck.work_byte_offset)?;
                Some(ResumeState {
                    ckpt: ck,
                    completed_chunks: std::collections::HashSet::new(),
                })
            }
            None => None,
        }
    } else {
        None
    };

    // ------------------------------------------------------------------
    let is_transitive = config.transitive || config.transitive_dfs;
    let user_window_size = window_size;
    let window_size = window_size.unwrap_or(DEFAULT_WINDOW_SIZE);
    if is_transitive {
        info!(
            "Computing depth (global mode, transitive BFS, {}MB query chunks)",
            TRANSITIVE_CHUNK_SIZE / 1_000_000
        );
    } else {
        info!("Computing depth (global mode, pre-scan + sweep-line)");
    }
    if min_interval_len > 0 {
        info!("Min interval length: {} bp (shorter intervals absorbed into neighbors)", min_interval_len);
    }

    // Build compact data structures
    let compact_lengths = CompactSequenceLengths::from_impg(impg, separator);
    let sample_index = compact_lengths.sample_index();
    let num_samples = sample_index.len();
    let num_sequences = impg.seq_index().len();

    debug!("Found {} samples, {} sequences", num_samples, num_sequences);

    // Build sequence inclusion filter (for --min-seq-length)
    let seq_included: Vec<bool> = (0..num_sequences as u32)
        .map(|id| compact_lengths.get_length(id) >= min_seq_length)
        .collect();
    if min_seq_length > 0 {
        let included = seq_included.iter().filter(|&&v| v).count();
        let excluded = num_sequences - included;
        info!(
            "Sequence length filter (>= {} bp): {} included, {} excluded",
            min_seq_length, included, excluded
        );
    }

    // Resolve ref sample ID
    let ref_sample_id: Option<u16> = if let Some(ref_name) = ref_sample {
        let id = sample_index.get_id(ref_name).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "Reference sample '{}' not found in alignment index",
                    ref_name
                ),
            )
        })?;
        if ref_only {
            info!(
                "Ref-only mode: '{}' (output filtered to this sample's anchored regions)",
                ref_name
            );
        } else {
            info!(
                "Ref-anchored mode: '{}' (Phase 1 anchor, coordinate system priority)",
                ref_name
            );
        }
        Some(id)
    } else {
        if ref_only {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--ref-only requires --ref to be specified",
            ));
        }
        None
    };

    // Load FAI sequences if provided
    let fai_seq_lengths: Option<SequenceLengths> = if let Some(fai_path) = fai_list {
        let fai = SequenceLengths::from_fai_list(fai_path, separator)?;
        debug!(
            "FAI: {} samples, {} sequences",
            fai.sample_to_seqs.len(),
            fai.lengths.len()
        );
        Some(fai)
    } else {
        None
    };

    // Pre-scan: compute alignment degrees for all included sequences.
    // Used for automatic hub detection (when --ref is not specified) and
    // for sorting sequences by connectivity (hubs first) in all modes.
    //
    // For MultiImpg, `compute_sample_degrees` internally iterates files in parallel
    // and uses transient sub-index loads (NOT the shared `sub_indices` cache), so
    // retained memory is bounded to `num_threads` sub-indices regardless of the
    // total number of per-file indices. The old tree-cache toggle is no longer
    // needed for the pre-scan — leave the main-phase cache setting untouched.
    info!("Pre-scanning alignment degrees...");
    let degrees =
        compute_alignment_degrees(impg, &compact_lengths, &seq_included, min_seq_length);
    // Pre-scan uses load_sub_index_uncached (does NOT populate
    // transient_header_cache), so the cache should already be empty here. Clear
    // both caches defensively in case a future code change reintroduces caching
    // upstream — keeping mmap pressure off vm.max_map_count is critical when
    // running with hundreds of thousands of per-file indices.
    impg.clear_sub_index_cache();
    impg.clear_transient_header_cache();
    let max_degree = degrees.iter().copied().max().unwrap_or(0);
    let included_count = seq_included.iter().filter(|&&v| v).count();
    info!(
        "Degree scan complete: {} sequences, max degree = {}",
        included_count, max_degree
    );

    // Build sequence processing order (degree descending, length descending, ref sample first)
    let sequence_order = build_sequence_order(
        impg,
        &compact_lengths,
        ref_sample_id,
        &seq_included,
        &degrees,
    );
    debug!(
        "Sequence processing order: {} sequences",
        sequence_order.len()
    );

    // Prepare output: TSV writer for normal mode, stats accumulators for --stats add-on.
    //
    // Uses a dedicated writer thread + bounded MPSC channel rather than the
    // previous `Mutex<BufWriter>`; at high thread counts the global mutex was
    // serialising every TSV chunk write, even though each worker had already
    // formatted its 4 MB buffer locally. The writer thread also writes the TSV
    // header before draining the channel.
    let writer: Option<DepthWriter> = if !stats_mode {
        // Resume mode: open both files in append mode and skip the TSV
        // header (already present from the previous run, possibly truncated
        // mid-row by the resume entry path — but the truncation is to the
        // last committed offset which is past the header by definition).
        // Non-resume but `--resume` set: open the work-log fresh and write
        // both headers. Plain run: legacy behaviour (TSV-only, fresh).
        let opts = WriterOpts {
            output_prefix,
            append_tsv: resume_state.is_some(),
            append_worklog: resume_state.is_some(),
            with_work_log: resume,
            skip_header: resume_state.is_some(),
        };
        Some(DepthWriter::open(opts)?)
    } else {
        None
    };

    // Stats accumulators (only allocated in --stats mode)
    let stats_accumulator: Option<Mutex<DepthStats>> = if stats_mode && !stats_combined {
        Some(Mutex::new(DepthStats::new()))
    } else {
        None
    };
    let stats_combined_acc: Option<Mutex<DepthStatsWithSamples>> = if stats_mode && stats_combined {
        Some(Mutex::new(DepthStatsWithSamples::new()))
    } else {
        None
    };

    let total_sequences = sequence_order.len();

    // Enable tree caching: trees loaded by get_or_load_tree are cached in sub-indices.
    // With jemalloc, this no longer causes RSS fragmentation. Trees are reused across
    // BFS nodes within the same sequence (avoiding repeated disk I/O), then freed
    // when clear_sub_index_cache drops the sub-index Arc (and all its cached trees).
    impg.set_tree_cache_enabled(true);

    let tracker = ConcurrentProcessedTracker::new(num_sequences);
    let global_used = std::sync::Arc::new(ConcurrentProcessedTracker::new(num_sequences));

    // Replay the work-log into the trackers if we're resuming. We do this
    // *after* num_sequences is known (so the IntervalSet vectors are sized)
    // but *before* Phase 1 starts — both `tracker` and `global_used` need
    // the union of every previously-committed `discovered_regions`.
    //
    // Replay also populates `completed_chunks`, the per-chunk_id skip set
    // consulted by Phase 1 / Phase 2 workers.
    let mut completed_chunks: std::collections::HashSet<u64> =
        std::collections::HashSet::new();
    if let Some(rs) = resume_state.as_ref() {
        let prefix = output_prefix.unwrap();
        let work_path: std::path::PathBuf = format!("{}{}", prefix, WORKLOG_SUFFIX).into();
        let mut replayed_records: u64 = 0;
        let mut replayed_regions: u64 = 0;
        replay_work_log(&work_path, rs.ckpt.work_byte_offset, |rec| {
            // Idempotent: IntervalSet::add is a no-op for already-covered
            // ranges, so re-applying a record yields the same end state. We
            // still de-dup chunk_ids in case the work-log somehow contains
            // two entries for the same id (it shouldn't, but better to be
            // strict on the skip set).
            tracker.mark_processed_batch(&rec.regions);
            global_used.mark_processed_batch(&rec.regions);
            replayed_records += 1;
            replayed_regions += rec.regions.len() as u64;
            completed_chunks.insert(rec.chunk_id);
            Ok(())
        })?;
        info!(
            "Resume replay: {} chunks, {} regions reconstructed from {}",
            replayed_records,
            replayed_regions,
            work_path.display()
        );
        // row_counter / intervals_counter are seeded from rs.ckpt at their
        // declaration further down (we can't borrow them yet here).
    }

    // =========================================================================
    // Two-phase parallel processing with hub-first guarantee.
    //
    // Phase 1 ensures high-connectivity (hub) sequences are processed first,
    // so their coordinate systems become the anchors for all homologous regions.
    //
    // Hub selection:
    //   --ref specified: ref sample's sequences are the hubs.
    //   --ref not specified: auto-detect hubs via alignment degree pre-scan.
    //     Sequences with degree >= max_degree/2 are classified as hubs.
    //     For star topology (hub degree ~800, leaf degree ~1), this cleanly
    //     separates the central sequence from leaves.
    //
    // Phase 1 parallelism:
    //   Transitive (-x): chunk-level (5MB chunks) for full thread utilization.
    //   Non-transitive: sequence-level (fast per-sequence, acceptable overhead).
    //
    // Phase 2: remaining sequences, most regions already claimed by Phase 1.
    // =========================================================================

    // Partition sequences into Phase 1 (hubs) and Phase 2 (rest)
    let (phase1_seqs, phase2_seqs): (Vec<u32>, Vec<u32>) = if let Some(ref_id) = ref_sample_id {
        // --ref specified: ref sequences go to Phase 1
        sequence_order
            .iter()
            .partition(|&&seq_id| compact_lengths.get_sample_id(seq_id) == ref_id)
    } else {
        // No --ref: auto-detect hubs by alignment degree
        let hub_threshold = (max_degree + 1) / 2; // ceiling of max_degree/2
        if hub_threshold > 1 {
            let (p1, p2): (Vec<u32>, Vec<u32>) = sequence_order.iter().partition(|&&seq_id| {
                degrees.get(seq_id as usize).copied().unwrap_or(0) >= hub_threshold
            });
            if !p1.is_empty() {
                info!(
                    "Auto-detected {} hub sequences (degree >= {}) for Phase 1",
                    p1.len(),
                    hub_threshold
                );
            }
            (p1, p2)
        } else {
            // No clear hubs (max_degree <= 1), skip Phase 1
            (Vec::new(), sequence_order.clone())
        }
    };

    let processed_count = AtomicUsize::new(0);
    let row_counter = AtomicUsize::new(
        resume_state
            .as_ref()
            .map(|rs| rs.ckpt.row_counter as usize)
            .unwrap_or(0),
    );
    let intervals_counter = AtomicUsize::new(
        resume_state
            .as_ref()
            .map(|rs| rs.ckpt.intervals_counter as usize)
            .unwrap_or(0),
    );

    // Build the checkpoint controller. `resume == true && writer.is_some()`
    // is the only configuration that produces a live ckpt-saving path; in
    // every other case the controller is a no-op and worker code stays
    // branch-free.
    let checkpoint_ctrl: CheckpointController = match (resume, writer.as_ref()) {
        (true, Some(w)) => CheckpointController::new(
            output_prefix.expect("--resume requires output_prefix"),
            invalidation_hash,
            w,
            &row_counter,
            &intervals_counter,
        ),
        _ => CheckpointController::disabled(),
    };
    let work_log_active = checkpoint_ctrl.enabled();

    // Helper closure: format and write results for a batch of AnchorRegionResults.
    // In stats mode: accumulates depth distribution + intervals into stats accumulators.
    // In normal mode: formats TSV rows into buf for the caller to write.
    let write_results = |region_results: Vec<AnchorRegionResult>,
                         buf: &mut Vec<u8>|
     -> io::Result<()> {
        // Thread-local accumulators (merged into global at end of batch)
        let mut local_stats: Option<DepthStats> = if stats_accumulator.is_some() {
            Some(DepthStats::new())
        } else {
            None
        };
        let mut local_combined: Option<DepthStatsWithSamples> = if stats_combined_acc.is_some() {
            Some(DepthStatsWithSamples::new())
        } else {
            None
        };

        for result in region_results {
            let seq_name = impg
                .seq_index()
                .get_name(result.anchor_seq_id)
                .unwrap_or("?");
            let anchor_sample_id = result.anchor_sample_id;

            let should_output = if ref_only {
                ref_sample_id == Some(anchor_sample_id)
            } else {
                true
            };

            if !should_output {
                continue;
            }

            let windowed_intervals = if user_window_size.is_some() {
                split_intervals_by_window(result.intervals, window_size)
            } else {
                result.intervals
            };

            let final_intervals = merge_short_intervals(windowed_intervals, min_interval_len);

            for interval in &final_intervals {
                if stats_mode {
                    // Stats add-on: accumulate depth distribution and intervals
                    let depth = interval.depth();
                    if let Some(ref mut ls) = local_stats {
                        ls.add_interval(
                            seq_name,
                            interval.start,
                            interval.end,
                            depth,
                            interval.pangenome_bases,
                        );
                    }
                    if let Some(ref mut lc) = local_combined {
                        let sample_names: Vec<String> = interval
                            .samples
                            .iter()
                            .filter_map(|&(sid, _, _, _)| {
                                sample_index.get_name(sid).map(|s| s.to_string())
                            })
                            .collect();
                        lc.add_interval(
                            seq_name,
                            interval.start,
                            interval.end,
                            depth,
                            sample_names,
                            interval.pangenome_bases,
                        );
                    }
                } else {
                    // Normal mode: write TSV row to buf
                    let rid = row_counter.fetch_add(1, Ordering::Relaxed) + 1;

                    write!(
                        buf,
                        "{}\t{}\t{}",
                        rid,
                        interval.end - interval.start,
                        interval.depth()
                    )?;

                    // Anchor position first
                    write!(buf, "\t{}:{}-{}", seq_name, interval.start, interval.end)?;
                    // Then other samples (non-anchor) sorted by sample_id
                    for &(sid, query_id, q_start, q_end) in &interval.samples {
                        if sid != anchor_sample_id {
                            if let Some(name) = impg.seq_index().get_name(query_id) {
                                write!(buf, ";{}:{}-{}", name, q_start, q_end)?;
                            }
                        }
                    }
                    writeln!(buf)?;

                    intervals_counter.fetch_add(1, Ordering::Relaxed);
                }
            }
        }

        // Merge thread-local stats into global accumulators (one lock per batch)
        if let Some(ls) = local_stats {
            if let Some(ref acc) = stats_accumulator {
                acc.lock().merge(&ls);
            }
        }
        if let Some(lc) = local_combined {
            if let Some(ref acc) = stats_combined_acc {
                acc.lock().merge(lc);
            }
        }

        Ok(())
    };

    // =========================================================================
    // Phase 1: Process hub sequences first (guaranteed to complete before Phase 2)
    // =========================================================================
    if !phase1_seqs.is_empty() {
        if is_transitive {
            // Transitive mode (raw BFS path): two-phase approach to bound memory.
            //
            // Old approach: par_iter over all chunks → each thread independently loads
            // sub-index files → T threads × sub_index_size peak memory (120 GB+).
            //
            // New approach (batch BFS):
            //   Phase A) batch_depth_bfs drives all chunk frontiers together.
            //            Each alignment file is loaded once per BFS round, shared
            //            across every chunk that needs it, then immediately freed.
            //            Peak memory = one sub-index at a time (sequential files).
            //            For hub workloads where N chunks all query the same sample
            //            files, this also reduces total I/O N-fold.
            //   Phase B) Parallel sweep-line: no sub-index loading, CPU-only.
            //            128 threads process the pre-computed hits concurrently.
            //
            // CIGAR BFS (--use-BFS): kept on the old par_iter path — CIGAR BFS
            // needs CIGAR-precise projection and doesn't benefit from the batch
            // optimisation; it is also rarely used in practice.

            // Split hub sequences into 5MB chunks
            let phase1_chunks: Vec<(u32, u16, i64, i64)> = phase1_seqs
                .iter()
                .flat_map(|&seq_id| {
                    let seq_len = compact_lengths.get_length(seq_id);
                    let sample_id = compact_lengths.get_sample_id(seq_id);
                    let mut chunks = Vec::new();
                    let mut pos = 0i64;
                    while pos < seq_len {
                        let chunk_end = (pos + TRANSITIVE_CHUNK_SIZE).min(seq_len);
                        chunks.push((seq_id, sample_id, pos, chunk_end));
                        pos = chunk_end;
                    }
                    chunks
                })
                .collect();

            let num_chunks = phase1_chunks.len();
            info!(
                "Phase 1: {} hub sequences -> {} chunks ({}MB each), batch BFS",
                phase1_seqs.len(),
                num_chunks,
                TRANSITIVE_CHUNK_SIZE / 1_000_000
            );

            let pb_phase1 = ProgressBar::new(num_chunks as u64);
            pb_phase1.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} chunks ({eta}) | Phase 1: hub sequences")
                    .unwrap()
                    .progress_chars("#>-")
            );

            let phase1_count = AtomicUsize::new(0);

            if config.use_cigar_bfs {
                // CIGAR path: keep existing per-chunk parallel approach.
                //
                // When checkpointing is active (`work_log_active`):
                //   - Each chunk's chunk_id is the deterministic
                //     `encode_chunk_id_phase1(idx)` of its position in the
                //     sorted `phase1_chunks` Vec.
                //   - Chunks already present in `completed_chunks` (recovered
                //     from a previous run's work-log) are short-circuited
                //     before any sub-index loading happens.
                //   - The per-chunk TSV bytes and a serialized work-log
                //     record are bundled atomically (`send_chunk_bundle`)
                //     so the writer thread advances both files in lockstep
                //     and the next ckpt commit captures consistent offsets.
                phase1_chunks
                    .par_iter()
                    .enumerate()
                    .try_for_each(
                        |(idx, &(seq_id, sample_id, chunk_start, chunk_end))| -> io::Result<()> {
                            let chunk_id = encode_chunk_id_phase1(idx);
                            if work_log_active && completed_chunks.contains(&chunk_id) {
                                let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                                pb_phase1.set_position(count as u64);
                                return Ok(());
                            }
                            let result = process_anchor_region(
                                impg,
                                config,
                                &compact_lengths,
                                num_samples,
                                seq_id,
                                sample_id,
                                chunk_start,
                                chunk_end,
                                &seq_included,
                                min_seq_length,
                                &global_used,
                            );
                            tracker.mark_processed_batch(&result.discovered_regions);
                            let mut tsv_buf: Vec<u8> = Vec::new();
                            // Capture discovered_regions before write_results
                            // moves the result out — needed for the work-log
                            // record so the next resume can rebuild the
                            // tracker state without recomputing the BFS.
                            let work_buf = if work_log_active {
                                encode_work_record(chunk_id, &result.discovered_regions)
                            } else {
                                Vec::new()
                            };
                            write_results(vec![result], &mut tsv_buf)?;
                            if let Some(ref w) = writer {
                                if work_log_active {
                                    w.send_chunk_bundle(
                                        std::mem::take(&mut tsv_buf),
                                        work_buf,
                                    )?;
                                } else if !tsv_buf.is_empty() {
                                    w.send(std::mem::take(&mut tsv_buf))?;
                                }
                            }
                            checkpoint_ctrl.note_chunk_done()?;
                            let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                            pb_phase1.set_position(count as u64);
                            Ok(())
                        },
                    )?;
            } else {
                // Raw BFS path (default): batch BFS — sequential file loading,
                // then parallel sweep-line.

                // Phase A: batch BFS (single-threaded file loading, all chunks coordinated)
                let chunk_coords: Vec<(u32, i64, i64)> = phase1_chunks
                    .iter()
                    .map(|&(seq_id, _, start, end)| (seq_id, start, end))
                    .collect();

                info!("Phase 1: running batch BFS across {} chunks...", num_chunks);
                let all_hits = batch_depth_bfs(
                    impg,
                    &chunk_coords,
                    config.max_depth,
                    config.min_transitive_len,
                    config.min_distance_between_ranges,
                );
                info!("Phase 1: batch BFS complete, starting parallel sweep-line...");

                // Phase B: parallel sweep-line (CPU-only, no sub-index loading)
                phase1_chunks
                    .par_iter()
                    .zip(all_hits.into_par_iter())
                    .try_for_each(
                        |(&(seq_id, sample_id, chunk_start, chunk_end), hits)| -> io::Result<()> {
                            let result = process_anchor_region_transitive_raw_with_hits(
                                hits,
                                &compact_lengths,
                                num_samples,
                                seq_id,
                                sample_id,
                                chunk_start,
                                chunk_end,
                                &seq_included,
                                min_seq_length,
                                &global_used,
                            );
                            tracker.mark_processed_batch(&result.discovered_regions);
                            let mut buf: Vec<u8> = Vec::new();
                            write_results(vec![result], &mut buf)?;
                            if !buf.is_empty() {
                                if let Some(ref w) = writer {
                                    w.send(std::mem::take(&mut buf))?;
                                }
                            }
                            let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                            pb_phase1.set_position(count as u64);
                            Ok(())
                        },
                    )?;
            }

            impg.clear_sub_index_cache();
            // Release the per-file `Arc<Impg>` headers that Phase 1's
            // batch BFS / chunked queries pinned in the transient cache.
            // At `--index-mode per-file` with hundreds of thousands of
            // alignment files this is what keeps Phase 2 from inheriting
            // a cache that has already crowded `vm.max_map_count`.
            impg.clear_transient_header_cache();
            pb_phase1.finish_and_clear();
        } else {
            // Non-transitive mode: chunk-level parallelism.
            //
            // Previous design ran one task per hub sequence and called
            // `query_raw_intervals_transient(seq_id)` to materialize *all* raw
            // alignments for the chromosome in a single Vec. For CHM13-scale
            // (~250 Mb chr1, 580 samples, 580 per-file indices) that single
            // raw_alns Vec can reach 0.3–1 GB, the sweep-line `events` Vec is 2×
            // that, and the streaming emitter's `seq_intervals` (when min_interval_len
            // = 0 it was never streamed) added 10s of GB more. With par_iter at
            // worker-count concurrency, the process OOMed before Phase 1 finished.
            //
            // New design: split each chromosome into 5 MB chunks (matching the
            // transitive path's `TRANSITIVE_CHUNK_SIZE`) and parallelize chunks
            // instead of sequences. Each chunk uses `query_raw_overlapping_transient`
            // for a bounded range query, so per-task peak is O(chunk_alignments)
            // and total memory peak across workers is `threads × O(chunk)`.
            //
            // Chunk boundaries are independent of alignment boundaries — claim_unprocessed
            // on the query side is atomic across chunks, so no double-counting.
            // The per-chunk emitter resets `pending`, so when min_interval_len > 0
            // a leading-short at a chunk seam is not merged across the seam; this
            // is acceptable at 5 MB granularity for transient-dropout absorption.
            let phase1_chunks: Vec<(u32, i64, i64)> = phase1_seqs
                .iter()
                .flat_map(|&seq_id| {
                    let seq_len = compact_lengths.get_length(seq_id);
                    let mut chunks = Vec::new();
                    if seq_len > 0 {
                        let mut start = 0i64;
                        while start < seq_len {
                            let end = (start + TRANSITIVE_CHUNK_SIZE).min(seq_len);
                            chunks.push((seq_id, start, end));
                            start = end;
                        }
                    }
                    chunks
                })
                .collect();

            info!(
                "Phase 1: {} hub sequences split into {} {}-MB chunks, chunk-level parallelism",
                phase1_seqs.len(),
                phase1_chunks.len(),
                TRANSITIVE_CHUNK_SIZE / 1_000_000,
            );

            let pb_phase1 = ProgressBar::new(phase1_chunks.len() as u64);
            pb_phase1.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} chunks ({eta}) | Phase 1: hub sequences")
                    .unwrap()
                    .progress_chars("#>-")
            );

            let phase1_count = AtomicUsize::new(0);

            phase1_chunks
                .par_iter()
                .try_for_each(|&(seq_id, chunk_start, chunk_end)| -> io::Result<()> {
                    let sample_id = compact_lengths.get_sample_id(seq_id);

                    // Range query bounded by the chunk: peak per-task memory is
                    // O(chunk_alignments) instead of O(chromosome_alignments).
                    let mut raw_alns = impg.query_raw_overlapping_transient(
                        seq_id,
                        chunk_start,
                        chunk_end,
                    );
                    if min_seq_length > 0 {
                        raw_alns.retain(|aln| {
                            seq_included
                                .get(aln.query_id as usize)
                                .copied()
                                .unwrap_or(false)
                        });
                    }
                    raw_alns.sort_unstable_by_key(|aln| aln.target_start);

                    let seq_name = impg.seq_index().get_name(seq_id).unwrap_or("?");
                    let should_output = if ref_only {
                        ref_sample_id == Some(sample_id)
                    } else {
                        true
                    };

                    let mut emitter = StreamingDepthEmitter {
                        buf: Vec::new(),
                        writer: &writer,
                        seq_name,
                        anchor_sample_id: sample_id,
                        num_samples,
                        seq_index: impg.seq_index(),
                        row_counter: &row_counter,
                        intervals_counter: &intervals_counter,
                        should_output,
                        stats_mode,
                        local_stats: if stats_accumulator.is_some() {
                            Some(DepthStats::new())
                        } else {
                            None
                        },
                        local_combined: if stats_combined_acc.is_some() {
                            Some(DepthStatsWithSamples::new())
                        } else {
                            None
                        },
                        sample_index: &sample_index,
                        min_interval_len,
                        seq_intervals: Vec::new(),
                        pending: None,
                        window_size: if user_window_size.is_some() {
                            Some(window_size)
                        } else {
                            None
                        },
                    };

                    let discovered = process_anchor_region_raw_streaming(
                        &raw_alns,
                        &compact_lengths,
                        num_samples,
                        seq_id,
                        sample_id,
                        chunk_start,
                        chunk_end,
                        &global_used,
                        &mut |interval| emitter.emit(interval),
                    );
                    emitter.flush()?;
                    tracker.mark_processed_batch(&discovered);
                    emitter.merge_stats_into(&stats_accumulator, &stats_combined_acc);

                    let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                    pb_phase1.set_position(count as u64);

                    Ok(())
                })?;

            // Same rationale as the transitive branch above: drop both the
            // sub-index BFS cache and the transient per-file header cache
            // before Phase 2 starts widening the working set further.
            impg.clear_sub_index_cache();
            impg.clear_transient_header_cache();
            pb_phase1.finish_and_clear();
        }

        // Mark all Phase 1 sequences as fully processed in both trackers.
        // global_used prevents Phase 2 from double-counting bases already counted as anchor.
        for &seq_id in &phase1_seqs {
            let seq_len = compact_lengths.get_length(seq_id);
            tracker.mark_processed(seq_id, 0, seq_len);
            global_used.mark_processed(seq_id, 0, seq_len);
        }

        // When checkpointing is active emit a synthetic work-log record at
        // the Phase 1 → Phase 2 boundary so a resume that loses Phase 1 in
        // mid-stream still ends up with the boundary mark applied. The
        // boundary record uses `CHUNK_ID_BOUNDARY` (= 0); workers consult
        // `completed_chunks` at chunk granularity, so we can safely re-emit
        // here only when the boundary has not yet been committed in a
        // previous run.
        if work_log_active && !completed_chunks.contains(&CHUNK_ID_BOUNDARY) {
            let regions: Vec<(u32, i64, i64)> = phase1_seqs
                .iter()
                .map(|&seq_id| (seq_id, 0i64, compact_lengths.get_length(seq_id)))
                .collect();
            let work_buf = encode_work_record(CHUNK_ID_BOUNDARY, &regions);
            if let Some(ref w) = writer {
                w.send_chunk_bundle(Vec::new(), work_buf)?;
            }
            // Force-commit immediately so the boundary is durably crossed
            // before Phase 2 starts producing chunks against it.
            checkpoint_ctrl.force_commit()?;
        }

        processed_count.store(phase1_seqs.len(), Ordering::Relaxed);
        info!(
            "Phase 1 complete: {} hub sequences processed",
            phase1_seqs.len()
        );
    }

    // =========================================================================
    // Phase 2: Process remaining sequences.
    //
    // Two distinct execution shapes:
    //
    //   Non-transitive (default): chunk-level parallelism, mirroring Phase 1.
    //     Per-task peak memory is bounded by O(chunk_alignments) instead of
    //     O(scaffold_alignments). This matters for pangenomes whose Phase 2
    //     pool contains long contigs/scaffolds (50–200 Mb is normal at
    //     CHM13/HPRC scale); the previous per-seq par_iter at ~80 threads
    //     could pin tens of GB of `raw_alns_cached` Vecs simultaneously.
    //
    //   Transitive (`-x`): seq-level parallelism is retained because the
    //     transitive path internally chunks at TRANSITIVE_CHUNK_SIZE and
    //     manages its own sub-index cache. Restructuring it would require
    //     also rewriting the BFS state plumbing — out of scope here.
    // =========================================================================
    let phase2_label = if phase1_seqs.is_empty() {
        ""
    } else {
        "Phase 2: "
    };

    if !is_transitive {
        // Precompute (seq_id, chunk_start, chunk_end) tuples by scanning the
        // current unprocessed gaps for every Phase 2 sequence and slicing each
        // gap at TRANSITIVE_CHUNK_SIZE. Phase 1 has already finished, so the
        // tracker state is stable when this runs.
        //
        // Two-pass build (count then fill) so the Vec is allocated exactly
        // once at its final capacity. With pangenome-scale workloads this
        // collect can return tens of millions of chunks (≥56M observed at
        // 581-sample VGP scale); the doubling growth pattern from a default
        // `flat_map().collect()` would create a transient ~2× peak right at
        // the Phase-1 → Phase-2 boundary, which on hosts already pressed up
        // against `vm.max_map_count` from the transient header cache was
        // enough to trip the allocator.
        let mut phase2_chunk_count_pre: usize = 0;
        for &seq_id in &phase2_seqs {
            let seq_len = compact_lengths.get_length(seq_id);
            if seq_len <= 0 {
                continue;
            }
            for (region_start, region_end) in tracker.get_unprocessed(seq_id, 0, seq_len) {
                let gap = region_end - region_start;
                if gap > 0 {
                    // ceil(gap / chunk_size)
                    phase2_chunk_count_pre +=
                        ((gap + TRANSITIVE_CHUNK_SIZE - 1) / TRANSITIVE_CHUNK_SIZE) as usize;
                }
            }
        }

        let mut phase2_chunks: Vec<(u32, i64, i64)> =
            Vec::with_capacity(phase2_chunk_count_pre);
        for &seq_id in &phase2_seqs {
            let seq_len = compact_lengths.get_length(seq_id);
            if seq_len <= 0 {
                continue;
            }
            for (region_start, region_end) in tracker.get_unprocessed(seq_id, 0, seq_len) {
                let mut start = region_start;
                while start < region_end {
                    let end = (start + TRANSITIVE_CHUNK_SIZE).min(region_end);
                    phase2_chunks.push((seq_id, start, end));
                    start = end;
                }
            }
        }
        debug_assert_eq!(phase2_chunks.len(), phase2_chunk_count_pre);

        info!(
            "Phase 2: {} sequences split into {} {}-MB chunks, chunk-level parallelism",
            phase2_seqs.len(),
            phase2_chunks.len(),
            TRANSITIVE_CHUNK_SIZE / 1_000_000,
        );

        let pb_depth = ProgressBar::new(phase2_chunks.len() as u64);
        pb_depth.set_style(
            ProgressStyle::default_bar()
                .template(&format!(
                    "{{spinner:.green}} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} chunks ({{eta}}) | {}remaining sequences",
                    phase2_label
                ))
                .unwrap()
                .progress_chars("#>-")
        );
        let phase2_chunk_count = AtomicUsize::new(0);

        // Two-phase batched architecture mirroring Phase 1 transitive
        // (`batch_depth_bfs`):
        //
        //   Phase A) `batch_query_raw_overlapping` groups the batch's queries by
        //            sub-index file and loads each file *once* per batch to serve
        //            every query that references it. Replaces the previous
        //            per-chunk `query_raw_overlapping_transient` which thrashed
        //            the transient header cache when consecutive rayon workers
        //            anchored on different leaf samples (Phase 1's pattern of
        //            "all chunks share the same hub seq's file set" does not
        //            hold for Phase 2 leaves).
        //
        //   Phase B) Parallel sweep + streaming emit. CPU-only, no sub-index
        //            loading. The hits Vec is already in original chunk order so
        //            zip() preserves the per-chunk closure semantics.
        //
        // Memory bound: the batch's hit Vec holds every alignment overlapping
        // every chunk in the batch. At PHASE2_BATCH_SIZE=65536 chunks × ~hundreds
        // of alignments per leaf 5MB chunk × ~64 bytes/RawAlignmentInterval the
        // peak is in the low-GB range. Larger batches further amortize header
        // loads but raise this peak; smaller batches lose amortization. 65k is a
        // pragmatic operating point for HPRC/VGP scale.
        const PHASE2_BATCH_SIZE: usize = 65_536;

        for batch in phase2_chunks.chunks(PHASE2_BATCH_SIZE) {
            // Phase A: file-grouped batch query. `batch_query_raw_overlapping`
            // is internally rayon-parallel across files; each file is loaded
            // transiently, answers all of its queries, then is freed.
            let queries: Vec<(u32, i64, i64)> =
                batch.iter().map(|&(s, cs, ce)| (s, cs, ce)).collect();
            let all_alns = impg.batch_query_raw_overlapping(&queries);

            // Phase B: parallel sweep + emit.
            batch
                .par_iter()
                .zip(all_alns.into_par_iter())
                .try_for_each(
                    |(&(seq_id, chunk_start, chunk_end), mut raw_alns)| -> io::Result<()> {
                        let sample_id = compact_lengths.get_sample_id(seq_id);

                        if min_seq_length > 0 {
                            raw_alns.retain(|aln| {
                                seq_included
                                    .get(aln.query_id as usize)
                                    .copied()
                                    .unwrap_or(false)
                            });
                        }
                        raw_alns.sort_unstable_by_key(|aln| aln.target_start);

                        let seq_name = impg.seq_index().get_name(seq_id).unwrap_or("?");
                        let should_output = if ref_only {
                            ref_sample_id == Some(sample_id)
                        } else {
                            true
                        };

                        let mut emitter = StreamingDepthEmitter {
                            buf: Vec::new(),
                            writer: &writer,
                            seq_name,
                            anchor_sample_id: sample_id,
                            num_samples,
                            seq_index: impg.seq_index(),
                            row_counter: &row_counter,
                            intervals_counter: &intervals_counter,
                            should_output,
                            stats_mode,
                            local_stats: if stats_accumulator.is_some() {
                                Some(DepthStats::new())
                            } else {
                                None
                            },
                            local_combined: if stats_combined_acc.is_some() {
                                Some(DepthStatsWithSamples::new())
                            } else {
                                None
                            },
                            sample_index: &sample_index,
                            min_interval_len,
                            seq_intervals: Vec::new(),
                            pending: None,
                            window_size: if user_window_size.is_some() {
                                Some(window_size)
                            } else {
                                None
                            },
                        };

                        let discovered = process_anchor_region_raw_streaming(
                            &raw_alns,
                            &compact_lengths,
                            num_samples,
                            seq_id,
                            sample_id,
                            chunk_start,
                            chunk_end,
                            &global_used,
                            &mut |interval| emitter.emit(interval),
                        );
                        emitter.flush()?;
                        tracker.mark_processed_batch(&discovered);
                        emitter.merge_stats_into(&stats_accumulator, &stats_combined_acc);

                        let count = phase2_chunk_count.fetch_add(1, Ordering::Relaxed) + 1;
                        pb_depth.set_position(count as u64);
                        Ok(())
                    },
                )?;
        }

        pb_depth.finish_and_clear();
        // Update overall processed counter so any downstream readers see the
        // full Phase 1 + Phase 2 sequence count.
        processed_count.store(phase1_seqs.len() + phase2_seqs.len(), Ordering::Relaxed);
    } else if config.use_cigar_bfs {
        // Phase 2 transitive — CIGAR BFS path (--use-BFS).
        //
        // Kept on the per-seq par_iter for the same reason Phase 1 keeps the
        // CIGAR BFS path on per-chunk par_iter (depth.rs above): CIGAR-precise
        // BFS does not slot into the raw `batch_depth_bfs` shape (its inner
        // hop projects through CIGAR ops, not raw alignment extents) and the
        // flag is rarely used in practice. Re-architecting it for batch loads
        // is out of scope for this optimisation.
        let pb_depth = ProgressBar::new(total_sequences as u64);
        pb_depth.set_style(
            ProgressStyle::default_bar()
                .template(&format!(
                    "{{spinner:.green}} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} seqs ({{eta}}) | {}remaining sequences",
                    phase2_label
                ))
                .unwrap()
                .progress_chars("#>-")
        );
        pb_depth.set_position(processed_count.load(Ordering::Relaxed) as u64);

        // Helper: send one chunk's TSV bytes + work-log record (when
        // checkpointing is active) and tick the checkpoint controller.
        let emit_chunk =
            |chunk_id: u64, result: AnchorRegionResult| -> io::Result<()> {
                let mut tsv_buf: Vec<u8> = Vec::new();
                let work_buf = if work_log_active {
                    encode_work_record(chunk_id, &result.discovered_regions)
                } else {
                    Vec::new()
                };
                write_results(vec![result], &mut tsv_buf)?;
                if let Some(ref w) = writer {
                    if work_log_active {
                        w.send_chunk_bundle(std::mem::take(&mut tsv_buf), work_buf)?;
                    } else if !tsv_buf.is_empty() {
                        w.send(std::mem::take(&mut tsv_buf))?;
                    }
                }
                checkpoint_ctrl.note_chunk_done()?;
                Ok(())
            };

        phase2_seqs
            .par_iter()
            .try_for_each(|&seq_id| -> io::Result<()> {
                let seq_len = compact_lengths.get_length(seq_id);
                if seq_len <= 0 {
                    let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
                    pb_depth.set_position(count as u64);
                    return Ok(());
                }

                // tracker.get_unprocessed is the natural skip mechanism for
                // resume: any chunk previously committed against this seq
                // shows up as a processed range, so its bytes are not in
                // `unprocessed` here and we won't re-do the work.
                let unprocessed = tracker.get_unprocessed(seq_id, 0, seq_len);
                if unprocessed.is_empty() {
                    let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
                    pb_depth.set_position(count as u64);
                    return Ok(());
                }

                let sample_id = compact_lengths.get_sample_id(seq_id);

                for (region_start, region_end) in unprocessed {
                    let gap_len = region_end - region_start;

                    if gap_len < config.min_transitive_len {
                        // Single-chunk raw fallback: short gaps below the
                        // BFS gate would otherwise be dropped entirely.
                        let chunk_id = encode_chunk_id_phase2(seq_id, region_start);
                        let mut raw_alns = impg.query_raw_overlapping_transient(
                            seq_id,
                            region_start,
                            region_end,
                        );
                        if min_seq_length > 0 {
                            raw_alns.retain(|aln| {
                                seq_included
                                    .get(aln.query_id as usize)
                                    .copied()
                                    .unwrap_or(false)
                            });
                        }
                        raw_alns.sort_unstable_by_key(|aln| aln.target_start);
                        let result = process_anchor_region_raw(
                            &raw_alns,
                            &compact_lengths,
                            num_samples,
                            seq_id,
                            sample_id,
                            region_start,
                            region_end,
                            &global_used,
                        );
                        tracker.mark_processed_batch(&result.discovered_regions);
                        emit_chunk(chunk_id, result)?;
                    } else if gap_len > TRANSITIVE_CHUNK_SIZE {
                        // Multi-chunk: split the gap into TRANSITIVE_CHUNK_SIZE
                        // slices so each `process_anchor_region` is its own
                        // ckpt-able unit.
                        let mut pos = region_start;
                        while pos < region_end {
                            let chunk_end = (pos + TRANSITIVE_CHUNK_SIZE).min(region_end);
                            let chunk_id = encode_chunk_id_phase2(seq_id, pos);
                            let result = process_anchor_region(
                                impg,
                                config,
                                &compact_lengths,
                                num_samples,
                                seq_id,
                                sample_id,
                                pos,
                                chunk_end,
                                &seq_included,
                                min_seq_length,
                                &global_used,
                            );
                            tracker.mark_processed_batch(&result.discovered_regions);
                            emit_chunk(chunk_id, result)?;
                            pos = chunk_end;
                        }
                    } else {
                        let chunk_id = encode_chunk_id_phase2(seq_id, region_start);
                        let result = process_anchor_region(
                            impg,
                            config,
                            &compact_lengths,
                            num_samples,
                            seq_id,
                            sample_id,
                            region_start,
                            region_end,
                            &seq_included,
                            min_seq_length,
                            &global_used,
                        );
                        tracker.mark_processed_batch(&result.discovered_regions);
                        emit_chunk(chunk_id, result)?;
                    }
                }

                let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
                pb_depth.set_position(count as u64);

                Ok(())
            })?;

        pb_depth.finish_and_clear();
    } else {
        // Phase 2 transitive — raw BFS path (default).
        //
        // Mirrors Phase 1 transitive's batched architecture (depth.rs above,
        // `batch_depth_bfs` + `process_anchor_region_transitive_raw_with_hits`).
        // The previous design ran one `par_iter` task per Phase 2 sequence and
        // each task independently called `query_transitive_bfs` + per-hop file
        // loads. With the Phase 2 working set spread across hundreds of leaf
        // samples whose per-file index sets do not overlap, the
        // `transient_header_cache` thrashed: every worker pulled in a different
        // sub-index per BFS hop, blowing out the cache cap, evicting, then
        // immediately re-loading. Phase 1 hides this because all hub chunks
        // reference roughly the same sub-index set, so the cache is bounded
        // and fully warm.
        //
        // Two-pass design:
        //   1) Sort Phase 2 chunks into transitive (gap >= min_transitive_len,
        //      sliced at TRANSITIVE_CHUNK_SIZE) and non-transitive fallback
        //      (gap < min_transitive_len, processed as a 1-hop sweep — the BFS
        //      `min_transitive_len` gate would otherwise drop these regions
        //      entirely).
        //   2) Process each pool in batches:
        //      - Non-transitive fallback: `batch_query_raw_overlapping` for the
        //        batch, then parallel sweep + emit. Same pattern as Phase 2
        //        non-transitive above.
        //      - Transitive: `batch_depth_bfs` for the batch (sequential file
        //        loading per BFS round across all chunks), then parallel
        //        sweep + emit via `process_anchor_region_transitive_raw_with_hits`.
        //
        // Memory bound: `batch_depth_bfs` keeps every chunk's BFS state alive
        // for the duration of the batch (queue + visited_ranges + hits), so
        // Phase 2 — which can have orders of magnitude more chunks than
        // Phase 1 — must process in fixed-size batches. PHASE2_TRANS_BATCH_SIZE
        // is the BFS pool size; PHASE2_NONTRANS_BATCH_SIZE is the fallback pool
        // size (matches the non-transitive path above).
        const PHASE2_NONTRANS_BATCH_SIZE: usize = 65_536;
        const PHASE2_TRANS_BATCH_SIZE: usize = 8_192;

        let mut transitive_chunks: Vec<(u32, i64, i64)> = Vec::new();
        let mut nontrans_chunks: Vec<(u32, i64, i64)> = Vec::new();

        for &seq_id in &phase2_seqs {
            let seq_len = compact_lengths.get_length(seq_id);
            if seq_len <= 0 {
                continue;
            }
            for (region_start, region_end) in tracker.get_unprocessed(seq_id, 0, seq_len) {
                let gap_len = region_end - region_start;
                if gap_len <= 0 {
                    continue;
                }
                if gap_len < config.min_transitive_len {
                    nontrans_chunks.push((seq_id, region_start, region_end));
                } else {
                    let mut pos = region_start;
                    while pos < region_end {
                        let chunk_end = (pos + TRANSITIVE_CHUNK_SIZE).min(region_end);
                        transitive_chunks.push((seq_id, pos, chunk_end));
                        pos = chunk_end;
                    }
                }
            }
        }

        let total_chunks = transitive_chunks.len() + nontrans_chunks.len();
        info!(
            "Phase 2: {} sequences split into {} transitive chunks ({}MB each) + {} sub-{}bp fallback chunks, batch BFS",
            phase2_seqs.len(),
            transitive_chunks.len(),
            TRANSITIVE_CHUNK_SIZE / 1_000_000,
            nontrans_chunks.len(),
            config.min_transitive_len,
        );

        let pb_depth = ProgressBar::new(total_chunks as u64);
        pb_depth.set_style(
            ProgressStyle::default_bar()
                .template(&format!(
                    "{{spinner:.green}} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} chunks ({{eta}}) | {}remaining sequences",
                    phase2_label
                ))
                .unwrap()
                .progress_chars("#>-")
        );
        let phase2_chunk_count = AtomicUsize::new(0);

        // Pass A: non-transitive fallback chunks (sub-min_transitive_len gaps).
        // Same batched query pattern as the non-transitive Phase 2 branch.
        for batch in nontrans_chunks.chunks(PHASE2_NONTRANS_BATCH_SIZE) {
            let queries: Vec<(u32, i64, i64)> =
                batch.iter().map(|&(s, cs, ce)| (s, cs, ce)).collect();
            let all_alns = impg.batch_query_raw_overlapping(&queries);

            batch
                .par_iter()
                .zip(all_alns.into_par_iter())
                .try_for_each(
                    |(&(seq_id, region_start, region_end), mut raw_alns)| -> io::Result<()> {
                        let sample_id = compact_lengths.get_sample_id(seq_id);

                        if min_seq_length > 0 {
                            raw_alns.retain(|aln| {
                                seq_included
                                    .get(aln.query_id as usize)
                                    .copied()
                                    .unwrap_or(false)
                            });
                        }
                        raw_alns.sort_unstable_by_key(|aln| aln.target_start);

                        let result = process_anchor_region_raw(
                            &raw_alns,
                            &compact_lengths,
                            num_samples,
                            seq_id,
                            sample_id,
                            region_start,
                            region_end,
                            &global_used,
                        );
                        tracker.mark_processed_batch(&result.discovered_regions);

                        let mut buf: Vec<u8> = Vec::new();
                        write_results(vec![result], &mut buf)?;
                        if !buf.is_empty() {
                            if let Some(ref w) = writer {
                                w.send(std::mem::take(&mut buf))?;
                            }
                        }

                        let count = phase2_chunk_count.fetch_add(1, Ordering::Relaxed) + 1;
                        pb_depth.set_position(count as u64);
                        Ok(())
                    },
                )?;
        }

        // Pass B: transitive chunks via batch BFS.
        // Each batch:
        //   Phase A) `batch_depth_bfs` drives all chunks' frontiers together;
        //            files load sequentially (one at a time) and serve every
        //            query in the batch that references them. Single-threaded
        //            file I/O bounds peak retained sub-index memory at one.
        //   Phase B) Parallel sweep + emit via
        //            `process_anchor_region_transitive_raw_with_hits` —
        //            CPU-only, no further file I/O.
        for batch in transitive_chunks.chunks(PHASE2_TRANS_BATCH_SIZE) {
            let chunk_coords: Vec<(u32, i64, i64)> =
                batch.iter().map(|&(s, cs, ce)| (s, cs, ce)).collect();
            let all_hits = batch_depth_bfs(
                impg,
                &chunk_coords,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
            );

            batch
                .par_iter()
                .zip(all_hits.into_par_iter())
                .try_for_each(
                    |(&(seq_id, chunk_start, chunk_end), hits)| -> io::Result<()> {
                        let sample_id = compact_lengths.get_sample_id(seq_id);
                        let result = process_anchor_region_transitive_raw_with_hits(
                            hits,
                            &compact_lengths,
                            num_samples,
                            seq_id,
                            sample_id,
                            chunk_start,
                            chunk_end,
                            &seq_included,
                            min_seq_length,
                            &global_used,
                        );
                        tracker.mark_processed_batch(&result.discovered_regions);

                        let mut buf: Vec<u8> = Vec::new();
                        write_results(vec![result], &mut buf)?;
                        if !buf.is_empty() {
                            if let Some(ref w) = writer {
                                w.send(std::mem::take(&mut buf))?;
                            }
                        }

                        let count = phase2_chunk_count.fetch_add(1, Ordering::Relaxed) + 1;
                        pb_depth.set_position(count as u64);
                        Ok(())
                    },
                )?;
        }

        pb_depth.finish_and_clear();
        // Update overall processed counter so any downstream readers see the
        // full Phase 1 + Phase 2 sequence count.
        processed_count.store(phase1_seqs.len() + phase2_seqs.len(), Ordering::Relaxed);
    }

    // Clear sub-index cache once after Phase 2 transitive processing.
    // Non-transitive already uses transient queries and has nothing to clear.
    // (The progress bar was already finished within each Phase 2 branch.)
    if is_transitive {
        impg.clear_sub_index_cache();
    }

    // Handle FAI sequences not in alignment index (depth=1 for unaligned sequences)
    //
    // Resume safety: when checkpointing is active and the FAI loop committed
    // successfully on a previous run, the `CHUNK_ID_FAI_DONE` sentinel is in
    // `completed_chunks` (re-built from the work-log). Skipping the loop
    // here is what prevents duplicate FAI rows: the previously-flushed rows
    // are already on disk past the resume-truncation offset, and re-running
    // the loop would append them a second time.
    let fai_already_done = work_log_active && completed_chunks.contains(&CHUNK_ID_FAI_DONE);
    if fai_already_done && fai_seq_lengths.is_some() {
        debug!("FAI emission already committed by a previous run; skipping");
    }
    if let Some(ref fai) = fai_seq_lengths {
        if !fai_already_done {
        let indexed_seqs: FxHashSet<&str> = (0..impg.seq_index().len() as u32)
            .filter_map(|id| impg.seq_index().get_name(id))
            .collect();

        for (seq_name, &seq_len) in &fai.lengths {
            if indexed_seqs.contains(seq_name.as_str()) {
                continue;
            }

            let sample_name = extract_sample(seq_name, separator);

            // ref_only filter for FAI sequences
            if ref_only {
                if let Some(ref_name) = ref_sample {
                    if sample_name != ref_name {
                        continue;
                    }
                }
            }

            if stats_mode {
                // Add depth=1 to stats accumulators (pangenome_bases = seq_len for single sample)
                if let Some(ref acc) = stats_accumulator {
                    acc.lock().add_interval(seq_name, 0, seq_len, 1, seq_len);
                }
                if let Some(ref acc) = stats_combined_acc {
                    acc.lock().add_interval(
                        seq_name,
                        0,
                        seq_len,
                        1,
                        vec![sample_name.to_string()],
                        seq_len,
                    );
                }
            } else if let Some(ref w) = writer {
                let rid = row_counter.fetch_add(1, Ordering::Relaxed) + 1;

                // Output as depth=1 for entire sequence. Match the normal row
                // format: id, length, depth, positions (semicolon-separated);
                // for an unaligned FAI scaffold the positions list has exactly
                // one entry — the sequence itself.
                let _ = num_samples;
                let _ = sample_index;
                let _ = sample_name;
                let mut buf: Vec<u8> = Vec::with_capacity(64);
                writeln!(
                    &mut buf,
                    "{}\t{}\t1\t{}:0-{}",
                    rid, seq_len, seq_name, seq_len
                )?;
                w.send(buf)?;
            }
        }

        // Mark the FAI loop done in the work-log with a sentinel record + a
        // forced commit. Subsequent resumes will see the sentinel in
        // `completed_chunks` and skip the FAI loop, preventing duplicate
        // rows.
        if work_log_active {
            if let Some(ref w) = writer {
                let work_buf = encode_work_record(CHUNK_ID_FAI_DONE, &[]);
                w.send_chunk_bundle(Vec::new(), work_buf)?;
            }
            checkpoint_ctrl.force_commit()?;
        }
        } // end !fai_already_done
    }

    // Finalize output
    if stats_mode {
        let prefix = output_prefix.unwrap_or("depth_stats");

        if let Some(acc) = stats_accumulator {
            let stats = acc.into_inner();

            // Write summary
            let summary_path = format!("{}.summary.txt", prefix);
            let mut summary_writer = BufWriter::new(std::fs::File::create(&summary_path)?);
            stats.write_summary(&mut summary_writer)?;
            summary_writer.flush()?;
            info!("Wrote summary to {}", summary_path);

            // Print summary to stdout
            stats.write_summary(&mut std::io::stdout())?;

            // Write per-depth BED files
            stats.write_depth_bed_files(prefix)?;

            info!(
                "Stats complete: {} total bases, max depth = {}",
                stats.total_bases,
                stats.max_depth()
            );
        }

        if let Some(acc) = stats_combined_acc {
            let mut stats = acc.into_inner();

            // Write summary
            let summary_path = format!("{}.summary.txt", prefix);
            let mut summary_writer = BufWriter::new(std::fs::File::create(&summary_path)?);
            stats.write_summary(&mut summary_writer)?;
            summary_writer.flush()?;
            info!("Wrote summary to {}", summary_path);

            // Print summary to stdout
            stats.write_summary(&mut std::io::stdout())?;

            // Write combined output file (with optional merging)
            stats.write_combined_output(prefix, min_interval_len)?;

            info!(
                "Stats complete: {} total bases, {} intervals, max depth = {}",
                stats.total_bases,
                stats.intervals.len(),
                stats.max_depth()
            );
        }
    } else {
        // Normal mode: drop the channel sender and join the dedicated writer
        // thread. CRITICAL: without this join, the program may exit before
        // the writer drains its queue, producing a truncated TSV file.
        //
        // When checkpointing is active we issue one final ckpt commit before
        // the writer is consumed (`finish()` drops the sender). The final
        // ckpt is what tells a future `--resume` invocation that the run
        // ended successfully — but to be doubly safe we also remove the
        // ckpt + work-log on success so the next run starts from a clean
        // slate.
        if checkpoint_ctrl.enabled() {
            checkpoint_ctrl.force_commit()?;
        }
        // The controller borrows `&writer` (alongside `&row_counter` etc.);
        // we have to drop it *before* moving `writer` into `finish()`.
        drop(checkpoint_ctrl);
        if let Some(w) = writer {
            w.finish()?;
        }
        if resume {
            if let Some(prefix) = output_prefix {
                DepthCheckpoint::cleanup_on_success(prefix);
            }
        }

        let total_intervals = row_counter.load(Ordering::Relaxed);
        info!(
            "Depth complete: {} sequences, {} intervals output",
            total_sequences, total_intervals
        );
    }

    Ok(())
}

// ============================================================================
// Region query mode and sweep-line helpers
// ============================================================================

/// Sweep-line over `CompactAlignmentInfo` (region-mode variant).
///
/// Unifies the region-query sweep with the global `sweep_line_depth` path: both
/// consume the same compact alignment representation and reuse `proj_offset` /
/// `map_target_to_query_linear` for coordinate mapping.
///
/// The region variant differs from `sweep_line_depth` in two ways:
///   1. Per-sample column output requires emitting **all** active alignments
///      per sample in a window (joined with `;`), not just the highest-overlap
///      alignment. `RegionDepthResult::add_sample_position` dedupes coincident
///      tuples, matching the pre-Stage-4 behavior.
///   2. The output is `RegionDepthResult` (string-keyed for the writer's column
///      layout) instead of `SparseDepthInterval`. Sample/seq IDs are translated
///      to names exactly once per active sample per window via the supplied
///      `SampleIndex` and the impg `SequenceIndex`.
///
/// Position tuples within a row are sorted by sample name in the caller's
/// output formatter (cf. `write_region_depth_output`); inside each sample's
/// position list the order matches alignment-insertion order, identical to the
/// pre-Stage-4 String-based path.
fn compute_region_sweep_compact(
    anchor_seq: &str,
    alignments: &[CompactAlignmentInfo],
    num_samples: usize,
    sample_index: &SampleIndex,
    seq_index: &crate::seqidx::SequenceIndex,
    config: &DepthConfig,
) -> Vec<RegionDepthResult> {
    if alignments.is_empty() {
        return Vec::new();
    }

    // Build sweep-line events (compact: u16 sample IDs + alignment idx).
    let mut events: Vec<CompactDepthEvent> = Vec::with_capacity(alignments.len() * 2);
    for (idx, aln) in alignments.iter().enumerate() {
        events.push(CompactDepthEvent {
            position: aln.target_start,
            is_start: true,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
        events.push(CompactDepthEvent {
            position: aln.target_end,
            is_start: false,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
    }
    events.sort_by_key(|e| e.packed_sort_key());

    // Sweep-line: track active alignments per sample (Vec<idx> indexed by sample_id).
    let mut results: Vec<RegionDepthResult> = Vec::new();
    let mut active_bitmap = SampleBitmap::new(num_samples);
    let mut active_alns: Vec<Vec<usize>> = vec![Vec::new(); num_samples];
    let mut prev_pos: Option<i64> = None;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev && active_bitmap.depth() > 0 {
                let mut result =
                    RegionDepthResult::new(anchor_seq.to_string(), prev, event.position);

                // Iterate active samples in u16 order (deterministic).
                for sample_id in active_bitmap.active_samples() {
                    let alns = &active_alns[sample_id as usize];
                    if alns.is_empty() {
                        continue;
                    }
                    let sample_name = sample_index.get_name(sample_id).unwrap_or("?");
                    for &idx in alns {
                        let aln = &alignments[idx];
                        let (q_start, q_end) = map_target_to_query_linear(
                            &[],
                            aln.target_start,
                            aln.target_end,
                            aln.query_start,
                            aln.query_end,
                            prev,
                            event.position,
                            aln.is_reverse,
                        );
                        let seq_name = seq_index.get_name(aln.query_id).unwrap_or("?");
                        result.add_sample_position(sample_name, seq_name, q_start, q_end);
                    }
                }

                result.update_depth();
                if result.depth > 0 {
                    results.push(result);
                }
            }
        }

        // Update active alignments
        if event.is_start {
            active_bitmap.add(event.sample_id);
            active_alns[event.sample_id as usize].push(event.alignment_idx);
        } else {
            active_bitmap.remove(event.sample_id);
            let v = &mut active_alns[event.sample_id as usize];
            if let Some(pos) = v.iter().position(|&idx| idx == event.alignment_idx) {
                v.swap_remove(pos);
            }
        }

        prev_pos = Some(event.position);
    }

    // Merge adjacent windows with same depth if configured
    if config.merge_adjacent {
        results = merge_adjacent_results(results);
    }

    results
}

/// Merge adjacent results with same depth
fn merge_adjacent_results(mut results: Vec<RegionDepthResult>) -> Vec<RegionDepthResult> {
    if results.len() <= 1 {
        return results;
    }

    let mut merged: Vec<RegionDepthResult> = Vec::new();
    let mut current = results.remove(0);

    for next in results {
        // Check if adjacent and same depth
        if current.ref_end == next.ref_start
            && current.depth == next.depth
            && current.ref_seq == next.ref_seq
        {
            // Merge
            current.ref_end = next.ref_end;
            // Merge sample positions
            for (sample, positions) in next.sample_positions {
                current
                    .sample_positions
                    .entry(sample)
                    .or_default()
                    .extend(positions);
            }
        } else {
            merged.push(current);
            current = next;
        }
    }
    merged.push(current);

    merged
}

/// Query depth for a specific region with per-sample position tracking
///
/// Algorithm:
/// 1. Parse the target region (seq_name:start-end)
/// 2. Find all alignments overlapping this region
/// 3. For each overlapping alignment, track the sample and its position
/// 4. If sample has multiple alignments, track all of them
/// 5. Output in tabular format with per-sample columns
pub fn query_region_depth(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    target_seq: &str,
    target_start: i64,
    target_end: i64,
    separator: &str,
    sample_filter: Option<&SampleFilter>,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<Vec<RegionDepthResult>> {
    debug!(
        "Querying region depth: {}:{}-{}",
        target_seq, target_start, target_end
    );

    let target_id = impg.seq_index().get_id(target_seq).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Sequence '{}' not found in index", target_seq),
        )
    })?;

    // Build compact sample/seq lookup once for this region. Maps query_seq_id ->
    // sample_id (u16) at intake instead of per-alignment String clones.
    let compact_lengths = CompactSequenceLengths::from_impg(impg, separator);
    let sample_idx = compact_lengths.sample_index();
    let num_samples = sample_idx.len();

    // Anchor sample id (used for self-alignment filtering).
    let target_sample_str = extract_sample(target_seq, separator);
    let anchor_sample_id = sample_idx.get_id(&target_sample_str).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Anchor sample '{}' not found after compact-index build",
                target_sample_str
            ),
        )
    })?;

    // Resolve the SampleFilter to a BitVec mask once (filter -> u16-id set).
    // `include_all == true` short-circuits the per-alignment lookup.
    let (include_all, sample_mask): (bool, BitVec) = match sample_filter {
        Some(f) if f.is_active() => {
            let mut mask: BitVec = bitvec![0; num_samples];
            for name in f.get_samples() {
                if let Some(sid) = sample_idx.get_id(name) {
                    mask.set(sid as usize, true);
                }
            }
            (false, mask)
        }
        _ => (true, BitVec::new()),
    };
    let sample_allowed = |sid: u16| -> bool {
        include_all || sample_mask.get(sid as usize).map_or(false, |b| *b)
    };

    // Collect all alignments for this region (compact: u16 sample, u32 seq).
    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    let is_transitive = config.transitive || config.transitive_dfs;

    // Forward direction: where target_seq is TARGET
    if is_transitive && config.use_cigar_bfs {
        // CIGAR-precise BFS (--use-BFS): identical to query command's BFS
        let overlaps = if config.transitive_dfs {
            impg.query_transitive_dfs(
                target_id,
                target_start,
                target_end,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                false,
                None,
            )
        } else {
            impg.query_transitive_bfs(
                target_id,
                target_start,
                target_end,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                false,
                None,
            )
        };

        let region_start = target_start;
        let region_end = target_end;

        // Pass 1: Build a map of each sequence's anchor coverage from hop 0 results
        let mut seq_anchor_coverage: FxHashMap<u32, Vec<HopZeroSeg>> = FxHashMap::default();
        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;
            if target_interval.metadata == target_id {
                let query_id = query_interval.metadata;
                let q_start = query_interval.first.min(query_interval.last) as i64;
                let q_end = query_interval.first.max(query_interval.last) as i64;
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                seq_anchor_coverage
                    .entry(query_id)
                    .or_default()
                    .push((q_start, q_end, t_start, t_end));
            }
        }
        depth_trace!(
            "PASS1 stage=region_use_bfs target={} region={}-{} overlaps={} keys={}",
            target_id,
            region_start,
            region_end,
            overlaps.len(),
            seq_anchor_coverage.len()
        );

        // Pass 2: Process all results for depth
        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;

            let query_id = query_interval.metadata;
            let query_sample_id = compact_lengths.get_sample_id(query_id);

            if is_self_alignment(query_sample_id, anchor_sample_id, query_id, target_interval.metadata) {
                continue;
            }
            if !sample_allowed(query_sample_id) {
                continue;
            }

            let is_reverse = query_interval.first > query_interval.last;
            let q_start = query_interval.first.min(query_interval.last) as i64;
            let q_end = query_interval.first.max(query_interval.last) as i64;

            let (a_start, a_end) = if target_interval.metadata == target_id {
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                (t_start.max(region_start), t_end.min(region_end))
            } else {
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                let p = project_hop0_coords(
                    seq_anchor_coverage.get(&target_interval.metadata),
                    t_start,
                    t_end,
                    region_start,
                    region_end,
                );
                depth_trace!(
                    "HOP2 stage=region_use_bfs sample_id={} q_id={} t_id={} t={}-{} a={}-{}",
                    query_sample_id, query_id, target_interval.metadata,
                    t_start, t_end, p.0, p.1
                );
                p
            };

            if a_start >= a_end {
                continue;
            }

            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                query_id,
                q_start,
                q_end,
                a_start,
                a_end,
                is_reverse,
            ));
        }
    } else if is_transitive {
        // Default transitive: raw-interval BFS with linear interpolation
        let hits = depth_transitive_bfs(
            impg,
            target_id,
            target_start,
            target_end,
            config.max_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            config.transitive_dfs,
        );

        let region_start = target_start;
        let region_end = target_end;

        // Pass 1: Build anchor coverage map from hop 0 results
        let mut seq_anchor_coverage: FxHashMap<u32, Vec<HopZeroSeg>> = FxHashMap::default();
        for hit in &hits {
            if hit.target_id == target_id {
                let q_start = hit.query_start.min(hit.query_end) as i64;
                let q_end = hit.query_start.max(hit.query_end) as i64;
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                seq_anchor_coverage
                    .entry(hit.query_id)
                    .or_default()
                    .push((q_start, q_end, t_start, t_end));
            }
        }
        depth_trace!(
            "PASS1 stage=region_raw_bfs target={} region={}-{} hits={} keys={}",
            target_id,
            region_start,
            region_end,
            hits.len(),
            seq_anchor_coverage.len()
        );

        // Pass 2: Process all hits for depth
        for hit in &hits {
            let query_sample_id = compact_lengths.get_sample_id(hit.query_id);

            if is_self_alignment(query_sample_id, anchor_sample_id, hit.query_id, hit.target_id) {
                continue;
            }
            if !sample_allowed(query_sample_id) {
                continue;
            }

            let q_start = hit.query_start.min(hit.query_end) as i64;
            let q_end = hit.query_start.max(hit.query_end) as i64;

            let (a_start, a_end) = if hit.target_id == target_id {
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                (t_start.max(region_start), t_end.min(region_end))
            } else {
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                let p = project_hop0_coords(
                    seq_anchor_coverage.get(&hit.target_id),
                    t_start,
                    t_end,
                    region_start,
                    region_end,
                );
                depth_trace!(
                    "HOP2 stage=region_raw_bfs sample_id={} q_id={} t_id={} t={}-{} a={}-{}",
                    query_sample_id, hit.query_id, hit.target_id,
                    t_start, t_end, p.0, p.1
                );
                p
            };

            if a_start >= a_end {
                continue;
            }

            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                hit.query_id,
                q_start,
                q_end,
                a_start,
                a_end,
                hit.is_reverse,
            ));
        }
    } else if config.use_cigar_bfs {
        // Non-transitive with --use-BFS: CIGAR-precise query
        let overlaps = impg.query(
            target_id,
            target_start,
            target_end,
            false,
            None,
            sequence_index,
            false,
        );

        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;

            let query_id = query_interval.metadata;
            let query_sample_id = compact_lengths.get_sample_id(query_id);

            if !sample_allowed(query_sample_id) {
                continue;
            }

            let is_reverse = query_interval.first > query_interval.last;
            let query_start = query_interval.first.min(query_interval.last) as i64;
            let query_end = query_interval.first.max(query_interval.last) as i64;
            let t_start = target_interval.first.min(target_interval.last) as i64;
            let t_end = target_interval.first.max(target_interval.last) as i64;

            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                query_id,
                query_start,
                query_end,
                t_start,
                t_end,
                is_reverse,
            ));
        }
    } else {
        // Default non-transitive: raw intervals + linear interpolation
        let mut raw_alns = impg.query_raw_overlapping(target_id, target_start, target_end);
        raw_alns.sort_unstable_by_key(|a| a.target_start);

        let region_start = target_start;
        let region_end = target_end;

        for aln in &raw_alns {
            let query_sample_id = compact_lengths.get_sample_id(aln.query_id);

            if !sample_allowed(query_sample_id) {
                continue;
            }

            let aln_target_start = aln.target_start as i64;
            let aln_target_end = aln.target_end as i64;
            let aln_query_start = aln.query_start as i64;
            let aln_query_end = aln.query_end as i64;

            // Clip target to region
            let clipped_target_start = aln_target_start.max(region_start);
            let clipped_target_end = aln_target_end.min(region_end);

            // Proportional query coordinate clipping (linear interpolation)
            let target_len = aln_target_end - aln_target_start;
            let query_len = aln_query_end - aln_query_start;
            let (clipped_query_start, clipped_query_end) = if target_len > 0 {
                let off_s = clipped_target_start - aln_target_start;
                let off_e = clipped_target_end - aln_target_start;
                if aln.is_reverse {
                    let cqe = aln_query_end - proj_offset(off_s, query_len, target_len);
                    let cqs = aln_query_end - proj_offset(off_e, query_len, target_len);
                    (cqs, cqe)
                } else {
                    let cqs = aln_query_start + proj_offset(off_s, query_len, target_len);
                    let cqe = aln_query_start + proj_offset(off_e, query_len, target_len);
                    (cqs, cqe)
                }
            } else {
                (aln_query_start, aln_query_end)
            };

            alignments.push(CompactAlignmentInfo::new(
                query_sample_id,
                aln.query_id,
                clipped_query_start.min(clipped_query_end),
                clipped_query_start.max(clipped_query_end),
                clipped_target_start,
                clipped_target_end,
                aln.is_reverse,
            ));
        }
    }

    // Reverse direction: where target_seq is QUERY.
    // Only needed for V1 (unidirectional) indices — V2 bidirectional indices already
    // include reverse-direction entries in every tree, so query_raw_overlapping /
    // query_transitive_bfs above has already captured them.
    if !impg.is_bidirectional() {
        let reverse_alignments = impg.query_reverse_for_depth(target_id);
        for (ref_start, ref_end, other_t_start, other_t_end, other_id) in reverse_alignments {
            // Check if our sequence's interval overlaps with query region
            if ref_end <= target_start || ref_start >= target_end {
                continue;
            }

            let other_sample_id = compact_lengths.get_sample_id(other_id);

            // Skip self
            if is_self_alignment(other_sample_id, anchor_sample_id, other_id, target_id) {
                continue;
            }
            if !sample_allowed(other_sample_id) {
                continue;
            }

            // Clip our-sequence coordinates to the queried region
            let clipped_ref_start = ref_start.max(target_start);
            let clipped_ref_end = ref_end.min(target_end);

            // Proportionally clip the other sequence's coordinates
            let our_len = ref_end - ref_start;
            let other_len = other_t_end - other_t_start;
            let (clipped_other_start, clipped_other_end) = if our_len > 0 {
                let off_s = clipped_ref_start - ref_start;
                let off_e = clipped_ref_end - ref_start;
                let cs = other_t_start + proj_offset(off_s, other_len, our_len);
                let ce = other_t_start + proj_offset(off_e, other_len, our_len);
                (cs.min(ce), cs.max(ce))
            } else {
                (other_t_start, other_t_end)
            };

            alignments.push(CompactAlignmentInfo::new(
                other_sample_id,
                other_id,
                clipped_other_start,
                clipped_other_end,
                clipped_ref_start,
                clipped_ref_end,
                false,
            ));
        }
    }

    // Add self (target sample) if in filter
    if sample_allowed(anchor_sample_id) {
        alignments.push(CompactAlignmentInfo::new(
            anchor_sample_id,
            target_id,
            target_start,
            target_end,
            target_start,
            target_end,
            false,
        ));
    }

    // Deduplicate alignments before sweep-line. Two sources of duplicates:
    //   1. V2 bidirectional indices store every alignment twice (target→query
    //      and query→target). Querying one side returns both copies.
    //   2. All-vs-all PAF inputs (e.g. both A_vs_B.paf and B_vs_A.paf) further
    //      double each entry because each direction is indexed separately.
    //   3. The non-transitive paths above include the input range as
    //      `impg.query`'s result[0], and we also unconditionally push self at
    //      the end of this function — duplicating the anchor sample.
    // `impg query` collapses these via `merge_query_adjusted_intervals` at
    // output time; depth's sweep-line, by contrast, would emit `pos;pos` in
    // the per-sample column. Dedup the alignment vec in place to match.
    // Sort key uses u32 query_id (compact) instead of the previous String name —
    // ordering still uniquely keys identical alignments.
    alignments.sort_by(|a, b| {
        a.query_id
            .cmp(&b.query_id)
            .then(a.query_start.cmp(&b.query_start))
            .then(a.query_end.cmp(&b.query_end))
            .then(a.target_start.cmp(&b.target_start))
            .then(a.target_end.cmp(&b.target_end))
            .then(a.is_reverse.cmp(&b.is_reverse))
    });
    alignments.dedup_by(|a, b| {
        a.query_id == b.query_id
            && a.query_start == b.query_start
            && a.query_end == b.query_end
            && a.target_start == b.target_start
            && a.target_end == b.target_end
            && a.is_reverse == b.is_reverse
    });

    // Compute depth windows using the unified compact sweep-line.
    #[cfg(feature = "depth-trace")]
    {
        let unique: std::collections::BTreeSet<u16> =
            alignments.iter().map(|a| a.sample_id).collect();
        depth_trace!(
            "SWEEP stage=region target_seq={} target_range={}-{} alignments={} unique_samples={} sample_ids={:?}",
            target_seq,
            target_start,
            target_end,
            alignments.len(),
            unique.len(),
            unique
        );
    }
    let results = compute_region_sweep_compact(
        target_seq,
        &alignments,
        num_samples,
        sample_idx,
        impg.seq_index(),
        config,
    );

    // Filter results to only include the query region
    let filtered_results: Vec<RegionDepthResult> = results
        .into_iter()
        .filter(|r| r.ref_end > target_start as i64 && r.ref_start < target_end as i64)
        .map(|mut r| {
            // Clip to query region
            r.ref_start = r.ref_start.max(target_start as i64);
            r.ref_end = r.ref_end.min(target_end as i64);
            r
        })
        .collect();

    Ok(filtered_results)
}

/// Parse target range string in format "seq_name:start-end"
pub fn parse_target_range_depth(target_range: &str) -> io::Result<(String, i64, i64)> {
    // Handle format: seq_name:start-end
    // Note: seq_name may contain colons (e.g., sample#hap#chr)
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid target range format '{}'. Expected format: seq_name:start-end",
                target_range
            ),
        ));
    }

    let range_str = parts[0];
    let seq_name = parts[1].to_string();

    let range_parts: Vec<&str> = range_str.split('-').collect();
    if range_parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid range format '{}'. Expected format: start-end",
                range_str
            ),
        ));
    }

    let start: i64 = range_parts[0].parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Invalid start position: {}", range_parts[0]),
        )
    })?;

    let end: i64 = range_parts[1].parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Invalid end position: {}", range_parts[1]),
        )
    })?;

    if start >= end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start ({}) must be less than end ({})", start, end),
        ));
    }

    Ok((seq_name, start, end))
}

/// Parse BED file for region queries
pub fn parse_bed_file_depth(bed_path: &str) -> io::Result<Vec<(String, i64, i64)>> {
    let content = std::fs::read_to_string(bed_path)?;
    let mut regions = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let seq_name = parts[0].to_string();
        let start: i64 = parts[1].parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid start position in BED file: {}", parts[1]),
            )
        })?;
        let end: i64 = parts[2].parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid end position in BED file: {}", parts[2]),
            )
        })?;

        regions.push((seq_name, start, end));
    }

    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_unprocessed_subinterval_semantics() {
        let tracker = ConcurrentProcessedTracker::new(1);
        tracker.mark_processed(0, 0, 5_000_000);
        tracker.mark_processed(0, 6_000_000, 8_000_000);

        let unprocessed = tracker.get_unprocessed(0, 4_000_000, 5_500_000);
        assert_eq!(unprocessed, vec![(5_000_000, 5_500_000)]);

        let fully_covered = tracker.get_unprocessed(0, 0, 5_000_000);
        assert!(fully_covered.is_empty());

        let fully_open = tracker.get_unprocessed(0, 8_000_000, 10_000_000);
        assert_eq!(fully_open, vec![(8_000_000, 10_000_000)]);
    }

    #[test]
    fn test_interval_set_add_merge_touching_and_overlapping() {
        let mut s = IntervalSet::new();
        s.add(10, 20);
        s.add(30, 40);
        s.add(20, 30); // touches both neighbors -> merged
        assert_eq!(s.intervals().collect::<Vec<_>>(), vec![(10, 40)]);
        assert_eq!(s.total_length(), 30);

        // Overlapping add that subsumes multiple
        let mut t = IntervalSet::new();
        t.add(0, 10);
        t.add(20, 30);
        t.add(40, 50);
        t.add(5, 45);
        assert_eq!(t.intervals().collect::<Vec<_>>(), vec![(0, 50)]);
        assert_eq!(t.total_length(), 50);

        // Idempotent re-add
        let mut u = IntervalSet::new_single(0, 100);
        u.add(30, 70);
        assert_eq!(u.intervals().collect::<Vec<_>>(), vec![(0, 100)]);
        assert_eq!(u.total_length(), 100);
    }

    #[test]
    fn test_interval_set_subtract_splits_and_trims() {
        // Full split
        let mut s = IntervalSet::new_single(0, 100);
        s.subtract(30, 70);
        assert_eq!(s.intervals().collect::<Vec<_>>(), vec![(0, 30), (70, 100)]);
        assert_eq!(s.total_length(), 60);

        // Left trim
        let mut t = IntervalSet::new_single(0, 100);
        t.subtract(0, 40);
        assert_eq!(t.intervals().collect::<Vec<_>>(), vec![(40, 100)]);

        // Right trim
        let mut u = IntervalSet::new_single(0, 100);
        u.subtract(60, 100);
        assert_eq!(u.intervals().collect::<Vec<_>>(), vec![(0, 60)]);

        // Subtract spanning multiple
        let mut v = IntervalSet::new();
        v.add(0, 20);
        v.add(30, 50);
        v.add(60, 80);
        v.subtract(10, 70);
        assert_eq!(v.intervals().collect::<Vec<_>>(), vec![(0, 10), (70, 80)]);
        assert_eq!(v.total_length(), 20);

        // No-op (disjoint)
        let mut w = IntervalSet::new_single(0, 50);
        w.subtract(100, 200);
        assert_eq!(w.intervals().collect::<Vec<_>>(), vec![(0, 50)]);
        assert_eq!(w.total_length(), 50);
    }

    #[test]
    fn test_interval_set_fuzz_equivalence_to_boolean_bitmap() {
        // Cross-check IntervalSet against a naive bitmap reference over
        // 10,000 random add/subtract ops in a bounded universe.
        use std::collections::BTreeMap as BM;
        let _ = BM::<i64, i64>::new(); // silence unused import if any

        // Deterministic LCG to keep the test reproducible without adding deps
        let mut state: u64 = 0xdead_beef_cafe_babe;
        let mut rand_u32 = || -> u32 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 32) as u32
        };

        const U: usize = 4096;
        let mut bitmap = vec![false; U];
        let mut set = IntervalSet::new();

        for _ in 0..10_000 {
            let a = (rand_u32() as usize) % U;
            let b = (rand_u32() as usize) % U;
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            let is_add = rand_u32() & 1 == 0;
            if is_add {
                for i in lo..hi { bitmap[i] = true; }
                set.add(lo as i64, hi as i64);
            } else {
                for i in lo..hi { bitmap[i] = false; }
                set.subtract(lo as i64, hi as i64);
            }

            // Materialize bitmap to intervals
            let mut expected: Vec<(i64, i64)> = Vec::new();
            let mut i = 0;
            while i < U {
                if bitmap[i] {
                    let start = i;
                    while i < U && bitmap[i] { i += 1; }
                    expected.push((start as i64, i as i64));
                } else {
                    i += 1;
                }
            }

            let got: Vec<(i64, i64)> = set.intervals().collect();
            assert_eq!(got, expected,
                "mismatch after op lo={} hi={} add={}", lo, hi, is_add);

            let expected_total: i64 = expected.iter().map(|(s, e)| e - s).sum();
            assert_eq!(set.total_length(), expected_total,
                "total_length mismatch: expected={}, got={}", expected_total, set.total_length());
        }
    }

    // Old Vec-backed IntervalSet, inlined verbatim for side-by-side
    // benchmarking against the BTreeMap-backed implementation.
    struct OldIntervalSetVec {
        intervals: Vec<(i64, i64)>,
        total_length: i64,
    }

    impl OldIntervalSetVec {
        fn new() -> Self { Self { intervals: Vec::new(), total_length: 0 } }
        fn new_single(s: i64, e: i64) -> Self {
            if s >= e { Self::new() }
            else { Self { intervals: vec![(s, e)], total_length: e - s } }
        }
        fn add(&mut self, start: i64, end: i64) {
            if start >= end { return; }
            if self.intervals.is_empty() {
                self.intervals.push((start, end));
                self.total_length = end - start;
                return;
            }
            let first_overlap_idx = self.intervals.partition_point(|&(_, e)| e < start);
            let last_overlap_end = first_overlap_idx
                + self.intervals[first_overlap_idx..].partition_point(|&(s, _)| s <= end);
            if first_overlap_idx == last_overlap_end {
                self.intervals.insert(first_overlap_idx, (start, end));
                self.total_length += end - start;
            } else {
                let ms = start.min(self.intervals[first_overlap_idx].0);
                let me = end.max(self.intervals[last_overlap_end - 1].1);
                let removed: i64 = self.intervals[first_overlap_idx..last_overlap_end]
                    .iter().map(|&(s, e)| e - s).sum();
                self.intervals.drain(first_overlap_idx..last_overlap_end);
                self.intervals.insert(first_overlap_idx, (ms, me));
                self.total_length = self.total_length - removed + (me - ms);
            }
        }
        fn subtract(&mut self, ss: i64, se: i64) {
            if ss >= se || self.intervals.is_empty() { return; }
            let first_overlap_idx = self.intervals.partition_point(|&(_, e)| e <= ss);
            let last_overlap_end = first_overlap_idx
                + self.intervals[first_overlap_idx..].partition_point(|&(s, _)| s < se);
            if first_overlap_idx == last_overlap_end { return; }
            let mut repl = Vec::new();
            let mut removed: i64 = 0;
            for &(s, e) in &self.intervals[first_overlap_idx..last_overlap_end] {
                if s < ss { repl.push((s, ss)); }
                if e > se { repl.push((se, e)); }
                let os = s.max(ss);
                let oe = e.min(se);
                removed += oe - os;
            }
            self.intervals.splice(first_overlap_idx..last_overlap_end, repl);
            self.total_length -= removed;
        }
        fn total_length(&self) -> i64 { self.total_length }
        fn intervals(&self) -> &[(i64, i64)] { &self.intervals }
    }

    // Simulates the inner claim_unprocessed path using the old impl.
    fn old_claim(set: &mut OldIntervalSetVec, start: i64, end: i64) -> Vec<(i64, i64)> {
        let unprocessed = if set.intervals().is_empty() {
            vec![(start, end)]
        } else {
            let mut result = OldIntervalSetVec::new_single(start, end);
            for &(s, e) in set.intervals() {
                if s >= end { break; }
                if e <= start { continue; }
                result.subtract(s, e);
            }
            result.intervals().to_vec()
        };
        for &(s, e) in &unprocessed {
            set.add(s, e);
        }
        unprocessed
    }

    fn lcg_step(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state
    }

    /// Ignored side-by-side perf smoke: old Vec-backed vs new BTreeMap-backed
    /// IntervalSet under the VGP-style claim_unprocessed workload.
    #[test]
    #[ignore]
    fn bench_interval_set_oldvsnew_claim_scaling() {
        use std::time::Instant;
        const UNIVERSE: u64 = 1_000_000_000_000;
        const CLAIM_W: u64 = 10_000;
        const OPS: usize = 1000;

        println!();
        println!(
            "  {:>10}  {:>12}  {:>12}  {:>8}",
            "size_N", "old_us/op", "new_us/op", "speedup"
        );
        for &target_size in &[1_000usize, 10_000, 100_000, 1_000_000] {
            // Seed both impls with the identical sequence of intervals so
            // their states stay in sync; count only disjoint inserts.
            let mut state = 0xdeadbeefcafebabe_u64;
            let mut old_set = OldIntervalSetVec::new();
            let mut new_set = IntervalSet::new();
            let mut inserted = 0usize;
            while inserted < target_size {
                let s = lcg_step(&mut state) % (UNIVERSE - CLAIM_W);
                let e = s + CLAIM_W;
                let b_old = old_set.total_length();
                old_set.add(s as i64, e as i64);
                new_set.add(s as i64, e as i64);
                if old_set.total_length() == b_old + CLAIM_W as i64 {
                    inserted += 1;
                }
            }
            assert_eq!(old_set.total_length(), new_set.total_length());

            // Save RNG checkpoint so both impls see identical op sequences.
            let ops_state_init = state;

            state = ops_state_init;
            let t_old = Instant::now();
            for _ in 0..OPS {
                let s = lcg_step(&mut state) % (UNIVERSE - CLAIM_W);
                let e = s + CLAIM_W;
                let _ = old_claim(&mut old_set, s as i64, e as i64);
            }
            let old_elapsed = t_old.elapsed();

            state = ops_state_init;
            let tracker = ConcurrentProcessedTracker::new(1);
            *tracker.processed[0].lock() = new_set;
            let t_new = Instant::now();
            for _ in 0..OPS {
                let s = lcg_step(&mut state) % (UNIVERSE - CLAIM_W);
                let e = s + CLAIM_W;
                let _ = tracker.claim_unprocessed(0, s as i64, e as i64);
            }
            let new_elapsed = t_new.elapsed();

            let old_us = old_elapsed.as_secs_f64() * 1e6 / OPS as f64;
            let new_us = new_elapsed.as_secs_f64() * 1e6 / OPS as f64;
            println!(
                "  {:>10}  {:>12.3}  {:>12.3}  {:>7.1}x",
                target_size,
                old_us,
                new_us,
                old_us / new_us
            );
        }
    }

    /// Ignored perf-smoke that approximates the VGP depth workload:
    /// repeatedly claim unprocessed ranges at random query offsets in a
    /// single IntervalSet that grows into the millions. With the old
    /// sorted-Vec backend this is O(N²); the BTreeMap backend keeps each
    /// op at O(log N). Run with:
    ///     cargo test --release --lib \
    ///         commands::depth::tests::bench_interval_set_claim_scaling \
    ///         -- --ignored --nocapture
    #[test]
    #[ignore]
    fn bench_interval_set_claim_scaling() {
        use std::time::Instant;

        const UNIVERSE: u64 = 1_000_000_000_000;
        const CLAIM_W: u64 = 10_000;
        const OPS: usize = 1000;

        for &target_size in &[1_000usize, 10_000, 100_000, 1_000_000] {
            let mut state = 0x1234_5678_9abc_def0_u64;
            let mut set = IntervalSet::new();
            let mut inserted = 0usize;
            while inserted < target_size {
                let s = lcg_step(&mut state) % (UNIVERSE - CLAIM_W);
                let e = s + CLAIM_W;
                let before = set.total_length();
                set.add(s as i64, e as i64);
                if set.total_length() == before + CLAIM_W as i64 {
                    inserted += 1;
                }
            }

            let tracker = ConcurrentProcessedTracker::new(1);
            *tracker.processed[0].lock() = set;

            let t0 = Instant::now();
            for _ in 0..OPS {
                let s = lcg_step(&mut state) % (UNIVERSE - CLAIM_W);
                let e = s + CLAIM_W;
                let _ = tracker.claim_unprocessed(0, s as i64, e as i64);
            }
            let elapsed = t0.elapsed();

            let final_size = tracker.processed[0].lock().intervals().count();
            println!(
                "IntervalSet size={:>8} -> {:>8} | {} ops in {:>9.3} ms | {:>8.2} us/op",
                target_size,
                final_size,
                OPS,
                elapsed.as_secs_f64() * 1000.0,
                elapsed.as_secs_f64() * 1e6 / OPS as f64,
            );
        }
    }

    #[test]
    fn test_is_self_alignment_same_sample_different_contigs() {
        assert!(is_self_alignment(1, 1, 10, 10));
        assert!(is_self_alignment(1, 1, 10, 20));
        assert!(!is_self_alignment(1, 2, 10, 20));
        assert!(!is_self_alignment(1, 2, 10, 10));
    }

    #[test]
    fn test_cigar_op_new_run_splits_long_runs() {
        use crate::impg::{CigarOp, CIGAR_OP_MAX_LEN};

        assert!(CigarOp::new_run(0, '=').is_empty());

        let single = CigarOp::new_run(1_000_000, '=');
        assert_eq!(single.len(), 1);
        assert_eq!(single[0].len() as i64, 1_000_000);

        let at_max = CigarOp::new_run(CIGAR_OP_MAX_LEN, '=');
        assert_eq!(at_max.len(), 1);
        assert_eq!(at_max[0].len() as i64, CIGAR_OP_MAX_LEN);

        let over = CigarOp::new_run(CIGAR_OP_MAX_LEN + 1, '=');
        assert_eq!(over.len(), 2);
        assert_eq!(
            over.iter().map(|op| op.len() as i64).sum::<i64>(),
            CIGAR_OP_MAX_LEN + 1
        );

        let six_hundred_mb: i64 = 600_000_000;
        let run = CigarOp::new_run(six_hundred_mb, '=');
        assert_eq!(
            run.iter().map(|op| op.len() as i64).sum::<i64>(),
            six_hundred_mb
        );
        for op in &run {
            assert!((op.len() as i64) <= CIGAR_OP_MAX_LEN);
        }
    }

    #[test]
    fn test_project_hop0_coords_strand_aware_gap_fallback() {
        let region_start = 0;
        let region_end = 10_000;

        let segments: Vec<HopZeroSeg> = vec![(100, 200, 50, 150), (800, 900, 60, 140)];

        // In the gap between hop-0 segments: drop the hit (empty range so the
        // caller's `if a_start >= a_end { continue }` skips it). The previous
        // whole-region fallback caused raw-BFS over-reports at high transitive
        // depth — see PLAN_depth_100pct.md §3 for the correctness rationale.
        assert_eq!(
            project_hop0_coords(Some(&segments), 500, 600, region_start, region_end),
            (region_start, region_start)
        );

        // Inside first segment: linear projection into [50,150).
        assert_eq!(
            project_hop0_coords(Some(&segments), 120, 180, region_start, region_end),
            (70, 130)
        );

        // Inside second segment: linear projection into [60,140). anc_len=80.
        // 60 + round(0.2*80)=76, 60 + round(0.8*80)=124.
        assert_eq!(
            project_hop0_coords(Some(&segments), 820, 880, region_start, region_end),
            (76, 124)
        );

        // No segments at all: drop the hit. The empty range routes the caller
        // to its "skip" branch instead of crediting the whole anchor span.
        assert_eq!(
            project_hop0_coords(None, 500, 600, region_start, region_end),
            (region_start, region_start)
        );
    }

    // ====================================================================
    // Streaming `pending` state-machine equivalence to merge_short_intervals.
    //
    // The non-windowed Phase 1/2 streaming path replaced the buffered
    // `seq_intervals` Vec + `merge_short_intervals` post-pass with a
    // single-element `pending` slot updated incrementally. These tests assert
    // that the streaming reformulation produces output byte-identical to the
    // original two-pass algorithm for every interval shape that matters.
    // ====================================================================

    fn make_iv(start: i64, end: i64, samples: &[u16]) -> SparseDepthInterval {
        let samples_vec: Vec<SamplePosition> = samples
            .iter()
            .enumerate()
            .map(|(i, &sid)| (sid, i as u32, start, end))
            .collect();
        SparseDepthInterval {
            start,
            end,
            samples: samples_vec,
            pangenome_bases: end - start,
        }
    }

    /// Simulate `StreamingDepthEmitter::merge_into_pending` over a Vec, then
    /// flush. Mirrors the production path exactly so the comparison is honest.
    fn stream_merge(intervals: Vec<SparseDepthInterval>, min_len: i64) -> Vec<SparseDepthInterval> {
        if min_len <= 0 {
            // Streaming-no-merge path: emitter calls emit_final directly.
            return intervals;
        }
        let mut out: Vec<SparseDepthInterval> = Vec::new();
        let mut pending: Option<SparseDepthInterval> = None;
        for new in intervals {
            let new_is_short = (new.end - new.start) < min_len;
            match pending.take() {
                None => pending = Some(new),
                Some(mut p) => {
                    if new_is_short {
                        StreamingDepthEmitter::absorb_right(&mut p, new);
                        pending = Some(p);
                    } else if (p.end - p.start) >= min_len {
                        out.push(p);
                        pending = Some(new);
                    } else {
                        let mut absorbed = new;
                        StreamingDepthEmitter::absorb_left(p, &mut absorbed);
                        pending = Some(absorbed);
                    }
                }
            }
        }
        if let Some(p) = pending {
            out.push(p);
        }
        out
    }

    /// Compare two interval Vecs on the fields the streaming pipeline writes
    /// to output: start, end, sample IDs (set), and pangenome_bases.
    fn assert_intervals_equivalent(a: &[SparseDepthInterval], b: &[SparseDepthInterval]) {
        assert_eq!(a.len(), b.len(), "interval count mismatch");
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert_eq!(x.start, y.start, "start mismatch at index {}", i);
            assert_eq!(x.end, y.end, "end mismatch at index {}", i);
            assert_eq!(
                x.pangenome_bases, y.pangenome_bases,
                "pangenome_bases mismatch at index {}",
                i
            );
            let mut xs: Vec<u16> = x.samples.iter().map(|s| s.0).collect();
            let mut ys: Vec<u16> = y.samples.iter().map(|s| s.0).collect();
            xs.sort_unstable();
            ys.sort_unstable();
            assert_eq!(xs, ys, "sample-id set mismatch at index {}", i);
        }
    }

    // ====================================================================
    // claim_any_unprocessed bool fast-path equivalence to claim_unprocessed.
    //
    // Both methods must end in the identical IntervalSet post-state, and
    // the bool result must equal `!claim_unprocessed(...).is_empty()` for
    // every input. Tested over a deterministic random sequence of claims so
    // we exercise the empty-set, fully-covered, partially-covered, and
    // touching-neighbor branches.
    // ====================================================================
    #[test]
    fn test_claim_any_unprocessed_matches_claim_unprocessed() {
        let mut rng_state: u64 = 0xC1A1_BEEF_FACE_F00D;
        let mut next = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            rng_state
        };

        for _trial in 0..50 {
            let baseline = ConcurrentProcessedTracker::new(1);
            let probe = ConcurrentProcessedTracker::new(1);

            for _ in 0..200 {
                let s = (next() % 100_000) as i64;
                let len = (next() % 5_000 + 1) as i64;
                let e = s + len;

                let baseline_vec = baseline.claim_unprocessed(0, s, e);
                let probe_bool = probe.claim_any_unprocessed(0, s, e);

                assert_eq!(
                    !baseline_vec.is_empty(),
                    probe_bool,
                    "bool divergence at [{},{}): baseline returned {:?}",
                    s, e, baseline_vec
                );

                // Compare full post-state intervals to confirm both trackers
                // ended up storing the same set after the operation.
                let baseline_state: Vec<(i64, i64)> = baseline.processed[0]
                    .lock()
                    .intervals()
                    .collect();
                let probe_state: Vec<(i64, i64)> = probe.processed[0]
                    .lock()
                    .intervals()
                    .collect();
                assert_eq!(
                    baseline_state, probe_state,
                    "post-state divergence after claim [{},{})",
                    s, e
                );
            }
        }
    }

    #[test]
    fn test_stream_merge_min_len_zero_passthrough() {
        let ivs = vec![
            make_iv(0, 100, &[1, 2]),
            make_iv(100, 200, &[3]),
            make_iv(200, 250, &[4]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 0);
        let streamed = stream_merge(ivs, 0);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_leading_shorts_into_first_long() {
        // s1, s2, L3 — pass 2 absorbs leading shorts into the first long.
        let ivs = vec![
            make_iv(0, 10, &[1]),
            make_iv(10, 20, &[2]),
            make_iv(20, 200, &[3]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_short_after_long_absorbs_left() {
        // L1, s2, L3 — pass 1 absorbs s2 into L1.
        let ivs = vec![
            make_iv(0, 200, &[1]),
            make_iv(200, 210, &[2]),
            make_iv(210, 400, &[3]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_all_shorts_kept_as_is() {
        // No long ever arrives — pass-2 boundary is None, intervals collapse
        // into one merged blob in pass 1 (consecutive shorts absorb into prev).
        let ivs = vec![
            make_iv(0, 10, &[1]),
            make_iv(10, 20, &[2]),
            make_iv(20, 30, &[3]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_leading_shorts_overflow_min_len() {
        // s1=80 + s2=80 → pass-1 merges to length 160 ≥ min_len=100.
        // Pass 2 sees first "long" at index 0, no absorption. The streaming
        // reformulation must NOT over-absorb the first long that follows.
        let ivs = vec![
            make_iv(0, 80, &[1]),
            make_iv(80, 160, &[2]),
            make_iv(160, 360, &[3]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
        // Sanity: this case really hit the [merged-shorts, long] split.
        assert_eq!(baseline.len(), 2);
    }

    #[test]
    fn test_stream_merge_alternating_shorts_and_longs() {
        let ivs = vec![
            make_iv(0, 5, &[1]),
            make_iv(5, 200, &[2]),
            make_iv(200, 210, &[3]),
            make_iv(210, 400, &[4]),
            make_iv(400, 405, &[5]),
            make_iv(405, 600, &[6]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_trailing_shorts_absorb_into_last_long() {
        // L1, s2, s3 — pass 1: s2 → L1, then s3 → L1 (the running prev).
        let ivs = vec![
            make_iv(0, 300, &[1]),
            make_iv(300, 310, &[2]),
            make_iv(310, 320, &[3]),
        ];
        let baseline = merge_short_intervals(ivs.clone(), 100);
        let streamed = stream_merge(ivs, 100);
        assert_intervals_equivalent(&baseline, &streamed);
    }

    #[test]
    fn test_stream_merge_single_interval() {
        let short_only = vec![make_iv(0, 10, &[1])];
        assert_intervals_equivalent(
            &merge_short_intervals(short_only.clone(), 100),
            &stream_merge(short_only, 100),
        );

        let long_only = vec![make_iv(0, 200, &[1])];
        assert_intervals_equivalent(
            &merge_short_intervals(long_only.clone(), 100),
            &stream_merge(long_only, 100),
        );
    }

    #[test]
    fn test_stream_merge_fuzz_equivalence() {
        // Deterministic LCG so the test is reproducible without bringing in
        // `rand` as a dev-dep.
        let mut rng_state: u64 = 0xDEAD_BEEF_CAFE_F00D;
        let mut next = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            rng_state
        };

        for _trial in 0..200 {
            let n = (next() % 30 + 1) as usize;
            let mut ivs = Vec::with_capacity(n);
            let mut pos: i64 = 0;
            for _ in 0..n {
                // Mix of short (1-30 bp) and long (50-500 bp) to exercise
                // the leading-shorts overflow case at varying min_len.
                let len = if next() % 2 == 0 {
                    (next() % 30 + 1) as i64
                } else {
                    (next() % 451 + 50) as i64
                };
                let nsamples = (next() % 4 + 1) as u16;
                let samples: Vec<u16> = (0..nsamples).collect();
                ivs.push(make_iv(pos, pos + len, &samples));
                pos += len;
            }

            for &min_len in &[0i64, 1, 50, 100, 250] {
                let baseline = merge_short_intervals(ivs.clone(), min_len);
                let streamed = stream_merge(ivs.clone(), min_len);
                assert_intervals_equivalent(&baseline, &streamed);
            }
        }
    }
}
