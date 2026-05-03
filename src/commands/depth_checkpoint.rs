//! Checkpoint / resume support for `impg depth` global mode.
//!
//! Scope: only the global, transitive, CIGAR-precise BFS path
//! (`compute_depth_global` + `--use-BFS`). All other modes are out of scope and
//! the calling code rejects `--resume` for them.
//!
//! Three on-disk artifacts share `<prefix>.depth.` as a base:
//!
//! - `<prefix>.depth.tsv`       — the main TSV output, append-only.
//! - `<prefix>.depth.work.bin`  — append-only binary log of `discovered_regions`
//!                                per chunk; replayed on resume to rebuild
//!                                `tracker` and `global_used`.
//! - `<prefix>.depth.ckpt`      — small bincode blob holding the consistent
//!                                byte offsets into both files plus shared
//!                                counters and an invalidation hash.
//!
//! Commit protocol (per N chunks, see `depth.rs`):
//!
//!   1. Workers send TSV bytes + work-log bytes for each completed chunk to
//!      the writer thread (single-threaded → ordering source-of-truth).
//!   2. Writer flushes both buffered writers and `sync_data`s the underlying
//!      files; the resulting `(tsv_offset, work_offset)` is captured.
//!   3. The ckpt is written to `<prefix>.depth.ckpt.tmp`, `sync_data`d, then
//!      renamed over the live `<prefix>.depth.ckpt` (atomic on POSIX).
//!
//! On startup with `--resume`:
//!
//!   1. Load the ckpt; verify `invalidation_hash`. Mismatch ⇒ refuse.
//!   2. Truncate `*.depth.tsv` to `tsv_byte_offset` and `*.depth.work.bin` to
//!      `work_byte_offset` so any post-commit dirty bytes go away.
//!   3. Stream-read the work log; for each record call back into the trackers
//!      so `tracker == global_used == { committed discovered_regions }`.
//!
//! On clean completion the caller deletes both `*.depth.ckpt` and
//! `*.depth.work.bin` as a success marker.

use std::fs::{self, File, OpenOptions};
use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// Bumped on any breaking layout change to the ckpt blob or work-log format.
pub const CKPT_SCHEMA_VERSION: u32 = 1;

/// Bumped on any breaking layout change to the work-log binary format.
pub const WORKLOG_VERSION: u16 = 1;
/// Identifies a `.depth.work.bin` file. ASCII "IMPW".
pub const WORKLOG_MAGIC: u32 = 0x49_4D_50_57;
/// Per-record sentinel inside the work-log. ASCII "RECD".
pub const WORKLOG_RECORD_MARKER: u32 = 0x52_45_43_44;

/// Filename suffixes; centralised so callers can't drift.
pub const TSV_SUFFIX: &str = ".depth.tsv";
pub const WORKLOG_SUFFIX: &str = ".depth.work.bin";
pub const CKPT_SUFFIX: &str = ".depth.ckpt";
pub const CKPT_TMP_SUFFIX: &str = ".depth.ckpt.tmp";

/// Persistent state captured at every commit barrier.
///
/// The data here is intentionally tiny — a few hundred bytes — because it has
/// to be `sync_data`'d atomically per commit. All bulk state (the union of
/// `discovered_regions` ≡ `tracker` ≡ `global_used`) lives in the work-log,
/// which is reconstructed on resume.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DepthCheckpoint {
    pub schema_version: u32,
    /// Salted hash over alignment file identities and depth configuration.
    /// A mismatch on resume aborts loudly rather than silently producing
    /// garbage by stitching incompatible state together.
    pub invalidation_hash: u64,
    /// Length of `*.depth.tsv` at commit time. Resume truncates the live TSV
    /// to this offset so any post-commit dirty bytes (a half-flushed worker
    /// chunk) are physically discarded.
    pub tsv_byte_offset: u64,
    /// Length of `*.depth.work.bin` at commit time. Resume truncates and
    /// stream-replays from byte 0 up to this offset to rebuild trackers.
    pub work_byte_offset: u64,
    /// Snapshot of `row_counter` (the next `#id` value to allocate). Resuming
    /// with this value keeps post-resume row IDs strictly larger than any
    /// already-flushed row, so the `#id` column stays unique within a run.
    pub row_counter: u64,
    pub intervals_counter: u64,
}

impl DepthCheckpoint {
    pub fn new(invalidation_hash: u64) -> Self {
        Self {
            schema_version: CKPT_SCHEMA_VERSION,
            invalidation_hash,
            tsv_byte_offset: 0,
            work_byte_offset: 0,
            row_counter: 0,
            intervals_counter: 0,
        }
    }

    /// Atomically write the ckpt: `tmp` + `sync_data` + `rename`.
    ///
    /// The directory is _not_ explicitly fsync'd here because the cost on
    /// shared HPC filesystems is non-trivial and the worst-case loss is the
    /// most recent ckpt update — the previous ckpt is still on disk. Callers
    /// that need strict crash-consistency for the *latest* commit should call
    /// [`fsync_dir`] after `save_atomic`.
    pub fn save_atomic(&self, prefix: &str) -> io::Result<()> {
        let final_path = format!("{}{}", prefix, CKPT_SUFFIX);
        let tmp_path = format!("{}{}", prefix, CKPT_TMP_SUFFIX);

        let bytes = bincode::serde::encode_to_vec(self, bincode::config::standard())
            .map_err(|e| io::Error::other(format!("ckpt encode: {}", e)))?;

        {
            let mut f = OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(&tmp_path)?;
            f.write_all(&bytes)?;
            f.sync_data()?;
        }
        fs::rename(&tmp_path, &final_path)?;
        Ok(())
    }

    /// Try to load `<prefix>.depth.ckpt`. Returns `Ok(None)` if the file is
    /// absent (a fresh run), `Err` for any other failure (corrupt bytes,
    /// schema mismatch, IO).
    pub fn try_load(prefix: &str) -> io::Result<Option<Self>> {
        let path = format!("{}{}", prefix, CKPT_SUFFIX);
        let bytes = match fs::read(&path) {
            Ok(b) => b,
            Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(e),
        };
        let (ckpt, _len): (DepthCheckpoint, usize) =
            bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
                .map_err(|e| io::Error::other(format!("ckpt decode: {}", e)))?;
        if ckpt.schema_version != CKPT_SCHEMA_VERSION {
            return Err(io::Error::other(format!(
                "ckpt schema version mismatch: file has {}, binary expects {}",
                ckpt.schema_version, CKPT_SCHEMA_VERSION
            )));
        }
        Ok(Some(ckpt))
    }

    /// Remove ckpt + work-log on successful completion. Both are best-effort:
    /// missing files are silently ignored.
    pub fn cleanup_on_success(prefix: &str) {
        let _ = fs::remove_file(format!("{}{}", prefix, CKPT_SUFFIX));
        let _ = fs::remove_file(format!("{}{}", prefix, CKPT_TMP_SUFFIX));
        let _ = fs::remove_file(format!("{}{}", prefix, WORKLOG_SUFFIX));
    }
}

/// Strict O_DIRECTORY-style fsync. Returns `Ok(())` on platforms where the
/// concept doesn't apply.
pub fn fsync_dir(path: &Path) -> io::Result<()> {
    let dir = match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let f = File::open(&dir)?;
    f.sync_all()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Hash / invalidation
// ---------------------------------------------------------------------------

/// Extra parameters that flow into the invalidation hash beyond the alignment
/// file fingerprints and `DepthConfig`. Everything that influences which
/// chunks Phase 1/Phase 2 produce or what `discovered_regions` end up in the
/// work-log must be hashed; otherwise resuming after a flag flip would mix
/// incompatible state.
#[derive(Debug, Clone, Copy)]
pub struct HashInputs<'a> {
    pub ref_sample: Option<&'a str>,
    pub ref_only: bool,
    pub min_seq_length: i64,
    pub min_interval_len: i64,
    pub window_size: Option<i64>,
    pub merge_adjacent: bool,
    pub stats_mode: bool,
    pub stats_combined: bool,
    pub separator: &'a str,
    pub fai_list: Option<&'a str>,
    pub use_cigar_bfs: bool,
    pub transitive: bool,
    pub transitive_dfs: bool,
    pub max_depth: u16,
    pub min_transitive_len: i64,
    pub min_distance_between_ranges: i64,
}

/// Build a stable u64 fingerprint of "the work the user asked for". Uses the
/// std `DefaultHasher` (Wyhash on stable since 1.66 — adequate for a sentinel
/// check; not used for security).
pub fn compute_invalidation_hash(alignment_files: &[String], inputs: &HashInputs<'_>) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();

    // Salt the hash with the schema version so a binary upgrade that changes
    // the ckpt or work-log layout invalidates older artifacts automatically.
    CKPT_SCHEMA_VERSION.hash(&mut hasher);
    WORKLOG_VERSION.hash(&mut hasher);

    // Alignment files: path + size + mtime (seconds + nanoseconds). mtime is
    // brittle under `cp` without `-p`; if a caller hits a false-positive
    // they can delete the ckpt to force a fresh run.
    let mut sorted_files = alignment_files.to_vec();
    sorted_files.sort();
    (sorted_files.len() as u64).hash(&mut hasher);
    for path in &sorted_files {
        path.hash(&mut hasher);
        if let Ok(meta) = fs::metadata(path) {
            meta.len().hash(&mut hasher);
            if let Ok(mtime) = meta.modified() {
                if let Ok(d) = mtime.duration_since(std::time::UNIX_EPOCH) {
                    d.as_secs().hash(&mut hasher);
                    d.subsec_nanos().hash(&mut hasher);
                }
            }
        } else {
            // File missing or unreadable: hash the path alone so the hash
            // stays deterministic; the resume path will fail loudly when it
            // tries to actually open the file.
            0u64.hash(&mut hasher);
        }
    }

    // Config that influences chunk identity / discovered_regions content.
    inputs.ref_sample.unwrap_or("").hash(&mut hasher);
    inputs.ref_only.hash(&mut hasher);
    inputs.min_seq_length.hash(&mut hasher);
    inputs.min_interval_len.hash(&mut hasher);
    inputs.window_size.unwrap_or(-1).hash(&mut hasher);
    inputs.merge_adjacent.hash(&mut hasher);
    inputs.stats_mode.hash(&mut hasher);
    inputs.stats_combined.hash(&mut hasher);
    inputs.separator.hash(&mut hasher);
    inputs.fai_list.unwrap_or("").hash(&mut hasher);
    inputs.use_cigar_bfs.hash(&mut hasher);
    inputs.transitive.hash(&mut hasher);
    inputs.transitive_dfs.hash(&mut hasher);
    inputs.max_depth.hash(&mut hasher);
    inputs.min_transitive_len.hash(&mut hasher);
    inputs.min_distance_between_ranges.hash(&mut hasher);

    hasher.finish()
}

// ---------------------------------------------------------------------------
// Work-log binary format
// ---------------------------------------------------------------------------
//
// Layout:
//   Header (16 bytes):
//     magic            : u32 LE = WORKLOG_MAGIC
//     version          : u16 LE = WORKLOG_VERSION
//     reserved         : u16 LE = 0
//     reserved2        : u64 LE = 0   (room for future record-count summary)
//
//   Record (variable):
//     record_marker    : u32 LE = WORKLOG_RECORD_MARKER
//     chunk_id         : u64 LE
//     num_regions      : u32 LE
//     payload          : num_regions × (u32 LE seq_id, i64 LE start, i64 LE end)
//
// On resume the entire file is streamed front-to-back; partial trailing bytes
// past the truncation offset cannot exist because the offset captured in the
// ckpt is only ever the post-flush+sync stream position.

pub const WORKLOG_HEADER_LEN: u64 = 16;

/// Encode one work-log record into a fresh `Vec<u8>` ready to be handed to
/// the writer thread. The hot path: called once per chunk by every worker,
/// so it deliberately does no extra allocations beyond the result Vec.
pub fn encode_work_record(chunk_id: u64, regions: &[(u32, i64, i64)]) -> Vec<u8> {
    // Record header: 4 + 8 + 4 = 16 bytes; payload: 20 bytes per region.
    let mut out = Vec::with_capacity(16 + 20 * regions.len());
    out.extend_from_slice(&WORKLOG_RECORD_MARKER.to_le_bytes());
    out.extend_from_slice(&chunk_id.to_le_bytes());
    out.extend_from_slice(&(regions.len() as u32).to_le_bytes());
    for &(seq_id, start, end) in regions {
        out.extend_from_slice(&seq_id.to_le_bytes());
        out.extend_from_slice(&start.to_le_bytes());
        out.extend_from_slice(&end.to_le_bytes());
    }
    out
}

/// Encode the worklog header. Called once when the file is first created.
pub fn encode_worklog_header() -> [u8; 16] {
    let mut buf = [0u8; 16];
    buf[0..4].copy_from_slice(&WORKLOG_MAGIC.to_le_bytes());
    buf[4..6].copy_from_slice(&WORKLOG_VERSION.to_le_bytes());
    // bytes 6..16 stay zero (reserved fields)
    buf
}

/// One decoded work-log record.
#[derive(Debug, Clone)]
pub struct WorkRecord {
    pub chunk_id: u64,
    pub regions: Vec<(u32, i64, i64)>,
}

/// Stream-decode `path`, calling `cb(record)` for each entry up to but not
/// past `expected_end_offset`. Errors:
/// - file missing while ckpt says we should have one ⇒ `InvalidData`
/// - magic / version mismatch ⇒ `InvalidData`
/// - truncated record ⇒ `UnexpectedEof`
///
/// The expected_end_offset is the `work_byte_offset` from the ckpt; bytes
/// past it (already truncated by the resume entry path) must not exist.
pub fn replay_work_log<F>(
    path: &Path,
    expected_end_offset: u64,
    mut cb: F,
) -> io::Result<u64>
where
    F: FnMut(WorkRecord) -> io::Result<()>,
{
    let mut file = match File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            if expected_end_offset == 0 {
                return Ok(0);
            }
            return Err(io::Error::other(format!(
                "work-log {:?} missing but ckpt expects {} bytes",
                path, expected_end_offset
            )));
        }
        Err(e) => return Err(e),
    };

    let actual_len = file.seek(SeekFrom::End(0))?;
    if actual_len < expected_end_offset {
        return Err(io::Error::other(format!(
            "work-log {:?} truncated: ckpt says {} bytes, file has {}",
            path, expected_end_offset, actual_len
        )));
    }
    file.seek(SeekFrom::Start(0))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    if expected_end_offset == 0 {
        return Ok(0);
    }
    if expected_end_offset < WORKLOG_HEADER_LEN {
        return Err(io::Error::other(format!(
            "work-log {:?} ckpt offset {} is shorter than header ({})",
            path, expected_end_offset, WORKLOG_HEADER_LEN
        )));
    }

    // Header
    let mut hdr = [0u8; 16];
    reader.read_exact(&mut hdr)?;
    let magic = u32::from_le_bytes(hdr[0..4].try_into().unwrap());
    let version = u16::from_le_bytes(hdr[4..6].try_into().unwrap());
    if magic != WORKLOG_MAGIC {
        return Err(io::Error::other(format!(
            "work-log {:?} bad magic: got 0x{:08x}, expected 0x{:08x}",
            path, magic, WORKLOG_MAGIC
        )));
    }
    if version != WORKLOG_VERSION {
        return Err(io::Error::other(format!(
            "work-log {:?} version mismatch: got {}, expected {}",
            path, version, WORKLOG_VERSION
        )));
    }

    let mut consumed = WORKLOG_HEADER_LEN;
    while consumed < expected_end_offset {
        let mut head = [0u8; 16]; // marker (4) + chunk_id (8) + num_regions (4)
        reader.read_exact(&mut head)?;
        let marker = u32::from_le_bytes(head[0..4].try_into().unwrap());
        if marker != WORKLOG_RECORD_MARKER {
            return Err(io::Error::other(format!(
                "work-log {:?} bad record marker at offset {}: 0x{:08x}",
                path, consumed, marker
            )));
        }
        let chunk_id = u64::from_le_bytes(head[4..12].try_into().unwrap());
        let num_regions = u32::from_le_bytes(head[12..16].try_into().unwrap()) as usize;

        let mut regions: Vec<(u32, i64, i64)> = Vec::with_capacity(num_regions);
        let mut buf = [0u8; 20];
        for _ in 0..num_regions {
            reader.read_exact(&mut buf)?;
            let seq_id = u32::from_le_bytes(buf[0..4].try_into().unwrap());
            let start = i64::from_le_bytes(buf[4..12].try_into().unwrap());
            let end = i64::from_le_bytes(buf[12..20].try_into().unwrap());
            regions.push((seq_id, start, end));
        }
        consumed += 16 + (20 * num_regions as u64);
        cb(WorkRecord { chunk_id, regions })?;
    }
    if consumed != expected_end_offset {
        return Err(io::Error::other(format!(
            "work-log {:?} ended at {} bytes, expected {} (record framing drift?)",
            path, consumed, expected_end_offset
        )));
    }
    Ok(consumed)
}

// ---------------------------------------------------------------------------
// chunk_id encoding
// ---------------------------------------------------------------------------
//
// `chunk_id` is a u64 identifier for one unit of work. It is used both as
// the worker's "skip this chunk on resume" key and as the work-log record
// header. The encoding must be:
//   - deterministic across resumes (same inputs → same id)
//   - free of cross-phase collisions (Phase 1 id != Phase 2 id != boundary)
//
// Layout:
//   bits 63..62 : phase tag (00 = boundary / sentinel, 01 = Phase 1,
//                 10 = Phase 2, 11 = reserved sentinels)
//   Phase 1 (tag = 01): bits 61..0 = chunk index in the deterministic
//                       phase1_chunks order.
//   Phase 2 (tag = 10): bits 61..36 = seq_id (low 26 bits; supports up to
//                       ~67 M unique sequences), bits 35..0 = chunk_start
//                       in bp (36 bits → up to 64 Gbp per seq, comfortably
//                       above any plausible chromosome).
//   Boundary  (tag = 00): always 0.
//   FAI-done  (tag = 11): u64::MAX. Marks "the post-Phase-2 unaligned-FAI
//                       loop has been committed"; consulted on resume to
//                       skip re-emission of those rows.
//
// Mismatched encodings between resumes are caught upstream by the schema
// version check, so we don't bother adding an extra sentinel here.

pub const CHUNK_ID_BOUNDARY: u64 = 0;
/// Sentinel: emitted once after the post-Phase-2 unaligned-FAI emission
/// loop completes and a barrier+commit has flushed it. Resume consults this
/// sentinel to skip re-emitting FAI rows that are already durable in TSV.
pub const CHUNK_ID_FAI_DONE: u64 = u64::MAX;
const PHASE1_TAG: u64 = 1u64 << 62;
const PHASE2_TAG: u64 = 2u64 << 62;

pub fn encode_chunk_id_phase1(idx: usize) -> u64 {
    PHASE1_TAG | ((idx as u64) & ((1u64 << 62) - 1))
}

pub fn encode_chunk_id_phase2(seq_id: u32, chunk_start: i64) -> u64 {
    let seq = ((seq_id as u64) & ((1u64 << 26) - 1)) << 36;
    let start = (chunk_start as u64) & ((1u64 << 36) - 1);
    PHASE2_TAG | seq | start
}

// ---------------------------------------------------------------------------
// Resume bundle
// ---------------------------------------------------------------------------

/// Aggregate state recovered on `--resume`. Constructed by replaying the
/// work-log; consumed by the depth pipeline before Phase 1 starts.
pub struct ResumeState {
    pub ckpt: DepthCheckpoint,
    /// All `chunk_id` values that have been committed before the most recent
    /// ckpt. Workers consult this set to short-circuit re-processing.
    pub completed_chunks: std::collections::HashSet<u64>,
}

/// Truncate a file to `len` bytes; create+sync if missing and `len == 0`.
/// Used by the resume path to throw away bytes written past the last commit.
///
/// Refuses to extend a file: if the current length is *less* than `len`, that
/// means the ckpt got durably written but the data file did not — almost
/// certainly a broken `sync_data` on the underlying filesystem. Calling
/// `set_len` here would zero-fill the missing bytes and silently corrupt the
/// resumed run. Bail loudly instead.
pub fn truncate_to(path: &Path, len: u64) -> io::Result<()> {
    match OpenOptions::new().write(true).open(path) {
        Ok(f) => {
            let actual = f.metadata()?.len();
            if actual < len {
                return Err(io::Error::other(format!(
                    "{}: ckpt expects at least {} bytes but file is {} bytes — \
                     refusing to extend (would zero-fill and corrupt). Delete \
                     the .ckpt and .work.bin to start fresh.",
                    path.display(),
                    len,
                    actual
                )));
            }
            f.set_len(len)?;
            f.sync_data()?;
            Ok(())
        }
        Err(e) if e.kind() == io::ErrorKind::NotFound && len == 0 => Ok(()),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn ckpt_round_trip() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("foo").to_string_lossy().into_owned();
        let ck = DepthCheckpoint {
            schema_version: CKPT_SCHEMA_VERSION,
            invalidation_hash: 0xDEAD_BEEF_CAFE_F00D,
            tsv_byte_offset: 1024,
            work_byte_offset: 2048,
            row_counter: 17,
            intervals_counter: 99,
        };
        ck.save_atomic(&prefix).unwrap();
        let loaded = DepthCheckpoint::try_load(&prefix).unwrap().unwrap();
        assert_eq!(loaded.tsv_byte_offset, 1024);
        assert_eq!(loaded.work_byte_offset, 2048);
        assert_eq!(loaded.row_counter, 17);
        assert_eq!(loaded.intervals_counter, 99);
        assert_eq!(loaded.invalidation_hash, 0xDEAD_BEEF_CAFE_F00D);
    }

    #[test]
    fn ckpt_missing_returns_none() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("nope").to_string_lossy().into_owned();
        assert!(DepthCheckpoint::try_load(&prefix).unwrap().is_none());
    }

    #[test]
    fn worklog_round_trip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("w.bin");

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_worklog_header());
        bytes.extend(encode_work_record(
            1,
            &[(0, 0, 100), (1, 50, 150)],
        ));
        bytes.extend(encode_work_record(2, &[(2, 200, 300)]));
        bytes.extend(encode_work_record(3, &[]));
        std::fs::write(&path, &bytes).unwrap();

        let mut got: Vec<WorkRecord> = Vec::new();
        let consumed = replay_work_log(&path, bytes.len() as u64, |rec| {
            got.push(rec);
            Ok(())
        })
        .unwrap();
        assert_eq!(consumed, bytes.len() as u64);
        assert_eq!(got.len(), 3);
        assert_eq!(got[0].chunk_id, 1);
        assert_eq!(got[0].regions, vec![(0, 0, 100), (1, 50, 150)]);
        assert_eq!(got[1].chunk_id, 2);
        assert_eq!(got[1].regions, vec![(2, 200, 300)]);
        assert_eq!(got[2].chunk_id, 3);
        assert!(got[2].regions.is_empty());
    }

    #[test]
    fn worklog_partial_replay_stops_at_offset() {
        // Mirror the resume path: the ckpt offset may be smaller than the
        // file. (After resume the file is physically truncated, but the
        // replay code must still respect the offset rather than reading past
        // it.)
        let dir = tempdir().unwrap();
        let path = dir.path().join("w.bin");
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_worklog_header());
        bytes.extend(encode_work_record(1, &[(0, 0, 10)]));
        let stop_at = bytes.len() as u64;
        bytes.extend(encode_work_record(2, &[(0, 10, 20)]));
        std::fs::write(&path, &bytes).unwrap();

        let mut chunk_ids = Vec::new();
        let consumed =
            replay_work_log(&path, stop_at, |rec| {
                chunk_ids.push(rec.chunk_id);
                Ok(())
            })
            .unwrap();
        assert_eq!(consumed, stop_at);
        assert_eq!(chunk_ids, vec![1u64]);
    }

    #[test]
    fn worklog_bad_magic_rejected() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("w.bin");
        let mut bytes = vec![0u8; 16];
        bytes[0..4].copy_from_slice(&0xDEADBEEFu32.to_le_bytes());
        std::fs::write(&path, &bytes).unwrap();
        let err = replay_work_log(&path, bytes.len() as u64, |_| Ok(())).unwrap_err();
        assert!(format!("{err}").contains("bad magic"));
    }

    #[test]
    fn chunk_id_encodings_disjoint() {
        let p1 = encode_chunk_id_phase1(42);
        let p2 = encode_chunk_id_phase2(7, 5_000_000);
        let bd = CHUNK_ID_BOUNDARY;
        let fai = CHUNK_ID_FAI_DONE;
        assert_ne!(p1, p2);
        assert_ne!(p1, bd);
        assert_ne!(p2, bd);
        assert_ne!(p1, fai);
        assert_ne!(p2, fai);
        assert_ne!(bd, fai);
        // Phase tag bits
        assert_eq!(p1 >> 62, 1);
        assert_eq!(p2 >> 62, 2);
        assert_eq!(bd >> 62, 0);
        assert_eq!(fai >> 62, 3);
        // Phase 2 round-trip
        let p2b = encode_chunk_id_phase2(7, 5_000_000);
        assert_eq!(p2, p2b);
        // Different start bp → different id
        assert_ne!(p2, encode_chunk_id_phase2(7, 10_000_000));
        // Different seq → different id
        assert_ne!(p2, encode_chunk_id_phase2(8, 5_000_000));
        // 36-bit chunk_start covers >4Gbp without truncation collisions.
        let big1 = encode_chunk_id_phase2(7, 4_500_000_000);
        let big2 = encode_chunk_id_phase2(7, 4_500_000_000 + 5_000_000);
        assert_ne!(big1, big2);
        assert_ne!(big1, encode_chunk_id_phase2(7, 5_000_000));
    }

    #[test]
    fn invalidation_hash_changes_with_config() {
        let inputs1 = HashInputs {
            ref_sample: None,
            ref_only: false,
            min_seq_length: 0,
            min_interval_len: 0,
            window_size: None,
            merge_adjacent: false,
            stats_mode: false,
            stats_combined: false,
            separator: "#",
            fai_list: None,
            use_cigar_bfs: true,
            transitive: true,
            transitive_dfs: false,
            max_depth: 0,
            min_transitive_len: 0,
            min_distance_between_ranges: 0,
        };
        let mut inputs2 = inputs1.clone();
        inputs2.ref_only = true;
        let h1 = compute_invalidation_hash(&[], &inputs1);
        let h2 = compute_invalidation_hash(&[], &inputs2);
        assert_ne!(h1, h2, "ref_only flip must change hash");
    }
}

