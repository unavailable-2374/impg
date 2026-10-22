//! Write-ahead journal that lets `impg partition` survive being killed.
//!
//! Why this exists: partition is a serial greedy loop that keeps every
//! partition in memory and writes the BED only at the very end. On a cluster
//! whose large-memory queues cap at 48 h, a run longer than that is killed
//! with nothing on disk -- and there is no resume.
//!
//! What is journaled: for every window, the overlap set that goes INTO
//! `mask_and_update_regions` (after merge and boundary extension). That
//! function, the final re-merge, the stall guard and the next window selection
//! are all deterministic given this input and the current masks, so replaying
//! the journal through the unchanged loop rebuilds exactly the state -- masks,
//! missing regions, partition numbers, collected partitions -- that the killed
//! run had. Only the expensive part, the transitive query, is skipped.
//!
//! What would break it, and what is done about it:
//! - Resuming with different parameters or a different alignment set would
//!   replay one run's masks into another's. The header carries every
//!   parameter that shapes the partition plus a hash of the sequence index
//!   (names and lengths, in id order, since records store ids) and a caller
//!   tag; any difference is a hard error, never a warning.
//! - Window selection that is not reproducible would silently diverge. Every
//!   record carries its window, and replay demands the loop reselect that exact
//!   window; a mismatch is a hard error.
//! - A kill mid-write leaves a torn tail record. Each record ends in a
//!   checksum; replay stops at the first short or corrupt record, the file is
//!   truncated back to the last good record, and those windows are simply
//!   queried again.
//!
//! Enabled only by IMPG_PARTITION_JOURNAL=<path>. Without it, partition runs
//! the upstream code path unchanged.
//!
//! Coordinate width is per file and fixed by its magic: IMPGPJ01 stores 4-byte
//! coordinates (the original format, still written whenever every sequence
//! fits in i32, which keeps multi-GB journals small), IMPGPJ02 stores 8-byte
//! coordinates and is used only when some sequence is longer than i32::MAX.
//! A resumed journal keeps the width it was created with.

use coitrees::Interval;
use log::info;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::time::{Duration, Instant};

const FILE_MAGIC_NARROW: &[u8; 8] = b"IMPGPJ01";
const FILE_MAGIC_WIDE: &[u8; 8] = b"IMPGPJ02";
const REC_MAGIC_NARROW: u32 = 0x5057_4A31; // "PWJ1"
const REC_MAGIC_WIDE: u32 = 0x5057_4A32; // "PWJ2"
/// Seconds between flush+fsync of the journal: a hard kill loses at most this
/// much work. IMPG_PARTITION_JOURNAL_SYNC_SECS overrides it -- the resume test
/// needs it short, because with kills shorter than the interval nothing is ever
/// flushed and every restart silently begins from zero (which is what the
/// first version of that test measured).
fn sync_every() -> Duration {
    let secs = std::env::var("IMPG_PARTITION_JOURNAL_SYNC_SECS")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(60);
    Duration::from_secs(secs)
}

fn fnv1a(bytes: &[u8], mut h: u64) -> u64 {
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}
const FNV_INIT: u64 = 0xcbf2_9ce4_8422_2325;

/// Hash of the sequence index in id order: records store numeric ids, so an
/// index whose ids map to other sequences must not be replayed into.
pub fn seq_index_fingerprint(names_lens: impl Iterator<Item = (String, usize)>) -> u64 {
    let mut h = FNV_INIT;
    for (name, len) in names_lens {
        h = fnv1a(name.as_bytes(), h);
        h = fnv1a(&[0], h);
        h = fnv1a(&(len as u64).to_le_bytes(), h);
    }
    h
}

pub type Replayed = Vec<Interval<u32>>;

/// On-disk coordinate width of one journal file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Width {
    Narrow,
    Wide,
}

impl Width {
    fn file_magic(self) -> &'static [u8; 8] {
        match self {
            Width::Narrow => FILE_MAGIC_NARROW,
            Width::Wide => FILE_MAGIC_WIDE,
        }
    }
    fn rec_magic(self) -> u32 {
        match self {
            Width::Narrow => REC_MAGIC_NARROW,
            Width::Wide => REC_MAGIC_WIDE,
        }
    }
    fn coord_bytes(self) -> usize {
        match self {
            Width::Narrow => 4,
            Width::Wide => 8,
        }
    }
    /// Record header: magic, seq id, start, end, interval count.
    fn head_bytes(self) -> usize {
        12 + 2 * self.coord_bytes()
    }
    /// Per interval: metadata, first, last.
    fn interval_bytes(self) -> usize {
        4 + 2 * self.coord_bytes()
    }
    fn put(self, buf: &mut Vec<u8>, v: i64) -> io::Result<()> {
        match self {
            Width::Narrow => {
                let v = i32::try_from(v).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("coordinate {v} does not fit a 4-byte partition journal"),
                    )
                })?;
                buf.extend_from_slice(&v.to_le_bytes());
            }
            Width::Wide => buf.extend_from_slice(&v.to_le_bytes()),
        }
        Ok(())
    }
    fn get(self, b: &[u8], o: usize) -> i64 {
        match self {
            Width::Narrow => i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as i64,
            Width::Wide => i64::from_le_bytes(b[o..o + 8].try_into().unwrap()),
        }
    }
}

pub struct Journal {
    path: String,
    width: Width,
    reader: Option<BufReader<File>>,
    good_offset: u64,
    writer: Option<BufWriter<File>>,
    last_sync: Instant,
    sync_every: Duration,
    pub replayed_windows: u64,
    pub appended_windows: u64,
}

impl Journal {
    /// Open (or create) the journal named by IMPG_PARTITION_JOURNAL.
    /// Returns Ok(None) when the variable is unset. `max_seq_len` picks the
    /// coordinate width of a new journal; an existing one keeps its own.
    pub fn from_env(header: &str, max_seq_len: i64) -> io::Result<Option<Journal>> {
        match std::env::var("IMPG_PARTITION_JOURNAL") {
            Ok(p) if !p.is_empty() => Journal::open(p, header, max_seq_len).map(Some),
            _ => Ok(None),
        }
    }

    fn open(path: String, header: &str, max_seq_len: i64) -> io::Result<Journal> {
        let exists = std::fs::metadata(&path).map(|m| m.len() > 0).unwrap_or(false);
        if !exists {
            let width = if max_seq_len > i32::MAX as i64 {
                Width::Wide
            } else {
                Width::Narrow
            };
            let mut f = File::create(&path)?;
            f.write_all(width.file_magic())?;
            f.write_all(&(header.len() as u32).to_le_bytes())?;
            f.write_all(header.as_bytes())?;
            f.sync_all()?;
            info!("partition journal: NEW {path} ({}-byte coordinates)", width.coord_bytes());
            info!("partition journal header: {header}");
            let good = (8 + 4 + header.len()) as u64;
            let writer = OpenOptions::new().append(true).open(&path)?;
            return Ok(Journal {
                path,
                width,
                reader: None,
                good_offset: good,
                writer: Some(BufWriter::with_capacity(8 << 20, writer)),
                last_sync: Instant::now(),
                sync_every: sync_every(),
                replayed_windows: 0,
                appended_windows: 0,
            });
        }
        let mut r = BufReader::with_capacity(8 << 20, File::open(&path)?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        let width = if &magic == FILE_MAGIC_NARROW {
            Width::Narrow
        } else if &magic == FILE_MAGIC_WIDE {
            Width::Wide
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{path}: not a partition journal (magic {magic:?})"),
            ));
        };
        if width == Width::Narrow && max_seq_len > i32::MAX as i64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{path}: journal stores 4-byte coordinates but a sequence is {max_seq_len} bp \
                     long; it cannot be resumed"
                ),
            ));
        }
        let mut len = [0u8; 4];
        r.read_exact(&mut len)?;
        let mut old = vec![0u8; u32::from_le_bytes(len) as usize];
        r.read_exact(&mut old)?;
        let old = String::from_utf8_lossy(&old).into_owned();
        if old != header {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{path}: journal was written by a run with different parameters or inputs; \
                     refusing to replay it.\n  journal: {old}\n  this run: {header}"
                ),
            ));
        }
        info!("partition journal: RESUME {path}");
        info!("partition journal header: {header}");
        let good = (8 + 4 + old.len()) as u64;
        Ok(Journal {
            path,
            width,
            reader: Some(r),
            good_offset: good,
            writer: None,
            last_sync: Instant::now(),
            sync_every: sync_every(),
            replayed_windows: 0,
            appended_windows: 0,
        })
    }

    /// Next journaled window, if replay is still in progress. The loop's own
    /// choice of window is passed in and must match the journal's.
    pub fn take_replay(&mut self, seq_id: u32, start: i64, end: i64) -> io::Result<Option<Replayed>> {
        let Some(r) = self.reader.as_mut() else {
            return Ok(None);
        };
        match read_record(r, self.width)? {
            Some((w, ivs, nbytes)) => {
                if w != (seq_id, start, end) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "{}: replay diverged at record {}: journal has window {:?}, the loop \
                             selected {:?}. Window selection is not reproducible; the journal \
                             cannot be trusted and must not be resumed.",
                            self.path,
                            self.replayed_windows,
                            w,
                            (seq_id, start, end)
                        ),
                    ));
                }
                self.good_offset += nbytes;
                self.replayed_windows += 1;
                Ok(Some(ivs))
            }
            None => {
                self.switch_to_append()?;
                Ok(None)
            }
        }
    }

    fn switch_to_append(&mut self) -> io::Result<()> {
        self.reader = None;
        let f = OpenOptions::new().write(true).open(&self.path)?;
        let len = f.metadata()?.len();
        if len > self.good_offset {
            info!(
                "partition journal: dropping {} bytes of torn tail after record {}",
                len - self.good_offset,
                self.replayed_windows
            );
        }
        f.set_len(self.good_offset)?;
        let mut f = f;
        f.seek(SeekFrom::Start(self.good_offset))?;
        info!(
            "partition journal: replayed {} windows, now querying",
            self.replayed_windows
        );
        self.writer = Some(BufWriter::with_capacity(8 << 20, f));
        Ok(())
    }

    pub fn append(&mut self, seq_id: u32, start: i64, end: i64, ivs: &[Interval<u32>]) -> io::Result<()> {
        if self.reader.is_some() {
            // Replay must be exhausted (take_replay returned None) first.
            return Err(io::Error::new(io::ErrorKind::Other, "journal append during replay"));
        }
        let width = self.width;
        let mut buf =
            Vec::with_capacity(width.head_bytes() + ivs.len() * width.interval_bytes() + 8);
        buf.extend_from_slice(&width.rec_magic().to_le_bytes());
        buf.extend_from_slice(&seq_id.to_le_bytes());
        width.put(&mut buf, start)?;
        width.put(&mut buf, end)?;
        buf.extend_from_slice(&(ivs.len() as u32).to_le_bytes());
        for iv in ivs {
            buf.extend_from_slice(&iv.metadata.to_le_bytes());
            width.put(&mut buf, iv.first)?;
            width.put(&mut buf, iv.last)?;
        }
        let h = fnv1a(&buf, FNV_INIT);
        buf.extend_from_slice(&h.to_le_bytes());
        let w = self.writer.as_mut().expect("journal writer");
        w.write_all(&buf)?;
        self.appended_windows += 1;
        if self.last_sync.elapsed() >= self.sync_every {
            self.sync()?;
        }
        Ok(())
    }

    pub fn sync(&mut self) -> io::Result<()> {
        if let Some(w) = self.writer.as_mut() {
            w.flush()?;
            w.get_ref().sync_data()?;
        }
        self.last_sync = Instant::now();
        Ok(())
    }

    /// Replay must have consumed the whole journal before the loop can end;
    /// otherwise the journal holds windows the loop never reached, which means
    /// replay diverged in a way the per-window check could not see.
    pub fn finish(&mut self) -> io::Result<()> {
        if let Some(r) = self.reader.as_mut() {
            if read_record(r, self.width)?.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}: partitioning finished with journal records left unreplayed",
                        self.path
                    ),
                ));
            }
            self.switch_to_append()?;
        }
        self.sync()?;
        info!(
            "partition journal: {} windows replayed, {} queried this run",
            self.replayed_windows, self.appended_windows
        );
        Ok(())
    }
}

/// One record, or None at EOF / a torn or corrupt record (both end replay).
fn read_record(
    r: &mut BufReader<File>,
    width: Width,
) -> io::Result<Option<((u32, i64, i64), Replayed, u64)>> {
    let (cb, hb, ib) = (width.coord_bytes(), width.head_bytes(), width.interval_bytes());
    let mut head = [0u8; 28];
    let head = &mut head[..hb];
    if !read_full(r, head)? {
        return Ok(None);
    }
    let head = &*head;
    let u32_at = |b: &[u8], o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    if u32_at(head, 0) != width.rec_magic() {
        return Ok(None);
    }
    let n = u32_at(head, 8 + 2 * cb) as usize;
    let mut body = vec![0u8; n * ib];
    if !read_full(r, &mut body)? {
        return Ok(None);
    }
    let mut tail = [0u8; 8];
    if !read_full(r, &mut tail)? {
        return Ok(None);
    }
    let h = fnv1a(&body, fnv1a(head, FNV_INIT));
    if h != u64::from_le_bytes(tail) {
        return Ok(None);
    }
    let ivs = (0..n)
        .map(|k| Interval {
            metadata: u32_at(&body, k * ib),
            first: width.get(&body, k * ib + 4),
            last: width.get(&body, k * ib + 4 + cb),
        })
        .collect();
    Ok(Some((
        (u32_at(head, 4), width.get(head, 8), width.get(head, 8 + cb)),
        ivs,
        (hb + n * ib + 8) as u64,
    )))
}

/// read_exact that reports a clean or torn EOF as `false` instead of an error.
fn read_full(r: &mut BufReader<File>, buf: &mut [u8]) -> io::Result<bool> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..]) {
            Ok(0) => return Ok(false),
            Ok(k) => got += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(metadata: u32, first: i64, last: i64) -> Interval<u32> {
        Interval { metadata, first, last }
    }

    fn signature(ivs: &[Interval<u32>]) -> Vec<(u32, i64, i64)> {
        ivs.iter().map(|x| (x.metadata, x.first, x.last)).collect()
    }

    fn temp_path(name: &str) -> String {
        let dir = std::env::temp_dir().join(format!("impg-journal-{}-{name}", std::process::id()));
        let _ = std::fs::remove_file(&dir);
        dir.to_string_lossy().into_owned()
    }

    /// Write two windows, reopen, and replay them; returns the file's magic.
    fn round_trip(name: &str, max_seq_len: i64, windows: &[(u32, i64, i64, Vec<Interval<u32>>)]) -> [u8; 8] {
        let path = temp_path(name);
        let mut j = Journal::open(path.clone(), "h", max_seq_len).unwrap();
        for (seq, s, e, ivs) in windows {
            j.append(*seq, *s, *e, ivs).unwrap();
        }
        j.finish().unwrap();
        drop(j);

        let mut j = Journal::open(path.clone(), "h", max_seq_len).unwrap();
        for (seq, s, e, ivs) in windows {
            let got = j.take_replay(*seq, *s, *e).unwrap().expect("record");
            assert_eq!(signature(&got), signature(ivs));
        }
        assert!(j.take_replay(0, 0, 0).unwrap().is_none());
        let mut magic = [0u8; 8];
        File::open(&path).unwrap().read_exact(&mut magic).unwrap();
        std::fs::remove_file(&path).unwrap();
        magic
    }

    #[test]
    fn small_genomes_keep_the_four_byte_format() {
        let windows = vec![
            (1, 0, 5000, vec![iv(1, 0, 5000), iv(2, 900, 100)]),
            (2, 5000, i32::MAX as i64, vec![iv(2, i32::MAX as i64 - 1, i32::MAX as i64)]),
        ];
        assert_eq!(&round_trip("narrow", i32::MAX as i64, &windows), FILE_MAGIC_NARROW);
    }

    #[test]
    fn long_sequences_switch_to_eight_byte_coordinates() {
        let big = 3_000_000_000i64;
        let windows = vec![
            (7, big - 20_000, big, vec![iv(7, big - 20_000, big), iv(8, big + 5, big - 5)]),
            (8, 0, 10, vec![iv(8, 0, 10)]),
        ];
        assert_eq!(&round_trip("wide", big + 10, &windows), FILE_MAGIC_WIDE);
    }

    #[test]
    fn four_byte_journal_rejects_out_of_range_coordinates() {
        let path = temp_path("overflow");
        let mut j = Journal::open(path.clone(), "h", 1000).unwrap();
        assert!(j.append(1, 0, i32::MAX as i64 + 1, &[]).is_err());
        drop(j);
        // Reopening it for a run that now has a longer sequence must refuse.
        assert!(Journal::open(path.clone(), "h", i32::MAX as i64 + 1).is_err());
        std::fs::remove_file(&path).unwrap();
    }

    /// Byte layout of a record written by the original i32-only journal.
    #[test]
    fn replays_records_written_by_the_original_format() {
        let path = temp_path("legacy");
        let header = "h";
        let mut bytes = Vec::new();
        bytes.extend_from_slice(FILE_MAGIC_NARROW);
        bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
        bytes.extend_from_slice(header.as_bytes());
        let mut rec = Vec::new();
        rec.extend_from_slice(&REC_MAGIC_NARROW.to_le_bytes());
        rec.extend_from_slice(&3u32.to_le_bytes());
        rec.extend_from_slice(&100i32.to_le_bytes());
        rec.extend_from_slice(&200i32.to_le_bytes());
        rec.extend_from_slice(&1u32.to_le_bytes());
        rec.extend_from_slice(&4u32.to_le_bytes());
        rec.extend_from_slice(&250i32.to_le_bytes());
        rec.extend_from_slice(&150i32.to_le_bytes());
        let h = fnv1a(&rec, FNV_INIT);
        rec.extend_from_slice(&h.to_le_bytes());
        bytes.extend_from_slice(&rec);
        std::fs::write(&path, &bytes).unwrap();

        let mut j = Journal::open(path.clone(), header, 1000).unwrap();
        let got = j.take_replay(3, 100, 200).unwrap().expect("record");
        assert_eq!(signature(&got), vec![(4, 250, 150)]);
        assert!(j.take_replay(3, 200, 300).unwrap().is_none());
        std::fs::remove_file(&path).unwrap();
    }
}
