//! Immutable coordinate pages and sparse manifests. Epochs append only changed
//! blocks; committing a new manifest never rewrites all discovered boundaries.
use super::bitmap::{Bits, SHIFT, WIDTH};
use super::relations::Diagonal;
use super::{invalid, Sequence};
use crate::commands::depth_storage::{RecordSet, Records, RowAccess, Sorter};
use parking_lot::Mutex;
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::{BTreeMap, VecDeque};
use std::fs::{self, File, OpenOptions};
use std::hash::{Hash, Hasher};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

pub(super) type Key = (u64, u64);
pub(super) fn key(seq: u64, pos: u64) -> Key {
    (seq, pos >> SHIFT)
}
pub(super) fn offset(pos: u64) -> u32 {
    (pos & (WIDTH - 1)) as u32
}
fn digest(bytes: &[u8]) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bytes.hash(&mut h);
    h.finish()
}
fn read_at(file: &File, mut offset: u64, mut out: &mut [u8]) -> io::Result<()> {
    while !out.is_empty() {
        let n = file.read_at(out, offset)?;
        if n == 0 {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        offset += n as u64;
        out = &mut out[n..];
    }
    Ok(())
}
fn write_at(file: &File, mut offset: u64, mut bytes: &[u8]) -> io::Result<()> {
    while !bytes.is_empty() {
        let n = file.write_at(bytes, offset)?;
        if n == 0 {
            return Err(io::Error::from(io::ErrorKind::WriteZero));
        }
        offset += n as u64;
        bytes = &bytes[n..];
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub(super) struct Page {
    pub bucket: usize,
    pub at: u64,
    pub len: usize,
    pub checksum: u64,
}
struct Container {
    file: File,
    next: AtomicU64,
}
struct Cached {
    values: HashMap<Page, Arc<Bits>>,
    clock: VecDeque<Page>,
    bytes: usize,
}
thread_local! { static LOCAL_BITS:RefCell<(HashMap<(u64,Page),Arc<Bits>>,VecDeque<(u64,Page)>)> = RefCell::new((HashMap::default(),VecDeque::new())); }
static NEXT_PAGES: AtomicU64 = AtomicU64::new(1);
pub(super) struct Pages {
    token: u64,
    local_limit: usize,
    files: Vec<Container>,
    cache: Mutex<Cached>,
    budget: usize,
}
impl Pages {
    pub fn open(dir: &Path, budget: usize) -> io::Result<Self> {
        fs::create_dir_all(dir)?;
        let mut files = Vec::new();
        for i in 0..16 {
            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .create(true)
                .truncate(false)
                .open(dir.join(format!("pages-{i:02}.bin")))?;
            let n = file.metadata()?.len();
            files.push(Container {
                file,
                next: AtomicU64::new(n),
            });
        }
        Ok(Self {
            token: NEXT_PAGES.fetch_add(1, Ordering::Relaxed),
            local_limit: (budget
                / rayon::current_num_threads().max(1)
                / (WIDTH as usize / 8 + 8192)
                / 2)
            .clamp(2, 1024),
            files,
            cache: Mutex::new(Cached {
                values: HashMap::default(),
                clock: VecDeque::new(),
                bytes: 0,
            }),
            budget: budget.max(1 << 20),
        })
    }
    pub fn put(&self, k: Key, bits: &Bits) -> io::Result<Page> {
        let mut bytes = Vec::with_capacity(bits.bytes().min(WIDTH as usize / 8 + 16));
        bits.encode(&mut bytes)?;
        self.put_raw(k, &bytes)
    }
    pub fn put_raw(&self, k: Key, bytes: &[u8]) -> io::Result<Page> {
        let bucket = ((k.0.wrapping_mul(0x9e3779b185ebca87) ^ k.1) as usize) % self.files.len();
        let file = &self.files[bucket];
        let at = file.next.fetch_add(bytes.len() as u64, Ordering::Relaxed);
        write_at(&file.file, at, bytes)?;
        Ok(Page {
            bucket,
            at,
            len: bytes.len(),
            checksum: digest(bytes),
        })
    }
    pub fn read_raw(&self, page: Page) -> io::Result<Vec<u8>> {
        let mut bytes = vec![0; page.len];
        read_at(
            &self
                .files
                .get(page.bucket)
                .ok_or_else(|| invalid("invalid page container"))?
                .file,
            page.at,
            &mut bytes,
        )?;
        if digest(&bytes) != page.checksum {
            return Err(invalid("page checksum mismatch"));
        }
        Ok(bytes)
    }
    pub fn get(&self, page: Page) -> io::Result<Arc<Bits>> {
        if let Some(bits) =
            LOCAL_BITS.with(|cache| cache.borrow().0.get(&(self.token, page)).cloned())
        {
            return Ok(bits);
        }
        let bits = self.get_shared(page)?;
        LOCAL_BITS.with(|cache| {
            let mut cache = cache.borrow_mut();
            while cache.0.len() >= self.local_limit {
                if let Some(old) = cache.1.pop_front() {
                    cache.0.remove(&old);
                } else {
                    break;
                }
            }
            cache.0.insert((self.token, page), bits.clone());
            cache.1.push_back((self.token, page));
        });
        Ok(bits)
    }
    fn get_shared(&self, page: Page) -> io::Result<Arc<Bits>> {
        if let Some(bits) = self.cache.lock().values.get(&page).cloned() {
            return Ok(bits);
        }
        let mut bytes = vec![0; page.len];
        read_at(
            &self
                .files
                .get(page.bucket)
                .ok_or_else(|| invalid("invalid boundary container"))?
                .file,
            page.at,
            &mut bytes,
        )?;
        if digest(&bytes) != page.checksum {
            return Err(invalid("boundary page checksum mismatch"));
        }
        let bits = Arc::new(Bits::decode(&mut &bytes[..])?);
        let size = bits.bytes();
        let mut cache = self.cache.lock();
        if let Some(old) = cache.values.get(&page) {
            return Ok(old.clone());
        }
        let mut tries = cache.clock.len();
        while cache.bytes + size > self.budget && tries > 0 {
            tries -= 1;
            if let Some(k) = cache.clock.pop_front() {
                if cache
                    .values
                    .get(&k)
                    .is_some_and(|v| Arc::strong_count(v) == 1)
                {
                    let b = cache.values.remove(&k).unwrap();
                    cache.bytes -= b.bytes();
                } else {
                    cache.clock.push_back(k);
                }
            }
        }
        // One page per active worker can be pinned. The allowance is bounded,
        // rather than switching to a global scan when the cache is full.
        if cache.bytes + size
            > self.budget + rayon::current_num_threads() * (WIDTH as usize / 8 + 8192)
        {
            return Err(invalid("too many pinned boundary pages"));
        }
        cache.bytes += size;
        cache.values.insert(page, bits.clone());
        cache.clock.push_back(page);
        Ok(bits)
    }
    pub fn sync(&self) -> io::Result<()> {
        self.files.par_iter().try_for_each(|c| c.file.sync_data())
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub(super) struct Change {
    pub key: Key,
    pub full: Page,
    pub delta: Page,
}
pub(super) type Snapshot = BTreeMap<Key, Page>;
pub(super) fn write_checked<T: Serialize + ?Sized>(path: &Path, value: &T) -> io::Result<()> {
    let bytes = bincode::serde::encode_to_vec(value, bincode::config::standard())
        .map_err(|e| invalid(e.to_string()))?;
    let temporary = path.with_extension("tmp");
    let mut out = BufWriter::new(File::create(&temporary)?);
    out.write_all(b"IMPDAT07")?;
    out.write_all(&(bytes.len() as u64).to_le_bytes())?;
    out.write_all(&digest(&bytes).to_le_bytes())?;
    out.write_all(&bytes)?;
    out.flush()?;
    out.get_ref().sync_all()?;
    fs::rename(temporary, path)
}
pub(super) fn read_checked<T: DeserializeOwned>(path: &Path) -> io::Result<T> {
    let file = File::open(path)?;
    let total = file.metadata()?.len();
    let mut reader = BufReader::new(file);
    let mut header = [0; 24];
    reader.read_exact(&mut header)?;
    let len = u64::from_le_bytes(header[8..16].try_into().unwrap());
    let sum = u64::from_le_bytes(header[16..24].try_into().unwrap());
    if &header[..8] != b"IMPDAT07" || len.checked_add(24) != Some(total) || len > 8u64 << 30 {
        return Err(invalid("invalid depth metadata frame"));
    }
    let mut bytes = vec![0; len as usize];
    reader.read_exact(&mut bytes)?;
    if digest(&bytes) != sum {
        return Err(invalid("depth metadata checksum mismatch"));
    }
    let (value, used) = bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
        .map_err(|e| invalid(e.to_string()))?;
    if used != bytes.len() {
        return Err(invalid("depth metadata trailing bytes"));
    }
    Ok(value)
}
pub(super) fn load_changes(path: &Path) -> io::Result<Vec<Change>> {
    read_checked(path)
}
pub(super) fn save_changes(dir: &Path, epoch: usize, changes: &[Change]) -> io::Result<PathBuf> {
    let path = dir.join(format!("epoch-{epoch:08}.bin"));
    write_checked(&path, changes)?;
    Ok(path)
}

/// Producer-local bitmaps compress duplicate candidates before they reach disk.
/// At most 8 MiB of point offsets are accumulated by each producer batch.
pub(super) struct Candidates<'a> {
    pages: &'a Pages,
    sort: Sorter<7>,
    pending: HashMap<Key, Vec<u32>>,
    bytes: usize,
}
impl<'a> Candidates<'a> {
    pub fn new(pages: &'a Pages, dir: &Path, budget: usize) -> Self {
        Self {
            pages,
            sort: Sorter::new(dir, budget.max(4096)),
            pending: HashMap::default(),
            bytes: 0,
        }
    }
    pub fn push(&mut self, seq: u64, pos: u64) -> io::Result<()> {
        self.pending
            .entry(key(seq, pos))
            .or_default()
            .push(offset(pos));
        self.bytes += 4;
        if self.bytes >= 8 << 20 {
            self.flush()?;
        }
        Ok(())
    }
    pub fn push_bits(&mut self, k: Key, bits: &Bits) -> io::Result<()> {
        if bits.len() == 0 {
            return Ok(());
        }
        let p = self.pages.put(k, bits)?;
        self.sort
            .push([k.0, k.1, p.bucket as u64, p.at, p.len as u64, p.checksum, 0])
    }
    fn flush(&mut self) -> io::Result<()> {
        let pending = std::mem::take(&mut self.pending);
        for (k, v) in pending {
            self.push_bits(k, &Bits::from_sorted(v))?;
        }
        self.bytes = 0;
        Ok(())
    }
    pub fn finish(mut self) -> io::Result<Records<7>> {
        self.flush()?;
        self.sort.finish()
    }
    pub fn finish_nonempty(mut self) -> io::Result<Option<RecordSet<7>>> {
        if self.pending.is_empty() && self.sort.is_empty() {
            Ok(None)
        } else {
            self.flush()?;
            self.sort.finish_adaptive().map(Some)
        }
    }
}
fn record_page(r: [u64; 7]) -> Page {
    Page {
        bucket: r[2] as usize,
        at: r[3],
        len: r[4] as usize,
        checksum: r[5],
    }
}

/// Parallel reductions split at key boundaries. No point-sized global sorting
/// or full-boundary copy is needed when a tiny frontier remains.
pub(super) fn reduce(
    pages: &Pages,
    records: &impl RowAccess<7>,
    old: &Snapshot,
) -> io::Result<Vec<Change>> {
    let step = (records.len() / rayon::current_num_threads().max(1) / 4).max(1024);
    let mut ranges = Vec::new();
    let mut begin = 0;
    while begin < records.len() {
        let mut end = (begin + step).min(records.len());
        if end < records.len() {
            let r = records.get(end - 1);
            end =
                records.lower_bound([r[0], r[1], u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        }
        ranges.push((begin, end));
        begin = end;
    }
    let groups: Vec<Vec<Change>> = ranges
        .into_par_iter()
        .map(|(mut i, end)| -> io::Result<Vec<Change>> {
            let mut out = Vec::new();
            while i < end {
                let r = records.get(i);
                let k = (r[0], r[1]);
                let before = if let Some(p) = old.get(&k) {
                    pages.get(*p)?
                } else {
                    Arc::new(Bits::default())
                };
                let mut combined = Bits::default();
                while i < end {
                    let r = records.get(i);
                    if (r[0], r[1]) != k {
                        break;
                    }
                    let proposed = pages.get(record_page(r))?;
                    combined = combined.union(&proposed);
                    i += 1;
                }
                let delta = combined.difference(&before);
                if delta.len() > 0 {
                    let full = before.union(&delta);
                    let full_page = pages.put(k, &full)?;
                    let delta_page = if before.len() == 0 {
                        full_page
                    } else {
                        pages.put(k, &delta)?
                    };
                    out.push(Change {
                        key: k,
                        full: full_page,
                        delta: delta_page,
                    });
                }
            }
            Ok(out)
        })
        .collect::<io::Result<_>>()?;
    let mut out: Vec<_> = groups.into_iter().flatten().collect();
    out.sort_unstable_by_key(|c| c.key);
    Ok(out)
}

/// A native alignment has one spatial bin per endpoint, even when it spans
/// many leaf blocks. Querying a leaf visits that bin and all its ancestors.
pub(super) fn bin(start: u64, end: u64) -> u64 {
    let (mut shift, mut level) = (SHIFT, 0u64);
    while shift < 64 && start >> shift != end >> shift {
        shift = (shift + 3).min(64);
        level += 1;
    }
    (level << 56) | if shift == 64 { 0 } else { start >> shift }
}
pub(super) fn ancestor_bins(tile: u64) -> impl Iterator<Item = u64> {
    (0..=15u64).map(move |level| {
        let shift = level * 3;
        (level << 56) | if shift >= 64 { 0 } else { tile >> shift }
    })
}
pub(super) fn incident(
    index: &Records<4>,
    k: Key,
    mut f: impl FnMut(u64, usize) -> io::Result<()>,
) -> io::Result<()> {
    for bin in ancestor_bins(k.1) {
        let mut i = index.lower_bound([k.0, bin, 0, 0]);
        while i < index.len {
            let r = index.get(i);
            if r[0] != k.0 || r[1] != bin {
                break;
            }
            f(r[2], r[3] as usize)?;
            i += 1;
        }
    }
    Ok(())
}

pub(super) struct AtomBlock {
    pub key: Key,
    pub page: Page,
    pub base: u64,
    pub starts: u64,
}
pub(super) struct Atoms {
    pub blocks: Vec<AtomBlock>,
    lookup: HashMap<Key, usize>,
    pub count: u64,
    lengths: Vec<u64>,
    pages: Arc<Pages>,
}
impl Atoms {
    pub fn new(snapshot: &Snapshot, pages: Arc<Pages>, sequences: &[Sequence]) -> io::Result<Self> {
        let mut blocks = Vec::with_capacity(snapshot.len());
        let mut lookup = HashMap::default();
        let mut count = 0;
        let entries: Vec<_> = snapshot.iter().map(|(&k, &p)| (k, p)).collect();
        let counts = entries
            .par_iter()
            .map(|&(k, page)| -> io::Result<_> {
                let bits = pages.get(page)?;
                let length = sequences[k.0 as usize].length;
                Ok(bits.len() as u64
                    - u64::from(key(k.0, length) == k && bits.contains(offset(length))))
            })
            .collect::<io::Result<Vec<_>>>()?;
        // Only the small prefix of block counts is serial; loading and
        // validating all boundary pages remains parallel on a cold resume.
        for ((k, page), starts) in entries.into_iter().zip(counts) {
            lookup.insert(k, blocks.len());
            blocks.push(AtomBlock {
                key: k,
                page,
                base: count,
                starts,
            });
            count += starts;
        }
        Ok(Self {
            blocks,
            lookup,
            count,
            lengths: sequences.iter().map(|s| s.length).collect(),
            pages,
        })
    }
    pub fn bits(&self, block: usize) -> io::Result<Arc<Bits>> {
        self.pages.get(self.blocks[block].page)
    }
    pub fn lookup(&self, seq: u64, pos: u64) -> io::Result<u64> {
        if pos >= self.lengths[seq as usize] {
            return Err(invalid("terminal boundary cannot start an atom"));
        }
        let idx = *self
            .lookup
            .get(&key(seq, pos))
            .ok_or_else(|| invalid("missing propagated atom boundary block"))?;
        let b = &self.blocks[idx];
        let bits = self.bits(idx)?;
        if !bits.contains(offset(pos)) {
            return Err(invalid("missing propagated atom boundary"));
        }
        Ok(b.base + bits.rank(offset(pos)) as u64)
    }
    pub fn successor(&self, seq: u64, pos: u64) -> io::Result<u64> {
        let i = *self
            .lookup
            .get(&key(seq, pos))
            .ok_or_else(|| invalid("missing atom block"))?;
        let b = &self.blocks[i];
        let bits = self.bits(i)?;
        if let Some(p) = bits.successor(offset(pos)) {
            return Ok((b.key.1 << SHIFT) + p as u64);
        }
        let next = self
            .blocks
            .get(i + 1)
            .filter(|n| n.key.0 == seq)
            .ok_or_else(|| invalid("missing sequence terminal boundary"))?;
        Ok((next.key.1 << SHIFT) + self.bits(i + 1)?.select(0) as u64)
    }
    pub fn get(&self, id: u64) -> io::Result<(u64, u64, u64)> {
        if id >= self.count {
            return Err(invalid("invalid initial atom ID"));
        }
        let i = self.blocks.partition_point(|b| b.base <= id) - 1;
        let b = &self.blocks[i];
        let bits = self.bits(i)?;
        let rank = (id - b.base) as usize;
        let p = bits.select(rank);
        let start = (b.key.1 << SHIFT) + p as u64;
        let end = if rank + 1 < bits.len() {
            (b.key.1 << SHIFT) + bits.select(rank + 1) as u64
        } else {
            self.successor(b.key.0, start)?
        };
        Ok((b.key.0, start, end))
    }
}

pub(super) fn mapped(d: Diagonal, side: usize, pos: u64) -> (u64, u64) {
    let (src, dst, other) = if side == 0 {
        (d.s, d.b, d.t)
    } else {
        (d.t, d.a, d.s)
    };
    (
        dst,
        other
            + if d.reverse {
                d.n - (pos - src)
            } else {
                pos - src
            },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn closed_test_endpoints_are_found_across_spatial_bins() {
        for (s, e) in [
            (0, WIDTH),
            (WIDTH - 1, WIDTH),
            (17, WIDTH * 90),
            (0, i64::MAX as u64),
        ] {
            let b = bin(s, e);
            for p in [s, e, s + (e - s) / 2] {
                assert!(
                    ancestor_bins(p >> SHIFT).any(|a| a == b),
                    "TEST {s}-{e} at {p}"
                );
            }
        }
    }
}
