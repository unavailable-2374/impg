//! Persistent normalized CIGARs, not expanded per-diagonal coordinate tuples.
//! Cache entries use local sequence dictionaries, so an unchanged source can
//! be reused when more samples are added to a run.
use super::{invalid, CoordinateMode, DepthAlignment, ImpgIndex, Sequence};
#[cfg(test)]
use crate::impg::CigarOp;
use memmap2::{Mmap, MmapOptions};
use parking_lot::Mutex;
use rustc_hash::FxHashMap as HashMap;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::VecDeque;
use std::fs::{self, File, OpenOptions};
use std::hash::{Hash, Hasher};
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

const BUCKETS: usize = 64;
const MAX_OP: u64 = (1 << 30) - 1;
const HEADER: usize = 80;
const STRIDE: usize = 256;
const CACHE_MAGIC: u64 = 0x3748434143504d49;

fn hash<T: Hash>(value: &T) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}
fn mix(state: u64, value: u64) -> u64 {
    (state.rotate_left(7) ^ value).wrapping_mul(0x9e3779b185ebca87)
}
fn crc_start(row: &[u64]) -> u64 {
    row.iter().fold(0xcbf29ce484222325, |s, &x| mix(s, x))
}
fn u64_at(bytes: &[u8], i: usize) -> u64 {
    u64::from_le_bytes(bytes[i * 8..i * 8 + 8].try_into().unwrap())
}
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
struct Signature {
    version: u32,
    files: Vec<(PathBuf, u64, u128)>,
    selected: Vec<(String, u64)>,
    mode: u8,
}
#[derive(Clone, Serialize, Deserialize)]
pub(super) struct Source {
    bucket: usize,
    signature: Signature,
    pub count: u64,
    offsets: u64,
    names: Vec<(String, u64)>,
}

struct Lock(File);
impl Lock {
    fn acquire(path: &Path) -> io::Result<Self> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(path)?;
        if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX) } != 0 {
            return Err(io::Error::last_os_error());
        }
        Ok(Self(file))
    }
}
impl Drop for Lock {
    fn drop(&mut self) {
        unsafe {
            libc::flock(self.0.as_raw_fd(), libc::LOCK_UN);
        }
    }
}

#[derive(Default)]
struct Journal {
    entries: HashMap<u64, Source>,
    offset: u64,
}
impl Journal {
    fn refresh(&mut self, path: &Path) -> io::Result<()> {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(path)?;
        file.seek(SeekFrom::Start(self.offset))?;
        let mut reader = BufReader::new(file);
        loop {
            let begin = reader.stream_position()?;
            let mut header = [0u8; 24];
            match reader.read_exact(&mut header) {
                Ok(()) => {}
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                    reader.get_mut().set_len(begin)?;
                    break;
                }
                Err(e) => return Err(e),
            }
            let n = u64_at(&header, 1);
            if u64_at(&header, 0) != CACHE_MAGIC || n > 1 << 30 {
                return Err(invalid("invalid normalized CIGAR cache journal"));
            }
            let mut bytes = vec![0; n as usize];
            match reader.read_exact(&mut bytes) {
                Ok(()) => {}
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                    reader.get_mut().set_len(begin)?;
                    break;
                }
                Err(e) => return Err(e),
            }
            if hash(&bytes) != u64_at(&header, 2) {
                return Err(invalid("normalized CIGAR cache journal checksum mismatch"));
            }
            let (entry, used): (Source, usize) =
                bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
                    .map_err(|e| invalid(format!("decode CIGAR cache: {e}")))?;
            if used != bytes.len() {
                return Err(invalid("CIGAR cache journal trailing data"));
            }
            self.entries.insert(hash(&entry.signature), entry);
            self.offset = reader.stream_position()?;
        }
        Ok(())
    }
    fn commit(&mut self, path: &Path, entry: &Source) -> io::Result<()> {
        let bytes = bincode::serde::encode_to_vec(entry, bincode::config::standard())
            .map_err(|e| invalid(e.to_string()))?;
        let mut file = OpenOptions::new().append(true).open(path)?;
        for x in [CACHE_MAGIC, bytes.len() as u64, hash(&bytes)] {
            file.write_all(&x.to_le_bytes())?;
        }
        file.write_all(&bytes)?;
        file.sync_all()?;
        self.offset = file.stream_position()?;
        self.entries.insert(hash(&entry.signature), entry.clone());
        Ok(())
    }
}

struct Builder {
    out: BufWriter<File>,
    position: u64,
    rows: Vec<u64>,
    names: Vec<(String, u64)>,
    ids: HashMap<u32, u64>,
}
impl Builder {
    fn append(
        &mut self,
        a: DepthAlignment,
        ops: &[u32],
        prefix: &[[u64; 3]],
        impg: &impl ImpgIndex,
    ) -> io::Result<()> {
        let mut local = |id: u32| -> io::Result<u64> {
            if let Some(&id) = self.ids.get(&id) {
                return Ok(id);
            }
            let name = impg
                .seq_index()
                .get_name(id)
                .ok_or_else(|| invalid("missing source sequence name"))?;
            let length = impg
                .seq_index()
                .get_len_from_id(id)
                .ok_or_else(|| invalid("missing source sequence length"))?
                as u64;
            let new = self.names.len() as u64;
            self.names.push((name.to_owned(), length));
            self.ids.insert(id, new);
            Ok(new)
        };
        let ta = local(a.target_id)?;
        let qb = local(a.query_id)?;
        let mut row = [
            ta,
            qb,
            a.target_start as u64,
            a.target_end as u64,
            a.query_start as u64,
            a.query_end as u64,
            u64::from(a.is_reverse),
            ops.len() as u64,
            prefix.len() as u64,
            0,
        ];
        row[9] = ops
            .iter()
            .fold(crc_start(&row[..9]), |s, &op| mix(s, op as u64));
        self.rows.push(self.position);
        for word in row {
            self.out.write_all(&word.to_le_bytes())?;
        }
        for &word in ops {
            self.out.write_all(&word.to_le_bytes())?;
        }
        for checkpoint in prefix {
            for word in checkpoint {
                self.out.write_all(&word.to_le_bytes())?;
            }
        }
        self.position += HEADER as u64 + ops.len() as u64 * 4 + prefix.len() as u64 * 24;
        Ok(())
    }
}

pub(super) struct Cache {
    pub dir: PathBuf,
    journals: Vec<Mutex<Journal>>,
}
impl Cache {
    pub fn new(dir: PathBuf) -> io::Result<Self> {
        fs::create_dir_all(&dir)?;
        Ok(Self {
            dir,
            journals: (0..BUCKETS)
                .map(|_| Mutex::new(Journal::default()))
                .collect(),
        })
    }
    pub fn prepare(
        &self,
        impg: &impl ImpgIndex,
        unit: usize,
        mode: CoordinateMode,
        mapping: &[Option<u64>],
        sequences: &[Sequence],
        priorities: &[usize],
        sequence_index: Option<&crate::sequence_index::UnifiedSequenceIndex>,
    ) -> io::Result<Source> {
        let ids = impg.depth_unit_sequence_ids(unit)?;
        let mut selected: Vec<_> = ids
            .iter()
            .filter_map(|&id| mapping.get(id as usize).copied().flatten())
            .map(|id| {
                let s = &sequences[id as usize];
                (s.name.clone(), s.length)
            })
            .collect();
        selected.sort();
        selected.dedup();
        let paths: Vec<_> = if impg.depth_scan_units() == impg.alignment_files().len() {
            vec![&impg.alignment_files()[unit]]
        } else {
            impg.alignment_files().iter().collect()
        };
        let files = paths
            .into_iter()
            .map(|p| {
                let p = fs::canonicalize(p)?;
                let m = fs::metadata(&p)?;
                Ok((
                    p,
                    m.len(),
                    m.modified()?
                        .duration_since(std::time::UNIX_EPOCH)
                        .map_err(|e| invalid(e.to_string()))?
                        .as_nanos(),
                ))
            })
            .collect::<io::Result<Vec<_>>>()?;
        let signature = Signature {
            version: 7,
            files,
            selected,
            mode: u8::from(mode == CoordinateMode::Approx),
        };
        let key = hash(&signature);
        let bucket = key as usize % BUCKETS;
        let journal_path = self.dir.join(format!("{bucket:02}.index"));
        let data_path = self.dir.join(format!("{bucket:02}.data"));
        let _lock = Lock::acquire(&self.dir.join(format!("{bucket:02}.lock")))?;
        let mut journal = self.journals[bucket].lock();
        journal.refresh(&journal_path)?;
        if let Some(source) = journal
            .entries
            .get(&key)
            .filter(|s| s.signature == signature)
        {
            if fs::metadata(&data_path)?.len() < source.offsets + source.count * 8 {
                return Err(invalid("truncated normalized CIGAR cache"));
            }
            return Ok(source.clone());
        }
        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(false)
            .open(data_path)?;
        let position = file.metadata()?.len();
        let mut out = BufWriter::with_capacity(1 << 20, file);
        out.seek(SeekFrom::Start(position))?;
        let builder = Mutex::new(Builder {
            out,
            position,
            rows: Vec::new(),
            names: Vec::new(),
            ids: HashMap::default(),
        });
        impg.scan_depth_unit_parallel(
            unit,
            mode == CoordinateMode::Exact,
            sequence_index,
            priorities,
            &|a| {
                let selected = |id: u32| mapping.get(id as usize).copied().flatten();
                let (Some(t), Some(q)) = (selected(a.target_id), selected(a.query_id)) else {
                    return Ok(());
                };
                if a.target_end < 0
                    || a.query_end < 0
                    || a.target_end as u64 > sequences[t as usize].length
                    || a.query_end as u64 > sequences[q as usize].length
                {
                    return Err(invalid("alignment endpoints exceed sequence length"));
                }
                let (ops, prefix) = normalize(&a, mode)?;
                builder.lock().append(a, &ops, &prefix, impg)
            },
        )?;
        let mut builder = builder.into_inner();
        let offsets = builder.position;
        for i in 0..builder.rows.len() {
            let value = builder.rows[i];
            builder.out.write_all(&value.to_le_bytes())?;
        }
        builder.out.flush()?;
        builder.out.get_ref().sync_all()?;
        let source = Source {
            bucket,
            signature,
            count: builder.rows.len() as u64,
            offsets,
            names: builder.names,
        };
        journal.commit(&journal_path, &source)?;
        Ok(source)
    }
}

fn push_op(ops: &mut Vec<u32>, kind: u32, mut n: u64) {
    while n > 0 {
        if let Some(last) = ops.last_mut().filter(|v| **v & 3 == kind) {
            let old = (*last >> 2) as u64;
            let add = (MAX_OP - old).min(n);
            if add > 0 {
                *last = ((old + add) as u32) << 2 | kind;
                n -= add;
                continue;
            }
        }
        let add = n.min(MAX_OP);
        ops.push((add as u32) << 2 | kind);
        n -= add;
    }
}
fn normalize(a: &DepthAlignment, mode: CoordinateMode) -> io::Result<(Vec<u32>, Vec<[u64; 3]>)> {
    if a.target_start < 0
        || a.query_start < 0
        || a.target_start >= a.target_end
        || a.query_start >= a.query_end
    {
        return Err(invalid("invalid normalized alignment endpoints"));
    }
    let tl = (a.target_end - a.target_start) as u64;
    let ql = (a.query_end - a.query_start) as u64;
    let mut ops = Vec::new();
    if mode == CoordinateMode::Exact {
        if a.cigar.is_empty() {
            return Err(invalid("exact global depth requires CIGAR"));
        }
        for op in &a.cigar {
            let k = match op.op() {
                'M' | '=' | 'X' => 0,
                'I' => 1,
                'D' => 2,
                _ => return Err(invalid("unsupported normalized CIGAR operation")),
            };
            if op.len() <= 0 {
                return Err(invalid("nonpositive CIGAR operation"));
            }
            push_op(&mut ops, k, op.len() as u64);
        }
    } else {
        let short = tl.min(ql);
        let long = tl.max(ql);
        let gap_kind = if tl <= ql { 1 } else { 2 };
        let bias = u128::from(tl > ql && a.is_reverse);
        let mut i = 0;
        let mut consumed = 0;
        while i < short {
            let mapped =
                ((((2 * i as u128 + 1) * long as u128) - bias) / (2 * short as u128)) as u64;
            let gap = mapped - i;
            let next = if long == short {
                short
            } else {
                let num =
                    (2 * short as u128 * (gap as u128 + 1) + bias).saturating_sub(long as u128);
                short
                    .min(num.div_ceil(2 * (long - short) as u128) as u64)
                    .max(i + 1)
            };
            push_op(&mut ops, gap_kind, mapped - consumed);
            push_op(&mut ops, 0, next - i);
            consumed = mapped + next - i;
            i = next;
        }
        push_op(&mut ops, gap_kind, long - consumed);
    }
    let (mut t, mut q) = (0u64, 0u64);
    let mut prefix = Vec::new();
    for (i, &word) in ops.iter().enumerate() {
        if i % STRIDE == 0 {
            prefix.push([i as u64, t, q]);
        }
        let n = (word >> 2) as u64;
        if word & 3 != 1 {
            t = t
                .checked_add(n)
                .ok_or_else(|| invalid("CIGAR target overflow"))?;
        }
        if word & 3 != 2 {
            q = q
                .checked_add(n)
                .ok_or_else(|| invalid("CIGAR query overflow"))?;
        }
    }
    if (t, q) != (tl, ql) {
        return Err(invalid(
            "CIGAR consumption does not match alignment endpoints",
        ));
    }
    Ok((ops, prefix))
}

struct Window {
    map: Mmap,
    base: u64,
}
struct Windows {
    maps: HashMap<(usize, u64), Arc<Window>>,
    lru: VecDeque<(usize, u64)>,
    bytes: usize,
}
thread_local! {
    static LOCAL_WINDOWS: RefCell<Vec<((u64,usize,u64),Arc<Window>)>> = const { RefCell::new(Vec::new()) };
}
static NEXT_STORE: AtomicU64 = AtomicU64::new(1);
pub(super) struct Store {
    token: u64,
    headers: Vec<Mutex<HashMap<u64, Alignment>>>,
    header_limit: usize,
    files: Vec<Option<Arc<File>>>,
    sources: Vec<Source>,
    prefix: Vec<u64>,
    maps: Vec<Vec<Option<u64>>>,
    lengths: Vec<Vec<u64>>,
    windows: Mutex<Windows>,
    window_size: u64,
    budget: usize,
}
#[derive(Clone, Copy, Debug)]
pub(super) struct Alignment {
    pub bucket: usize,
    pub offset: u64,
    pub a: u64,
    pub b: u64,
    pub s: u64,
    pub e: u64,
    pub t: u64,
    pub u: u64,
    pub reverse: bool,
    ops: u64,
    checkpoints: u64,
    checksum: u64,
    crc_prefix: u64,
}
#[derive(Clone, Copy, Debug)]
pub(super) struct Diagonal {
    pub a: u64,
    pub s: u64,
    pub b: u64,
    pub t: u64,
    pub n: u64,
    pub reverse: bool,
}
impl Store {
    pub fn read_catalog(path: &Path) -> io::Result<Vec<Source>> {
        if path
            .file_name()
            .is_some_and(|n| n == "sources.manifest.bin")
        {
            let chunks: Vec<PathBuf> = super::blocks::read_checked(path)?;
            let mut sources = Vec::new();
            for chunk in chunks {
                sources.extend(super::blocks::read_checked::<Vec<Source>>(&chunk)?);
            }
            Ok(sources)
        } else {
            // Earlier schema-7 checkpoints stored one source vector.
            super::blocks::read_checked(path)
        }
    }
    pub fn release_thread_windows(&self) {
        rayon::broadcast(|_| LOCAL_WINDOWS.with(|local| local.borrow_mut().clear()));
    }
    #[cfg(test)]
    pub fn from_test_alignments(
        dir: &Path,
        sequences: &[Sequence],
        alignments: &[DepthAlignment],
        mode: CoordinateMode,
    ) -> io::Result<Self> {
        let mut file = BufWriter::new(File::create(dir.join("00.data"))?);
        let mut positions = Vec::new();
        let mut position = 0u64;
        for a in alignments {
            let (ops, prefix) = normalize(a, mode)?;
            let mut row = [
                a.target_id as u64,
                a.query_id as u64,
                a.target_start as u64,
                a.target_end as u64,
                a.query_start as u64,
                a.query_end as u64,
                a.is_reverse as u64,
                ops.len() as u64,
                prefix.len() as u64,
                0,
            ];
            row[9] = ops
                .iter()
                .fold(crc_start(&row[..9]), |s, &op| mix(s, op as u64));
            positions.push(position);
            for v in row {
                file.write_all(&v.to_le_bytes())?;
            }
            for &v in &ops {
                file.write_all(&v.to_le_bytes())?;
            }
            for p in &prefix {
                for v in p {
                    file.write_all(&v.to_le_bytes())?;
                }
            }
            position += HEADER as u64 + ops.len() as u64 * 4 + prefix.len() as u64 * 24;
        }
        for p in &positions {
            file.write_all(&p.to_le_bytes())?;
        }
        file.flush()?;
        let names: Vec<_> = sequences
            .iter()
            .map(|s| (s.name.clone(), s.length))
            .collect();
        let source = Source {
            bucket: 0,
            signature: Signature {
                version: 7,
                files: Vec::new(),
                selected: names.clone(),
                mode: mode as u8,
            },
            count: positions.len() as u64,
            offsets: position,
            names,
        };
        Self::open(dir, vec![source], sequences, 16 << 20)
    }
    pub fn open(
        dir: &Path,
        mut sources: Vec<Source>,
        sequences: &[Sequence],
        budget: usize,
    ) -> io::Result<Self> {
        let names: HashMap<_, _> = sequences
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.as_str(), (i as u64, s.length)))
            .collect();
        let mut prefix = vec![0];
        let mut maps = Vec::new();
        let mut lengths = Vec::new();
        let mut files: Vec<Option<Arc<File>>> = (0..BUCKETS).map(|_| None).collect();
        for source in &mut sources {
            if files[source.bucket].is_none() {
                files[source.bucket] = Some(Arc::new(File::open(
                    dir.join(format!("{:02}.data", source.bucket)),
                )?));
            }
            prefix.push(prefix.last().copied().unwrap() + source.count);
            maps.push(
                source
                    .names
                    .iter()
                    .map(|(name, length)| match names.get(name.as_str()) {
                        Some(&(id, n)) if n == *length => Ok(Some(id)),
                        Some(_) => Err(invalid(format!(
                            "cached source/FAI length conflict for '{name}'"
                        ))),
                        None => Ok(None),
                    })
                    .collect::<io::Result<_>>()?,
            );
            lengths.push(source.names.iter().map(|(_, length)| *length).collect());
            // Names/signatures remain in the durable cache. Runtime lookup
            // needs only dense local-to-canonical IDs and declared lengths.
            source.names = Vec::new();
            source.signature.selected = Vec::new();
            source.signature.files = Vec::new();
        }
        let workers = rayon::current_num_threads().max(1);
        let cache_budget = (budget - budget / 16).max(workers * 4096 * 8);
        // Four TLS windows, an active CIGAR reader and its replacement can
        // coexist. Leave room for all of them before shared-cache eviction.
        let window_size = ((cache_budget / workers / 8).min(64 << 20).max(4096) & !4095) as u64;
        Ok(Self {
            token: NEXT_STORE.fetch_add(1, Ordering::Relaxed),
            headers: (0..64).map(|_| Mutex::new(HashMap::default())).collect(),
            header_limit: (budget / 32 / 192 / 64).max(1),
            files,
            sources,
            prefix,
            maps,
            lengths,
            windows: Mutex::new(Windows {
                maps: HashMap::default(),
                lru: VecDeque::new(),
                bytes: 0,
            }),
            window_size,
            budget: cache_budget,
        })
    }
    pub fn len(&self) -> u64 {
        *self.prefix.last().unwrap()
    }
    fn window(&self, bucket: usize, offset: u64) -> io::Result<Arc<Window>> {
        let base = offset / self.window_size * self.window_size;
        let key = (self.token, bucket, base);
        if let Some(window) = LOCAL_WINDOWS.with(|local| {
            local
                .borrow()
                .iter()
                .find(|(k, _)| *k == key)
                .map(|(_, w)| w.clone())
        }) {
            return Ok(window);
        }
        let window = self.shared_window(bucket, offset)?;
        LOCAL_WINDOWS.with(|local| {
            let mut local = local.borrow_mut();
            if local.len() >= 4 {
                local.remove(0);
            }
            local.push((key, window.clone()));
        });
        Ok(window)
    }
    fn shared_window(&self, bucket: usize, offset: u64) -> io::Result<Arc<Window>> {
        let base = offset / self.window_size * self.window_size;
        let key = (bucket, base);
        let mut w = self.windows.lock();
        if let Some(map) = w.maps.get(&key) {
            return Ok(map.clone());
        }
        let file = self.files[bucket]
            .as_ref()
            .ok_or_else(|| invalid("missing CIGAR container"))?;
        let total = file.metadata()?.len();
        if offset >= total {
            return Err(invalid("truncated normalized CIGAR record"));
        }
        let len = (total - base).min(self.window_size) as usize;
        let mut attempts = w.lru.len();
        while w.bytes + len > self.budget && attempts > 0 {
            attempts -= 1;
            let Some(old) = w.lru.pop_front() else {
                break;
            };
            if w.maps.get(&old).is_some_and(|p| Arc::strong_count(p) == 1) {
                let old = w.maps.remove(&old).unwrap();
                w.bytes -= old.map.len();
            } else {
                w.lru.push_back(old);
            }
        }
        if w.bytes + len > self.budget {
            return Err(invalid(
                "CIGAR window cache has too many pinned readers for its memory budget",
            ));
        }
        let map = Arc::new(Window {
            map: unsafe { MmapOptions::new().offset(base).len(len).map(&**file)? },
            base,
        });
        w.bytes += len;
        w.maps.insert(key, map.clone());
        w.lru.push_back(key);
        Ok(map)
    }
    fn bytes(&self, bucket: usize, mut offset: u64, mut out: &mut [u8]) -> io::Result<()> {
        while !out.is_empty() {
            let w = self.window(bucket, offset)?;
            let at = (offset - w.base) as usize;
            let n = out.len().min(w.map.len() - at);
            out[..n].copy_from_slice(&w.map[at..at + n]);
            offset += n as u64;
            out = &mut out[n..];
        }
        Ok(())
    }
    fn value(&self, bucket: usize, offset: u64) -> io::Result<u64> {
        let mut bytes = [0; 8];
        self.bytes(bucket, offset, &mut bytes)?;
        Ok(u64::from_le_bytes(bytes))
    }
    pub fn alignment(&self, id: u64) -> io::Result<Alignment> {
        let shard = id as usize % self.headers.len();
        if let Some(a) = self.headers[shard].lock().get(&id).copied() {
            return Ok(a);
        }
        let a = self.read_alignment(id)?;
        let mut headers = self.headers[shard].lock();
        if headers.len() >= self.header_limit {
            headers.clear();
        }
        headers.insert(id, a);
        Ok(a)
    }
    fn read_alignment(&self, id: u64) -> io::Result<Alignment> {
        if id >= self.len() {
            return Err(invalid("invalid normalized alignment ID"));
        }
        let unit = self.prefix.partition_point(|&p| p <= id) - 1;
        let source = &self.sources[unit];
        let offset = self.value(source.bucket, source.offsets + (id - self.prefix[unit]) * 8)?;
        let mut bytes = [0; HEADER];
        self.bytes(source.bucket, offset, &mut bytes)?;
        let get = |i| u64_at(&bytes, i);
        let extent = get(7)
            .checked_mul(4)
            .and_then(|n| get(8).checked_mul(24).and_then(|p| n.checked_add(p)))
            .and_then(|n| offset.checked_add(HEADER as u64)?.checked_add(n));
        if get(2) >= get(3)
            || get(4) >= get(5)
            || get(6) > 1
            || get(7) == 0
            || get(8) != get(7).div_ceil(STRIDE as u64)
            || extent.is_none_or(|end| end > source.offsets)
            || self.lengths[unit]
                .get(get(0) as usize)
                .is_none_or(|n| get(3) > *n)
            || self.lengths[unit]
                .get(get(1) as usize)
                .is_none_or(|n| get(5) > *n)
        {
            return Err(invalid("invalid normalized CIGAR record header"));
        }
        let a = self.maps[unit]
            .get(get(0) as usize)
            .copied()
            .flatten()
            .ok_or_else(|| invalid("filtered cached target ID"))?;
        let b = self.maps[unit]
            .get(get(1) as usize)
            .copied()
            .flatten()
            .ok_or_else(|| invalid("filtered cached query ID"))?;
        Ok(Alignment {
            bucket: source.bucket,
            offset,
            a,
            b,
            s: get(2),
            e: get(3),
            t: get(4),
            u: get(5),
            reverse: get(6) != 0,
            ops: get(7),
            checkpoints: get(8),
            checksum: get(9),
            crc_prefix: crc_start(&(0..9).map(get).collect::<Vec<_>>()),
        })
    }
    pub fn source_ranges(&self) -> Vec<(u64, u64)> {
        let mut out = Vec::new();
        for pair in self.prefix.windows(2) {
            let mut begin = pair[0];
            while begin < pair[1] {
                let end = (begin + 8192).min(pair[1]);
                out.push((begin, end));
                begin = end;
            }
        }
        out
    }
    pub fn diagonals(
        &self,
        a: Alignment,
        side: usize,
        begin: u64,
        end: u64,
        visit: impl FnMut(Diagonal) -> io::Result<()>,
    ) -> io::Result<()> {
        self.walk(a, side, begin, end, false, visit)
    }
    pub fn seed_diagonals(
        &self,
        a: Alignment,
        visit: impl FnMut(Diagonal) -> io::Result<()>,
    ) -> io::Result<()> {
        self.walk(a, 0, a.s, a.e, true, visit)
    }
    fn walk(
        &self,
        a: Alignment,
        side: usize,
        begin: u64,
        end: u64,
        validate: bool,
        mut visit: impl FnMut(Diagonal) -> io::Result<()>,
    ) -> io::Result<()> {
        let desired = if side == 0 {
            begin.saturating_sub(a.s)
        } else if a.reverse {
            a.u.saturating_sub(end)
        } else {
            begin.saturating_sub(a.t)
        };
        let limit = if side == 0 {
            end.saturating_sub(a.s)
        } else if a.reverse {
            a.u.saturating_sub(begin)
        } else {
            end.saturating_sub(a.t)
        };
        let prefix_base = a.offset + HEADER as u64 + a.ops * 4;
        let mut lo = 0;
        let mut hi = a.checkpoints;
        if !validate {
            while lo < hi {
                let mid = lo + (hi - lo) / 2;
                let at = self.value(a.bucket, prefix_base + mid * 24 + 8 * (side as u64 + 1))?;
                if at < desired {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
        }
        let checkpoint = lo.saturating_sub(1);
        let (mut i, mut t, mut q) = if validate {
            (0, 0, 0)
        } else {
            (
                self.value(a.bucket, prefix_base + checkpoint * 24)?,
                self.value(a.bucket, prefix_base + checkpoint * 24 + 8)?,
                self.value(a.bucket, prefix_base + checkpoint * 24 + 16)?,
            )
        };
        let mut pending: Option<Diagonal> = None;
        let mut window: Option<Arc<Window>> = None;
        let mut crc = a.crc_prefix;
        let mut emit = |d: Diagonal| -> io::Result<()> {
            let s = if side == 0 { d.s } else { d.t };
            if s <= end && s + d.n >= begin {
                visit(d)?;
            }
            Ok(())
        };
        while i < a.ops {
            if !validate && (if side == 0 { t } else { q }) > limit {
                break;
            }
            if validate && i % STRIDE as u64 == 0 {
                let base = prefix_base + i / STRIDE as u64 * 24;
                if [
                    self.value(a.bucket, base)?,
                    self.value(a.bucket, base + 8)?,
                    self.value(a.bucket, base + 16)?,
                ] != [i, t, q]
                {
                    return Err(invalid("normalized CIGAR seek index mismatch"));
                }
            }
            let at = a.offset + HEADER as u64 + i * 4;
            if window
                .as_ref()
                .is_none_or(|w| at < w.base || at + 4 > w.base + w.map.len() as u64)
            {
                window = Some(self.window(a.bucket, at)?);
            }
            let w = window.as_ref().unwrap();
            let j = (at - w.base) as usize;
            let op = u32::from_le_bytes(w.map[j..j + 4].try_into().unwrap());
            let n = (op >> 2) as u64;
            if n == 0 || op & 3 == 3 {
                return Err(invalid("invalid normalized CIGAR operation"));
            }
            if validate {
                crc = mix(crc, op as u64);
            }
            if op & 3 == 0 {
                let next = Diagonal {
                    a: a.a,
                    s: a.s + t,
                    b: a.b,
                    t: if a.reverse {
                        a.u.checked_sub(q + n)
                            .ok_or_else(|| invalid("cached CIGAR query underflow"))?
                    } else {
                        a.t + q
                    },
                    n,
                    reverse: a.reverse,
                };
                if let Some(old) = pending.as_mut() {
                    old.n += n;
                    if a.reverse {
                        old.t = next.t;
                    }
                } else {
                    pending = Some(next);
                }
            } else if let Some(d) = pending.take() {
                emit(d)?;
            }
            if op & 3 != 1 {
                t += n;
            }
            if op & 3 != 2 {
                q += n;
            }
            i += 1;
        }
        if let Some(d) = pending {
            emit(d)?;
        }
        if validate && (crc != a.checksum || (t, q) != (a.e - a.s, a.u - a.t)) {
            return Err(invalid(
                "normalized CIGAR cache checksum/consumption mismatch",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_catalog_reads_test_shards_and_legacy_schema_seven_vectors() {
        let dir = tempfile::tempdir().unwrap();
        let names = vec![("TEST_A#0#x".to_owned(), 100)];
        let source = Source {
            bucket: 0,
            signature: Signature {
                version: 7,
                files: Vec::new(),
                selected: names.clone(),
                mode: 0,
            },
            count: 7,
            offsets: 23,
            names: names.clone(),
        };
        let legacy = dir.path().join("TEST_legacy.bin");
        super::super::blocks::write_checked(&legacy, &vec![source.clone()]).unwrap();
        let chunk = dir.path().join("TEST_chunk.bin");
        super::super::blocks::write_checked(&chunk, &vec![source]).unwrap();
        let manifest = dir.path().join("sources.manifest.bin");
        super::super::blocks::write_checked(&manifest, &vec![chunk]).unwrap();
        for p in [legacy, manifest] {
            let entries = Store::read_catalog(&p).unwrap();
            assert_eq!(entries.len(), 1);
            assert_eq!(entries[0].names, names);
            assert_eq!(entries[0].count, 7);
            assert_eq!(entries[0].offsets, 23);
        }
    }
    #[test]
    fn prefix_seek_matches_complete_test_cigar_in_both_axes_and_strands() {
        for reverse in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let mut cigar = Vec::new();
            for i in 0..1200 {
                cigar.push(CigarOp::new(
                    1 + (i % 5),
                    match i % 3 {
                        0 => '=',
                        1 => 'I',
                        _ => 'D',
                    },
                ));
            }
            let tl = cigar.iter().map(|o| o.target_delta() as i64).sum::<i64>();
            let ql = cigar
                .iter()
                .filter(|o| o.op() != 'D')
                .map(|o| o.len() as i64)
                .sum::<i64>();
            let a = DepthAlignment {
                target_id: 0,
                query_id: 1,
                target_start: 1048400,
                target_end: 1048400 + tl,
                query_start: 1048500,
                query_end: 1048500 + ql,
                is_reverse: reverse,
                cigar,
            };
            let sequences = vec![
                Sequence {
                    name: "TEST_A#0#x".into(),
                    sample: "TEST_A".into(),
                    length: 2 << 20,
                },
                Sequence {
                    name: "TEST_B#0#y".into(),
                    sample: "TEST_B".into(),
                    length: 2 << 20,
                },
            ];
            let store =
                Store::from_test_alignments(dir.path(), &sequences, &[a], CoordinateMode::Exact)
                    .unwrap();
            let a = store.alignment(0).unwrap();
            let mut complete = Vec::new();
            store
                .seed_diagonals(a, |d| {
                    complete.push(d);
                    Ok(())
                })
                .unwrap();
            for side in 0..2 {
                let (s, e) = if side == 0 { (a.s, a.e) } else { (a.t, a.u) };
                for begin in (s..=e).step_by(13) {
                    for length in [0, 1, 37, 256] {
                        let end = (begin + length).min(e);
                        let encode = |d: &Diagonal| (d.a, d.s, d.b, d.t, d.n, d.reverse);
                        let expected = complete
                            .iter()
                            .filter(|d| {
                                let s = if side == 0 { d.s } else { d.t };
                                s <= end && s + d.n >= begin
                            })
                            .map(encode)
                            .collect::<Vec<_>>();
                        let mut actual = Vec::new();
                        store
                            .diagonals(a, side, begin, end, |d| {
                                actual.push(encode(&d));
                                Ok(())
                            })
                            .unwrap();
                        assert_eq!(
                            actual, expected,
                            "TEST reverse={reverse} side={side} begin={begin} end={end}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn normalized_test_cigars_preserve_consumption_and_reverse_center_ties() {
        // TEST-only coordinate fixtures.
        for reverse in [false, true] {
            for t in 1..30 {
                for q in 1..30 {
                    let a = DepthAlignment {
                        target_id: 0,
                        query_id: 1,
                        target_start: 0,
                        target_end: t,
                        query_start: 0,
                        query_end: q,
                        is_reverse: reverse,
                        cigar: Vec::new(),
                    };
                    let (ops, _) = normalize(&a, CoordinateMode::Approx).unwrap();
                    let (mut x, mut y) = (0, 0);
                    let mut edges = Vec::new();
                    for word in ops {
                        let n = (word >> 2) as i64;
                        match word & 3 {
                            0 => {
                                for i in 0..n {
                                    edges
                                        .push((x + i, if reverse { q - 1 - y - i } else { y + i }));
                                }
                                x += n;
                                y += n;
                            }
                            1 => y += n,
                            2 => x += n,
                            _ => unreachable!(),
                        }
                    }
                    assert_eq!((x, y), (t, q));
                    let mut expected = Vec::new();
                    let short = t.min(q);
                    let long = t.max(q);
                    for i in 0..short {
                        let p = ((2 * i + 1) * long) / (2 * short);
                        expected.push(if t <= q {
                            (i, if reverse { q - 1 - p } else { p })
                        } else {
                            (if reverse { t - 1 - p } else { p }, i)
                        });
                    }
                    edges.sort();
                    expected.sort();
                    assert_eq!(edges, expected, "TEST lengths={t},{q},reverse={reverse}");
                }
            }
        }
    }
}
