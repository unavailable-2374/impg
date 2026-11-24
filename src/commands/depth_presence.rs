//! Direct-alignment presence layer for `impg depth`.
//!
//! The global layer answers, for every base of every genome, which other samples
//! align to it directly. It never builds a transitive closure: each sample is
//! processed against the alignments it takes part in, so memory scales with one
//! genome's alignments rather than with the whole collection.
//!
//! `build` reads every PAF once, splits each alignment at indels of at least
//! `gap` bases into blocks and appends each block to one shard per sample (both
//! sides of the alignment). A sweep then writes, per sample X:
//!
//! * `presence/<X>.bin`: `{seq, start, end, other}` as little-endian u32, merged per
//!   `(seq, other)` and sorted by `seq, other, start`. The sample set of any piece of
//!   X can be rebuilt from this table, so nothing below is lossy.
//! * `depth/<X>.bin`: `{seq, start, end, depth}`, adjacent equal depth merged. Depth
//!   counts X itself, so a base no other sample aligns to has depth 1.
//! * `stats/<X>.tsv`: `STAT`, `HIST` (bases per depth), `COV` (bases of X covered by
//!   each other sample) and `SEQ` (depth sum per sequence) rows.
//!
//! Coordinates are 0-based half-open. The store also holds `seqs.tsv` and
//! `samples.txt`, which `query` uses to read it back.
//!
//! A sample with a stats file is complete, and a finished scan leaves `scan.done`,
//! so an interrupted build continues where it stopped.

use flate2::read::MultiGzDecoder;
use log::{info, warn};
use rustc_hash::FxHashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::Instant;

#[cfg(not(target_endian = "little"))]
compile_error!("depth_presence stores records as little-endian u32 and assumes a little-endian host");

/// One interval record. On disk it is four little-endian u32 values.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Rec {
    pub seq: u32,
    pub s: u32,
    pub e: u32,
    pub other: u32,
}

const REC: usize = std::mem::size_of::<Rec>();
const SHARD_BUFFER: usize = 16384;

fn as_bytes(r: &[Rec]) -> &[u8] {
    // Rec is repr(C) with four u32 fields and no padding.
    unsafe { std::slice::from_raw_parts(r.as_ptr() as *const u8, std::mem::size_of_val(r)) }
}

fn read_recs(path: &Path) -> io::Result<Vec<Rec>> {
    let bytes = fs::metadata(path)?.len() as usize;
    if bytes % REC != 0 {
        return Err(io::Error::other(format!("{} is not a whole number of records", path.display())));
    }
    let mut v = vec![Rec::default(); bytes / REC];
    let buf = unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, bytes) };
    File::open(path)?.read_exact(buf)?;
    Ok(v)
}

/// Sequences and samples of a store. Samples are sorted by name; sequences are
/// grouped by sample in that order and keep their FAI order within a sample.
pub struct Table {
    pub name: Vec<String>,
    pub len: Vec<u64>,
    pub own: Vec<u32>,
    pub id: FxHashMap<String, u32>,
    pub samples: Vec<String>,
}

impl Table {
    /// Build from a file listing FAI paths. The sample is the part of a sequence
    /// name before the first `separator`.
    pub fn from_fai_list(list: &str, separator: &str) -> io::Result<Table> {
        let mut per_sample: FxHashMap<String, Vec<(String, u64)>> = FxHashMap::default();
        let mut seen: FxHashMap<String, ()> = FxHashMap::default();
        for path in BufReader::new(File::open(list)?).lines() {
            let path = path?;
            let path = path.trim();
            if path.is_empty() || path.starts_with('#') {
                continue;
            }
            for line in BufReader::new(File::open(path)?).lines() {
                let line = line?;
                let mut f = line.split('\t');
                let (Some(name), Some(len)) = (f.next(), f.next()) else {
                    return Err(io::Error::other(format!("bad FAI line in {path}: {line}")));
                };
                let len: u64 = len
                    .parse()
                    .map_err(|_| io::Error::other(format!("bad length in {path}: {line}")))?;
                if len > u32::MAX as u64 {
                    return Err(io::Error::other(format!("{name} is longer than 2^32-1 bases")));
                }
                if seen.insert(name.to_string(), ()).is_some() {
                    return Err(io::Error::other(format!("sequence {name} appears twice in the FAI files")));
                }
                let sample = name.split(separator).next().unwrap_or(name).to_string();
                per_sample.entry(sample).or_default().push((name.to_string(), len));
            }
        }
        let mut samples: Vec<String> = per_sample.keys().cloned().collect();
        samples.sort();
        let mut t = Table { name: Vec::new(), len: Vec::new(), own: Vec::new(), id: FxHashMap::default(), samples };
        for (si, s) in t.samples.iter().enumerate() {
            for (n, l) in &per_sample[s] {
                t.id.insert(n.clone(), t.name.len() as u32);
                t.name.push(n.clone());
                t.len.push(*l);
                t.own.push(si as u32);
            }
        }
        if t.samples.is_empty() {
            return Err(io::Error::other(format!("no sequences found through {list}")));
        }
        Ok(t)
    }

    /// Read `seqs.tsv` and `samples.txt` from a store.
    pub fn read(dir: &Path) -> io::Result<Table> {
        let mut t = Table { name: Vec::new(), len: Vec::new(), own: Vec::new(), id: FxHashMap::default(), samples: Vec::new() };
        for line in BufReader::new(File::open(dir.join("samples.txt"))?).lines() {
            let line = line?;
            if !line.is_empty() {
                t.samples.push(line);
            }
        }
        for line in BufReader::new(File::open(dir.join("seqs.tsv"))?).lines() {
            let line = line?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 4 || f[0].parse::<usize>().ok() != Some(t.name.len()) {
                return Err(io::Error::other(format!("bad seqs.tsv line: {line}")));
            }
            t.id.insert(f[1].to_string(), t.name.len() as u32);
            t.name.push(f[1].to_string());
            t.len.push(f[2].parse().map_err(io::Error::other)?);
            t.own.push(f[3].parse().map_err(io::Error::other)?);
        }
        Ok(t)
    }

    fn write(&self, dir: &Path) -> io::Result<()> {
        let mut w = BufWriter::new(File::create(dir.join("seqs.tsv"))?);
        for i in 0..self.name.len() {
            writeln!(w, "{}\t{}\t{}\t{}", i, self.name[i], self.len[i], self.own[i])?;
        }
        w.flush()?;
        let mut w = BufWriter::new(File::create(dir.join("samples.txt"))?);
        for s in &self.samples {
            writeln!(w, "{s}")?;
        }
        w.flush()
    }
}

#[derive(Clone, Debug, Default)]
struct FileResult {
    status: &'static str,
    lines: u64,
    self_lines: u64,
    blocks: u64,
    qbases: u64,
    tbases: u64,
}

struct Shards {
    files: Vec<Mutex<File>>,
    failed: AtomicBool,
}

impl Shards {
    fn write(&self, sample: usize, recs: &[Rec]) {
        let mut f = self.files[sample].lock().unwrap_or_else(|e| e.into_inner());
        if f.write_all(as_bytes(recs)).is_err() {
            self.failed.store(true, Ordering::SeqCst);
        }
    }
}

/// Per-thread buffers, one per sample, flushed when full and on drop.
struct ShardBuffers<'a> {
    shards: &'a Shards,
    buf: Vec<Vec<Rec>>,
}

impl ShardBuffers<'_> {
    fn put(&mut self, sample: u32, r: Rec) {
        let b = &mut self.buf[sample as usize];
        if b.capacity() == 0 {
            b.reserve(SHARD_BUFFER);
        }
        b.push(r);
        if b.len() >= SHARD_BUFFER {
            self.shards.write(sample as usize, b);
            b.clear();
        }
    }
}

impl Drop for ShardBuffers<'_> {
    fn drop(&mut self) {
        for (i, b) in self.buf.iter_mut().enumerate() {
            if !b.is_empty() {
                self.shards.write(i, b);
                b.clear();
            }
        }
    }
}

fn find_cigar(rest: &[u8]) -> Option<usize> {
    let pat = b"cg:Z:";
    (0..rest.len().saturating_sub(pat.len() - 1)).find(|&i| &rest[i..i + pat.len()] == pat).map(|i| i + pat.len())
}

fn parse_u64(b: &[u8]) -> Option<u64> {
    if b.is_empty() {
        return None;
    }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u64)?;
    }
    Some(v)
}

/// Split one PAF record into blocks and hand both sides to the shard buffers.
fn scan_line(line: &[u8], t: &Table, gap: u64, out: &mut ShardBuffers, r: &mut FileResult) -> Result<(), &'static str> {
    let mut f: [&[u8]; 9] = [&[]; 9];
    let mut p = 0usize;
    for slot in f.iter_mut() {
        let Some(tab) = line[p..].iter().position(|&c| c == b'\t') else {
            return Err("bad_line");
        };
        *slot = &line[p..p + tab];
        p += tab + 1;
    }
    let rest = &line[p..];
    let qn = std::str::from_utf8(f[0]).map_err(|_| "bad_line")?;
    let tn = std::str::from_utf8(f[5]).map_err(|_| "bad_line")?;
    let (Some(&qid), Some(&tid)) = (t.id.get(qn), t.id.get(tn)) else {
        return Err("unknown_seq");
    };
    let (qsample, tsample) = (t.own[qid as usize], t.own[tid as usize]);
    if qsample == tsample {
        r.self_lines += 1;
        return Ok(());
    }
    let (Some(qs), Some(qe), Some(ts), Some(te)) = (parse_u64(f[2]), parse_u64(f[3]), parse_u64(f[7]), parse_u64(f[8])) else {
        return Err("bad_line");
    };
    let rev = f[4].first() == Some(&b'-');
    let Some(c0) = find_cigar(rest) else {
        return Err("no_cigar");
    };
    let cg = &rest[c0..];
    let (mut on, mut tb, mut qb, mut tn_, mut qn_) = (false, 0u64, 0u64, 0u64, 0u64);
    let (mut tp, mut qp, mut v) = (0u64, 0u64, 0u64);
    let emit = |on: &mut bool, tb: u64, qb: u64, tn_: u64, qn_: u64, out: &mut ShardBuffers, r: &mut FileResult| {
        if !*on {
            return;
        }
        out.put(tsample, Rec { seq: tid, s: (ts + tb) as u32, e: (ts + tn_) as u32, other: qsample });
        let (lo, hi) = if rev { (qe - qn_, qe - qb) } else { (qs + qb, qs + qn_) };
        out.put(qsample, Rec { seq: qid, s: lo as u32, e: hi as u32, other: tsample });
        r.blocks += 1;
        r.tbases += tn_ - tb;
        r.qbases += qn_ - qb;
        *on = false;
    };
    for &h in cg {
        if h.is_ascii_digit() {
            v = v * 10 + (h - b'0') as u64;
            continue;
        }
        match h {
            b'=' | b'M' | b'X' => {
                if !on {
                    on = true;
                    tb = tp;
                    qb = qp;
                }
                tn_ = tp + v;
                qn_ = qp + v;
                tp += v;
                qp += v;
            }
            b'I' | b'D' => {
                if v >= gap {
                    emit(&mut on, tb, qb, tn_, qn_, out, r);
                }
                if h == b'I' {
                    qp += v;
                } else {
                    tp += v;
                }
            }
            _ => break,
        }
        v = 0;
    }
    emit(&mut on, tb, qb, tn_, qn_, out, r);
    if tp != te.wrapping_sub(ts) || qp != qe.wrapping_sub(qs) {
        return Err("cigar_mismatch");
    }
    Ok(())
}

fn scan_file(path: &str, t: &Table, gap: u64, out: &mut ShardBuffers) -> FileResult {
    let mut r = FileResult { status: "ok", ..Default::default() };
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => {
            r.status = "open_failed";
            return r;
        }
    };
    let mut rd = BufReader::with_capacity(1 << 22, MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file)));
    let mut line = Vec::with_capacity(1 << 16);
    loop {
        line.clear();
        match rd.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {
                if line.last() != Some(&b'\n') {
                    r.status = "no_final_newline";
                    break;
                }
                r.lines += 1;
                if let Err(s) = scan_line(&line[..line.len() - 1], t, gap, out, &mut r) {
                    r.status = s;
                }
            }
            Err(_) => {
                r.status = "gz_error";
                break;
            }
        }
    }
    r
}

fn scan(files: &[String], t: &Table, gap: u64, out: &Path) -> io::Result<()> {
    let rec = out.join("rec");
    fs::create_dir_all(&rec)?;
    let _ = fs::remove_file(out.join("scan.done"));
    let s = t.samples.len();
    let mut handles = Vec::with_capacity(s);
    for i in 0..s {
        handles.push(Mutex::new(
            OpenOptions::new().create(true).truncate(true).write(true).open(rec.join(format!("{i}.bin")))?,
        ));
    }
    let shards = Shards { files: handles, failed: AtomicBool::new(false) };
    let results: Vec<Mutex<FileResult>> = files.iter().map(|_| Mutex::new(FileResult::default())).collect();
    let done = AtomicUsize::new(0);
    let start = Instant::now();
    let next = AtomicUsize::new(0);
    let threads = rayon::current_num_threads().max(1);
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| {
                let mut bufs = ShardBuffers { shards: &shards, buf: vec![Vec::new(); s] };
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= files.len() {
                        break;
                    }
                    let r = scan_file(&files[i], t, gap, &mut bufs);
                    *results[i].lock().unwrap() = r;
                    let k = done.fetch_add(1, Ordering::Relaxed) + 1;
                    if k % 5000 == 0 {
                        info!("presence scan {}/{} files, {:.0}s", k, files.len(), start.elapsed().as_secs_f64());
                    }
                }
            });
        }
    });
    if shards.failed.load(Ordering::SeqCst) {
        return Err(io::Error::other("writing a sample shard failed"));
    }
    for f in &shards.files {
        f.lock().unwrap().sync_all()?;
    }
    let mut bad = 0usize;
    let mut self_lines = 0u64;
    let mut m = BufWriter::new(File::create(out.join("manifest.tsv"))?);
    writeln!(m, "path\tstatus\tlines\tblocks\tqbases\ttbases\tself_lines")?;
    for (p, r) in files.iter().zip(&results) {
        let r = r.lock().unwrap();
        writeln!(m, "{}\t{}\t{}\t{}\t{}\t{}\t{}", p, r.status, r.lines, r.blocks, r.qbases, r.tbases, r.self_lines)?;
        if r.status != "ok" {
            bad += 1;
        }
        self_lines += r.self_lines;
    }
    m.flush()?;
    info!("presence scan finished: {} files, {} with problems, {:.0}s", files.len(), bad, start.elapsed().as_secs_f64());
    if self_lines > 0 {
        warn!("{self_lines} alignment records join two sequences of the same sample and were not counted");
    }
    if bad > 0 {
        return Err(io::Error::other(format!(
            "{bad} alignment files could not be read completely; see {}",
            out.join("manifest.tsv").display()
        )));
    }
    fs::write(out.join("scan.done"), "ok\n")
}

struct SampleResult {
    stats: String,
}

/// Merge, segment and write one sample. Returns the stats text.
fn sweep_sample(x: usize, t: &Table, seqs_of: &[u32], out: &Path) -> io::Result<SampleResult> {
    let s = t.samples.len();
    let name = &t.samples[x];
    let mut rs = read_recs(&out.join("rec").join(format!("{x}.bin")))?;
    rs.retain(|r| r.e > r.s);
    let kept = rs.len() as u64;
    rs.sort_unstable_by(|a, b| (a.seq, a.other, a.s).cmp(&(b.seq, b.other, b.s)));
    let mut cov = vec![0u64; s];
    {
        let mut w = 0usize;
        let mut i = 0usize;
        while i < rs.len() {
            let mut cur = rs[i];
            let mut j = i + 1;
            while j < rs.len() && rs[j].seq == cur.seq && rs[j].other == cur.other {
                if rs[j].s <= cur.e {
                    if rs[j].e > cur.e {
                        cur.e = rs[j].e;
                    }
                } else {
                    cov[cur.other as usize] += (cur.e - cur.s) as u64;
                    rs[w] = cur;
                    w += 1;
                    cur = rs[j];
                }
                j += 1;
            }
            cov[cur.other as usize] += (cur.e - cur.s) as u64;
            rs[w] = cur;
            w += 1;
            i = j;
        }
        rs.truncate(w);
    }
    let iv = rs.len() as u64;
    let ptmp = out.join("presence").join(format!("{name}.bin.tmp"));
    {
        let mut f = BufWriter::with_capacity(1 << 22, File::create(&ptmp)?);
        f.write_all(as_bytes(&rs))?;
        f.flush()?;
    }
    // events: (seq, pos, delta)
    let mut ev: Vec<(u32, u32, i32)> = Vec::with_capacity(rs.len() * 2);
    for r in &rs {
        ev.push((r.seq, r.s, 1));
        ev.push((r.seq, r.e, -1));
    }
    drop(rs);
    ev.sort_unstable_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    let dtmp = out.join("depth").join(format!("{name}.bin.tmp"));
    let mut dw = BufWriter::with_capacity(1 << 22, File::create(&dtmp)?);
    let (mut segs, mut bases, mut dsum, mut toggles, mut overflow) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut lh = vec![0u64; 40];
    let mut lb = vec![0u64; 40];
    let mut hist = vec![0u64; s + 2];
    let mut seq_text = String::new();
    let mut j = 0usize;
    for &sq in seqs_of {
        let ln = t.len[sq as usize];
        let mut ss = 0u64;
        let mut sqd = 0u64;
        let mut depth: i64 = 1;
        let mut run = Rec { seq: sq, s: 0, e: 0, other: 0 };
        let mut seg = |s0: u64, e0: u64, d: i64, run: &mut Rec, dw: &mut BufWriter<File>| -> io::Result<()> {
            if e0 <= s0 {
                return Ok(());
            }
            let len = e0 - s0;
            segs += 1;
            bases += len;
            dsum += len * d as u64;
            sqd += len * d as u64;
            let b = 64 - len.leading_zeros() as usize;
            lh[b] += 1;
            lb[b] += len;
            hist[d as usize] += len;
            if run.e == s0 as u32 && run.other == d as u32 && run.e > run.s {
                run.e = e0 as u32;
            } else {
                if run.e > run.s {
                    dw.write_all(as_bytes(std::slice::from_ref(run)))?;
                }
                *run = Rec { seq: sq, s: s0 as u32, e: e0 as u32, other: d as u32 };
            }
            Ok(())
        };
        while j < ev.len() && ev[j].0 == sq {
            let p = ev[j].1;
            let mut nd = depth;
            while j < ev.len() && ev[j].0 == sq && ev[j].1 == p {
                nd += ev[j].2 as i64;
                toggles += 1;
                j += 1;
            }
            if p as u64 > ln {
                overflow += 1;
            }
            seg(ss, p as u64, depth, &mut run, &mut dw)?;
            ss = p as u64;
            depth = nd;
        }
        seg(ss, ln, depth, &mut run, &mut dw)?;
        if run.e > run.s {
            dw.write_all(as_bytes(std::slice::from_ref(&run)))?;
        }
        seq_text.push_str(&format!("SEQ\t{}\t{}\t{}\t{}\n", name, sq, ln, sqd));
    }
    dw.flush()?;
    drop(dw);
    if j != ev.len() {
        return Err(io::Error::other(format!("{name}: a record names a sequence this sample does not own")));
    }
    let mut st = format!("STAT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t", name, x, kept, iv, segs, bases, dsum, toggles, overflow);
    for v in &lh[1..34] {
        st.push_str(&format!("{v},"));
    }
    st.push('\t');
    for v in &lb[1..34] {
        st.push_str(&format!("{v},"));
    }
    st.push('\n');
    for (d, v) in hist.iter().enumerate() {
        if *v > 0 {
            st.push_str(&format!("HIST\t{}\t{}\t{}\n", name, d, v));
        }
    }
    for (o, v) in cov.iter().enumerate() {
        if o != x {
            st.push_str(&format!("COV\t{}\t{}\t{}\n", name, t.samples[o], v));
        }
    }
    st.push_str(&seq_text);
    Ok(SampleResult { stats: st })
}

fn sweep(t: &Table, out: &Path, mem_budget: u64) -> io::Result<()> {
    if !out.join("scan.done").exists() {
        return Err(io::Error::other("scan.done is missing; the scan did not finish"));
    }
    for d in ["presence", "depth", "stats"] {
        fs::create_dir_all(out.join(d))?;
    }
    let s = t.samples.len();
    let mut seqs_of: Vec<Vec<u32>> = vec![Vec::new(); s];
    for (i, &o) in t.own.iter().enumerate() {
        seqs_of[o as usize].push(i as u32);
    }
    let mut order: Vec<(u64, usize)> = Vec::new();
    for x in 0..s {
        if out.join("stats").join(format!("{}.tsv", t.samples[x])).exists() {
            continue;
        }
        order.push((fs::metadata(out.join("rec").join(format!("{x}.bin")))?.len(), x));
    }
    order.sort_unstable_by(|a, b| b.cmp(a));
    info!("presence sweep: {} samples to do, {} already complete", order.len(), s - order.len());
    let used = Mutex::new(0u64);
    let cv = Condvar::new();
    let next = AtomicUsize::new(0);
    let err: Mutex<Option<io::Error>> = Mutex::new(None);
    let start = Instant::now();
    let done = AtomicUsize::new(0);
    std::thread::scope(|sc| {
        for _ in 0..rayon::current_num_threads().max(1) {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= order.len() || err.lock().unwrap().is_some() {
                    break;
                }
                let (bytes, x) = order[i];
                // records, events (1.5x the record bytes) and write buffers
                let need = bytes * 5 / 2 + (1 << 24);
                {
                    let mut u = used.lock().unwrap();
                    while !(*u == 0 || *u + need <= mem_budget) {
                        u = cv.wait(u).unwrap();
                    }
                    *u += need;
                }
                let r = sweep_sample(x, t, &seqs_of[x], out).and_then(|r| {
                    let name = &t.samples[x];
                    let stmp = out.join("stats").join(format!("{name}.tsv.tmp"));
                    fs::write(&stmp, r.stats)?;
                    fs::rename(out.join("presence").join(format!("{name}.bin.tmp")), out.join("presence").join(format!("{name}.bin")))?;
                    fs::rename(out.join("depth").join(format!("{name}.bin.tmp")), out.join("depth").join(format!("{name}.bin")))?;
                    fs::rename(stmp, out.join("stats").join(format!("{name}.tsv")))
                });
                {
                    let mut u = used.lock().unwrap();
                    *u -= need;
                }
                cv.notify_all();
                if let Err(e) = r {
                    *err.lock().unwrap() = Some(io::Error::other(format!("sample {}: {e}", t.samples[x])));
                    break;
                }
                let k = done.fetch_add(1, Ordering::Relaxed) + 1;
                info!("presence sweep {}/{} {} ({:.0}s)", k, order.len(), t.samples[x], start.elapsed().as_secs_f64());
            });
        }
    });
    if let Some(e) = err.into_inner().unwrap() {
        return Err(e);
    }
    Ok(())
}

/// Concatenate the per-sample stats into one file per row type.
fn summarize(t: &Table, out: &Path) -> io::Result<()> {
    let mut w: Vec<BufWriter<File>> = Vec::new();
    let kinds = ["STAT", "HIST", "COV", "SEQ"];
    for k in kinds {
        w.push(BufWriter::new(File::create(out.join(format!("all.{k}.tsv")))?));
    }
    for s in &t.samples {
        for line in BufReader::new(File::open(out.join("stats").join(format!("{s}.tsv")))?).lines() {
            let line = line?;
            if let Some(i) = kinds.iter().position(|k| line.starts_with(k) && line.as_bytes().get(k.len()) == Some(&b'\t')) {
                writeln!(w[i], "{line}")?;
            }
        }
    }
    for mut x in w {
        x.flush()?;
    }
    Ok(())
}

/// Build the global presence store in `out`.
pub fn build(files: &[String], fai_list: &str, separator: &str, gap: u64, mem_budget: u64, keep_shards: bool, out: &Path) -> io::Result<()> {
    fs::create_dir_all(out)?;
    let t = Table::from_fai_list(fai_list, separator)?;
    info!("presence: {} samples, {} sequences, {} alignment files, gap {}", t.samples.len(), t.name.len(), files.len(), gap);
    if out.join("seqs.tsv").exists() {
        let old = Table::read(out)?;
        if old.name != t.name || old.len != t.len || old.samples != t.samples {
            return Err(io::Error::other(format!(
                "{} already holds a store built from different sequences; use a new directory",
                out.display()
            )));
        }
    } else {
        t.write(out)?;
    }
    let settings = format!("gap\t{gap}\nfiles\t{}\n", files.len());
    let settings_path = out.join("settings.tsv");
    if out.join("scan.done").exists() && fs::read_to_string(&settings_path).ok().as_deref() == Some(settings.as_str()) {
        info!("presence: scan.done found with the same settings, skipping the scan");
    } else {
        // a new scan invalidates every sample result
        for d in ["presence", "depth", "stats"] {
            let _ = fs::remove_dir_all(out.join(d));
        }
        fs::write(&settings_path, &settings)?;
        scan(files, &t, gap, out)?;
    }
    sweep(&t, out, mem_budget)?;
    summarize(&t, out)?;
    if !keep_shards {
        fs::remove_dir_all(out.join("rec"))?;
    }
    info!("presence store complete: {}", out.display());
    Ok(())
}

/// Output of [`query`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueryMode {
    /// One row per other sample: covered bases and fraction of the region.
    Samples,
    /// Pieces of the region where the sample set changes, with present and absent samples.
    Segments,
}

/// Parse `name:start-end` (0-based half-open).
pub fn parse_region(s: &str) -> io::Result<(String, u64, u64)> {
    let colon = s.rfind(':').ok_or_else(|| io::Error::other(format!("bad region {s}")))?;
    let dash = s[colon..].rfind('-').map(|d| d + colon).ok_or_else(|| io::Error::other(format!("bad region {s}")))?;
    let a = s[colon + 1..dash].parse().map_err(|_| io::Error::other(format!("bad region {s}")))?;
    let b = s[dash + 1..].parse().map_err(|_| io::Error::other(format!("bad region {s}")))?;
    Ok((s[..colon].to_string(), a, b))
}

/// Call `f(start, end)` for every record of `(seq, other)` overlapping `[qs, qe)`, clipped to it.
/// `rec` is sorted by `(seq, other, start)` and disjoint within `(seq, other)`.
fn for_overlaps(rec: &[Rec], seq: u32, other: u32, qs: u32, qe: u32, mut f: impl FnMut(u32, u32)) {
    let mut i = rec.partition_point(|r| (r.seq, r.other, r.s) < (seq, other, qs));
    if i > 0 && rec[i - 1].seq == seq && rec[i - 1].other == other && rec[i - 1].e > qs {
        i -= 1;
    }
    while i < rec.len() && rec[i].seq == seq && rec[i].other == other && rec[i].s < qe {
        let (a, b) = (rec[i].s.max(qs), rec[i].e.min(qe));
        if b > a {
            f(a, b);
        }
        i += 1;
    }
}

/// Load the hop-only records of sample `x` and check that they were built on the same table.
fn read_hop(dir: &Path, t: &Table, x: u32) -> io::Result<Vec<Rec>> {
    let path = dir.join(format!("{}.bin", t.samples[x as usize]));
    let rec = read_recs(&path).map_err(|e| io::Error::other(format!("{}: {e}", path.display())))?;
    let ok = rec.iter().all(|r| {
        r.s < r.e
            && t.own.get(r.seq as usize) == Some(&x)
            && r.e as u64 <= t.len[r.seq as usize]
            && (r.other as usize) < t.samples.len()
            && r.other != x
    });
    let sorted = rec
        .windows(2)
        .all(|p| (p[0].seq, p[0].other, p[0].e) <= (p[1].seq, p[1].other, p[1].s));
    if !ok || !sorted {
        return Err(io::Error::other(format!(
            "{} does not match the presence store (different sequence table, or not sorted)",
            path.display()
        )));
    }
    Ok(rec)
}

/// Answer region queries from a presence store.
///
/// `hop` is an optional directory of two-hop supplementary files (`<X>.bin`, same record
/// format, bases without direct presence only, built on the same sequence table). It adds
/// columns and never changes the direct ones, except that a sample present only through
/// the hop layer is listed in `hop_only` instead of `absent`.
pub fn query(
    dir: &Path,
    hop: Option<&Path>,
    regions: &[(String, u64, u64)],
    mode: QueryMode,
    w: &mut impl Write,
) -> io::Result<()> {
    let t = Table::read(dir)?;
    let s = t.samples.len() as u32;
    match mode {
        QueryMode::Samples => write!(w, "#region\tsample\tcovered_bp\tfraction")?,
        QueryMode::Segments => write!(w, "#region\tstart\tend\tdepth\tpresent\tabsent")?,
    }
    match (hop, mode) {
        (None, _) => writeln!(w)?,
        (Some(_), QueryMode::Samples) => writeln!(w, "\thop_only_bp\thop_only_fraction")?,
        (Some(_), QueryMode::Segments) => writeln!(w, "\thop_only_count\thop_only")?,
    }
    let mut hop_rec: Vec<Rec> = Vec::new();
    let mut loaded = u32::MAX;
    let mut rec: Vec<Rec> = Vec::new();
    for (name, qs0, qe0) in regions {
        let &sq = t.id.get(name).ok_or_else(|| io::Error::other(format!("unknown sequence {name}")))?;
        let x = t.own[sq as usize];
        let ln = t.len[sq as usize];
        let (qs, mut qe) = (*qs0, *qe0);
        if qe > ln {
            warn!("{name}: end {qe} clipped to length {ln}");
            qe = ln;
        }
        if qs >= qe {
            return Err(io::Error::other(format!("empty region {name}:{qs}-{qe}")));
        }
        if x != loaded {
            rec = read_recs(&dir.join("presence").join(format!("{}.bin", t.samples[x as usize])))?;
            if let Some(h) = hop {
                hop_rec = read_hop(h, &t, x)?;
            }
            loaded = x;
        }
        let label = format!("{name}:{qs}-{qe}");
        let (qs32, qe32) = (qs as u32, qe as u32);
        // (pos, sample, layer, starts); layer 0 = direct, 1 = hop-only
        let mut ev: Vec<(u32, u32, u8, bool)> = Vec::new();
        let width = (qe - qs) as f64;
        for o in 0..s {
            if o == x {
                continue;
            }
            let mut cov = [0u64; 2];
            for (layer, r) in [&rec, &hop_rec].into_iter().enumerate() {
                for_overlaps(r, sq, o, qs32, qe32, |a, b| {
                    cov[layer] += (b - a) as u64;
                    ev.push((a, o, layer as u8, true));
                    ev.push((b, o, layer as u8, false));
                });
            }
            if mode == QueryMode::Samples {
                write!(w, "{}\t{}\t{}\t{:.6}", label, t.samples[o as usize], cov[0], cov[0] as f64 / width)?;
                if hop.is_some() {
                    write!(w, "\t{}\t{:.6}", cov[1], cov[1] as f64 / width)?;
                }
                writeln!(w)?;
            }
        }
        if mode == QueryMode::Samples {
            continue;
        }
        ev.sort_unstable();
        let mut on = [vec![false; s as usize], vec![false; s as usize]];
        let mut ss = qs32;
        let join = |v: &[&str]| if v.is_empty() { "-".to_string() } else { v.join(",") };
        let emit = |a: u32, b: u32, on: &[Vec<bool>; 2], w: &mut dyn Write| -> io::Result<()> {
            if b <= a {
                return Ok(());
            }
            let (mut present, mut absent, mut hops) = (Vec::new(), Vec::new(), Vec::new());
            for o in (0..s as usize).filter(|&o| o != x as usize) {
                let name = t.samples[o].as_str();
                if on[0][o] {
                    present.push(name);
                } else if on[1][o] {
                    hops.push(name);
                } else {
                    absent.push(name);
                }
            }
            write!(w, "{label}\t{a}\t{b}\t{}\t{}\t{}", present.len() + 1, join(&present), join(&absent))?;
            if hop.is_some() {
                write!(w, "\t{}\t{}", hops.len(), join(&hops))?;
            }
            writeln!(w)
        };
        let mut j = 0;
        while j < ev.len() {
            let p = ev[j].0;
            emit(ss, p, &on, w)?;
            // sorted so that an end precedes a start of the same sample, layer and position
            while j < ev.len() && ev[j].0 == p {
                on[ev[j].2 as usize][ev[j].1 as usize] = ev[j].3;
                j += 1;
            }
            ss = p;
        }
        emit(ss, qe32, &on, w)?;
    }
    Ok(())
}

/// Read a BED file of query regions.
pub fn read_bed(path: &str) -> io::Result<Vec<(String, u64, u64)>> {
    let mut v = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            return Err(io::Error::other(format!("bad BED line: {line}")));
        }
        v.push((
            f[0].to_string(),
            f[1].parse().map_err(|_| io::Error::other(format!("bad BED line: {line}")))?,
            f[2].parse().map_err(|_| io::Error::other(format!("bad BED line: {line}")))?,
        ));
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn gz(path: &Path, text: &str) {
        let mut e = flate2::write::GzEncoder::new(File::create(path).unwrap(), flate2::Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    /// Three samples: A aligns to B over [100,300) with a 60 bp deletion inside,
    /// A aligns to C over [200,250) on the reverse strand.
    fn fixture(dir: &Path) -> (Vec<String>, String) {
        let fa = dir.join("a.fai");
        fs::write(&fa, "A#0#c1\t1000\t0\t0\t0\n").unwrap();
        let fb = dir.join("b.fai");
        fs::write(&fb, "B#0#c1\t1000\t0\t0\t0\n").unwrap();
        let fc = dir.join("c.fai");
        fs::write(&fc, "C#0#c1\t500\t0\t0\t0\n").unwrap();
        let list = dir.join("fai.list");
        fs::write(&list, format!("{}\n{}\n{}\n", fa.display(), fb.display(), fc.display())).unwrap();
        // query A [100,300) = 200 bases: 100= 60D(target) 100=  -> target B [400,660)
        let ab = dir.join("A_vs_B.paf.gz");
        gz(&ab, "A#0#c1\t1000\t100\t300\t+\tB#0#c1\t1000\t400\t660\t200\t260\t60\tcg:Z:100=60D100=\n");
        // query C [10,60) reverse onto A [200,250)
        let ca = dir.join("C_vs_A.paf.gz");
        gz(&ca, "C#0#c1\t500\t10\t60\t-\tA#0#c1\t1000\t200\t250\t50\t50\t60\tcg:Z:50=\n");
        (vec![ab.display().to_string(), ca.display().to_string()], list.display().to_string())
    }

    #[test]
    fn builds_and_queries() {
        let tmp = tempfile::tempdir().unwrap();
        let (files, list) = fixture(tmp.path());
        let out = tmp.path().join("store");
        build(&files, &list, "#", 50, 1 << 30, false, &out).unwrap();
        assert!(!out.join("rec").exists());
        let a = read_recs(&out.join("presence/A.bin")).unwrap();
        // B over A: [100,200) and [200,300) merge (abutting); C over A: [200,250)
        assert_eq!(a, vec![Rec { seq: 0, s: 100, e: 300, other: 1 }, Rec { seq: 0, s: 200, e: 250, other: 2 }]);
        let b = read_recs(&out.join("presence/B.bin")).unwrap();
        // the deletion splits B into two blocks
        assert_eq!(b, vec![Rec { seq: 1, s: 400, e: 500, other: 0 }, Rec { seq: 1, s: 560, e: 660, other: 0 }]);
        let c = read_recs(&out.join("presence/C.bin")).unwrap();
        assert_eq!(c, vec![Rec { seq: 2, s: 10, e: 60, other: 0 }]);
        let d = read_recs(&out.join("depth/A.bin")).unwrap();
        assert_eq!(
            d,
            vec![
                Rec { seq: 0, s: 0, e: 100, other: 1 },
                Rec { seq: 0, s: 100, e: 200, other: 2 },
                Rec { seq: 0, s: 200, e: 250, other: 3 },
                Rec { seq: 0, s: 250, e: 300, other: 2 },
                Rec { seq: 0, s: 300, e: 1000, other: 1 },
            ]
        );
        let stat = fs::read_to_string(out.join("stats/A.tsv")).unwrap();
        assert!(stat.contains("COV\tA\tB\t200\n") && stat.contains("COV\tA\tC\t50\n"), "{stat}");
        assert!(stat.contains("HIST\tA\t3\t50\n"), "{stat}");

        let mut buf = Vec::new();
        query(&out, None, &[("A#0#c1".into(), 150, 260)], QueryMode::Samples, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("A#0#c1:150-260\tB\t110\t1.000000"), "{text}");
        assert!(text.contains("A#0#c1:150-260\tC\t50\t0.454545"), "{text}");
        let mut buf = Vec::new();
        query(&out, None, &[("A#0#c1".into(), 150, 260)], QueryMode::Segments, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let rows: Vec<&str> = text.lines().skip(1).collect();
        assert_eq!(rows, vec!["A#0#c1:150-260\t150\t200\t2\tB\tC", "A#0#c1:150-260\t200\t250\t3\tB,C\t-", "A#0#c1:150-260\t250\t260\t2\tB\tC"]);
    }

    #[test]
    fn query_reports_hop_layer_separately() {
        let tmp = tempfile::tempdir().unwrap();
        let (files, list) = fixture(tmp.path());
        let out = tmp.path().join("store");
        build(&files, &list, "#", 50, 1 << 30, false, &out).unwrap();
        let hop = tmp.path().join("hop");
        fs::create_dir(&hop).unwrap();
        // C is direct over [200,250); hop-only continues it over [250,280)
        fs::write(hop.join("A.bin"), as_bytes(&[Rec { seq: 0, s: 250, e: 280, other: 2 }])).unwrap();
        let region = [("A#0#c1".to_string(), 190, 290)];
        let mut buf = Vec::new();
        query(&out, Some(hop.as_path()), &region, QueryMode::Segments, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let rows: Vec<&str> = text.lines().skip(1).map(|l| l.split_once('\t').unwrap().1).collect();
        assert_eq!(
            rows,
            vec!["190\t200\t2\tB\tC\t0\t-", "200\t250\t3\tB,C\t-\t0\t-", "250\t280\t2\tB\t-\t1\tC", "280\t290\t2\tB\tC\t0\t-"]
        );
        let mut buf = Vec::new();
        query(&out, Some(hop.as_path()), &region, QueryMode::Samples, &mut buf).unwrap();
        assert!(String::from_utf8(buf).unwrap().contains("\tC\t50\t0.500000\t30\t0.300000\n"));
        // a record on another sample's sequence means a different table
        fs::write(hop.join("A.bin"), as_bytes(&[Rec { seq: 1, s: 0, e: 10, other: 2 }])).unwrap();
        assert!(query(&out, Some(hop.as_path()), &region, QueryMode::Samples, &mut Vec::new()).is_err());
    }

    #[test]
    fn truncated_file_fails_the_build() {
        let tmp = tempfile::tempdir().unwrap();
        let (files, list) = fixture(tmp.path());
        let bytes = fs::read(&files[0]).unwrap();
        fs::write(&files[0], &bytes[..bytes.len() - 12]).unwrap();
        let out = tmp.path().join("store");
        let e = build(&files, &list, "#", 50, 1 << 30, false, &out).unwrap_err();
        assert!(e.to_string().contains("could not be read completely"), "{e}");
        let m = fs::read_to_string(out.join("manifest.tsv")).unwrap();
        assert!(m.lines().nth(1).unwrap().contains("\tgz_error\t") || m.lines().nth(1).unwrap().contains("\tno_final_newline\t"), "{m}");
        assert!(!out.join("scan.done").exists());
    }

    #[test]
    fn region_parsing() {
        assert_eq!(parse_region("A#0#c1:5-10").unwrap(), ("A#0#c1".to_string(), 5, 10));
        assert!(parse_region("nope").is_err());
    }
}
