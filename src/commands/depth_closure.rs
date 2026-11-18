//! Global homology depth: coordinate constraints -> interval refinement ->
//! oriented equivalence classes -> deterministic, complete output partition.
//! Unlike the legacy anchor sweep, discovering a position never claims it.
use super::depth::{extract_sample, SampleFilter};
#[cfg(test)]
use super::depth_storage::{extend_boundaries, UnionFind};
use super::depth_storage::{
    merge_adaptive, scratch_file, write_record, RecordSet, Records, Sorter,
};
use crate::impg_index::{DepthAlignment, ImpgIndex};
use clap::ValueEnum;
#[cfg(test)]
use coitrees::{BasicCOITree, Interval, IntervalTree};
use log::info;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{hash_map::DefaultHasher, BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
#[cfg(test)]
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

mod bitmap;
mod blocks;
mod engine;
mod parallel_union;
mod relations;

pub fn compute(impg: &impl ImpgIndex, options: &Options<'_>) -> io::Result<()> {
    engine::compute(impg, options)
}

const VERSION: u32 = 7;
const MAGIC: &[u8] = b"IMPDCL07\n";

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, ValueEnum, Serialize, Deserialize)]
pub enum CopyMode {
    #[default]
    All,
    Representative,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, ValueEnum, Serialize, Deserialize)]
pub enum CoordinateMode {
    #[default]
    Exact,
    Approx,
}

pub struct Options<'a> {
    pub copies: CopyMode,
    pub coordinates: CoordinateMode,
    pub separator: &'a str,
    pub output_prefix: Option<&'a str>,
    pub fai_list: Option<&'a str>,
    pub samples: Option<&'a SampleFilter>,
    pub sequence_index: Option<&'a crate::sequence_index::UnifiedSequenceIndex>,
    pub min_seq_length: u64,
    pub reference: Option<&'a str>,
    pub ref_only: bool,
    pub window_size: Option<i64>,
    pub stats: bool,
    pub combined: bool,
    pub resume: bool,
    pub explicit_order: Option<&'a [String]>,
    pub dump_order: Option<&'a str>,
}

#[derive(Clone, Debug, Hash)]
struct Sequence {
    name: String,
    sample: String,
    length: u64,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn universe(
    impg: &impl ImpgIndex,
    options: &Options<'_>,
) -> io::Result<(Vec<Sequence>, Vec<Option<u64>>)> {
    let mut lengths = BTreeMap::new();
    for id in 0..impg.seq_index().len() as u32 {
        let name = impg
            .seq_index()
            .get_name(id)
            .ok_or_else(|| invalid("missing indexed sequence name"))?;
        let length =
            impg.seq_index()
                .get_len_from_id(id)
                .ok_or_else(|| invalid(format!("missing length for '{name}'")))? as u64;
        lengths.insert(name.to_string(), length);
    }
    if let Some(list) = options.fai_list {
        for path in fs::read_to_string(list)?
            .lines()
            .map(str::trim)
            .filter(|s| !s.is_empty() && !s.starts_with('#'))
        {
            for (line, row) in BufReader::new(File::open(path)?).lines().enumerate() {
                let row = row?;
                if row.trim().is_empty() {
                    continue;
                }
                let columns: Vec<_> = row.split('\t').collect();
                if columns.len() < 2 || columns[0].is_empty() {
                    return Err(invalid(format!("invalid FAI row {} in '{path}'", line + 1)));
                }
                let length: u64 = columns[1]
                    .parse()
                    .map_err(|_| invalid(format!("invalid FAI length at '{path}':{}", line + 1)))?;
                if let Some(old) = lengths.insert(columns[0].to_string(), length) {
                    if old != length {
                        return Err(invalid(format!(
                            "length conflict for '{}': index/FAI {old} != {length}",
                            columns[0]
                        )));
                    }
                }
            }
        }
    }
    if let Some(order) = options.explicit_order {
        for name in order {
            if !lengths.contains_key(name) {
                return Err(invalid(format!(
                    "sequence '{name}' in requested order is absent from index/FAI"
                )));
            }
        }
    }
    let sequences: Vec<_> = lengths
        .into_iter()
        .filter_map(|(name, length)| {
            let sample = extract_sample(&name, options.separator);
            (length > 0
                && length >= options.min_seq_length
                && options.samples.is_none_or(|s| s.includes(&sample)))
            .then_some(Sequence {
                name,
                sample,
                length,
            })
        })
        .collect();
    if sequences.iter().any(|s| s.length > i64::MAX as u64) {
        return Err(invalid("sequence length exceeds signed coordinate range"));
    }
    if let Some(reference) = options.reference {
        if !sequences.iter().any(|s| s.sample == reference) {
            return Err(invalid(format!(
                "reference sample '{reference}' is absent from participating sequences"
            )));
        }
    }
    if options.ref_only && options.reference.is_none() {
        return Err(invalid("--ref-only requires --ref"));
    }
    let names: BTreeMap<_, _> = sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i as u64))
        .collect();
    let mapping = (0..impg.seq_index().len() as u32)
        .map(|id| {
            impg.seq_index()
                .get_name(id)
                .and_then(|n| names.get(n).copied())
        })
        .collect();
    Ok((sequences, mapping))
}

fn memory_budget() -> io::Result<usize> {
    if let Ok(value) = std::env::var("IMPG_DEPTH_MEMORY_BYTES") {
        let bytes = value
            .parse::<usize>()
            .map_err(|_| invalid("IMPG_DEPTH_MEMORY_BYTES must be an integer byte count"))?;
        if bytes < 1 << 20 {
            return Err(invalid("IMPG_DEPTH_MEMORY_BYTES must be at least 1048576"));
        }
        return Ok(bytes);
    }
    let mut available = 64u64 << 30;
    if let Ok(text) = fs::read_to_string("/proc/meminfo") {
        if let Some(line) = text.lines().find(|s| s.starts_with("MemAvailable:")) {
            if let Some(kb) = line
                .split_whitespace()
                .nth(1)
                .and_then(|s| s.parse::<u64>().ok())
            {
                available = kb.saturating_mul(1024);
            }
        }
    }
    if let Ok(limit) = fs::read_to_string("/sys/fs/cgroup/memory.max") {
        if let Ok(limit) = limit.trim().parse::<u64>() {
            available = available.min(limit);
        }
    }
    Ok((available.saturating_mul(4) / 5)
        .min(128u64 << 30)
        .max(1 << 20) as usize)
}

// Largest-first order is only work scheduling. Equivalence classes are
// independent of this ordering. Coordinates use canonical sequence-name IDs.
#[cfg(test)]
type Runs = Records<7>;
#[cfg(test)]
fn push_run(
    out: &mut Sorter<7>,
    mut a: u64,
    mut s: u64,
    mut b: u64,
    mut t: u64,
    length: u64,
    reverse: bool,
) -> io::Result<()> {
    if length == 0 {
        return Ok(());
    }
    if (a, s) > (b, t) {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut s, &mut t);
    }
    if a == b && s == t && !reverse {
        return Ok(());
    }
    out.push([u64::MAX - length, a, s, b, t, length, reverse as u64])
}

#[cfg(test)]
fn alignment_runs(
    alignment: DepthAlignment,
    sequences: &[Sequence],
    mapping: &[Option<u64>],
    mode: CoordinateMode,
    out: &mut Sorter<7>,
) -> io::Result<()> {
    let lookup = |id: u32| {
        mapping
            .get(id as usize)
            .copied()
            .ok_or_else(|| invalid(format!("unmapped alignment sequence ID {id}")))
    };
    let (Some(a), Some(b)) = (lookup(alignment.target_id)?, lookup(alignment.query_id)?) else {
        return Ok(());
    };
    let (ts, te, qs, qe) = (
        alignment.target_start,
        alignment.target_end,
        alignment.query_start,
        alignment.query_end,
    );
    if ts < 0
        || qs < 0
        || ts >= te
        || qs >= qe
        || te as u64 > sequences[a as usize].length
        || qe as u64 > sequences[b as usize].length
    {
        return Err(invalid(format!(
            "invalid depth alignment coordinates: {}:{ts}-{te} / {}:{qs}-{qe}",
            sequences[a as usize].name, sequences[b as usize].name
        )));
    }
    if mode == CoordinateMode::Approx {
        let (tl, ql) = ((te - ts) as u64, (qe - qs) as u64);
        let (short, long) = (tl.min(ql), tl.max(ql));
        // Estimate a monotone one-to-one relation, not overlapping rounded
        // bounding boxes. Jump to the next rounding change in O(number of runs).
        let (mut i, mut previous_gap, mut begin) = (0u64, 0u64, 0u64);
        while i < short {
            let mapped = (((2 * i as u128 + 1) * long as u128) / (2 * short as u128)) as u64;
            let gap = mapped - i;
            if i == 0 {
                previous_gap = gap;
            }
            if gap != previous_gap {
                emit_approx_run(
                    out,
                    a,
                    ts as u64,
                    tl,
                    b,
                    qs as u64,
                    ql,
                    alignment.is_reverse,
                    begin,
                    i - begin,
                    previous_gap,
                )?;
                begin = i;
                previous_gap = gap;
            }
            let difference = long - short;
            i = if difference == 0 {
                short
            } else {
                let numerator = 2 * short as u128 * (gap as u128 + 1);
                let next = numerator
                    .saturating_sub(long as u128)
                    .div_ceil(2 * difference as u128);
                short.min((next as u64).max(i + 1))
            };
        }
        return emit_approx_run(
            out,
            a,
            ts as u64,
            tl,
            b,
            qs as u64,
            ql,
            alignment.is_reverse,
            begin,
            short - begin,
            previous_gap,
        );
    }
    if alignment.cigar.is_empty() {
        return Err(invalid("exact global depth requires a nonempty CIGAR"));
    }
    let (mut t, mut q) = (ts, if alignment.is_reverse { qe } else { qs });
    let direction = if alignment.is_reverse { -1 } else { 1 };
    // Coalesce M/= /X before emitting, without confusing mismatch with lack of homology.
    let (mut run_t, mut run_q, mut run_length) = (t, q, 0u64);
    for op in &alignment.cigar {
        let n = op.len() as i64;
        if n <= 0 {
            return Err(invalid("nonpositive CIGAR operation"));
        }
        match op.op() {
            'M' | '=' | 'X' => {
                if run_length == 0 {
                    run_t = t;
                    run_q = q;
                }
                run_length += n as u64;
                t = t
                    .checked_add(n)
                    .ok_or_else(|| invalid("CIGAR target overflow"))?;
                q = q
                    .checked_add(direction * n)
                    .ok_or_else(|| invalid("CIGAR query overflow"))?;
            }
            'I' | 'D' => {
                if run_length > 0 {
                    let low = if alignment.is_reverse {
                        run_q - run_length as i64
                    } else {
                        run_q
                    };
                    push_run(
                        out,
                        a,
                        run_t as u64,
                        b,
                        low as u64,
                        run_length,
                        alignment.is_reverse,
                    )?;
                    run_length = 0;
                }
                if op.op() == 'I' {
                    q = q
                        .checked_add(direction * n)
                        .ok_or_else(|| invalid("CIGAR query overflow"))?;
                } else {
                    t = t
                        .checked_add(n)
                        .ok_or_else(|| invalid("CIGAR target overflow"))?;
                }
            }
            _ => return Err(invalid("unsupported CIGAR operation in global depth")),
        }
        if t > te || t < ts || q < qs || q > qe {
            return Err(invalid("CIGAR consumption exceeds alignment endpoints"));
        }
    }
    if t != te || q != if alignment.is_reverse { qs } else { qe } {
        return Err(invalid(
            "CIGAR consumption does not match alignment endpoints",
        ));
    }
    if run_length > 0 {
        let low = if alignment.is_reverse {
            run_q - run_length as i64
        } else {
            run_q
        };
        push_run(
            out,
            a,
            run_t as u64,
            b,
            low as u64,
            run_length,
            alignment.is_reverse,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn emit_approx_run(
    out: &mut Sorter<7>,
    a: u64,
    ts: u64,
    tl: u64,
    b: u64,
    qs: u64,
    ql: u64,
    reverse: bool,
    i: u64,
    n: u64,
    gap: u64,
) -> io::Result<()> {
    // The shorter side is oriented forward in physical coordinates, making
    // reciprocal records produce the same estimated relation.
    let (t, q) = if tl <= ql {
        (
            ts + i,
            if reverse {
                qs + ql - (i + gap + n)
            } else {
                qs + i + gap
            },
        )
    } else {
        (
            if reverse {
                ts + tl - (i + gap + n)
            } else {
                ts + i + gap
            },
            qs + i,
        )
    };
    push_run(out, a, t, b, q, n, reverse)
}

#[cfg(test)]
fn initial_boundaries(
    runs: &Runs,
    sequences: &[Sequence],
    dir: &Path,
    budget: usize,
) -> io::Result<Records<2>> {
    let mut sorter = Sorter::new(dir, budget / 4);
    for (id, sequence) in sequences.iter().enumerate() {
        sorter.push([id as u64, 0])?;
        sorter.push([id as u64, sequence.length])?;
    }
    for i in 0..runs.len {
        let [_, a, s, b, t, n, _] = runs.get(i);
        for point in [[a, s], [a, s + n], [b, t], [b, t + n]] {
            sorter.push(point)?;
        }
    }
    sorter.finish()
}

#[cfg(test)]
struct BoundaryWork {
    sets: BTreeMap<u64, BTreeSet<u64>>,
    pending: std::collections::VecDeque<[u64; 2]>,
}

/// Incremental in-memory incident-run index, with a byte-budget fallback to
/// the external delta sweep. Components share only read-only relation data.
#[cfg(test)]
fn propagate_cached(
    runs: &Runs,
    bounds: &Records<2>,
    delta: &Records<2>,
    dir: &Path,
    budget: usize,
) -> io::Result<Option<(Records<2>, Records<2>)>> {
    let sequence_count = if bounds.len == 0 {
        0
    } else {
        bounds.get(bounds.len - 1)[0] as usize + 1
    };
    let base_bytes = runs
        .len
        .saturating_mul(192)
        .saturating_add(sequence_count.saturating_mul(128));
    if base_bytes.saturating_add(bounds.len.saturating_mul(104)) > budget / 2
        || runs.len > u32::MAX as usize / 2
    {
        return Ok(None);
    }
    let mut offsets = vec![0u64; sequence_count + 1];
    for i in 0..bounds.len {
        let [seq, pos] = bounds.get(i);
        offsets[seq as usize + 1] = offsets[seq as usize + 1].max(pos);
    }
    for i in 1..offsets.len() {
        let Some(sum) = offsets[i].checked_add(offsets[i - 1]) else {
            return Ok(None);
        };
        if sum > i64::MAX as u64 {
            return Ok(None);
        }
        offsets[i] = sum;
    }
    let mut parent: Vec<_> = (0..sequence_count).collect();
    fn root(parent: &[usize], mut id: usize) -> usize {
        while parent[id] != id {
            id = parent[id];
        }
        id
    }
    let mut intervals = Vec::with_capacity(runs.len * 2);
    for i in 0..runs.len {
        let [_, a, s, b, t, n, _] = runs.get(i);
        let (x, y) = (root(&parent, a as usize), root(&parent, b as usize));
        parent[x.max(y)] = x.min(y);
        for (side, seq, start) in [(0, a, s), (1, b, t)] {
            intervals.push(Interval::new(
                (offsets[seq as usize] + start) as i64,
                (offsets[seq as usize] + start + n) as i64,
                i as u64 * 2 + side,
            ));
        }
    }
    let tree = BasicCOITree::<u64, u32>::new(&intervals);
    drop(intervals);
    let mut work: BTreeMap<usize, BoundaryWork> = BTreeMap::new();
    for seq in 0..sequence_count {
        work.entry(root(&parent, seq))
            .or_insert_with(|| BoundaryWork {
                sets: BTreeMap::new(),
                pending: std::collections::VecDeque::new(),
            })
            .sets
            .insert(seq as u64, BTreeSet::new());
    }
    for i in 0..bounds.len {
        let point = bounds.get(i);
        work.get_mut(&root(&parent, point[0] as usize))
            .unwrap()
            .sets
            .get_mut(&point[0])
            .unwrap()
            .insert(point[1]);
    }
    for i in 0..delta.len {
        let point = delta.get(i);
        work.get_mut(&root(&parent, point[0] as usize))
            .unwrap()
            .pending
            .push_back(point);
    }
    let components = work.len();
    let count = AtomicUsize::new(bounds.len);
    let spill = AtomicBool::new(false);
    let results: Vec<_> = work
        .into_values()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|mut work| {
            while !spill.load(Ordering::Relaxed) {
                let Some([seq, pos]) = work.pending.pop_front() else {
                    break;
                };
                let coordinate = (offsets[seq as usize] + pos) as i64;
                let mut unfinished = false;
                tree.query(coordinate, coordinate, |hit| {
                    if spill.load(Ordering::Relaxed) {
                        unfinished = true;
                        return;
                    }
                    let [_, a, s, b, t, n, reverse] = runs.get((hit.metadata / 2) as usize);
                    let (src, start, dst, other) = if hit.metadata % 2 == 0 {
                        (a, s, b, t)
                    } else {
                        (b, t, a, s)
                    };
                    if src != seq || pos < start || pos > start + n {
                        return;
                    }
                    let next = other
                        + if reverse == 0 {
                            pos - start
                        } else {
                            n - (pos - start)
                        };
                    if work
                        .sets
                        .get_mut(&dst)
                        .expect("incident run stays in its sequence component")
                        .insert(next)
                    {
                        work.pending.push_back([dst, next]);
                        let used = count.fetch_add(1, Ordering::Relaxed) + 1;
                        if base_bytes.saturating_add(used.saturating_mul(104)) > budget / 2 {
                            spill.store(true, Ordering::Relaxed);
                        }
                    }
                });
                if unfinished {
                    work.pending.push_front([seq, pos]);
                }
            }
            work
        })
        .collect();
    let mut next = Sorter::new(dir, budget / 8);
    let mut pending = Sorter::new(dir, budget / 8);
    for work in results {
        for (seq, points) in work.sets {
            for pos in points {
                next.push([seq, pos])?;
            }
        }
        for point in work.pending {
            pending.push(point)?;
        }
    }
    info!("Incremental homology refinement: {components} independent components; {} boundaries; external continuation={}", count.load(Ordering::Relaxed), spill.load(Ordering::Relaxed));
    Ok(Some((next.finish()?, pending.finish()?)))
}

#[cfg(test)]
fn propagate(
    runs: &Runs,
    mut bounds: Records<2>,
    mut delta: Records<2>,
    dir: &Path,
    budget: usize,
    mut save: impl FnMut(&Records<2>, &Records<2>) -> io::Result<()>,
) -> io::Result<Records<2>> {
    if delta.len > 0 {
        if let Some((next, pending)) = propagate_cached(runs, &bounds, &delta, dir, budget)? {
            save(&next, &pending)?;
            for old in [&bounds.path, &delta.path] {
                if old.exists() {
                    fs::remove_file(old)?;
                }
            }
            bounds = next;
            delta = pending;
        }
    }
    let mut round = 0;
    while delta.len > 0 {
        let mut proposals = Sorter::new(dir, budget / 4);
        for i in 0..runs.len {
            let [_, a, s, b, t, n, reverse] = runs.get(i);
            for (src, begin, dst, other) in [(a, s, b, t), (b, t, a, s)] {
                let mut j = delta.lower_bound([src, begin]);
                while j < delta.len {
                    let [seq, pos] = delta.get(j);
                    if seq != src || pos > begin + n {
                        break;
                    }
                    let offset = pos - begin;
                    proposals
                        .push([dst, other + if reverse == 0 { offset } else { n - offset }])?;
                    j += 1;
                }
            }
        }
        let proposed = proposals.finish()?;
        let (next, added) = extend_boundaries(&bounds, &proposed, dir)?;
        round += 1;
        info!(
            "Homology refinement round {round}: {} boundaries, {} new",
            next.len, added.len
        );
        save(&next, &added)?;
        // Checkpoints now refer only to the new generation.
        for old in [&bounds.path, &delta.path, &proposed.path] {
            if old != &next.path && old != &added.path && old.exists() {
                fs::remove_file(old)?;
            }
        }
        bounds = next;
        delta = added;
    }
    Ok(bounds)
}

#[cfg(test)]
fn atoms_from_bounds(bounds: &Records<2>, dir: &Path) -> io::Result<Records<3>> {
    let (path, file) = scratch_file(dir, "atoms")?;
    let mut out = BufWriter::new(file);
    for i in 1..bounds.len {
        let ([a, s], [b, e]) = (bounds.get(i - 1), bounds.get(i));
        if a == b && s < e {
            write_record(&mut out, [a, s, e])?;
        }
    }
    out.flush()?;
    Records::open(path)
}

#[cfg(test)]
fn union_atoms(
    runs: &Runs,
    atoms: &Records<3>,
    dir: &Path,
    budget: usize,
) -> io::Result<(UnionFind, Records<2>)> {
    let mut uf = UnionFind::create(dir, atoms.len)?;
    let mut conflicts = Sorter::<1>::new(dir, budget / 8);
    for i in 0..runs.len {
        let [_, a, s, b, t, n, reverse] = runs.get(i);
        let mut j = atoms.lower_bound([a, s, 0]);
        while j < atoms.len {
            let [seq, start, end] = atoms.get(j);
            if seq != a || start >= s + n {
                break;
            }
            let other = t + if reverse == 0 {
                start - s
            } else {
                n - (end - s)
            };
            let k = atoms.lower_bound([b, other, 0]);
            if k == atoms.len || atoms.get(k) != [b, other, other + end - start] {
                return Err(invalid(
                    "unpropagated coordinate boundary in homology union",
                ));
            }
            if !uf.join(j, k, reverse, end - start) {
                conflicts.push([j as u64])?;
            }
            j += 1;
        }
    }
    let conflicts = conflicts.finish()?;
    let mut bad_roots = Sorter::<1>::new(dir, budget / 8);
    for i in 0..conflicts.len {
        bad_roots.push([uf.find(conflicts.get(i)[0] as usize).0 as u64])?;
    }
    let bad_roots = bad_roots.finish()?;
    let mut splits = Sorter::<2>::new(dir, budget / 8);
    for i in 0..atoms.len {
        let root = uf.find(i).0 as u64;
        let k = bad_roots.lower_bound([root]);
        if k < bad_roots.len && bad_roots.get(k) == [root] {
            let [seq, start, end] = atoms.get(i);
            if end - start > 1 {
                splits.push([seq, start + (end - start) / 2])?;
                splits.push([seq, start + (end - start).div_ceil(2)])?;
            }
        }
    }
    uf.flush()?;
    fs::remove_file(&conflicts.path)?;
    fs::remove_file(&bad_roots.path)?;
    Ok((uf, splits.finish()?))
}

trait MemberAccess: Sync {
    fn len(&self) -> usize;
    fn get(&self, index: usize) -> io::Result<[u64; 6]>;
    fn root(&self, index: usize) -> u64;
}
impl MemberAccess for Records<6> {
    fn len(&self) -> usize {
        self.len
    }
    fn get(&self, index: usize) -> io::Result<[u64; 6]> {
        Ok(Records::get(self, index))
    }
    fn root(&self, index: usize) -> u64 {
        Records::get(self, index)[0]
    }
}

// Members: root, sequence, start, end, parity, ambiguous.
#[cfg(test)]
fn collect_members(
    atoms: &Records<3>,
    uf: &mut UnionFind,
    dir: &Path,
    budget: usize,
) -> io::Result<Records<6>> {
    let mut sorter = Sorter::new(dir, budget / 4);
    for i in 0..atoms.len {
        let [seq, start, end] = atoms.get(i);
        let (root, parity) = uf.find(i);
        sorter.push([
            root as u64,
            seq,
            start,
            end,
            parity,
            uf.ambiguous(root) as u64,
        ])?;
    }
    sorter.finish()
}

// Groups: canonical sequence/start/end, member offset/count, canonical parity.
fn collect_groups(
    members: &impl MemberAccess,
    dir: &Path,
    budget: usize,
) -> io::Result<RecordSet<6>> {
    let workers = rayon::current_num_threads().max(1);
    if members.len() >= 65536 && workers > 1 {
        let mut ranges = Vec::new();
        let step = (members.len() / workers / 4).max(16384);
        let mut start = 0;
        while start < members.len() {
            let mut end = (start + step).min(members.len());
            if end < members.len() {
                let root = members.root(end - 1);
                let (mut lo, mut hi) = (end, members.len());
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    if members.root(mid) <= root {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                end = lo;
            }
            ranges.push((start, end));
            start = end;
        }
        let rows = ranges
            .into_par_iter()
            .map(|(mut i, end)| -> io::Result<_> {
                let mut sorter = Sorter::new(dir, (budget / workers / 16).max(4096));
                while i < end {
                    let [root, seq, start, finish, parity, _] = members.get(i)?;
                    let mut j = i + 1;
                    while j < end && members.root(j) == root {
                        j += 1;
                    }
                    sorter.push([seq, start, finish, i as u64, (j - i) as u64, parity])?;
                    i = j;
                }
                sorter.finish_adaptive()
            })
            .collect::<io::Result<Vec<_>>>()?;
        return merge_adaptive(dir, budget / 8, rows);
    }
    let mut sorter = Sorter::new(dir, budget / 16);
    let mut i = 0;
    while i < members.len() {
        let [root, seq, start, end, parity, _] = members.get(i)?;
        let mut j = i + 1;
        while j < members.len() && members.root(j) == root {
            j += 1;
        }
        sorter.push([seq, start, end, i as u64, (j - i) as u64, parity])?;
        i = j;
    }
    sorter.finish_adaptive()
}

#[derive(Serialize, Deserialize)]
struct State {
    version: u32,
    fingerprint: u64,
    dir: PathBuf,
    units: Vec<Option<PathBuf>>,
    runs: Option<PathBuf>,
    bounds: Option<PathBuf>,
    delta: Option<PathBuf>,
    members: Option<PathBuf>,
    groups: Option<PathBuf>,
    output_group: usize,
    output_position: Option<u64>,
    output_offset: u64,
    output_rows: u64,
    #[serde(skip)]
    output_checkpoint: Option<Instant>,
    engine: Option<engine::Checkpoint>,
}

impl State {
    fn save(&self, prefix: Option<&str>) -> io::Result<()> {
        let Some(prefix) = prefix else {
            return Ok(());
        };
        // This descriptor is the durable index of the external work files.
        let work_path = format!("{prefix}.depth.work.bin.tmp");
        let mut work = File::create(&work_path)?;
        work.write_all(MAGIC)?;
        bincode::serde::encode_into_std_write(
            (self.version, self.fingerprint, &self.dir),
            &mut work,
            bincode::config::standard(),
        )
        .map_err(|e| invalid(format!("encode homology work descriptor: {e}")))?;
        work.sync_all()?;
        fs::rename(work_path, format!("{prefix}.depth.work.bin"))?;
        let path = format!("{prefix}.depth.ckpt.tmp");
        let mut file = File::create(&path)?;
        file.write_all(MAGIC)?;
        bincode::serde::encode_into_std_write(self, &mut file, bincode::config::standard())
            .map_err(|e| invalid(format!("encode homology checkpoint: {e}")))?;
        file.sync_all()?;
        fs::rename(path, format!("{prefix}.depth.ckpt"))
    }
    fn load(prefix: &str) -> io::Result<Option<Self>> {
        let bytes = match fs::read(format!("{prefix}.depth.ckpt")) {
            Ok(bytes) => bytes,
            Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(e),
        };
        if !bytes.starts_with(MAGIC) {
            return Err(invalid("legacy depth checkpoint cannot resume global homology closure; use a fresh output prefix"));
        }
        let (state, used): (Self, _) =
            bincode::serde::decode_from_slice(&bytes[MAGIC.len()..], bincode::config::standard())
                .map_err(|e| invalid(format!("decode homology checkpoint: {e}")))?;
        if used + MAGIC.len() != bytes.len() || state.version != VERSION {
            return Err(invalid("homology checkpoint version/length mismatch"));
        }
        let work_path = format!("{prefix}.depth.work.bin");
        let work = fs::read(&work_path)?;
        if !work.starts_with(MAGIC) {
            return Err(invalid("invalid homology work descriptor magic"));
        }
        let ((version, fingerprint, dir), used): ((u32, u64, PathBuf), _) =
            bincode::serde::decode_from_slice(&work[MAGIC.len()..], bincode::config::standard())
                .map_err(|e| invalid(format!("decode homology work descriptor: {e}")))?;
        if version != state.version || fingerprint != state.fingerprint || dir != state.dir {
            return Err(invalid("homology checkpoint/work descriptor mismatch"));
        }
        OpenOptions::new()
            .write(true)
            .open(work_path)?
            .set_len((MAGIC.len() + used) as u64)?;
        Ok(Some(state))
    }
}

fn fingerprint(impg: &impl ImpgIndex, sequences: &[Sequence], o: &Options<'_>) -> io::Result<u64> {
    let mut hash = DefaultHasher::new();
    VERSION.hash(&mut hash);
    env!("CARGO_PKG_VERSION").hash(&mut hash);
    sequences.hash(&mut hash);
    (
        o.copies,
        o.coordinates,
        o.separator,
        o.reference,
        o.ref_only,
        o.window_size,
        o.stats,
        o.combined,
    )
        .hash(&mut hash);
    for path in impg.alignment_files().iter().chain(impg.sequence_files()) {
        let m = fs::metadata(path)?;
        path.hash(&mut hash);
        m.len().hash(&mut hash);
        m.modified()?.hash(&mut hash);
    }
    impg.depth_scan_units().hash(&mut hash);
    Ok(hash.finish())
}

#[derive(Clone, Copy)]
struct Span {
    first: [u64; 6],
    last: [u64; 6],
    next: usize,
}

fn member(
    members: &impl MemberAccess,
    first: [u64; 6],
    last: [u64; 6],
    i: usize,
) -> io::Result<(usize, u64, u64, u64)> {
    let a = members.get(first[3] as usize + i)?;
    let b = members.get(last[3] as usize + i)?;
    let orientation = a[4] ^ first[5];
    if a[1] != b[1] || orientation != b[4] ^ last[5] {
        return Err(invalid("inconsistent member extension"));
    }
    Ok((
        a[1] as usize,
        if orientation == 0 { a[2] } else { b[2] },
        if orientation == 0 { b[3] } else { a[3] },
        orientation,
    ))
}

fn can_extend(members: &impl MemberAccess, left: [u64; 6], right: [u64; 6]) -> io::Result<bool> {
    if left[0] != right[0] || left[2] != right[1] || left[4] != right[4] {
        return Ok(false);
    }
    for i in 0..left[4] as usize {
        let a = members.get(left[3] as usize + i)?;
        let b = members.get(right[3] as usize + i)?;
        let orientation = a[4] ^ left[5];
        if a[5] != 0
            || b[5] != 0
            || a[1] != b[1]
            || orientation != b[4] ^ right[5]
            || (if orientation == 0 {
                a[3] != b[2]
            } else {
                b[3] != a[2]
            })
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn next_span(groups: &RecordSet<6>, members: &impl MemberAccess, start: usize) -> io::Result<Span> {
    let first = groups.get(start);
    let (mut last, mut next) = (first, start + 1);
    while next < groups.len && can_extend(members, last, groups.get(next))? {
        last = groups.get(next);
        next += 1;
    }
    Ok(Span { first, last, next })
}

struct StatsOutput {
    prefix: String,
    combined: Option<BufWriter<File>>,
    counts: BTreeMap<usize, (u128, u128)>,
    paths: BTreeSet<usize>,
    beds: BTreeMap<usize, BufWriter<File>>,
}
impl StatsOutput {
    fn new(prefix: &str, combined: bool) -> io::Result<Self> {
        let combined = if combined {
            let mut out = BufWriter::new(File::create(format!("{prefix}.combined.bed"))?);
            writeln!(out, "#sequence\tstart\tend\tdepth\tsamples")?;
            Some(out)
        } else {
            None
        };
        Ok(Self {
            prefix: prefix.to_string(),
            combined,
            counts: BTreeMap::new(),
            paths: BTreeSet::new(),
            beds: BTreeMap::new(),
        })
    }
    fn add(
        &mut self,
        depth: usize,
        width: u64,
        genome_bases: u128,
        names: &[&str],
        positions: impl Iterator<Item = io::Result<(String, u64, u64)>>,
    ) -> io::Result<()> {
        let totals = self.counts.entry(depth).or_default();
        totals.0 += width as u128;
        totals.1 += genome_bases;
        if self.combined.is_none() && !self.beds.contains_key(&depth) {
            if self.beds.len() >= 32 {
                let key = *self.beds.keys().next().unwrap();
                self.beds.remove(&key).unwrap().flush()?;
            }
            let fresh = self.paths.insert(depth);
            let file = OpenOptions::new()
                .write(true)
                .create(true)
                .append(!fresh)
                .truncate(fresh)
                .open(format!("{}.depth{depth}.bed", self.prefix))?;
            self.beds.insert(depth, BufWriter::new(file));
        }
        let samples = names.join(",");
        for position in positions {
            let (name, start, end) = position?;
            if let Some(out) = self.combined.as_mut() {
                writeln!(out, "{name}\t{start}\t{end}\t{depth}\t{samples}")?;
            } else {
                writeln!(self.beds.get_mut(&depth).unwrap(), "{name}\t{start}\t{end}")?;
            }
        }
        Ok(())
    }
    fn finish(mut self) -> io::Result<()> {
        if let Some(out) = self.combined.as_mut() {
            out.flush()?;
        }
        for out in self.beds.values_mut() {
            out.flush()?;
        }
        let mut out = BufWriter::new(File::create(format!("{}.summary.txt", self.prefix))?);
        writeln!(
            out,
            "Global homology depth\n#depth\tgroup_bases\tpangenome_bases"
        )?;
        for (depth, (groups, genome)) in self.counts {
            writeln!(out, "{depth}\t{groups}\t{genome}")?;
        }
        out.flush()
    }
}

fn emit_output(
    sequences: &[Sequence],
    members: &impl MemberAccess,
    groups: &RecordSet<6>,
    options: &Options<'_>,
    state: &mut State,
    checkpoint: Option<&str>,
) -> io::Result<()> {
    let path = options.output_prefix.map(|p| format!("{p}.depth.tsv"));
    let mut file = if let Some(path) = &path {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(state.output_offset == 0)
            .open(path)?;
        if state.output_offset > file.metadata()?.len() {
            return Err(invalid("checkpoint TSV offset exceeds file length"));
        }
        file.set_len(state.output_offset)?;
        file.seek(SeekFrom::Start(state.output_offset))?;
        Some(file)
    } else {
        None
    };
    let destination: Box<dyn Write> = if let Some(file) = &file {
        Box::new(file.try_clone()?)
    } else {
        Box::new(io::stdout())
    };
    let mut out = BufWriter::with_capacity(8 << 20, destination);
    if state.output_offset == 0 {
        write!(out, "#id\tlength\tdepth\tpositions")?;
        if options.copies == CopyMode::Representative {
            write!(out, "\trepresentatives")?;
        }
        writeln!(out)?;
        if options.coordinates == CoordinateMode::Approx {
            writeln!(
                out,
                "##coordinate_mode=approx; coordinates and homology assignments are estimated"
            )?;
        }
    }
    if !options.stats && options.window_size.is_none() && state.output_position.is_none() {
        return emit_parallel(
            sequences, members, groups, options, state, checkpoint, &mut out, &mut file,
        );
    }
    let mut stats = if options.stats {
        Some(StatsOutput::new(
            options.output_prefix.unwrap_or("depth_stats"),
            options.combined,
        )?)
    } else {
        None
    };
    let mut cursor = state.output_group;
    while cursor < groups.len {
        let span = next_span(groups, members, cursor)?;
        let count = span.first[4] as usize;
        let span_width = span.last[2] - span.first[1];
        let mut representatives: BTreeMap<&str, usize> = BTreeMap::new();
        let mut chosen = 0;
        let mut found_ref = false;
        for i in 0..count {
            let (seq, start, end, _) = member(members, span.first, span.last, i)?;
            if start >= end || end > sequences[seq].length || end - start != span_width {
                return Err(invalid("invalid group member coordinates"));
            }
            let sample = sequences[seq].sample.as_str();
            let replace = if let Some(&old) = representatives.get(sample) {
                let (old_seq, old_start, old_end, _) = member(members, span.first, span.last, old)?;
                (
                    std::cmp::Reverse(end - start),
                    &sequences[seq].name,
                    start,
                    end,
                ) < (
                    std::cmp::Reverse(old_end - old_start),
                    &sequences[old_seq].name,
                    old_start,
                    old_end,
                )
            } else {
                true
            };
            if replace {
                representatives.insert(sample, i);
            }
            if !found_ref && options.reference == Some(sample) {
                chosen = i;
                found_ref = true;
            }
        }
        if !options.ref_only || found_ref {
            let (_, ref_start, ref_end, ref_orientation) =
                member(members, span.first, span.last, chosen)?;
            let depth = representatives.len();
            let mut start = state.output_position.take().unwrap_or(ref_start);
            if start < ref_start || start > ref_end {
                return Err(invalid(
                    "checkpoint window position is outside the homology group",
                ));
            }
            while start < ref_end {
                let end = if let Some(window) = options.window_size {
                    ref_end.min(((start / window as u64) + 1).saturating_mul(window as u64))
                } else {
                    ref_end
                };
                let project = |i: usize| -> io::Result<(usize, u64, u64)> {
                    let (seq, s, e, orientation) = member(members, span.first, span.last, i)?;
                    let (a, b) = (start - ref_start, end - ref_start);
                    Ok((
                        seq,
                        if orientation ^ ref_orientation == 0 {
                            s + a
                        } else {
                            e - b
                        },
                        if orientation ^ ref_orientation == 0 {
                            s + b
                        } else {
                            e - a
                        },
                    ))
                };
                state.output_rows += 1;
                write!(out, "{}\t{}\t{}\t", state.output_rows, end - start, depth)?;
                let order = std::iter::once(chosen).chain((0..count).filter(|&i| i != chosen));
                let mut genome_bases = 0u128;
                for (k, i) in order.enumerate() {
                    let (seq, s, e) = project(i)?;
                    if k > 0 {
                        write!(out, ";")?;
                    }
                    write!(out, "{}:{s}-{e}", sequences[seq].name)?;
                    genome_bases += (e - s) as u128;
                }
                if options.copies == CopyMode::Representative {
                    write!(out, "\t")?;
                    for (k, &i) in representatives.values().enumerate() {
                        let (seq, s, e) = project(i)?;
                        if k > 0 {
                            write!(out, ";")?;
                        }
                        write!(out, "{}:{s}-{e}", sequences[seq].name)?;
                    }
                }
                writeln!(out)?;
                if let Some(stats) = stats.as_mut() {
                    let positions = (0..count)
                        .map(|i| project(i).map(|(seq, s, e)| (sequences[seq].name.clone(), s, e)));
                    stats.add(
                        depth,
                        end - start,
                        genome_bases,
                        &representatives.keys().copied().collect::<Vec<_>>(),
                        positions,
                    )?;
                }
                start = end;
                if state.output_rows % 1024 == 0 {
                    state.output_group = cursor;
                    state.output_position = Some(start);
                    commit_output(&mut out, &mut file, state, checkpoint, false)?;
                }
            }
        }
        cursor = span.next;
        state.output_position = None;
        if cursor % 1024 == 0 || cursor == groups.len {
            state.output_group = cursor;
            commit_output(&mut out, &mut file, state, checkpoint, cursor == groups.len)?;
        }
    }
    commit_output(&mut out, &mut file, state, checkpoint, true)?;
    if let Some(stats) = stats {
        stats.finish()?;
    }
    Ok(())
}

// Coordinates are decoded once per span by independent formatting tasks.
// Batches preserve canonical order and bound queued text to approximately
// 8 MiB; a single unusually large row is handled as its own batch.
fn emit_parallel(
    sequences: &[Sequence],
    members: &impl MemberAccess,
    groups: &RecordSet<6>,
    options: &Options<'_>,
    state: &mut State,
    checkpoint: Option<&str>,
    out: &mut impl Write,
    file: &mut Option<File>,
) -> io::Result<()> {
    let render = |span: Span| -> io::Result<Option<Vec<u8>>> {
        let count = span.first[4] as usize;
        let positions = (0..count)
            .map(|i| member(members, span.first, span.last, i))
            .collect::<io::Result<Vec<_>>>()?;
        let span_width = span.last[2] - span.first[1];
        let mut reps = BTreeMap::<&str, usize>::new();
        let mut chosen = 0;
        let mut found = false;
        for (i, &(seq, s, e, _)) in positions.iter().enumerate() {
            if s >= e || e > sequences[seq].length || e - s != span_width {
                return Err(invalid("invalid group member coordinates"));
            }
            let sample = sequences[seq].sample.as_str();
            if !found && options.reference == Some(sample) {
                chosen = i;
                found = true;
            }
            let key = |j: usize| {
                let (seq, s, e, _) = positions[j];
                (std::cmp::Reverse(e - s), &sequences[seq].name, s, e)
            };
            if reps.get(sample).is_none_or(|&old| key(i) < key(old)) {
                reps.insert(sample, i);
            }
        }
        if options.ref_only && !found {
            return Ok(None);
        }
        let width = positions[chosen].2 - positions[chosen].1;
        let mut text = Vec::new();
        write!(text, "{width}\t{}\t", reps.len())?;
        for (k, i) in std::iter::once(chosen)
            .chain((0..count).filter(|&i| i != chosen))
            .enumerate()
        {
            let (seq, s, e, _) = positions[i];
            if k > 0 {
                text.push(b';');
            }
            write!(text, "{}:{s}-{e}", sequences[seq].name)?;
        }
        if options.copies == CopyMode::Representative {
            text.push(b'\t');
            for (k, &i) in reps.values().enumerate() {
                let (seq, s, e, _) = positions[i];
                if k > 0 {
                    text.push(b';');
                }
                write!(text, "{}:{s}-{e}", sequences[seq].name)?;
            }
        }
        text.push(b'\n');
        Ok(Some(text))
    };
    let flush = |spans: &mut Vec<Span>,
                 state: &mut State,
                 out: &mut dyn Write,
                 file: &mut Option<File>|
     -> io::Result<()> {
        if spans.is_empty() {
            return Ok(());
        }
        if spans.iter().any(|s| s.first[4] >= 65536) {
            for span in spans.iter().copied() {
                if span.first[4] >= 65536 {
                    if emit_large_span(
                        sequences,
                        members,
                        span,
                        options,
                        state.output_rows + 1,
                        out,
                    )? {
                        state.output_rows += 1;
                    }
                } else if let Some(text) = render(span)? {
                    state.output_rows += 1;
                    write!(out, "{}\t", state.output_rows)?;
                    out.write_all(&text)?;
                }
                state.output_group = span.next;
                commit_output(out, file, state, checkpoint, false)?;
            }
            spans.clear();
            return Ok(());
        }
        let rows = spans
            .par_iter()
            .map(|&s| render(s))
            .collect::<io::Result<Vec<_>>>()?;
        for text in rows.into_iter().flatten() {
            state.output_rows += 1;
            write!(out, "{}\t", state.output_rows)?;
            out.write_all(&text)?;
        }
        state.output_group = spans.last().unwrap().next;
        spans.clear();
        commit_output(out, file, state, checkpoint, false)
    };
    let initial = state.output_group;
    let mut carry: Option<Span> = None;
    let mut pending = Vec::new();
    let mut size = 0usize;
    let longest = sequences.iter().map(|s| s.name.len()).max().unwrap_or(0) + 48;
    for begin in (initial..groups.len).step_by(32768) {
        let end = (begin + 32768).min(groups.len);
        let joins = (begin..end)
            .into_par_iter()
            .map(|i| {
                if i == initial {
                    Ok(false)
                } else {
                    can_extend(members, groups.get(i - 1), groups.get(i))
                }
            })
            .collect::<io::Result<Vec<_>>>()?;
        for (i, join) in (begin..end).zip(joins) {
            let group = groups.get(i);
            if join {
                let s = carry
                    .as_mut()
                    .ok_or_else(|| invalid("missing output span carry"))?;
                s.last = group;
                s.next = i + 1;
            } else {
                if let Some(span) = carry.take() {
                    size = size.saturating_add(
                        (span.first[4] as usize)
                            .saturating_mul(longest)
                            .saturating_mul(if options.copies == CopyMode::Representative {
                                2
                            } else {
                                1
                            }),
                    );
                    pending.push(span);
                    if size >= 8 << 20 {
                        flush(&mut pending, state, out, file)?;
                        size = 0;
                    }
                }
                carry = Some(Span {
                    first: group,
                    last: group,
                    next: i + 1,
                });
            }
        }
        // Flush only completed spans. The final span can extend into the next
        // bounded comparison batch and remains uncommitted until complete.
        flush(&mut pending, state, out, file)?;
        size = 0;
    }
    if let Some(span) = carry {
        pending.push(span);
    }
    flush(&mut pending, state, out, file)?;
    commit_output(out, file, state, checkpoint, true)
}

/// A single homology class can contain millions of copies. Scan its samples
/// in parallel, then stream bounded coordinate chunks in canonical order;
/// neither all decoded coordinates nor the entire TSV row is buffered.
fn emit_large_span(
    sequences: &[Sequence],
    members: &impl MemberAccess,
    span: Span,
    options: &Options<'_>,
    row: u64,
    out: &mut dyn Write,
) -> io::Result<bool> {
    type Choice = (usize, u64, u64, usize);
    let count = span.first[4] as usize;
    let width = span.last[2] - span.first[1];
    let (reps, reference) = (0..count)
        .into_par_iter()
        .try_fold(
            || (BTreeMap::<&str, Choice>::new(), None::<usize>),
            |(mut reps, mut reference), i| -> io::Result<_> {
                let (seq, s, e, _) = member(members, span.first, span.last, i)?;
                if s >= e || e > sequences[seq].length || e - s != width {
                    return Err(invalid("invalid group member coordinates"));
                }
                let sample = sequences[seq].sample.as_str();
                let choice = (seq, s, e, i);
                reps.entry(sample)
                    .and_modify(|old| *old = (*old).min(choice))
                    .or_insert(choice);
                if options.reference == Some(sample) {
                    reference = Some(reference.map_or(i, |old| old.min(i)));
                }
                Ok((reps, reference))
            },
        )
        .try_reduce(
            || (BTreeMap::new(), None),
            |(mut a, ar), (b, br)| -> io::Result<_> {
                for (sample, choice) in b {
                    a.entry(sample)
                        .and_modify(|old| *old = (*old).min(choice))
                        .or_insert(choice);
                }
                Ok((
                    a,
                    match (ar, br) {
                        (Some(a), Some(b)) => Some(a.min(b)),
                        (a, b) => a.or(b),
                    },
                ))
            },
        )?;
    if options.ref_only && reference.is_none() {
        return Ok(false);
    }
    let chosen = reference.unwrap_or(0);
    let (seq, s, e, _) = member(members, span.first, span.last, chosen)?;
    write!(
        out,
        "{row}\t{width}\t{}\t{}:{s}-{e}",
        reps.len(),
        sequences[seq].name
    )?;
    let longest = sequences.iter().map(|s| s.name.len()).max().unwrap_or(0) + 48;
    let chunk = (4usize * 1024 * 1024 / longest).clamp(1, 65536);
    let workers = rayon::current_num_threads().max(1);
    let batch = chunk.saturating_mul(workers);
    for begin in (0..count).step_by(batch) {
        let end = (begin + batch).min(count);
        let jobs = (end - begin).div_ceil(chunk);
        let text = (0..jobs)
            .into_par_iter()
            .map(|j| -> io::Result<_> {
                let mut text = Vec::new();
                for i in begin + j * chunk..(begin + (j + 1) * chunk).min(end) {
                    if i == chosen {
                        continue;
                    }
                    let (seq, s, e, _) = member(members, span.first, span.last, i)?;
                    write!(text, ";{}:{s}-{e}", sequences[seq].name)?;
                }
                Ok(text)
            })
            .collect::<io::Result<Vec<_>>>()?;
        for text in text {
            out.write_all(&text)?;
        }
    }
    if options.copies == CopyMode::Representative {
        write!(out, "\t")?;
        for (k, &(seq, s, e, _)) in reps.values().enumerate() {
            if k > 0 {
                write!(out, ";")?;
            }
            write!(out, "{}:{s}-{e}", sequences[seq].name)?;
        }
    }
    writeln!(out)?;
    Ok(true)
}

fn commit_output(
    out: &mut (impl Write + ?Sized),
    file: &mut Option<File>,
    state: &mut State,
    checkpoint: Option<&str>,
    force: bool,
) -> io::Result<()> {
    let now = Instant::now();
    let last = state.output_checkpoint.get_or_insert(now);
    let hook = checkpoint.is_some()
        && (std::env::var("IMPG_TEST_EXIT_AFTER_STAGE").as_deref() == Ok("output")
            || std::env::var_os("IMPG_TEST_EXIT_AFTER_COMMIT").is_some());
    if !force && !hook && (checkpoint.is_none() || now.duration_since(*last).as_secs() < 300) {
        return Ok(());
    }
    out.flush()?;
    if let Some(file) = file.as_mut() {
        file.sync_data()?;
        state.output_offset = file.stream_position()?;
    }
    state.save(checkpoint)?;
    state.output_checkpoint = Some(Instant::now());
    test_exit_after_stage(checkpoint, "output");
    if checkpoint.is_some() {
        test_exit_after_commit();
    }
    Ok(())
}

fn test_exit_after_commit() {
    if std::env::var_os("IMPG_TEST_EXIT_AFTER_COMMIT").is_some() {
        eprintln!("[impg test hook] stopping after global homology checkpoint");
        std::process::exit(99);
    }
}

fn test_exit_after_stage(checkpoint: Option<&str>, stage: &str) {
    if checkpoint.is_some() && std::env::var("IMPG_TEST_EXIT_AFTER_STAGE").as_deref() == Ok(stage) {
        eprintln!("[impg test hook] stopping after {stage} checkpoint");
        std::process::exit(99);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impg::CigarOp;

    fn test_sequences(lengths: &[u64]) -> Vec<Sequence> {
        lengths
            .iter()
            .enumerate()
            .map(|(i, &length)| Sequence {
                name: format!("TEST_{}#0#chr{i}", i % 3),
                sample: format!("TEST_{}", i % 3),
                length,
            })
            .collect()
    }

    fn test_alignment(
        a: u32,
        s: i64,
        b: u32,
        t: i64,
        reverse: bool,
        ops: &[(i32, char)],
    ) -> DepthAlignment {
        let cigar: Vec<_> = ops.iter().map(|&(n, op)| CigarOp::new(n, op)).collect();
        let te = s + cigar.iter().map(|op| op.target_delta() as i64).sum::<i64>();
        let qe = t + cigar
            .iter()
            .map(|op| op.query_delta(crate::alignment_record::Strand::Forward) as i64)
            .sum::<i64>();
        DepthAlignment {
            target_id: a,
            query_id: b,
            target_start: s,
            target_end: te,
            query_start: t,
            query_end: qe,
            is_reverse: reverse,
            cigar,
        }
    }

    // Independent TEST-only oracle: explicitly enumerate base correspondences.
    // Production code never allocates a per-base union-find.
    fn oracle(
        sequences: &[Sequence],
        alignments: &[DepthAlignment],
        mode: CoordinateMode,
    ) -> Vec<Vec<usize>> {
        let mut offsets = vec![0usize];
        for s in sequences {
            offsets.push(offsets.last().unwrap() + s.length as usize);
        }
        let mut parent: Vec<_> = (0..*offsets.last().unwrap()).collect();
        fn root(parent: &[usize], mut i: usize) -> usize {
            while parent[i] != i {
                i = parent[i];
            }
            i
        }
        let mut join = |a: u32, p: i64, b: u32, q: i64| {
            let x = root(&parent, offsets[a as usize] + p as usize);
            let y = root(&parent, offsets[b as usize] + q as usize);
            parent[x] = y;
        };
        for r in alignments {
            if mode == CoordinateMode::Approx {
                let (tl, ql) = (
                    (r.target_end - r.target_start) as u64,
                    (r.query_end - r.query_start) as u64,
                );
                for i in 0..tl.min(ql) {
                    let j = (((2 * i + 1) as u128 * tl.max(ql) as u128) / (2 * tl.min(ql)) as u128)
                        as i64;
                    let (t, q) = if tl <= ql {
                        (
                            r.target_start + i as i64,
                            if r.is_reverse {
                                r.query_end - 1 - j
                            } else {
                                r.query_start + j
                            },
                        )
                    } else {
                        (
                            if r.is_reverse {
                                r.target_end - 1 - j
                            } else {
                                r.target_start + j
                            },
                            r.query_start + i as i64,
                        )
                    };
                    join(r.target_id, t, r.query_id, q);
                }
            } else {
                let (mut t, mut q) = (
                    r.target_start,
                    if r.is_reverse {
                        r.query_end - 1
                    } else {
                        r.query_start
                    },
                );
                let direction = if r.is_reverse { -1 } else { 1 };
                for op in &r.cigar {
                    for _ in 0..op.len() {
                        match op.op() {
                            'M' | '=' | 'X' => {
                                join(r.target_id, t, r.query_id, q);
                                t += 1;
                                q += direction;
                            }
                            'I' => q += direction,
                            'D' => t += 1,
                            _ => unreachable!(),
                        }
                    }
                }
            }
        }
        let mut groups: BTreeMap<_, Vec<_>> = BTreeMap::new();
        for i in 0..parent.len() {
            groups.entry(root(&parent, i)).or_default().push(i);
        }
        let mut groups: Vec<_> = groups.into_values().collect();
        groups.sort();
        groups
    }

    fn interval_result(
        sequences: &[Sequence],
        alignments: Vec<DepthAlignment>,
        mode: CoordinateMode,
        budget: usize,
    ) -> Vec<Vec<usize>> {
        let incremental = engine::test_partition(sequences, &alignments, mode).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let mapping: Vec<_> = (0..sequences.len()).map(|i| Some(i as u64)).collect();
        let mut sorter = Sorter::new(dir.path(), 1024);
        for alignment in alignments {
            alignment_runs(alignment, sequences, &mapping, mode, &mut sorter).unwrap();
        }
        let runs = sorter.finish().unwrap();
        let mut bounds = initial_boundaries(&runs, sequences, dir.path(), 4096).unwrap();
        let mut delta = Records::open(&bounds.path).unwrap();
        let (atoms, mut uf) = loop {
            bounds = propagate(&runs, bounds, delta, dir.path(), budget, |_, _| Ok(())).unwrap();
            let atoms = atoms_from_bounds(&bounds, dir.path()).unwrap();
            let (uf, splits) = union_atoms(&runs, &atoms, dir.path(), 4096).unwrap();
            if splits.len == 0 {
                break (atoms, uf);
            }
            (bounds, delta) = extend_boundaries(&bounds, &splits, dir.path()).unwrap();
            assert!(delta.len > 0);
        };
        let members = collect_members(&atoms, &mut uf, dir.path(), 4096).unwrap();
        let groups = collect_groups(&members, dir.path(), 4096).unwrap();
        let mut offsets = vec![0usize];
        for s in sequences {
            offsets.push(offsets.last().unwrap() + s.length as usize);
        }
        let mut seen = vec![false; *offsets.last().unwrap()];
        let mut actual = Vec::new();
        let mut cursor = 0;
        while cursor < groups.len {
            let span = next_span(&groups, &members, cursor).unwrap();
            let (_, start, end, _) = member(&members, span.first, span.last, 0).unwrap();
            for offset in 0..end - start {
                let mut column = Vec::new();
                for i in 0..span.first[4] as usize {
                    let (seq, s, e, reverse) = member(&members, span.first, span.last, i).unwrap();
                    let base = offsets[seq]
                        + if reverse == 0 {
                            s + offset
                        } else {
                            e - 1 - offset
                        } as usize;
                    assert!(!seen[base], "TEST duplicate base {base}");
                    seen[base] = true;
                    column.push(base);
                }
                column.sort();
                actual.push(column);
            }
            cursor = span.next;
        }
        assert!(seen.into_iter().all(|b| b), "TEST missing bases");
        actual.sort();
        assert_eq!(
            incremental, actual,
            "TEST incremental/reference partition mismatch"
        );
        actual
    }

    fn check(lengths: &[u64], make: impl Fn() -> Vec<DepthAlignment>) {
        let sequences = test_sequences(lengths);
        for mode in [CoordinateMode::Exact, CoordinateMode::Approx] {
            let expected = oracle(&sequences, &make(), mode);
            assert_eq!(
                interval_result(&sequences, make(), mode, 1 << 20),
                expected,
                "TEST mode {mode:?}"
            );
            let mut reversed = make();
            reversed.reverse();
            assert_eq!(
                interval_result(&sequences, reversed, mode, 4096),
                expected,
                "TEST source order {mode:?}"
            );
        }
    }

    #[test]
    fn depth_closure_chain_and_multiple_copies_match_base_oracle() {
        check(&[12, 12, 12, 12], || {
            vec![
                test_alignment(0, 0, 1, 0, false, &[(12, '=')]),
                test_alignment(1, 0, 2, 0, false, &[(12, 'X')]),
                test_alignment(3, 0, 1, 0, false, &[(12, 'M')]),
            ]
        });
    }
    #[test]
    fn depth_closure_indels_reverse_and_partial_copies_match_base_oracle() {
        check(&[20, 20, 20, 20], || {
            vec![
                test_alignment(
                    0,
                    1,
                    1,
                    2,
                    false,
                    &[(3, '='), (2, 'I'), (4, 'X'), (3, 'D'), (3, 'M')],
                ),
                test_alignment(1, 3, 2, 1, true, &[(2, 'M'), (1, 'D'), (6, 'M')]),
                test_alignment(3, 7, 0, 3, true, &[(6, 'M')]),
            ]
        });
    }
    #[test]
    fn depth_closure_folded_and_shifted_cycles_match_base_oracle() {
        check(&[11], || {
            vec![test_alignment(0, 0, 0, 0, true, &[(11, 'M')])]
        });
        check(&[11, 10], || {
            vec![
                test_alignment(0, 0, 1, 0, false, &[(10, 'M')]),
                test_alignment(0, 1, 1, 0, false, &[(10, 'M')]),
            ]
        });
    }
    #[test]
    fn depth_closure_cached_refinement_spills_without_losing_boundaries() {
        let sequences = test_sequences(&[120, 120]);
        let alignments = vec![
            test_alignment(0, 0, 1, 0, false, &[(119, 'M')]),
            test_alignment(0, 1, 1, 0, false, &[(119, 'M')]),
        ];
        let expected = oracle(&sequences, &alignments, CoordinateMode::Exact);
        assert_eq!(
            interval_result(&sequences, alignments, CoordinateMode::Exact, 4096),
            expected
        );
    }
    #[test]
    fn depth_closure_linear_estimate_matches_definition_at_varied_ratios() {
        for t in 1..12 {
            for q in 1..12 {
                check(&[12, 12], || {
                    vec![test_alignment(
                        0,
                        0,
                        1,
                        0,
                        true,
                        &[(t.min(q), 'M')]
                            .into_iter()
                            .chain(if t > q {
                                vec![(t - q, 'D')]
                            } else {
                                vec![(q - t, 'I')]
                            })
                            .filter(|(n, _)| *n > 0)
                            .collect::<Vec<_>>(),
                    )]
                });
            }
        }
    }
    #[test]
    fn depth_closure_fuzz_matches_independent_base_components() {
        let mut state = 0x84ab_6789_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            state
        };
        for _ in 0..60 {
            let specifications: Vec<_> = (0..8)
                .map(|_| {
                    let a = (next() % 4) as u32;
                    let b = (next() % 4) as u32;
                    let s = (next() % 8) as i64;
                    let t = (next() % 8) as i64;
                    let n = (next() % (12 - s.max(t) as u64) + 1) as i32;
                    (a, s, b, t, next() % 3 == 0, n)
                })
                .collect();
            check(&[12, 12, 12, 12], || {
                specifications
                    .iter()
                    .map(|&(a, s, b, t, r, n)| test_alignment(a, s, b, t, r, &[(n, 'M')]))
                    .collect()
            });
        }
    }
}
