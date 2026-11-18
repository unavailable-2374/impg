//! Fixed-width, externally sorted depth intermediates. Scratch files contain
//! derived coordinates only; no input alignment or sequence file is modified.
use memmap2::Mmap;
#[cfg(test)]
use memmap2::MmapMut;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

pub(super) fn scratch_file(dir: &Path, label: &str) -> io::Result<(PathBuf, File)> {
    loop {
        let path = dir.join(format!(
            "{label}-{}.bin",
            NEXT_FILE.fetch_add(1, Ordering::Relaxed)
        ));
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(file) => return Ok((path, file)),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
}

pub(super) fn write_record<const N: usize>(
    out: &mut impl Write,
    record: [u64; N],
) -> io::Result<()> {
    #[cfg(target_endian = "little")]
    {
        // Arrays contain only initialized u64 values and have no padding.
        return out
            .write_all(unsafe { std::slice::from_raw_parts(record.as_ptr().cast::<u8>(), N * 8) });
    }
    #[cfg(target_endian = "big")]
    {
        for value in record {
            out.write_all(&value.to_le_bytes())?;
        }
        Ok(())
    }
}

pub(super) fn write_rows_parallel<const N: usize>(
    dir: &Path,
    label: &str,
    rows: &[[u64; N]],
) -> io::Result<PathBuf> {
    use std::os::unix::fs::FileExt;
    let (path, file) = scratch_file(dir, label)?;
    if rows.len() < 65536 {
        let mut out = BufWriter::new(file);
        for &row in rows {
            write_record(&mut out, row)?;
        }
        out.flush()?;
        return Ok(path);
    }
    let length = rows
        .len()
        .checked_mul(N * 8)
        .ok_or_else(|| io::Error::other("depth record file size overflow"))?;
    file.set_len(length as u64)?;
    let chunk = (4usize * 1024 * 1024 / (N * 8)).max(1);
    rows.par_chunks(chunk)
        .enumerate()
        .try_for_each(|(i, rows)| -> io::Result<()> {
            #[cfg(target_endian = "little")]
            let bytes = unsafe {
                std::slice::from_raw_parts(rows.as_ptr().cast::<u8>(), rows.len() * N * 8)
            };
            #[cfg(target_endian = "big")]
            let encoded = {
                let mut bytes = Vec::with_capacity(rows.len() * N * 8);
                for &row in rows {
                    write_record(&mut bytes, row)?;
                }
                bytes
            };
            #[cfg(target_endian = "big")]
            let bytes = &encoded[..];
            file.write_all_at(bytes, (i * chunk * N * 8) as u64)
        })?;
    Ok(path)
}

fn read_record<const N: usize>(input: &mut impl Read) -> io::Result<Option<[u64; N]>> {
    let mut record = [0; N];
    let mut first = [0; 8];
    match input.read(&mut first[..1])? {
        0 => return Ok(None),
        _ => input.read_exact(&mut first[1..])?,
    }
    record[0] = u64::from_le_bytes(first);
    for value in &mut record[1..] {
        let mut bytes = [0; 8];
        input.read_exact(&mut bytes)?;
        *value = u64::from_le_bytes(bytes);
    }
    Ok(Some(record))
}

pub(super) struct Records<const N: usize> {
    pub path: PathBuf,
    map: Option<Mmap>,
    pub len: usize,
}

/// Bounded adaptive intermediates. Small/medium results remain in RAM;
/// oversized results use the same fixed-width external files as before.
pub(super) struct RecordSet<const N: usize> {
    pub path: Option<PathBuf>,
    disk: Option<Records<N>>,
    memory: Vec<[u64; N]>,
    pub len: usize,
}
impl<const N: usize> RecordSet<N> {
    pub fn memory(memory: Vec<[u64; N]>) -> Self {
        let len = memory.len();
        Self {
            path: None,
            disk: None,
            memory,
            len,
        }
    }
    pub fn disk(disk: Records<N>) -> Self {
        let len = disk.len;
        Self {
            path: Some(disk.path.clone()),
            disk: Some(disk),
            memory: Vec::new(),
            len,
        }
    }
    pub fn get(&self, index: usize) -> [u64; N] {
        self.disk
            .as_ref()
            .map_or_else(|| self.memory[index], |d| d.get(index))
    }
}

pub(super) fn merge_adaptive<const N: usize>(
    dir: &Path,
    budget: usize,
    rows: Vec<RecordSet<N>>,
) -> io::Result<RecordSet<N>> {
    let count = rows
        .iter()
        .try_fold(0usize, |n, r| n.checked_add(r.len))
        .ok_or_else(|| io::Error::other("depth row count overflow"))?;
    if rows.iter().all(|r| r.disk.is_none()) && (count as u128) * (N as u128) * 8 <= budget as u128
    {
        let mut values = Vec::with_capacity(count);
        for row in rows {
            values.extend(row.memory);
        }
        values.par_sort_unstable();
        values.dedup();
        return Ok(RecordSet::memory(values));
    }
    let mut sorter = Sorter::new(dir, budget.max(4096));
    for row in rows {
        if let Some(disk) = row.disk {
            if disk.len > 0 {
                sorter.add_sorted_file(&disk.path)?;
            }
            fs::remove_file(disk.path)?;
        } else {
            for value in row.memory {
                sorter.push(value)?;
            }
        }
    }
    sorter.finish().map(RecordSet::disk)
}

pub(super) trait RowAccess<const N: usize>: Sync {
    fn len(&self) -> usize;
    fn get(&self, index: usize) -> [u64; N];
    fn lower_bound(&self, key: [u64; N]) -> usize {
        let (mut lo, mut hi) = (0, self.len());
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.get(mid) < key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }
}
impl<const N: usize> RowAccess<N> for Records<N> {
    fn len(&self) -> usize {
        self.len
    }
    fn get(&self, i: usize) -> [u64; N] {
        Records::get(self, i)
    }
}
impl<const N: usize> RowAccess<N> for RecordSet<N> {
    fn len(&self) -> usize {
        self.len
    }
    fn get(&self, i: usize) -> [u64; N] {
        RecordSet::get(self, i)
    }
}

impl<const N: usize> Records<N> {
    pub fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let bytes = file.metadata()?.len();
        if bytes % (N as u64 * 8) != 0 {
            return Err(io::Error::other(format!(
                "truncated depth intermediate: {}",
                path.display()
            )));
        }
        let len = usize::try_from(bytes / (N as u64 * 8))
            .map_err(|_| io::Error::other("depth intermediate exceeds address space"))?;
        // The backing file is immutable for this Records lifetime.
        let map = if bytes == 0 {
            None
        } else {
            Some(unsafe { Mmap::map(&file)? })
        };
        Ok(Self { path, map, len })
    }

    pub fn get(&self, index: usize) -> [u64; N] {
        assert!(index < self.len);
        let bytes = &self.map.as_ref().unwrap()[index * N * 8..(index + 1) * N * 8];
        std::array::from_fn(|i| u64::from_le_bytes(bytes[i * 8..i * 8 + 8].try_into().unwrap()))
    }

    pub fn lower_bound(&self, key: [u64; N]) -> usize {
        let (mut lo, mut hi) = (0, self.len);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.get(mid) < key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }
}

pub(super) struct Sorter<const N: usize> {
    dir: PathBuf,
    limit: usize,
    buffer: Vec<[u64; N]>,
    runs: Vec<PathBuf>,
}

impl<const N: usize> Sorter<N> {
    /// Transfer derived sorted runs to a later merge without first merging
    /// them locally. The caller owns these scratch files, unlike borrowed
    /// checkpoint inputs passed through add_sorted_file.
    pub fn finish_runs(mut self) -> io::Result<Vec<PathBuf>> {
        self.flush()?;
        Ok(self.runs)
    }

    pub fn add_owned_runs(&mut self, runs: impl IntoIterator<Item = PathBuf>) {
        self.runs.extend(runs);
    }

    pub fn finish_adaptive(mut self) -> io::Result<RecordSet<N>> {
        if self.runs.is_empty() {
            self.buffer.par_sort_unstable();
            self.buffer.dedup();
            Ok(RecordSet::memory(self.buffer))
        } else {
            self.finish().map(RecordSet::disk)
        }
    }
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty() && self.runs.is_empty()
    }
    pub fn new(dir: &Path, bytes: usize) -> Self {
        Self {
            dir: dir.to_path_buf(),
            limit: (bytes / (N * 8)).max(1),
            buffer: Vec::new(),
            runs: Vec::new(),
        }
    }

    pub fn push(&mut self, record: [u64; N]) -> io::Result<()> {
        self.buffer.push(record);
        if self.buffer.len() >= self.limit {
            self.flush()?;
        }
        Ok(())
    }

    pub fn add_sorted_file(&mut self, path: &Path) -> io::Result<()> {
        loop {
            let target = self.dir.join(format!(
                "input-{}.bin",
                NEXT_FILE.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::hard_link(path, &target) {
                Ok(()) => {
                    self.runs.push(target);
                    return Ok(());
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        if self.buffer.len() >= 32768 {
            self.buffer.par_sort_unstable();
        } else {
            self.buffer.sort_unstable();
        }
        self.buffer.dedup();
        let path = write_rows_parallel(&self.dir, "sort", &self.buffer)?;
        self.buffer.clear();
        self.runs.push(path);
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<Records<N>> {
        self.flush()?;
        if self.runs.is_empty() {
            let (path, _) = scratch_file(&self.dir, "empty")?;
            return Records::open(path);
        }
        // Bound mappings and descriptors before the range-parallel final
        // merge. Large sample sets may produce far more runs than the
        // process VMA/descriptor limits permit opening at once.
        while self.runs.len() > 32 {
            let mut next = Vec::new();
            for batch in self.runs.chunks(32 * 8) {
                next.extend(
                    batch
                        .par_chunks(32)
                        .map(|paths| merge_serial::<N>(&self.dir, paths))
                        .collect::<io::Result<Vec<_>>>()?,
                );
            }
            self.runs = next;
        }
        if self.runs.len() > 1 && rayon::current_num_threads() > 1 {
            let input = self
                .runs
                .iter()
                .map(Records::<N>::open)
                .collect::<io::Result<Vec<_>>>()?;
            let count: usize = input.iter().map(|r| r.len).sum();
            if count >= 65536 {
                return parallel_merge(&self.dir, input);
            }
        }
        while self.runs.len() > 1 {
            let mut next = Vec::new();
            for paths in self.runs.chunks(32) {
                let (path, file) = scratch_file(&self.dir, "merge")?;
                let mut out = BufWriter::new(file);
                let mut inputs: Vec<_> = paths
                    .iter()
                    .map(|p| File::open(p).map(BufReader::new))
                    .collect::<io::Result<_>>()?;
                let mut heap = BinaryHeap::new();
                for (i, reader) in inputs.iter_mut().enumerate() {
                    if let Some(record) = read_record::<N>(reader)? {
                        heap.push(Reverse((record, i)));
                    }
                }
                let mut previous = None;
                while let Some(Reverse((record, i))) = heap.pop() {
                    if previous != Some(record) {
                        write_record(&mut out, record)?;
                        previous = Some(record);
                    }
                    if let Some(record) = read_record::<N>(&mut inputs[i])? {
                        heap.push(Reverse((record, i)));
                    }
                }
                out.flush()?;
                drop(inputs);
                for old in paths {
                    fs::remove_file(old)?;
                }
                next.push(path);
            }
            self.runs = next;
        }
        Records::open(&self.runs[0])
    }
}

/// Partition by complete record keys, so duplicates can never straddle output
/// ranges. Every range merges independently; the final copy uses disjoint
/// positional writes and preserves the exact global lexicographic order.
fn parallel_merge<const N: usize>(dir: &Path, inputs: Vec<Records<N>>) -> io::Result<Records<N>> {
    use std::os::unix::fs::FileExt;
    let partitions = rayon::current_num_threads().clamp(2, 112);
    let mut samples = Vec::new();
    for input in &inputs {
        for i in 1..partitions * 4 {
            if input.len > 0 {
                samples.push(input.get(input.len * i / (partitions * 4)));
            }
        }
    }
    samples.par_sort_unstable();
    let mut pivots = Vec::new();
    for i in 1..partitions {
        if !samples.is_empty() {
            pivots.push(samples[samples.len() * i / partitions]);
        }
    }
    pivots.dedup();
    let parts = (0..=pivots.len())
        .into_par_iter()
        .map(|part| -> io::Result<_> {
            let mut cursors = Vec::with_capacity(inputs.len());
            let mut heap = BinaryHeap::new();
            for (i, input) in inputs.iter().enumerate() {
                let start = if part == 0 {
                    0
                } else {
                    input.lower_bound(pivots[part - 1])
                };
                let end = if part == pivots.len() {
                    input.len
                } else {
                    input.lower_bound(pivots[part])
                };
                cursors.push((start, end));
                if start < end {
                    heap.push(Reverse((input.get(start), i)));
                }
            }
            let (path, file) = scratch_file(dir, "parallel-merge")?;
            let mut out = BufWriter::with_capacity(1 << 20, file);
            let mut previous = None;
            let mut count = 0u64;
            while let Some(Reverse((record, i))) = heap.pop() {
                if previous != Some(record) {
                    write_record(&mut out, record)?;
                    previous = Some(record);
                    count += 1;
                }
                cursors[i].0 += 1;
                if cursors[i].0 < cursors[i].1 {
                    heap.push(Reverse((inputs[i].get(cursors[i].0), i)));
                }
            }
            out.flush()?;
            Ok((path, count * N as u64 * 8))
        })
        .collect::<io::Result<Vec<_>>>()?;
    let (path, file) = scratch_file(dir, "merged")?;
    let mut position = 0;
    let parts: Vec<_> = parts
        .into_iter()
        .map(|(p, n)| {
            let at = position;
            position += n;
            (p, at, n)
        })
        .collect();
    file.set_len(position)?;
    parts
        .par_iter()
        .try_for_each(|(path, at, n)| -> io::Result<()> {
            let input = File::open(path)?;
            let mut buffer = vec![0u8; 1 << 20];
            let mut offset = 0;
            while offset < *n {
                let size = buffer.len().min((n - offset) as usize);
                input.read_exact_at(&mut buffer[..size], offset)?;
                file.write_all_at(&buffer[..size], at + offset)?;
                offset += size as u64;
            }
            fs::remove_file(path)
        })?;
    for input in inputs {
        fs::remove_file(input.path)?;
    }
    Records::open(path)
}

fn merge_serial<const N: usize>(dir: &Path, paths: &[PathBuf]) -> io::Result<PathBuf> {
    let (path, file) = scratch_file(dir, "merge-level")?;
    let mut out = BufWriter::with_capacity(1 << 20, file);
    let mut readers = paths
        .iter()
        .map(|p| File::open(p).map(BufReader::new))
        .collect::<io::Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (i, reader) in readers.iter_mut().enumerate() {
        if let Some(row) = read_record::<N>(reader)? {
            heap.push(Reverse((row, i)));
        }
    }
    let mut previous = None;
    while let Some(Reverse((row, i))) = heap.pop() {
        if previous != Some(row) {
            write_record(&mut out, row)?;
            previous = Some(row);
        }
        if let Some(row) = read_record::<N>(&mut readers[i])? {
            heap.push(Reverse((row, i)));
        }
    }
    out.flush()?;
    drop(readers);
    for old in paths {
        fs::remove_file(old)?;
    }
    Ok(path)
}

/// Merge new sorted boundary proposals, returning the union and only the new
/// boundaries. This lets each boundary propagate through incident runs once.
#[cfg(test)]
pub(super) fn extend_boundaries(
    old: &Records<2>,
    proposed: &Records<2>,
    dir: &Path,
) -> io::Result<(Records<2>, Records<2>)> {
    let (path, file) = scratch_file(dir, "boundaries")?;
    let (delta_path, delta_file) = scratch_file(dir, "delta")?;
    let (mut out, mut delta) = (BufWriter::new(file), BufWriter::new(delta_file));
    let (mut i, mut j) = (0, 0);
    while i < old.len || j < proposed.len {
        if j == proposed.len || (i < old.len && old.get(i) < proposed.get(j)) {
            write_record(&mut out, old.get(i))?;
            i += 1;
        } else if i < old.len && old.get(i) == proposed.get(j) {
            write_record(&mut out, old.get(i))?;
            i += 1;
            j += 1;
        } else {
            let record = proposed.get(j);
            write_record(&mut out, record)?;
            write_record(&mut delta, record)?;
            j += 1;
        }
    }
    out.flush()?;
    delta.flush()?;
    Ok((Records::open(path)?, Records::open(delta_path)?))
}

/// File-backed oriented union-find. Nodes refer to equal-length intervals;
/// parity relates their offsets, not their literal nucleotide identities.
#[cfg(test)]
pub(super) struct UnionFind {
    pub path: PathBuf,
    map: Option<MmapMut>,
}

#[cfg(test)]
impl Drop for UnionFind {
    fn drop(&mut self) {
        drop(self.map.take());
        let _ = fs::remove_file(&self.path);
    }
}

#[cfg(test)]
impl UnionFind {
    pub fn create(dir: &Path, count: usize) -> io::Result<Self> {
        let (path, file) = scratch_file(dir, "union")?;
        let mut out = BufWriter::new(file);
        for id in 0..count {
            write_record(&mut out, [id as u64, 0, 0, 0])?;
        }
        out.flush()?;
        let file = out.into_inner().map_err(|e| e.into_error())?;
        let map = if count == 0 {
            None
        } else {
            Some(unsafe { MmapMut::map_mut(&file)? })
        };
        Ok(Self { path, map })
    }

    fn get(&self, id: usize, field: usize) -> u64 {
        let p = id * 32 + field * 8;
        u64::from_le_bytes(self.map.as_ref().unwrap()[p..p + 8].try_into().unwrap())
    }
    fn set(&mut self, id: usize, field: usize, value: u64) {
        let p = id * 32 + field * 8;
        self.map.as_mut().unwrap()[p..p + 8].copy_from_slice(&value.to_le_bytes());
    }
    pub fn find(&mut self, id: usize) -> (usize, u64) {
        let (mut root, mut parity) = (id, 0);
        while self.get(root, 0) as usize != root {
            parity ^= self.get(root, 2);
            root = self.get(root, 0) as usize;
        }
        let (mut node, mut prefix) = (id, 0);
        while self.get(node, 0) as usize != node {
            let parent = self.get(node, 0) as usize;
            let edge = self.get(node, 2);
            self.set(node, 0, root as u64);
            self.set(node, 2, parity ^ prefix);
            prefix ^= edge;
            node = parent;
        }
        (root, parity)
    }
    /// False means an orientation conflict requiring further interval splits.
    pub fn join(&mut self, a: usize, b: usize, reverse: u64, width: u64) -> bool {
        let ((mut a, pa), (mut b, pb)) = (self.find(a), self.find(b));
        let parity = pa ^ pb ^ reverse;
        if a == b {
            if parity != 0 {
                if width > 1 {
                    return false;
                }
                self.set(a, 3, 1); // Orientation of a single base is ambiguous.
            }
            return true;
        }
        let (ra, rb) = (self.get(a, 1), self.get(b, 1));
        if ra < rb || (ra == rb && a > b) {
            std::mem::swap(&mut a, &mut b);
        }
        self.set(b, 0, a as u64);
        self.set(b, 2, parity);
        if ra == rb {
            self.set(a, 1, self.get(a, 1) + 1);
        }
        self.set(a, 3, self.get(a, 3) | self.get(b, 3));
        true
    }
    pub fn ambiguous(&self, root: usize) -> bool {
        self.get(root, 3) != 0
    }
    pub fn flush(&self) -> io::Result<()> {
        if let Some(map) = &self.map {
            map.flush()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn transferred_test_runs_merge_without_local_consolidation() {
        let dir = tempfile::tempdir().unwrap();
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| {
                // TEST keys overlap across producers and external merge levels.
                let mut expected = Vec::new();
                let mut merge = Sorter::<2>::new(dir.path(), 1 << 15);
                let mut files = Vec::new();
                for producer in 0..8 {
                    let mut sort = Sorter::<2>::new(dir.path(), 1 << 15);
                    for i in (0..20000u64).rev() {
                        let row = [(i * 17 + producer * 37) % 70001, i % 3];
                        sort.push(row).unwrap();
                        expected.push(row);
                    }
                    let runs = sort.finish_runs().unwrap();
                    assert!(runs.len() > 1);
                    files.extend(runs.clone());
                    merge.add_owned_runs(runs);
                }
                assert!(files.len() > 32);
                expected.sort_unstable();
                expected.dedup();
                let actual = merge.finish().unwrap();
                assert_eq!(
                    (0..actual.len).map(|i| actual.get(i)).collect::<Vec<_>>(),
                    expected
                );
                assert!(files.iter().all(|p| !p.exists()));
            });
    }
    #[test]
    fn external_sort_spills_and_deduplicates() {
        let dir = tempfile::tempdir().unwrap();
        let mut sort = Sorter::<2>::new(dir.path(), 16);
        for i in (0..100).rev() {
            sort.push([i % 17, i % 3]).unwrap();
        }
        let rows = sort.finish().unwrap();
        let values: Vec<_> = (0..rows.len).map(|i| rows.get(i)).collect();
        assert!(values.windows(2).all(|p| p[0] < p[1]));
        assert_eq!(values.len(), 51);
    }
    #[test]
    fn parallel_external_merge_preserves_test_keys_and_duplicates_at_pivots() {
        let dir = tempfile::tempdir().unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            let mut sort = Sorter::<2>::new(dir.path(), 1 << 15);
            let mut expected = Vec::new();
            for i in (0..160000).rev() {
                let row = [i % 20003, (i / 7) % 3];
                sort.push(row).unwrap();
                expected.push(row);
            }
            expected.sort_unstable();
            expected.dedup();
            let actual = sort.finish().unwrap();
            assert_eq!(
                (0..actual.len).map(|i| actual.get(i)).collect::<Vec<_>>(),
                expected
            );
        });
    }
    #[test]
    fn parallel_record_writes_preserve_every_test_chunk_boundary() {
        let dir = tempfile::tempdir().unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            let expected: Vec<_> = (0..600000u64).map(|i| [i, i.wrapping_mul(7919)]).collect();
            let path = write_rows_parallel(dir.path(), "TEST-parallel-rows", &expected).unwrap();
            let rows = Records::<2>::open(path).unwrap();
            assert_eq!(rows.len, expected.len());
            for (i, row) in expected.iter().enumerate() {
                assert_eq!(rows.get(i), *row);
            }
        });
    }
    #[test]
    fn oriented_union_detects_folds() {
        let dir = tempfile::tempdir().unwrap();
        let mut uf = UnionFind::create(dir.path(), 3).unwrap();
        assert!(uf.join(0, 1, 1, 10));
        assert!(uf.join(1, 2, 0, 10));
        assert!(!uf.join(0, 2, 0, 10));
        assert!(uf.join(0, 2, 0, 1));
        let root = uf.find(0).0;
        assert!(uf.ambiguous(root));
    }
}
