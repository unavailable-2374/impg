//! Sparse/dense coordinate blocks. A bit denotes a discovered boundary,
//! never ownership of a genomic base or permission to discard a relation.
use std::io;

pub(super) const SHIFT: u32 = 20;
pub(super) const WIDTH: u64 = 1 << SHIFT;
const WORDS: usize = WIDTH as usize / 64;
const DENSE_THRESHOLD: usize = WIDTH as usize / 32;

#[derive(Clone)]
pub(super) enum Bits {
    Sparse(Vec<u32>),
    Dense {
        words: Vec<u64>,
        ranks: Vec<u32>,
        count: u32,
    },
}

impl Default for Bits {
    fn default() -> Self {
        Self::Sparse(Vec::new())
    }
}

impl Bits {
    /// First and last set offsets in [begin,end), without scanning the points
    /// between them. Sparse frontiers can seek directly into long CIGARs.
    pub fn range_bounds(&self, begin: u32, end: u32) -> Option<(u32, u32)> {
        if begin >= end {
            return None;
        }
        match self {
            Self::Sparse(points) => {
                let a = points.partition_point(|&p| p < begin);
                let b = points.partition_point(|&p| p < end);
                (a < b).then(|| (points[a], points[b - 1]))
            }
            Self::Dense { words, .. } => {
                let first = self.range(begin, end).next()?;
                let mut i = (end as usize - 1) / 64;
                let bit = (end - 1) % 64;
                let mut word = words[i]
                    & if bit == 63 {
                        u64::MAX
                    } else {
                        (1u64 << (bit + 1)) - 1
                    };
                while word == 0 {
                    i -= 1;
                    word = words[i];
                }
                Some((first, (i * 64 + 63 - word.leading_zeros() as usize) as u32))
            }
        }
    }
    pub fn union(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Dense { words: a, .. }, Self::Dense { words: b, .. }) => {
                Self::dense(a.iter().zip(b).map(|(a, b)| a | b).collect())
            }
            (Self::Dense { words, .. }, Self::Sparse(points))
            | (Self::Sparse(points), Self::Dense { words, .. }) => {
                let mut words = words.clone();
                for &p in points {
                    words[p as usize / 64] |= 1u64 << (p % 64);
                }
                Self::dense(words)
            }
            (Self::Sparse(a), Self::Sparse(b)) => {
                let (mut i, mut j) = (0, 0);
                let mut out = Vec::with_capacity(a.len() + b.len());
                while i < a.len() && j < b.len() {
                    let p = a[i].min(b[j]);
                    out.push(p);
                    if a[i] == p {
                        i += 1;
                    }
                    if b[j] == p {
                        j += 1;
                    }
                }
                out.extend_from_slice(&a[i..]);
                out.extend_from_slice(&b[j..]);
                // Both inputs and their merge are already sorted and unique.
                if out.len() < DENSE_THRESHOLD {
                    Self::Sparse(out)
                } else {
                    let mut words = vec![0u64; WORDS];
                    for p in out {
                        words[p as usize / 64] |= 1u64 << (p % 64);
                    }
                    Self::dense(words)
                }
            }
        }
    }
    pub fn difference(&self, old: &Self) -> Self {
        match (self, old) {
            (Self::Dense { words: a, .. }, Self::Dense { words: b, .. }) => {
                Self::dense(a.iter().zip(b).map(|(a, b)| a & !b).collect())
            }
            (Self::Dense { words, .. }, Self::Sparse(points)) => {
                let mut words = words.clone();
                for &p in points {
                    words[p as usize / 64] &= !(1u64 << (p % 64));
                }
                Self::dense(words)
            }
            _ => Self::from_sorted(self.values().filter(|&p| !old.contains(p)).collect()),
        }
    }
    pub fn successor(&self, pos: u32) -> Option<u32> {
        match self {
            Self::Sparse(v) => v.get(v.partition_point(|&p| p <= pos)).copied(),
            Self::Dense { words, .. } => {
                let mut i = pos as usize / 64;
                let mut word = if pos % 64 == 63 {
                    0
                } else {
                    words[i] & (!0u64 << (pos % 64 + 1))
                };
                loop {
                    if word != 0 {
                        return Some((i * 64 + word.trailing_zeros() as usize) as u32);
                    }
                    i += 1;
                    word = *words.get(i)?;
                }
            }
        }
    }
    pub fn range(&self, begin: u32, end: u32) -> impl Iterator<Item = u32> + '_ {
        let mut iter = self.values();
        match self {
            Self::Sparse(v) => iter.index = v.partition_point(|&p| p < begin),
            Self::Dense { words, .. } => {
                let i = begin as usize / 64;
                iter.word = i + 1;
                iter.pending = words.get(i).copied().unwrap_or(0) & (!0u64 << (begin % 64));
            }
        }
        iter.take_while(move |&p| p < end)
    }
    pub fn from_sorted(mut values: Vec<u32>) -> Self {
        values.sort_unstable();
        values.dedup();
        assert!(values.last().is_none_or(|&v| (v as u64) < WIDTH));
        if values.len() < DENSE_THRESHOLD {
            return Self::Sparse(values);
        }
        let mut words = vec![0u64; WORDS];
        for &v in &values {
            words[v as usize / 64] |= 1 << (v % 64);
        }
        Self::dense(words)
    }
    fn dense(words: Vec<u64>) -> Self {
        let mut ranks = Vec::with_capacity(WORDS / 8 + 1);
        let mut count = 0;
        for chunk in words.chunks(8) {
            ranks.push(count);
            count += chunk.iter().map(|w| w.count_ones()).sum::<u32>();
        }
        ranks.push(count);
        Self::Dense {
            words,
            ranks,
            count,
        }
    }
    pub fn len(&self) -> usize {
        match self {
            Self::Sparse(v) => v.len(),
            Self::Dense { count, .. } => *count as usize,
        }
    }
    pub fn contains(&self, offset: u32) -> bool {
        if offset as u64 >= WIDTH {
            return false;
        }
        match self {
            Self::Sparse(v) => v.binary_search(&offset).is_ok(),
            Self::Dense { words, .. } => words[offset as usize / 64] & (1 << (offset % 64)) != 0,
        }
    }
    /// Number of boundaries strictly before `offset`.
    pub fn rank(&self, offset: u32) -> usize {
        if offset as u64 >= WIDTH {
            return self.len();
        }
        match self {
            Self::Sparse(v) => v.partition_point(|&v| v < offset),
            Self::Dense { words, ranks, .. } => {
                let word = offset as usize / 64;
                let block = word / 8;
                ranks[block] as usize
                    + words[block * 8..word]
                        .iter()
                        .map(|w| w.count_ones() as usize)
                        .sum::<usize>()
                    + (words[word] & ((1u64 << (offset % 64)).wrapping_sub(1))).count_ones()
                        as usize
            }
        }
    }
    pub fn select(&self, index: usize) -> u32 {
        assert!(index < self.len());
        match self {
            Self::Sparse(v) => v[index],
            Self::Dense { words, ranks, .. } => {
                let block = ranks.partition_point(|&v| v as usize <= index) - 1;
                let mut remaining = index - ranks[block] as usize;
                for (i, &word) in words[block * 8..block * 8 + 8].iter().enumerate() {
                    let count = word.count_ones() as usize;
                    if remaining < count {
                        let mut bits = word;
                        for _ in 0..remaining {
                            bits &= bits - 1;
                        }
                        return ((block * 8 + i) * 64 + bits.trailing_zeros() as usize) as u32;
                    }
                    remaining -= count;
                }
                unreachable!("rank directory disagrees with dense words")
            }
        }
    }
    pub fn values(&self) -> BitIter<'_> {
        BitIter {
            bits: self,
            index: 0,
            word: 0,
            pending: 0,
        }
    }
    pub fn bytes(&self) -> usize {
        match self {
            Self::Sparse(v) => v.capacity() * 4,
            Self::Dense { words, ranks, .. } => words.capacity() * 8 + ranks.capacity() * 4,
        }
    }
    /// Return the union and only newly discovered positions. Existing points
    /// stay present, but no old frontier is enqueued for a second time.
    pub fn extend(&self, proposed: &[u32]) -> (Self, Self) {
        let new: Vec<_> = proposed
            .iter()
            .copied()
            .filter(|&p| !self.contains(p))
            .collect();
        if new.is_empty() {
            return (self.clone(), Self::default());
        }
        let delta = Self::from_sorted(new);
        let merged = match self {
            Self::Sparse(v) => Self::from_sorted(v.iter().copied().chain(delta.values()).collect()),
            Self::Dense { words, .. } => {
                let mut words = words.clone();
                for p in delta.values() {
                    words[p as usize / 64] |= 1 << (p % 64);
                }
                Self::dense(words)
            }
        };
        (merged, delta)
    }
    pub fn encode(&self, out: &mut impl io::Write) -> io::Result<()> {
        match self {
            Self::Sparse(v) => {
                out.write_all(&[0])?;
                out.write_all(&(v.len() as u32).to_le_bytes())?;
                for &p in v {
                    out.write_all(&p.to_le_bytes())?;
                }
            }
            Self::Dense { words, .. } => {
                out.write_all(&[1])?;
                for &word in words {
                    out.write_all(&word.to_le_bytes())?;
                }
            }
        }
        Ok(())
    }
    pub fn decode(input: &mut impl io::Read) -> io::Result<Self> {
        let mut kind = [0];
        input.read_exact(&mut kind)?;
        match kind[0] {
            0 => {
                let mut n = [0; 4];
                input.read_exact(&mut n)?;
                let n = u32::from_le_bytes(n) as usize;
                if n > WIDTH as usize {
                    return Err(io::Error::other("invalid sparse boundary count"));
                }
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    let mut p = [0; 4];
                    input.read_exact(&mut p)?;
                    v.push(u32::from_le_bytes(p));
                }
                if v.windows(2).any(|p| p[0] >= p[1])
                    || v.last().is_some_and(|&p| p as u64 >= WIDTH)
                {
                    return Err(io::Error::other("invalid sorted boundary block"));
                }
                Ok(Self::Sparse(v))
            }
            1 => {
                let mut words = vec![0; WORDS];
                for w in &mut words {
                    let mut bytes = [0; 8];
                    input.read_exact(&mut bytes)?;
                    *w = u64::from_le_bytes(bytes);
                }
                Ok(Self::dense(words))
            }
            _ => Err(io::Error::other("invalid boundary block encoding")),
        }
    }
}

pub(super) struct BitIter<'a> {
    bits: &'a Bits,
    index: usize,
    word: usize,
    pending: u64,
}
impl Iterator for BitIter<'_> {
    type Item = u32;
    #[inline]
    fn next(&mut self) -> Option<u32> {
        match self.bits {
            Bits::Sparse(v) => {
                let p = v.get(self.index).copied();
                self.index += usize::from(p.is_some());
                p
            }
            Bits::Dense { words, .. } => loop {
                if self.pending != 0 {
                    let p = (self.word - 1) * 64 + self.pending.trailing_zeros() as usize;
                    self.pending &= self.pending - 1;
                    return Some(p as u32);
                }
                self.pending = *words.get(self.word)?;
                self.word += 1;
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_dense_rank_select_and_frontier_match_test_positions() {
        // TEST data, not experimental genome coordinates.
        for step in [1, 7, 31, 32, 33, 1009, 100_000] {
            let values: Vec<_> = (0..WIDTH as u32).step_by(step).collect();
            let bits = Bits::from_sorted(values.clone());
            for (i, &p) in values.iter().enumerate() {
                assert_eq!(bits.select(i), p);
                assert_eq!(bits.rank(p), i);
                assert!(bits.contains(p));
            }
            assert_eq!(bits.rank(WIDTH as u32), values.len());
            let mut bytes = Vec::new();
            bits.encode(&mut bytes).unwrap();
            let decoded = Bits::decode(&mut &bytes[..]).unwrap();
            assert!(bits.values().eq(decoded.values()));
            let (next, delta) = bits.extend(&[0, 1, 2, WIDTH as u32 - 1]);
            let expected: Vec<_> = [0, 1, 2, WIDTH as u32 - 1]
                .into_iter()
                .filter(|&p| !bits.contains(p))
                .collect();
            assert_eq!(delta.values().collect::<Vec<_>>(), expected);
            assert_eq!(next.len(), bits.len() + expected.len());
        }
    }
}
