//! Concurrent oriented disjoint sets with immutable disk checkpoints.
//! Fixed, deterministic randomized priorities avoid a contended root-rank
//! update. Every link strictly increases that total order, preventing cycles.
//! One atomic word holds parent, parity and a persistent odd-cycle witness.
use super::blocks::{Page, Pages};
use super::invalid;
use parking_lot::{Condvar, Mutex};
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;
use std::collections::{BTreeMap, VecDeque};
use std::io;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;

const MAX_PAGE_BYTES: usize = 2 << 20;
const WITNESS: u64 = 1 << 63;
const PARENT_MASK: u64 = WITNESS - 1;
const MAX_ID: u64 = (1 << 62) - 2;
fn priority(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
fn parent(id: u64, word: u64) -> u64 {
    let p = (word & PARENT_MASK) >> 1;
    if p == 0 {
        id
    } else {
        p - 1
    }
}
fn parity(word: u64) -> u64 {
    word & 1
}
fn link(parent: u64, reverse: u64, old: u64) -> u64 {
    ((parent + 1) << 1) | reverse | (old & WITNESS)
}

struct NodePage {
    values: Box<[AtomicU64]>,
    dirty: AtomicBool,
}
impl NodePage {
    fn new(size: usize, bytes: Option<Vec<u8>>) -> io::Result<Self> {
        let values = if let Some(bytes) = bytes {
            if bytes.len() != size {
                return Err(invalid("invalid parent page size"));
            }
            bytes
                .chunks_exact(8)
                .map(|v| AtomicU64::new(u64::from_le_bytes(v.try_into().unwrap())))
                .collect::<Vec<_>>()
                .into_boxed_slice()
        } else {
            (0..size / 8)
                .map(|_| AtomicU64::new(0))
                .collect::<Vec<_>>()
                .into_boxed_slice()
        };
        Ok(Self {
            values,
            dirty: AtomicBool::new(false),
        })
    }
    fn bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.values.len() * 8);
        for v in &self.values {
            out.extend_from_slice(&v.load(Ordering::Relaxed).to_le_bytes());
        }
        out
    }
}
struct Cache {
    hot: HashMap<u64, Arc<NodePage>>,
    clock: VecDeque<u64>,
    current: BTreeMap<u64, Page>,
    changed: BTreeMap<u64, Page>,
    pending: HashSet<u64>,
    loading: usize,
}
pub(super) struct Parents {
    count: u64,
    pub page_bytes: usize,
    slots: usize,
    pages: Pages,
    cache: Mutex<Cache>,
    ready: Condvar,
    limit: usize,
    l0: usize,
}
impl Parents {
    pub fn open(
        dir: &Path,
        count: u64,
        budget: usize,
        previous: &[(u64, Page)],
    ) -> io::Result<Self> {
        if count > MAX_ID {
            return Err(invalid("too many homology atoms for packed parent IDs"));
        }
        let workers = rayon::current_num_threads().max(1);
        let desired = (budget / workers / 8).clamp(4096, MAX_PAGE_BYTES);
        let page_bytes = 1usize << (usize::BITS - 1 - desired.leading_zeros());
        Self::open_layout(dir, count, budget, previous, page_bytes)
    }
    pub fn open_layout(
        dir: &Path,
        count: u64,
        budget: usize,
        previous: &[(u64, Page)],
        page_bytes: usize,
    ) -> io::Result<Self> {
        if page_bytes < 4096 || page_bytes > MAX_PAGE_BYTES || !page_bytes.is_power_of_two() {
            return Err(invalid("invalid parent page layout"));
        }
        let workers = rayon::current_num_threads().max(1);
        let limit = budget / page_bytes;
        if limit < workers + 2 {
            return Err(invalid(
                "parent memory budget cannot hold one active page per worker",
            ));
        }
        let l0 = (limit / workers / 4).clamp(1, 64);
        Ok(Self {
            count,
            page_bytes,
            slots: page_bytes / 8,
            pages: Pages::open(dir, 1 << 20)?,
            cache: Mutex::new(Cache {
                hot: HashMap::default(),
                clock: VecDeque::new(),
                current: previous.iter().copied().collect(),
                changed: BTreeMap::new(),
                pending: HashSet::default(),
                loading: 0,
            }),
            ready: Condvar::new(),
            limit,
            l0,
        })
    }
    fn page(&self, id: u64) -> io::Result<Arc<NodePage>> {
        if id >= self.count {
            return Err(invalid("parent ID outside atom partition"));
        }
        let number = id / self.slots as u64;
        let mut c = self.cache.lock();
        while c.pending.contains(&number) {
            self.ready.wait(&mut c);
        }
        if let Some(p) = c.hot.get(&number) {
            return Ok(p.clone());
        }
        let mut attempts = c.clock.len();
        let mut evicted = None;
        while c.hot.len() + c.loading >= self.limit && attempts > 0 {
            attempts -= 1;
            let Some(old) = c.clock.pop_front() else {
                break;
            };
            if c.hot.get(&old).is_some_and(|p| Arc::strong_count(p) == 1) {
                let page = c.hot.remove(&old).unwrap();
                c.pending.insert(old);
                evicted = Some((old, page));
                break;
            } else {
                c.clock.push_back(old);
            }
        }
        if c.hot.len() + c.loading >= self.limit {
            return Err(invalid("parent page budget exhausted by pinned readers"));
        }
        let previous = c.current.get(&number).copied();
        c.pending.insert(number);
        c.loading += 1;
        drop(c);

        // Only callers for these particular pages wait for disk I/O. Keeping
        // an evicted page pending prevents another loader reading its stale
        // checkpoint while the latest atomics are being written back.
        // Each loader retains at most one evicted node page in addition to
        // its reserved cache slot; serialization buffers are page sized too.
        let written = evicted
            .as_ref()
            .filter(|(_, p)| p.dirty.load(Ordering::Relaxed))
            .map(|(old, p)| self.pages.put_raw((0, *old), &p.bytes()))
            .transpose();
        let loaded = match &written {
            Ok(_) => previous
                .map(|r| self.pages.read_raw(r))
                .transpose()
                .and_then(|raw| NodePage::new(self.page_bytes, raw))
                .map(Arc::new),
            Err(e) => Err(io::Error::new(e.kind(), e.to_string())),
        };

        let mut c = self.cache.lock();
        c.loading -= 1;
        c.pending.remove(&number);
        if let Some((old, page)) = evicted {
            c.pending.remove(&old);
            match written {
                Ok(Some(r)) => {
                    c.current.insert(old, r);
                    c.changed.insert(old, r);
                }
                Ok(None) => {}
                Err(_) => {
                    // A failed write must not discard the only current copy.
                    c.hot.insert(old, page);
                    c.clock.push_back(old);
                }
            }
        }
        if let Ok(page) = &loaded {
            c.hot.insert(number, page.clone());
            c.clock.push_back(number);
        }
        self.ready.notify_all();
        loaded
    }
    pub fn view(&self) -> View<'_> {
        View {
            store: self,
            pages: HashMap::default(),
            clock: VecDeque::new(),
            links: 0,
            probes: 0,
        }
    }
    /// Only called after workers have joined. All changes are appended as new
    /// immutable pages; a crash cannot mutate an earlier committed generation.
    pub fn checkpoint(&self) -> io::Result<Vec<(u64, Page)>> {
        let mut c = self.cache.lock();
        let dirty: Vec<_> = c
            .hot
            .iter()
            .filter_map(|(&id, page)| {
                page.dirty
                    .swap(false, Ordering::Relaxed)
                    .then(|| (id, page.clone()))
            })
            .collect();
        let updates = dirty
            .into_par_iter()
            .map(|(id, page)| self.pages.put_raw((0, id), &page.bytes()).map(|p| (id, p)))
            .collect::<io::Result<Vec<_>>>()?;
        for (id, r) in updates {
            c.current.insert(id, r);
            c.changed.insert(id, r);
        }
        self.pages.sync()?;
        Ok(std::mem::take(&mut c.changed).into_iter().collect())
    }
}
pub(super) struct View<'a> {
    store: &'a Parents,
    pages: HashMap<u64, Arc<NodePage>>,
    clock: VecDeque<u64>,
    pub links: u64,
    pub probes: u64,
}
impl View<'_> {
    #[inline]
    fn page(&mut self, id: u64) -> io::Result<&NodePage> {
        let n = id / self.store.slots as u64;
        if !self.pages.contains_key(&n) {
            if self.pages.len() >= self.store.l0 {
                if let Some(old) = self.clock.pop_front() {
                    self.pages.remove(&old);
                }
            }
            self.pages.insert(n, self.store.page(id)?);
            self.clock.push_back(n);
        }
        Ok(self.pages.get(&n).unwrap())
    }
    #[inline]
    fn load(&mut self, id: u64) -> io::Result<u64> {
        let offset = id as usize % self.store.slots;
        Ok(self.page(id)?.values[offset].load(Ordering::SeqCst))
    }
    fn cas(&mut self, id: u64, old: u64, new: u64) -> io::Result<bool> {
        let offset = id as usize % self.store.slots;
        let p = self.page(id)?;
        let yes = p.values[offset]
            .compare_exchange(old, new, Ordering::SeqCst, Ordering::SeqCst)
            .is_ok();
        if yes {
            if !p.dirty.load(Ordering::Relaxed) {
                p.dirty.store(true, Ordering::Relaxed);
            }
        }
        Ok(yes)
    }
    pub fn find(&mut self, id: u64) -> io::Result<(u64, u64)> {
        let (mut v, mut sign) = (id, 0);
        loop {
            let old = self.load(v)?;
            let p = parent(v, old);
            if p == v {
                return Ok((v, sign));
            }
            let p_word = self.load(p)?;
            let grand = parent(p, p_word);
            if grand != p {
                self.cas(v, old, link(grand, parity(old) ^ parity(p_word), old))?;
            }
            sign ^= parity(old);
            v = p;
        }
    }
    pub fn witness(&mut self, id: u64) -> io::Result<bool> {
        Ok(self.load(id)? & WITNESS != 0)
    }
    pub fn mark_witness(&mut self, id: u64) -> io::Result<()> {
        loop {
            let old = self.load(id)?;
            if old & WITNESS != 0 {
                return Ok(());
            }
            if self.cas(id, old, old | WITNESS)? {
                return Ok(());
            }
        }
    }
    pub fn join(&mut self, a: u64, b: u64, reverse: u64) -> io::Result<()> {
        self.probes += 1;
        loop {
            let ((mut x, p), (mut y, q)) = (self.find(a)?, self.find(b)?);
            let sign = p ^ q ^ reverse;
            if x == y {
                if sign != 0 {
                    self.mark_witness(x)?;
                }
                return Ok(());
            }
            if (priority(x), x) > (priority(y), y) {
                std::mem::swap(&mut x, &mut y);
            }
            let old = self.load(x)?;
            if parent(x, old) != x {
                continue;
            }
            if self.cas(x, old, link(y, sign, old))? {
                self.links += 1;
                return Ok(());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;
    #[test]
    fn concurrent_eviction_preserves_parities_and_checkpoint() {
        // TEST signed star larger than the eight-page cache. A scrambled
        // access order forces overlapping page reads and dirty writebacks.
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| {
                let dir = tempfile::tempdir().unwrap();
                let count = 8193;
                let p = Parents::open_layout(dir.path(), count, 8 * 4096, &[], 4096).unwrap();
                (1..count).into_par_iter().for_each_init(
                    || p.view(),
                    |v, i| {
                        let id = ((i - 1) * 4051) % (count - 1) + 1;
                        v.join(0, id, id % 2).unwrap();
                    },
                );
                let committed = p.checkpoint().unwrap();
                assert!(committed.len() > 8);
                drop(p);
                let restored =
                    Parents::open_layout(dir.path(), count, 8 * 4096, &committed, 4096).unwrap();
                let (root, sign) = restored.view().find(0).unwrap();
                (0..count).into_par_iter().for_each_init(
                    || restored.view(),
                    |v, id| {
                        let (r, s) = v.find(id).unwrap();
                        assert_eq!(r, root);
                        assert_eq!(s ^ sign, id % 2);
                        assert!(!v.witness(r).unwrap());
                    },
                );
            });
    }
    #[test]
    fn concurrent_test_parities_and_copy_on_write_recovery() {
        // TEST graph: one signed chain plus its reversing-cycle witness.
        let dir = tempfile::tempdir().unwrap();
        let p = Parents::open(dir.path(), 1000, 32 << 20, &[]).unwrap();
        (0..999u64)
            .into_par_iter()
            .for_each(|i| p.view().join(i, i + 1, i % 2).unwrap());
        let (mut v, mut expected) = (p.view(), 0);
        for i in 0..1000 {
            let (root, s) = v.find(i).unwrap();
            assert_eq!(root, v.find(0).unwrap().0);
            assert_eq!(s ^ v.find(0).unwrap().1, expected);
            expected ^= i % 2;
        }
        let committed = p.checkpoint().unwrap();
        p.view().join(0, 0, 1).unwrap();
        assert!(p.view().witness(p.view().find(0).unwrap().0).unwrap());
        drop(v);
        drop(p);
        let restored = Parents::open(dir.path(), 1000, 32 << 20, &committed).unwrap();
        assert!(!restored
            .view()
            .witness(restored.view().find(0).unwrap().0)
            .unwrap());
        for i in 0..1000 {
            assert_eq!(
                restored.view().find(i).unwrap().0,
                restored.view().find(0).unwrap().0
            );
        }
    }
}
