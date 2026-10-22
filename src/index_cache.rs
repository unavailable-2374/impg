//! Soft byte-budgeted LRU. Eviction discards only reloadable index data.
//! In-flight Arc references are pinned; the budget is not a process RSS limit.
use std::collections::BTreeSet;
use std::sync::Arc;

struct Entry<T> {
    value: Arc<T>,
    bytes: usize,
    stamp: u64,
}

pub(crate) struct IndexCache<T> {
    entries: Vec<Option<Entry<T>>>,
    order: BTreeSet<(u64, usize)>,
    clock: u64,
    pub(crate) budget: Option<usize>,
    pub(crate) bytes: usize,
    pub(crate) hits: u64,
    pub(crate) misses: u64,
    pub(crate) evictions: u64,
}

impl<T> IndexCache<T> {
    pub(crate) fn new(len: usize, budget: Option<usize>) -> Self {
        Self { entries: (0..len).map(|_| None).collect(), order: BTreeSet::new(),
            clock: 0, budget, bytes: 0, hits: 0, misses: 0, evictions: 0 }
    }

    pub(crate) fn get(&mut self, id: usize) -> Option<Arc<T>> {
        match self.entries[id].as_mut() {
            Some(entry) => {
                self.hits += 1;
                if self.budget.is_none() { return Some(Arc::clone(&entry.value)); }
                self.order.remove(&(entry.stamp, id));
                self.clock += 1;
                entry.stamp = self.clock;
                self.order.insert((entry.stamp, id));
                Some(Arc::clone(&entry.value))
            }
            None => { self.misses += 1; None }
        }
    }

    // Loading is outside the lock. If another thread won the race, reuse its
    // instance, never replace a live cached object with a duplicate.
    pub(crate) fn insert(&mut self, id: usize, value: Arc<T>, bytes: usize) -> Arc<T> {
        if self.entries[id].is_some() { return self.get(id).unwrap(); }
        self.clock += 1;
        self.bytes += bytes;
        self.order.insert((self.clock, id));
        self.entries[id] = Some(Entry { value: Arc::clone(&value), bytes, stamp: self.clock });
        value
    }

    pub(crate) fn update(&mut self, id: usize, value: &Arc<T>, bytes: usize) {
        if let Some(entry) = self.entries[id].as_mut() {
            if Arc::ptr_eq(&entry.value, value) {
                self.bytes = self.bytes - entry.bytes + bytes;
                entry.bytes = bytes;
            }
        }
    }

    pub(crate) fn trim(&mut self) {
        let Some(budget) = self.budget else { return; };
        while self.bytes > budget {
            let victim = self.order.iter().copied().find(|&(_, id)| {
                Arc::strong_count(&self.entries[id].as_ref().unwrap().value) == 1
            });
            let Some((stamp, id)) = victim else { break; };
            self.order.remove(&(stamp, id));
            let entry = self.entries[id].take().unwrap();
            self.bytes -= entry.bytes;
            self.evictions += 1;
        }
    }

    pub(crate) fn len(&self) -> usize { self.order.len() }

    pub(crate) fn clear_unpinned(&mut self) {
        let victims: Vec<_> = self.order.iter().copied().filter(|&(_, id)| {
            Arc::strong_count(&self.entries[id].as_ref().unwrap().value) == 1
        }).collect();
        for (stamp, id) in victims {
            self.order.remove(&(stamp, id));
            self.bytes -= self.entries[id].take().unwrap().bytes;
            self.evictions += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lru_budget_and_pinning() {
        let mut c = IndexCache::new(4, Some(20));
        drop(c.insert(0, Arc::new(0), 10));
        drop(c.insert(1, Arc::new(1), 10));
        let pinned = c.get(0).unwrap();
        drop(c.insert(2, Arc::new(2), 20));
        c.trim();
        assert_eq!(c.bytes, 10);
        assert_eq!(c.evictions, 2);
        assert_eq!(*pinned, 0);
        c.update(0, &pinned, 40);
        c.trim();
        assert_eq!(c.bytes, 40); // in-flight entry cannot be evicted
        drop(pinned);
        c.trim();
        assert_eq!(c.bytes, 0);
    }
    #[test]
    fn touch_changes_victim_and_duplicate_load_keeps_original() {
        let mut c = IndexCache::new(3, Some(20));
        drop(c.insert(0, Arc::new(7), 10));
        drop(c.insert(1, Arc::new(8), 10));
        drop(c.get(0));
        drop(c.insert(2, Arc::new(9), 10));
        c.trim();
        assert!(c.get(1).is_none());
        assert_eq!(*c.insert(0, Arc::new(99), 100), 7);
        assert_eq!(c.bytes, 20);
    }
    #[test]
    fn disabled_cache_and_unbounded_cache() {
        let mut c = IndexCache::new(1, Some(0));
        let value = c.insert(0, Arc::new(1), 50);
        c.trim();
        assert_eq!(c.len(), 1);
        drop(value);
        c.trim();
        assert_eq!(c.len(), 0);
        let mut c = IndexCache::new(1, None);
        drop(c.insert(0, Arc::new(1), 50));
        c.trim();
        assert_eq!(c.bytes, 50);
        c.clear_unpinned();
        assert_eq!(c.bytes, 0);
    }
}
