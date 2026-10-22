//! Incremental, exact longest-missing-range selection. No biological filtering.
//!
//! The legacy selector reduces by (length, sequence ID), taking the last range
//! on a within-sequence tie. SortedRanges is ordered by start, hence the exact
//! priority is (length, sequence ID, start, end). Old heap entries are validated
//! against the authoritative missing map, never trusted as current state.
use crate::impg::SortedRanges;
use rustc_hash::FxHashMap;
use std::collections::BinaryHeap;

#[derive(Default)]
pub(super) struct LongestMissing {
    heap: BinaryHeap<(i32, u32, i32, i32)>,
    live_ranges: usize,
}

impl LongestMissing {
    pub(super) fn new(missing: &FxHashMap<u32, SortedRanges>) -> Self {
        let entries: Vec<_> = missing.iter().flat_map(|(&id, ranges)| {
            ranges.iter().map(move |&(s, e)| (e - s, id, s, e))
        }).collect();
        Self { live_ranges: entries.len(), heap: BinaryHeap::from(entries) }
    }

    pub(super) fn replace(&mut self, id: u32, old: &[(i32, i32)], new: &[(i32, i32)]) {
        self.live_ranges = self.live_ranges.checked_sub(old.len()).expect("longest cache range count") + new.len();
        // Entries that remain present are already in the heap. Avoid duplicate
        // insertion when a large mask span includes untouched missing ranges.
        for &(s, e) in new {
            if old.binary_search_by_key(&s, |&(start, _)| start)
                .map(|i| old[i].1 == e).unwrap_or(false) { continue; }
            self.heap.push((e - s, id, s, e));
        }
    }

    pub(super) fn select(&mut self, missing: &FxHashMap<u32, SortedRanges>) -> Option<(u32, i32, i32)> {
        // Bound stale-entry memory without rescanning the map per query.
        if self.heap.len() > self.live_ranges.saturating_mul(4).saturating_add(1_000_000) {
            *self = Self::new(missing);
        }
        while let Some(&(_, id, s, e)) = self.heap.peek() {
            let valid = missing.get(&id).map(|r| {
                r.ranges.binary_search_by_key(&s, |&(start, _)| start)
                    .map(|i| r.ranges[i].1 == e).unwrap_or(false)
            }).unwrap_or(false);
            if valid { return Some((id, s, e)); }
            self.heap.pop();
        }
        debug_assert_eq!(self.live_ranges, 0);
        None
    }
}
