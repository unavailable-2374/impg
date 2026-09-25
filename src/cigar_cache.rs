//! Process-wide cache of parsed PAF CIGARs, keyed by alignment record.
//!
//! Why: every projection re-reads its record's CIGAR from the PAF -- for a
//! BGZF file that is open() + header check + seek + inflate + parse, per
//! record, per query. `impg partition` walks each seed sequence window by
//! window, so a 1 Mb alignment is fetched again for each of the ~200 5-kb
//! windows it spans, and on Lustre every open() is a metadata round trip. At
//! 167k per-file PAFs x 580 records per window that, not the interval trees,
//! is where a genome-scale partition spends its time.
//!
//! Output cannot change: a record's CIGAR is a pure function of (file, byte
//! offset, length, strand/reversal bits), and exactly that is the key, so a
//! hit returns the ops the uncached path would have parsed. The file is
//! identified by two independent 64-bit hashes of its path (a false hit needs
//! both to collide AND the same offset and length).
//!
//! Enabled by IMPG_CIGAR_CACHE_GB=<budget>; unset or 0 means off, and every
//! call takes the uncached path unchanged. Eviction is FIFO per shard: seeds
//! advance along a sequence, so a record stops being hit soon after the walk
//! passes it, and recency tracking would buy little for its locking cost.

use crate::impg::CigarOp;
use log::info;
use rustc_hash::FxHashMap;
use std::collections::VecDeque;
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

const SHARDS: usize = 512;
/// Per-entry bookkeeping on top of the ops themselves (key twice, Arc, Vec
/// header, map slot); counted so the budget is a budget of real memory.
const ENTRY_OVERHEAD: usize = 96;

type Key = (u64, u64, u64, u32);

struct Shard {
    map: FxHashMap<Key, Arc<Vec<CigarOp>>>,
    order: VecDeque<(Key, usize)>,
    bytes: usize,
}

pub struct CigarCache {
    shards: Vec<Mutex<Shard>>,
    shard_budget: usize,
    hits: AtomicU64,
    misses: AtomicU64,
    evictions: AtomicU64,
    bytes: AtomicU64,
}

static CACHE: OnceLock<Option<CigarCache>> = OnceLock::new();

pub fn global() -> Option<&'static CigarCache> {
    CACHE
        .get_or_init(|| {
            let gb: f64 = std::env::var("IMPG_CIGAR_CACHE_GB")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0);
            if gb <= 0.0 {
                return None;
            }
            let budget = (gb * 1024.0 * 1024.0 * 1024.0) as usize;
            info!("CIGAR cache: {gb} GB budget over {SHARDS} shards");
            Some(CigarCache {
                shards: (0..SHARDS)
                    .map(|_| {
                        Mutex::new(Shard {
                            map: FxHashMap::default(),
                            order: VecDeque::new(),
                            bytes: 0,
                        })
                    })
                    .collect(),
                shard_budget: budget / SHARDS,
                hits: AtomicU64::new(0),
                misses: AtomicU64::new(0),
                evictions: AtomicU64::new(0),
                bytes: AtomicU64::new(0),
            })
        })
        .as_ref()
}

fn path_hashes(path: &str) -> (u64, u64) {
    let mut fx = rustc_hash::FxHasher::default();
    path.hash(&mut fx);
    // FNV-1a: independent of FxHash, so the pair behaves as a 128-bit id.
    let mut fnv: u64 = 0xcbf2_9ce4_8422_2325;
    for b in path.as_bytes() {
        fnv ^= *b as u64;
        fnv = fnv.wrapping_mul(0x0000_0100_0000_01b3);
    }
    (fx.finish(), fnv)
}

impl CigarCache {
    /// `record` is QueryMetadata's raw strand/reversed/offset word: it carries
    /// the orientation bits, so the two directions of a bidirectional record
    /// (whose ops differ by inversion) never share an entry.
    pub fn get_or_insert_with<F>(&self, path: &str, record: u64, len: u32, parse: F) -> Arc<Vec<CigarOp>>
    where
        F: FnOnce() -> Vec<CigarOp>,
    {
        let (h1, h2) = path_hashes(path);
        let key: Key = (h1, h2, record, len);
        let shard_idx = {
            let mut h = rustc_hash::FxHasher::default();
            key.hash(&mut h);
            (h.finish() as usize) % SHARDS
        };
        if let Some(v) = self.shards[shard_idx].lock().unwrap().map.get(&key) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Arc::clone(v);
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        // Parse outside the lock: two threads missing the same key both parse
        // it and insert identical values -- wasteful but never wrong.
        let ops = Arc::new(parse());
        let size = ops.len() * std::mem::size_of::<CigarOp>() + ENTRY_OVERHEAD;
        if size > self.shard_budget {
            return ops; // larger than a whole shard: serve it, do not cache it
        }
        let mut shard = self.shards[shard_idx].lock().unwrap();
        if shard.map.contains_key(&key) {
            return ops;
        }
        while shard.bytes + size > self.shard_budget {
            let Some((old, old_size)) = shard.order.pop_front() else { break };
            shard.map.remove(&old);
            shard.bytes -= old_size;
            self.bytes.fetch_sub(old_size as u64, Ordering::Relaxed);
            self.evictions.fetch_add(1, Ordering::Relaxed);
        }
        shard.map.insert(key, Arc::clone(&ops));
        shard.order.push_back((key, size));
        shard.bytes += size;
        self.bytes.fetch_add(size as u64, Ordering::Relaxed);
        ops
    }

    pub fn stats_line(&self) -> String {
        let h = self.hits.load(Ordering::Relaxed);
        let m = self.misses.load(Ordering::Relaxed);
        format!(
            "CIGAR cache: {} hits, {} misses ({:.1}% hit), {} evictions, {:.2} GB held",
            h,
            m,
            if h + m > 0 { 100.0 * h as f64 / (h + m) as f64 } else { 0.0 },
            self.evictions.load(Ordering::Relaxed),
            self.bytes.load(Ordering::Relaxed) as f64 / (1u64 << 30) as f64
        )
    }
}
