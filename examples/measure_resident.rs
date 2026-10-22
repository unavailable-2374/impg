//! Empirically calibrate `SUB_INDEX_RESIDENT_EXPANSION` (the on-disk → resident
//! multiplier used by the byte-budgeted sub-index cache in `MultiImpg`).
//!
//! It replicates the binary's jemalloc allocator + config so RSS behaves exactly
//! as in a real `impg depth` run, then mirrors the Phase-1 holding pattern: load
//! each per-file `.impg` via the same `Impg::load_from_file` path `get_sub_index`
//! uses, enable tree caching, and query the *ref* target's full range (the depth
//! command queries only the anchored ref target per neighbour file, so only that
//! target's COITree becomes resident). The loaded `Arc<Impg>`s are held — exactly
//! what the cache pins — and RSS is sampled before/after to derive the real
//! resident/on-disk ratio and its breakdown (base header vs cached trees).
//!
//! Usage:
//!   cargo run --release --example measure_resident -- <impg-list> <ref-prefix> [N]
//!     <impg-list>   newline-delimited file of `.impg` paths
//!     <ref-prefix>  only targets whose sequence name starts with this are
//!                   queried (the depth anchor sample, e.g. GCA_009914755.4).
//!                   Pass "" to query every target (worst case: all trees).
//!     [N]           optional cap on number of files sampled

// Replicate the binary's global allocator so RSS matches production. (Examples
// don't inherit main.rs's allocator; dependencies of the crate are available.)
#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(feature = "jemalloc")]
const JEMALLOC_MALLOC_CONF: &[u8] =
    b"narenas:8,dirty_decay_ms:1000,muzzy_decay_ms:0,metadata_thp:auto,abort_conf:false\0";

#[cfg(feature = "jemalloc")]
#[allow(non_upper_case_globals)]
#[export_name = "_rjem_malloc_conf"]
pub static malloc_conf: Option<&'static std::ffi::c_char> = Some(unsafe {
    union U {
        bytes: &'static u8,
        c: &'static std::ffi::c_char,
    }
    U {
        bytes: &JEMALLOC_MALLOC_CONF[0],
    }
    .c
});

use impg::impg::Impg;
use impg::impg_index::ImpgIndex;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::Arc;

/// Resident set size in bytes from /proc/self/statm (field 2 = resident pages).
fn rss_bytes() -> u64 {
    let s = std::fs::read_to_string("/proc/self/statm").unwrap_or_default();
    let pages: u64 = s
        .split_whitespace()
        .nth(1)
        .and_then(|f| f.parse().ok())
        .unwrap_or(0);
    pages * 4096
}

/// Sample RSS after letting jemalloc's dirty-page decay settle, so transient
/// per-query CIGAR buffers (freed immediately) don't inflate the reading.
fn settled_rss() -> u64 {
    std::thread::sleep(std::time::Duration::from_millis(1500));
    let mut last = rss_bytes();
    for _ in 0..6 {
        std::thread::sleep(std::time::Duration::from_millis(500));
        let now = rss_bytes();
        if now >= last {
            return now;
        }
        last = now;
    }
    last
}

fn mb(b: u64) -> f64 {
    b as f64 / (1u64 << 20) as f64
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: {} <impg-list> <ref-prefix|\"\"> [N]",
            args.first()
                .map(String::as_str)
                .unwrap_or("measure_resident")
        );
        std::process::exit(2);
    }
    let list_path = &args[1];
    let ref_prefix = args[2].clone();
    let cap: usize = args
        .get(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);

    let mut paths: Vec<PathBuf> = std::fs::read_to_string(list_path)
        .expect("read impg list")
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(PathBuf::from)
        .collect();
    paths.truncate(cap);
    assert!(!paths.is_empty(), "no .impg paths in list");

    let on_disk: u64 = paths
        .iter()
        .map(|p| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0))
        .sum();

    println!(
        "measuring {} files, {:.1} MB on disk, ref-prefix={:?}",
        paths.len(),
        mb(on_disk),
        if ref_prefix.is_empty() {
            "<all targets>"
        } else {
            &ref_prefix
        }
    );

    let rss0 = settled_rss();

    // Phase A: load every per-file index (header only — no trees yet), hold Arcs.
    let mut held: Vec<Arc<Impg>> = Vec::with_capacity(paths.len());
    for p in &paths {
        let f = File::open(p).expect("open .impg");
        let af = vec![p.to_string_lossy().into_owned()];
        let impg = Impg::load_from_file(
            BufReader::new(f),
            &af,
            p.to_string_lossy().into_owned(),
            None,
        )
        .expect("load_from_file");
        impg.set_tree_cache_enabled(true);
        held.push(Arc::new(impg));
    }
    let rss_loaded = settled_rss();

    // Phase B: load the ref target(s)' COITrees into each Arc's tree cache via
    // `get_or_load_tree` — exactly what the depth hop-0 sweep pins per neighbour
    // file. This loads the persistent resident (the tree's QueryMetadata array)
    // WITHOUT the transient per-query CIGAR (read from the PAF and dropped, never
    // cached), so the reading isolates what the cache actually holds.
    let mut queried_targets = 0usize;
    for impg in &held {
        for tid in impg.target_ids() {
            let name = impg.seq_index.get_name(tid).unwrap_or("");
            if !ref_prefix.is_empty() && !name.starts_with(&ref_prefix) {
                continue;
            }
            let _ = impg.get_or_load_tree(tid);
            queried_targets += 1;
        }
    }
    let rss_queried = settled_rss();

    let base = rss_loaded.saturating_sub(rss0);
    let trees = rss_queried.saturating_sub(rss_loaded);
    let total = rss_queried.saturating_sub(rss0);

    println!("--- resident breakdown (jemalloc, settled) ---");
    println!("ref targets queried : {queried_targets}");
    println!(
        "header/base         : {:>9.1} MB  ({:.3}x on-disk, {:.3} MB/file)",
        mb(base),
        base as f64 / on_disk as f64,
        mb(base) / held.len() as f64
    );
    println!(
        "cached trees        : {:>9.1} MB  ({:.3}x on-disk, {:.3} MB/file)",
        mb(trees),
        trees as f64 / on_disk as f64,
        mb(trees) / held.len() as f64
    );
    println!(
        "TOTAL resident      : {:>9.1} MB  ({:.3}x on-disk, {:.3} MB/file)",
        mb(total),
        total as f64 / on_disk as f64,
        mb(total) / held.len() as f64
    );
    println!(
        ">>> measured resident/on-disk multiplier = {:.3}x  (SUB_INDEX_RESIDENT_EXPANSION currently 2)",
        total as f64 / on_disk as f64
    );

    // Keep the cache resident across the final measurement.
    std::hint::black_box(&held);
}
