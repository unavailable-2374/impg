//! Read-only reconstruction of cached tree payloads from real production logs.
use coitrees::IntervalNode;
use impg::impg::{Impg, QueryMetadata};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::Deserialize;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, SeekFrom};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Default)]
struct Stats {
    files: u64,
    loaded_files: u64,
    trees: u64,
    nodes: u64,
    metadata_sequences: u64,
    metadata_targets: u64,
    sampled_files: u64,
    sampled_nodes: u64,
    bins: [u64; 6],
}
impl Stats {
    fn add(mut self, b: Self) -> Self {
        self.files += b.files;
        self.loaded_files += b.loaded_files;
        self.trees += b.trees;
        self.nodes += b.nodes;
        self.metadata_sequences += b.metadata_sequences;
        self.metadata_targets += b.metadata_targets;
        self.sampled_files += b.sampled_files;
        self.sampled_nodes += b.sampled_nodes;
        for i in 0..6 {
            self.bins[i] += b.bins[i];
        }
        self
    }
}
#[derive(Deserialize)]
struct DiskInterval {
    first: i32,
    last: i32,
    #[allow(dead_code)]
    metadata: QueryMetadata,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err(
            "usage: profile_partition_cache production.log alignments.list replayed_window_count"
                .into(),
        );
    }
    let cutoff: u64 = args[3].parse()?;
    let start = Instant::now();
    let mut targets = FxHashSet::<String>::default();
    let mut reader = BufReader::with_capacity(4 << 20, File::open(&args[1])?);
    let mut line = String::new();
    let mut count = 0u64;
    let mut short_windows = [0u64; 6];
    let mut querying = false;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.contains("windows, now querying") {
            let expected = format!("replayed {cutoff} windows, now querying");
            if !line.contains(&expected) {
                return Err("unexpected production replay boundary".into());
            }
            querying = true;
        }
        if !querying {
            continue;
        }
        let Some(pos) = line.find("Computed partition") else {
            continue;
        };
        let rest = &line[pos + "Computed partition".len()..];
        let partition: u64 = rest
            .split_whitespace()
            .next()
            .ok_or("missing partition number")?
            .parse()?;
        if partition < cutoff {
            return Err("partition number precedes replay boundary".into());
        }
        let pos = line.find("(query ").ok_or("missing query window")? + 7;
        let end = line[pos..].find(", len: ").ok_or("missing query length")? + pos;
        let region = &line[pos..end];
        let colon = region.rfind(':').ok_or("missing coordinate separator")?;
        let name = &region[..colon];
        if !targets.contains(name) {
            targets.insert(name.to_owned());
        }
        let length: u64 = line[end + 7..].trim_end().trim_end_matches(')').parse()?;
        short_windows[bin(length)] += 1;
        count += 1;
    }
    if !querying || targets.is_empty() {
        return Err("no executed production queries found".into());
    }
    println!("log={} cutoff={} logged_query_windows={} distinct_queried_targets={} elapsed_s={:.1} window_bins={:?}",args[1],cutoff,count,targets.len(),start.elapsed().as_secs_f64(),short_windows);
    let files: Vec<_> = std::fs::read_to_string(&args[2])?
        .lines()
        .filter(|x| !x.is_empty())
        .map(str::to_owned)
        .collect();
    let done = AtomicUsize::new(0);
    let total: Stats = files
        .par_iter()
        .enumerate()
        .map(|(idx, paf)| -> io::Result<Stats> {
            let path = format!("{paf}.impg");
            let index = Impg::load_from_file(
                BufReader::new(File::open(&path)?),
                std::slice::from_ref(paf),
                path.clone(),
                None,
            )?;
            let mut input = BufReader::with_capacity(
                if idx % 1000 == 0 { 512 << 10 } else { 256 },
                File::open(&path)?,
            );
            let mut s = Stats {
                files: 1,
                ..Default::default()
            };
            let sampled = idx % 1000 == 0;
            for (&target, &offset) in &index.forest_map.entries {
                let name = index
                    .seq_index
                    .get_name(target)
                    .ok_or_else(|| io::Error::other("target missing from dictionary"))?;
                if !targets.contains(name) {
                    continue;
                }
                input.seek(SeekFrom::Start(offset))?;
                let (actual, n): (u32, u64) =
                    bincode::serde::decode_from_std_read(&mut input, bincode::config::standard())
                        .map_err(io::Error::other)?;
                if actual != target {
                    return Err(io::Error::other(format!(
                        "{path}: tree ID mismatch {actual} != {target}"
                    )));
                }
                s.trees += 1;
                s.nodes += n;
                if sampled {
                    for _ in 0..n {
                        let iv: DiskInterval = bincode::serde::decode_from_std_read(
                            &mut input,
                            bincode::config::standard(),
                        )
                        .map_err(io::Error::other)?;
                        let length = (i64::from(iv.last) - i64::from(iv.first)).unsigned_abs();
                        s.bins[bin(length)] += 1;
                        s.sampled_nodes += 1;
                    }
                }
            }
            if s.trees > 0 {
                s.loaded_files = 1;
                s.metadata_sequences = index.seq_index.len() as u64;
                s.metadata_targets = index.forest_map.entries.len() as u64;
                if sampled {
                    s.sampled_files = 1;
                }
            }
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            if n % 5000 == 0 {
                eprintln!(
                    "read_only_index_scan files={n}/{} elapsed_s={:.1}",
                    files.len(),
                    start.elapsed().as_secs_f64()
                );
            }
            Ok(s)
        })
        .try_reduce(Stats::default, |a, b| Ok(a.add(b)))?;
    let node_size = std::mem::size_of::<IntervalNode<QueryMetadata, u32>>();
    println!("files={} inferred_loaded_subindices={} inferred_cached_trees={} inferred_cached_nodes={} node_size={} node_payload_bytes={} node_payload_TiB={:.6} metadata_sequences={} metadata_targets={} elapsed_s={:.1}",total.files,total.loaded_files,total.trees,total.nodes,node_size,total.nodes*u64::try_from(node_size)?,(total.nodes as f64*node_size as f64)/(1024f64.powi(4)),total.metadata_sequences,total.metadata_targets,start.elapsed().as_secs_f64());
    println!("SYSTEMATIC_REAL_INDEX_SAMPLE stride=1000 sampled_files={} sampled_cached_nodes={} length_bins_le100_le500_le1000_le5000_le20000_gt20000={:?}",total.sampled_files,total.sampled_nodes,total.bins);
    Ok(())
}
fn bin(n: u64) -> usize {
    if n <= 100 {
        0
    } else if n <= 500 {
        1
    } else if n <= 1000 {
        2
    } else if n <= 5000 {
        3
    } else if n <= 20000 {
        4
    } else {
        5
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::IntervalTree;
    #[test]
    fn real_index_header_count_matches_loaded_tree() {
        let paf="/scratch/10779/shuocao2374/VGP/ILS/partition/verify/fresh/GCA_965226575.1_vs_GCF_017976325.1.paf.gz";
        let path = format!("{paf}.impg");
        let index = Impg::load_from_file(
            BufReader::new(File::open(&path).unwrap()),
            &[paf.to_owned()],
            path.clone(),
            None,
        )
        .unwrap();
        let target = *index.forest_map.entries.keys().min().unwrap();
        let mut reader = BufReader::new(File::open(path).unwrap());
        reader
            .seek(SeekFrom::Start(index.forest_map.entries[&target]))
            .unwrap();
        let (id, count): (u32, u64) =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard()).unwrap();
        assert_eq!(id, target);
        let tree = index.get_or_load_tree(target).unwrap();
        assert_eq!(count as usize, tree.len());
        let mut decoded = 0;
        for _ in 0..count {
            let iv: DiskInterval =
                bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                    .unwrap();
            assert!(iv.last >= iv.first);
            decoded += 1;
        }
        assert_eq!(decoded as usize, tree.len());
    }
}
