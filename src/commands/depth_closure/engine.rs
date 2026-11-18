//! Delta propagation on spatial boundary pages, followed by parallel signed
//! union. Native CIGARs remain compressed; physical atom coordinates are
//! recovered with rank/select instead of a global coordinate array.
use super::bitmap::{Bits, SHIFT, WIDTH};
use super::blocks::{self, Atoms, Candidates, Change, Key, Page, Pages, Snapshot};
use super::parallel_union::Parents;
use super::relations::{Cache, Store};
use super::*;
use rustc_hash::FxHashMap as HashMap;
use std::sync::Arc;
use std::time::Instant;

#[derive(Default, Serialize, Deserialize)]
pub(super) struct Checkpoint {
    cache: PathBuf,
    sources: Option<PathBuf>,
    incidence: Option<PathBuf>,
    epochs: Vec<PathBuf>,
    refined: bool,
    parent_epochs: Vec<PathBuf>,
    parent_page_bytes: Option<usize>,
    union_cursor: usize,
    union_complete: bool,
}

struct CompactMembers<'a> {
    rows: RecordSet<2>,
    atoms: &'a Atoms,
}
impl MemberAccess for CompactMembers<'_> {
    fn len(&self) -> usize {
        self.rows.len
    }
    fn root(&self, index: usize) -> u64 {
        self.rows.get(index)[0]
    }
    fn get(&self, index: usize) -> io::Result<[u64; 6]> {
        let [root, code] = self.rows.get(index);
        let (seq, mut s, mut e) = self.atoms.get(code >> 4)?;
        let h = (e - s) / 2;
        match (code >> 2) & 3 {
            0 => {}
            1 => e = s + h,
            2 => {
                s += h;
                e = s + 1;
            }
            3 => s = e - h,
            _ => unreachable!(),
        }
        Ok([root, seq, s, e, code & 1, (code >> 1) & 1])
    }
}

fn merge<const N: usize>(
    dir: &Path,
    budget: usize,
    rows: Vec<Records<N>>,
) -> io::Result<Records<N>> {
    let mut sorter = Sorter::new(dir, budget.max(4096));
    for row in rows {
        if row.len > 0 {
            sorter.add_sorted_file(&row.path)?;
        }
        fs::remove_file(&row.path)?;
    }
    sorter.finish()
}

fn commit_changes(
    state: &mut State,
    pages: &Pages,
    changes: &[Change],
    checkpoint: Option<&str>,
    pending: &mut BTreeMap<Key, Change>,
    last_commit: &mut Instant,
    force: bool,
) -> io::Result<()> {
    if checkpoint.is_none() {
        return Ok(());
    }
    for c in changes {
        pending.insert(c.key, *c);
    }
    if !force && last_commit.elapsed().as_secs() < 300 {
        return Ok(());
    }
    // Keep the latest full page for every modified tile, but enqueue only the
    // current frontier on resume. Earlier generations remain immutable.
    let current: HashMap<_, _> = changes.iter().map(|c| (c.key, c.delta)).collect();
    let empty = pages.put((0, 0), &Bits::default())?;
    let committed: Vec<_> = pending
        .values()
        .map(|c| Change {
            key: c.key,
            full: c.full,
            delta: current.get(&c.key).copied().unwrap_or(empty),
        })
        .collect();
    pages.sync()?;
    let epoch = state.engine.as_ref().unwrap().epochs.len();
    let path = blocks::save_changes(&state.dir, epoch, &committed)?;
    state.engine.as_mut().unwrap().epochs.push(path);
    state.save(checkpoint)?;
    pending.clear();
    *last_commit = Instant::now();
    Ok(())
}

fn apply(snapshot: &mut Snapshot, changes: &[Change]) {
    for c in changes {
        snapshot.insert(c.key, c.full);
    }
}

// Each task owns a bounded set of native alignment IDs. Dense tiles therefore
// still provide many independent tasks, even for one connected chromosome.
fn propagate(
    store: &Store,
    pages: &Pages,
    incidence: &Records<4>,
    snapshot: &Snapshot,
    frontier: &[Change],
    dir: &Path,
    budget: usize,
) -> io::Result<Vec<Change>> {
    let mut ordered = frontier.to_vec();
    ordered.par_sort_unstable_by_key(|c| (std::cmp::Reverse(c.delta.len), c.key));
    // Enumerate at most one tile ahead. Native IDs are shared by that tile's
    // tasks; neither all incident edges nor one scratch file per task is
    // materialized for the whole frontier.
    let jobs = ordered.iter().flat_map(|change| {
        let mut ids = Vec::new();
        let result = blocks::incident(incidence, change.key, |id, side| {
            ids.push((id, side));
            Ok(())
        });
        if let Err(e) = result {
            return vec![Err(e)].into_iter();
        }
        let ids = Arc::new(ids);
        (0..ids.len())
            .step_by(256)
            .map(|s| {
                Ok((
                    change.key,
                    change.delta,
                    ids.clone(),
                    s,
                    (s + 256).min(ids.len()),
                ))
            })
            .collect::<Vec<_>>()
            .into_iter()
    });
    let per_worker = (budget / rayon::current_num_threads().max(1) / 16).clamp(4096, 8 << 20);
    let rows = jobs
        .par_bridge()
        .try_fold(
            || Candidates::new(pages, dir, per_worker),
            |mut out, job| -> io::Result<_> {
                let (key, delta, ids, s, e) = job?;
                let bits = pages.get(delta)?;
                let lo = key.1 << SHIFT;
                let hi = lo + WIDTH - 1;
                let mut known: HashMap<Key, Arc<Bits>> = HashMap::default();
                for &(id, side) in &ids[s..e] {
                    let a = store.alignment(id)?;
                    let (s, e) = if side == 0 { (a.s, a.e) } else { (a.t, a.u) };
                    if s > hi || e < lo {
                        continue;
                    }
                    let begin = s.max(lo) - lo;
                    let end = e.min(hi) - lo + 1;
                    let Some((first, last)) = bits.range_bounds(begin as u32, end as u32) else {
                        continue;
                    };
                    store.diagonals(a, side, lo + first as u64, lo + last as u64, |d| {
                        let s = if side == 0 { d.s } else { d.t };
                        let begin = s.max(lo);
                        let end = (s + d.n).min(hi);
                        if begin > end {
                            return Ok(());
                        }
                        for p in bits.range((begin - lo) as u32, (end - lo + 1) as u32) {
                            let (seq, pos) = blocks::mapped(d, side, lo + p as u64);
                            let k = blocks::key(seq, pos);
                            if !known.contains_key(&k) {
                                if known.len() >= 64 {
                                    known.clear();
                                }
                                let b = snapshot
                                    .get(&k)
                                    .map(|&p| pages.get(p))
                                    .transpose()?
                                    .unwrap_or_default();
                                known.insert(k, b);
                            }
                            if !known[&k].contains(blocks::offset(pos)) {
                                out.push(seq, pos)?;
                            }
                        }
                        Ok(())
                    })?;
                }
                Ok(out)
            },
        )
        .map(|out| out.and_then(Candidates::finish_nonempty))
        .collect::<io::Result<Vec<_>>>()?;
    let rows = rows.into_iter().flatten().collect();
    let rows = merge_adaptive(dir, budget / 16, rows)?;
    let changes = blocks::reduce(pages, &rows, snapshot)?;
    if let Some(path) = rows.path {
        fs::remove_file(path)?;
    }
    Ok(changes)
}

fn union_block(
    store: &Store,
    incidence: &Records<4>,
    atoms: &Atoms,
    parents: &Parents,
    block: usize,
) -> io::Result<()> {
    let b = &atoms.blocks[block];
    if b.starts == 0 {
        return Ok(());
    }
    let lo = b.key.1 << SHIFT;
    let hi = lo + WIDTH - 1;
    let bits = atoms.bits(block)?;
    let mut ids = Vec::new();
    blocks::incident(incidence, b.key, |id, side| {
        if side == 0 {
            ids.push(id);
        }
        Ok(())
    })?;
    ids.par_chunks(128).try_for_each(|ids| -> io::Result<()> {
        let mut view = parents.view();
        for &id in ids {
            let a = store.alignment(id)?;
            if a.s > hi || a.e <= lo {
                continue;
            }
            store.diagonals(a, 0, lo, hi, |d| {
                let begin = lo.max(d.s);
                let end = (lo + WIDTH).min(d.s + d.n);
                if begin >= end {
                    return Ok(());
                }
                let mut points = bits
                    .range((begin - lo) as u32, (end - lo) as u32)
                    .peekable();
                let mut qid = None;
                let mut source_id = None;
                while let Some(p) = points.next() {
                    let start = lo + p as u64;
                    let finish = if let Some(&next) = points.peek() {
                        lo + next as u64
                    } else {
                        atoms.successor(b.key.0, start)?
                    };
                    if finish > d.s + d.n {
                        return Err(invalid("unpropagated diagonal end boundary"));
                    }
                    let target = d.t
                        + if d.reverse {
                            d.n - (finish - d.s)
                        } else {
                            start - d.s
                        };
                    let q = match qid {
                        Some(q) => q,
                        None => atoms.lookup(d.b, target)?,
                    };
                    let source = source_id.unwrap_or_else(|| b.base + bits.rank(p) as u64);
                    view.join(source, q, d.reverse as u64)?;
                    source_id = Some(source + 1);
                    qid = if points.peek().is_some() {
                        Some(if d.reverse {
                            q.checked_sub(1)
                                .ok_or_else(|| invalid("reverse atom cursor underflow"))?
                        } else {
                            q + 1
                        })
                    } else {
                        None
                    };
                }
                Ok(())
            })?;
        }
        Ok(())
    })
}

fn parent_commit(
    state: &mut State,
    parents: &Parents,
    cursor: usize,
    complete: bool,
    checkpoint: Option<&str>,
) -> io::Result<()> {
    if checkpoint.is_none() {
        let c = state.engine.as_mut().unwrap();
        c.union_cursor = cursor;
        c.union_complete = complete;
        return Ok(());
    }
    let refs = parents.checkpoint()?;
    let c = state.engine.as_mut().unwrap();
    let path = state
        .dir
        .join(format!("parents-{:08}.bin", c.parent_epochs.len()));
    blocks::write_checked(&path, &refs)?;
    c.parent_epochs.push(path);
    c.parent_page_bytes = Some(parents.page_bytes);
    c.union_cursor = cursor;
    c.union_complete = complete;
    state.save(checkpoint)
}

fn parent_memory_budget(budget: usize, count: u64) -> usize {
    let minimum = (budget / 4).max((rayon::current_num_threads() + 2) * 4096);
    // Fit the whole parent table when possible. Keep at least one eighth for
    // CIGAR windows, one eighth for boundary pages and one eighth headroom.
    let maximum = (budget - budget / 4 - budget / 8).max(minimum);
    let page = 2u128 << 20;
    let bytes = (u128::from(count) * 8).div_ceil(page) * page;
    bytes.clamp(minimum as u128, maximum as u128) as usize
}

fn merge_member_runs(
    dir: &Path,
    budget: usize,
    runs: Vec<Vec<PathBuf>>,
) -> io::Result<RecordSet<2>> {
    let mut sorter = Sorter::new(dir, budget.max(4096));
    sorter.add_owned_runs(runs.into_iter().flatten());
    sorter.finish().map(RecordSet::disk)
}

fn make_members(
    atoms: &Atoms,
    parents: Parents,
    dir: &Path,
    budget: usize,
) -> io::Result<RecordSet<2>> {
    let ranges: Vec<_> = (0..atoms.count)
        .step_by(65536)
        .map(|a| (a, (a + 65536).min(atoms.count)))
        .collect();
    // Witness bits remain on their original nodes while unions are concurrent.
    // This separate barrier resolves every witness to its final component root.
    ranges
        .par_iter()
        .try_for_each(|&(s, e)| -> io::Result<()> {
            let mut view = parents.view();
            for id in s..e {
                if view.witness(id)? {
                    let root = view.find(id)?.0;
                    view.mark_witness(root)?;
                }
            }
            Ok(())
        })?;
    let produce = |s, e| -> io::Result<Vec<[u64; 2]>> {
        let mut rows = Vec::with_capacity((e - s) as usize);
        let mut view = parents.view();
        for id in s..e {
            let (root, parity) = view.find(id)?;
            let bad = view.witness(root)?;
            if bad {
                let (_, s, e) = atoms.get(id)?;
                if e - s > 1 {
                    rows.push([root << 1, (id << 4) | (1 << 2)]);
                    if (e - s) % 2 != 0 {
                        rows.push([(root << 1) | 1, (id << 4) | (2 << 2) | 2]);
                    }
                    rows.push([root << 1, (id << 4) | (3 << 2) | 1]);
                    continue;
                }
            }
            rows.push([root << 1, (id << 4) | (u64::from(bad) << 1) | parity]);
        }
        Ok(rows)
    };
    // Two compact words per membership. Larger datasets spill sorted batches;
    // they never require a 6-word physical coordinate table in memory.
    if (atoms.count as u128) * 16 <= (budget as u128) / 4 {
        let chunks = ranges
            .par_iter()
            .map(|&(s, e)| produce(s, e))
            .collect::<io::Result<Vec<_>>>()?;
        let count = chunks.iter().map(Vec::len).sum();
        if (count as u128) * 16 > budget as u128 / 4 {
            let limit = (budget / rayon::current_num_threads().max(1) / 32).max(4096);
            let rows = chunks
                .into_par_iter()
                .try_fold(
                    || Sorter::new(dir, limit),
                    |mut sorter, chunk| -> io::Result<_> {
                        for row in chunk {
                            sorter.push(row)?;
                        }
                        Ok(sorter)
                    },
                )
                .map(|s| s.and_then(Sorter::finish_runs))
                .collect::<io::Result<Vec<_>>>()?;
            drop(parents);
            return merge_member_runs(dir, budget / 16, rows);
        }
        drop(parents);
        let mut rows = Vec::with_capacity(count);
        for chunk in chunks {
            rows.extend(chunk);
        }
        rows.par_sort_unstable();
        return Ok(RecordSet::memory(rows));
    }
    let limit = (budget / rayon::current_num_threads().max(1) / 16).max(4096);
    let rows = ranges
        .par_iter()
        .try_fold(
            || Sorter::new(dir, limit),
            |mut sorter, &(s, e)| -> io::Result<_> {
                for row in produce(s, e)? {
                    sorter.push(row)?;
                }
                Ok(sorter)
            },
        )
        .map(|sorter| sorter.and_then(Sorter::finish_runs))
        .collect::<io::Result<Vec<_>>>()?;
    // Parent lookups have finished. Retaining the entire atomic table while
    // reading the external membership runs only increases memory pressure.
    drop(parents);
    merge_member_runs(dir, budget / 8, rows)
}

pub(super) fn compute(impg: &impl ImpgIndex, options: &Options<'_>) -> io::Result<()> {
    if options.window_size.is_some_and(|n| n <= 0) {
        return Err(invalid("--window-size must be positive"));
    }
    if options.resume && (options.output_prefix.is_none() || options.stats) {
        return Err(invalid(
            "--resume requires --output-prefix and is not supported with --stats",
        ));
    }
    let (sequences, mapping) = universe(impg, options)?;
    let fingerprint = fingerprint(impg, &sequences, options)?;
    let budget = memory_budget()?;
    let checkpoint = options.output_prefix.filter(|_| options.resume);
    let mut temporary = None;
    let mut state = match checkpoint.map(State::load).transpose()?.flatten() {
        Some(state) => {
            if state.fingerprint != fingerprint {
                return Err(invalid(
                    "homology checkpoint input/configuration fingerprint mismatch",
                ));
            }
            if !state.dir.is_dir() || state.engine.is_none() {
                return Err(invalid(
                    "missing incremental homology checkpoint work files",
                ));
            }
            state
        }
        None => {
            let dir = tempfile::Builder::new()
                .prefix("impg-depth-closure-v7-")
                .tempdir()?;
            let path = if options.resume {
                dir.keep()
            } else {
                let path = dir.path().to_owned();
                temporary = Some(dir);
                path
            };
            let cache = if let Some(path) = std::env::var_os("IMPG_DEPTH_CACHE_DIR") {
                PathBuf::from(path)
            } else {
                options
                    .output_prefix
                    .and_then(|p| Path::new(p).parent())
                    .filter(|p| !p.as_os_str().is_empty())
                    .unwrap_or(Path::new("."))
                    .join(".impg-depth-cache-v7")
            };
            fs::create_dir_all(&cache)?;
            let state = State {
                version: VERSION,
                fingerprint,
                dir: path,
                units: Vec::new(),
                runs: None,
                bounds: None,
                delta: None,
                members: None,
                groups: None,
                output_group: 0,
                output_position: None,
                output_offset: 0,
                output_rows: 0,
                output_checkpoint: None,
                engine: Some(Checkpoint {
                    cache: fs::canonicalize(cache)?,
                    ..Default::default()
                }),
            };
            state.save(checkpoint)?;
            state
        }
    };
    info!(
        "Incremental homology depth: {} sequences; {} workers; working memory budget {:.1} GiB",
        sequences.len(),
        rayon::current_num_threads(),
        budget as f64 / (1u64 << 30) as f64
    );
    let explicit: BTreeMap<_, _> = options
        .explicit_order
        .unwrap_or(&[])
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let mut order: Vec<_> = sequences.iter().enumerate().collect();
    order.sort_by_key(|(_, s)| {
        (
            explicit.get(s.name.as_str()).copied().unwrap_or(usize::MAX),
            options.reference.is_some_and(|r| r != s.sample),
            std::cmp::Reverse(s.length),
            &s.name,
        )
    });
    let mut ranks = vec![0; sequences.len()];
    for (rank, &(id, _)) in order.iter().enumerate() {
        ranks[id] = rank;
    }
    let priorities: Vec<_> = mapping
        .iter()
        .map(|id| id.map(|i| ranks[i as usize]).unwrap_or(usize::MAX))
        .collect();
    if let Some(path) = options.dump_order {
        let mut out = BufWriter::new(File::create(path)?);
        for (_, s) in order {
            writeln!(out, "{}", s.name)?;
        }
        out.flush()?;
    }
    if state.engine.as_ref().unwrap().sources.is_none() {
        let start = Instant::now();
        let cache = Cache::new(state.engine.as_ref().unwrap().cache.clone())?;
        let mut sources = Vec::new();
        let mut source_chunks = Vec::new();
        // Bound simultaneous native shards. Within each shard, target/record
        // batches fan out on the same pool, without workers waiting on permits.
        for begin in (0..impg.depth_scan_units()).step_by(16) {
            // Blocking file/cache locks stay outside the Rayon pool. Otherwise
            // a nested target scan can deadlock if all workers steal sources
            // waiting on the bucket held by that scan.
            let next = std::thread::scope(|scope| {
                let handles = (begin..(begin + 16).min(impg.depth_scan_units()))
                    .map(|unit| {
                        let cache = &cache;
                        let mapping = &mapping;
                        let sequences = &sequences;
                        let priorities = &priorities;
                        scope.spawn(move || {
                            let s = cache.prepare(
                                impg,
                                unit,
                                options.coordinates,
                                &mapping,
                                &sequences,
                                &priorities,
                                options.sequence_index,
                            )?;
                            if checkpoint.is_some() {
                                test_exit_after_commit();
                            }
                            Ok(s)
                        })
                    })
                    .collect::<Vec<_>>();
                handles
                    .into_iter()
                    .map(|h| {
                        h.join()
                            .map_err(|_| invalid("native source preparation thread panicked"))?
                    })
                    .collect::<io::Result<Vec<_>>>()
            })?;
            sources.extend(next);
            if sources.len() >= 256 {
                let chunk = state
                    .dir
                    .join(format!("sources-{:08}.bin", source_chunks.len()));
                blocks::write_checked(&chunk, &sources)?;
                source_chunks.push(chunk);
                sources.clear();
            }
        }
        if !sources.is_empty() {
            let chunk = state
                .dir
                .join(format!("sources-{:08}.bin", source_chunks.len()));
            blocks::write_checked(&chunk, &sources)?;
            source_chunks.push(chunk);
        }
        let path = state.dir.join("sources.manifest.bin");
        blocks::write_checked(&path, &source_chunks)?;
        state.engine.as_mut().unwrap().sources = Some(path);
        state.save(checkpoint)?;
        info!(
            "Normalized native sources ready in {:.3}s",
            start.elapsed().as_secs_f64()
        );
        test_exit_after_stage(checkpoint, "sources");
    }
    let sources = Store::read_catalog(state.engine.as_ref().unwrap().sources.as_ref().unwrap())?;
    let mut store = Store::open(
        &state.engine.as_ref().unwrap().cache,
        sources,
        &sequences,
        budget / 2,
    )?;
    let pages = Arc::new(Pages::open(&state.dir.join("boundaries"), budget / 8)?);
    let mut snapshot = Snapshot::new();
    let mut frontier = Vec::new();
    let mut pending_changes = BTreeMap::new();
    let mut last_boundary_commit = Instant::now();
    for epoch in &state.engine.as_ref().unwrap().epochs {
        frontier = blocks::load_changes(epoch)?;
        apply(&mut snapshot, &frontier);
    }
    // An encoded empty sparse page is one tag byte and a four-byte count.
    frontier.retain(|c| c.delta.len > 5);
    if state.engine.as_ref().unwrap().incidence.is_none() {
        let start = Instant::now();
        let per_worker = (budget / rayon::current_num_threads().max(1) / 32).clamp(4096, 8 << 20);
        let rows = store
            .source_ranges()
            .into_par_iter()
            .try_fold(
                || {
                    (
                        Candidates::new(&pages, &state.dir, per_worker),
                        Sorter::new(&state.dir, per_worker),
                    )
                },
                |(mut candidates, mut index), (s, e)| -> io::Result<_> {
                    for id in s..e {
                        let a = store.alignment(id)?;
                        index.push([a.a, blocks::bin(a.s, a.e), id, 0])?;
                        index.push([a.b, blocks::bin(a.t, a.u), id, 1])?;
                        store.seed_diagonals(a, |d| {
                            candidates.push(d.a, d.s)?;
                            candidates.push(d.a, d.s + d.n)?;
                            candidates.push(d.b, d.t)?;
                            candidates.push(d.b, d.t + d.n)
                        })?;
                    }
                    Ok((candidates, index))
                },
            )
            .map(|result| {
                result.and_then(|(candidates, index)| Ok((candidates.finish()?, index.finish()?)))
            })
            .collect::<io::Result<Vec<_>>>()?;
        let (mut points, indices): (Vec<_>, Vec<_>) = rows.into_iter().unzip();
        let mut terminals = Candidates::new(&pages, &state.dir, per_worker);
        for (id, s) in sequences.iter().enumerate() {
            terminals.push(id as u64, 0)?;
            terminals.push(id as u64, s.length)?;
        }
        points.push(terminals.finish()?);
        let points = merge(&state.dir, budget / 16, points)?;
        let index = merge(&state.dir, budget / 16, indices)?;
        File::open(&index.path)?.sync_all()?;
        frontier = blocks::reduce(&pages, &points, &snapshot)?;
        state.engine.as_mut().unwrap().incidence = Some(index.path);
        commit_changes(
            &mut state,
            &pages,
            &frontier,
            checkpoint,
            &mut pending_changes,
            &mut last_boundary_commit,
            true,
        )?;
        apply(&mut snapshot, &frontier);
        fs::remove_file(points.path)?;
        info!(
            "Seeded {} spatial boundary blocks from {} native alignments in {:.3}s",
            frontier.len(),
            store.len(),
            start.elapsed().as_secs_f64()
        );
        test_exit_after_stage(checkpoint, "boundaries");
    }
    let incidence = Records::open(state.engine.as_ref().unwrap().incidence.as_ref().unwrap())?;
    let mut round = 0usize;
    while !state.engine.as_ref().unwrap().refined {
        let start = Instant::now();
        let changes = propagate(
            &store, &pages, &incidence, &snapshot, &frontier, &state.dir, budget,
        )?;
        let count = changes
            .iter()
            .try_fold(0usize, |n, c| pages.get(c.delta).map(|b| n + b.len()))?;
        state.engine.as_mut().unwrap().refined = changes.is_empty();
        let force = changes.is_empty()
            || std::env::var("IMPG_TEST_EXIT_AFTER_STAGE").as_deref() == Ok("refinement");
        commit_changes(
            &mut state,
            &pages,
            &changes,
            checkpoint,
            &mut pending_changes,
            &mut last_boundary_commit,
            force,
        )?;
        apply(&mut snapshot, &changes);
        round += 1;
        info!(
            "Delta refinement {}: {} new boundaries in {} blocks, {:.3}s",
            round,
            count,
            changes.len(),
            start.elapsed().as_secs_f64()
        );
        frontier = changes;
        test_exit_after_stage(checkpoint, "refinement");
    }
    let atoms = Atoms::new(&snapshot, pages.clone(), &sequences)?;
    if atoms.count >= 1 << 60 {
        return Err(invalid("too many initial atoms for compact membership IDs"));
    }
    let parent_budget = parent_memory_budget(budget, atoms.count);
    if state.members.is_none()
        && !state.engine.as_ref().unwrap().union_complete
        && parent_budget > budget / 4
    {
        // Frontier construction has finished. Transfer window-cache capacity
        // to parents before allocating them, without increasing total budget.
        store.release_thread_windows();
        drop(store);
        let sources =
            Store::read_catalog(state.engine.as_ref().unwrap().sources.as_ref().unwrap())?;
        store = Store::open(
            &state.engine.as_ref().unwrap().cache,
            sources,
            &sequences,
            budget - budget / 4 - parent_budget,
        )?;
    }
    let member_rows = if let Some(path) = state.members.as_ref() {
        RecordSet::disk(Records::open(path)?)
    } else {
        let start = Instant::now();
        let c = state.engine.as_ref().unwrap();
        let mut previous = BTreeMap::<u64, Page>::new();
        for path in &c.parent_epochs {
            for (id, p) in blocks::read_checked::<Vec<(u64, Page)>>(path)? {
                previous.insert(id, p);
            }
        }
        let previous: Vec<_> = previous.into_iter().collect();
        let parents = if let Some(layout) = c.parent_page_bytes {
            Parents::open_layout(
                &state.dir.join("parents"),
                atoms.count,
                parent_budget,
                &previous,
                layout,
            )?
        } else {
            Parents::open(
                &state.dir.join("parents"),
                atoms.count,
                parent_budget,
                &previous,
            )?
        };
        let mut order: Vec<_> = (0..atoms.blocks.len())
            .filter(|&i| atoms.blocks[i].starts > 0)
            .collect();
        order.par_sort_unstable_by_key(|&i| (std::cmp::Reverse(atoms.blocks[i].starts), i));
        info!(
            "Union partition: {} atoms in {} spatial blocks; parent cache budget {:.1} GiB",
            atoms.count,
            order.len(),
            parent_budget as f64 / (1u64 << 30) as f64
        );
        let mut committed = Instant::now();
        for begin in (c.union_cursor..order.len()).step_by(128) {
            let end = (begin + 128).min(order.len());
            order[begin..end]
                .par_iter()
                .try_for_each(|&i| union_block(&store, &incidence, &atoms, &parents, i))?;
            if committed.elapsed().as_secs() >= 300 || end == order.len() {
                parent_commit(&mut state, &parents, end, end == order.len(), checkpoint)?;
                info!("Union checkpoint: {end}/{} spatial blocks", order.len());
                committed = Instant::now();
            }
        }
        store.release_thread_windows();
        drop(store);
        let members = make_members(&atoms, parents, &state.dir, budget)?;
        if let Some(path) = members.path.as_ref() {
            File::open(path)?.sync_all()?;
        }
        // In-memory tables are reproducible from the durable final parent
        // generation. Resume reconstructs them before applying TSV offsets.
        state.members = members.path.clone();
        state.save(checkpoint)?;
        info!(
            "Parallel homology union: {} atoms, {:.3}s",
            atoms.count,
            start.elapsed().as_secs_f64()
        );
        test_exit_after_stage(checkpoint, "closure");
        members
    };
    let members = CompactMembers {
        rows: member_rows,
        atoms: &atoms,
    };
    let groups = if let Some(path) = state.groups.as_ref() {
        RecordSet::disk(Records::open(path)?)
    } else {
        let start = Instant::now();
        let groups = collect_groups(&members, &state.dir, budget)?;
        if let Some(path) = groups.path.as_ref() {
            File::open(path)?.sync_all()?;
        }
        state.groups = groups.path.clone();
        state.save(checkpoint)?;
        info!(
            "Canonical homology groups ready in {:.3}s",
            start.elapsed().as_secs_f64()
        );
        test_exit_after_stage(checkpoint, "groups");
        groups
    };
    info!("Writing {} canonical homology groups", groups.len);
    emit_output(
        &sequences, &members, &groups, options, &mut state, checkpoint,
    )?;
    info!(
        "Global homology depth complete: {} output records",
        state.output_rows
    );
    drop(groups);
    drop(members);
    drop(atoms);
    drop(pages);
    drop(incidence);
    if let Some(prefix) = checkpoint {
        fs::remove_file(format!("{prefix}.depth.ckpt"))?;
        fs::remove_file(format!("{prefix}.depth.work.bin"))?;
        fs::remove_dir_all(&state.dir)?;
    }
    drop(temporary);
    Ok(())
}

#[cfg(test)]
pub(super) fn test_partition(
    sequences: &[Sequence],
    alignments: &[DepthAlignment],
    mode: CoordinateMode,
) -> io::Result<Vec<Vec<usize>>> {
    let dir = tempfile::tempdir()?;
    let store = Store::from_test_alignments(dir.path(), sequences, alignments, mode)?;
    let pages = Arc::new(Pages::open(&dir.path().join("TEST-boundaries"), 4 << 20)?);
    let mut candidates = Candidates::new(&pages, dir.path(), 4096);
    let mut index = Sorter::new(dir.path(), 4096);
    for id in 0..store.len() {
        let a = store.alignment(id)?;
        index.push([a.a, blocks::bin(a.s, a.e), id, 0])?;
        index.push([a.b, blocks::bin(a.t, a.u), id, 1])?;
        store.seed_diagonals(a, |d| {
            candidates.push(d.a, d.s)?;
            candidates.push(d.a, d.s + d.n)?;
            candidates.push(d.b, d.t)?;
            candidates.push(d.b, d.t + d.n)
        })?;
    }
    for (id, s) in sequences.iter().enumerate() {
        candidates.push(id as u64, 0)?;
        candidates.push(id as u64, s.length)?;
    }
    let rows = candidates.finish()?;
    let index = index.finish()?;
    let mut snapshot = Snapshot::new();
    let mut frontier = blocks::reduce(&pages, &rows, &snapshot)?;
    apply(&mut snapshot, &frontier);
    while !frontier.is_empty() {
        frontier = propagate(
            &store,
            &pages,
            &index,
            &snapshot,
            &frontier,
            dir.path(),
            1 << 20,
        )?;
        apply(&mut snapshot, &frontier);
    }
    let atoms = Atoms::new(&snapshot, pages, sequences)?;
    let parents = Parents::open(
        &dir.path().join("TEST-parents"),
        atoms.count,
        (rayon::current_num_threads() + 4) * (64 << 10),
        &[],
    )?;
    (0..atoms.blocks.len())
        .into_par_iter()
        .try_for_each(|b| union_block(&store, &index, &atoms, &parents, b))?;
    let members = CompactMembers {
        rows: make_members(&atoms, parents, dir.path(), 1 << 20)?,
        atoms: &atoms,
    };
    let groups = collect_groups(&members, dir.path(), 4096)?;
    let mut offsets = vec![0usize];
    for s in sequences {
        offsets.push(offsets.last().unwrap() + s.length as usize);
    }
    let mut seen = vec![false; *offsets.last().unwrap()];
    let mut result = Vec::new();
    let mut cursor = 0;
    while cursor < groups.len {
        let span = next_span(&groups, &members, cursor)?;
        let (_, start, end, _) = member(&members, span.first, span.last, 0)?;
        for off in 0..end - start {
            let mut column = Vec::new();
            for i in 0..span.first[4] as usize {
                let (seq, s, e, rev) = member(&members, span.first, span.last, i)?;
                let base = offsets[seq] + (if rev == 0 { s + off } else { e - 1 - off }) as usize;
                assert!(!seen[base], "TEST new engine duplicate base {base}");
                seen[base] = true;
                column.push(base);
            }
            column.sort();
            result.push(column);
        }
        cursor = span.next;
    }
    assert!(seen.into_iter().all(|s| s), "TEST new engine omitted bases");
    result.sort();
    Ok(result)
}
