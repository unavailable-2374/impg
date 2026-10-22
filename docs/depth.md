# `impg depth` — what it computes, what depends on the order, and what does not

This is the reference document for the VGP pangenome-depth analysis. For
user-facing usage see `notes/DEPTH_COMMAND.md`; **that file's output-format
example is stale** — the current TSV has four columns, not one per sample.

## The output

```
#id	length	depth	positions
1	281	2	GCA_035083965.1#0#CM068849.1:30000000-30000281;GCA_019207075.1#0#CM033043.1:30047195-30047467
```

One row per interval. `depth` is the number of distinct samples covering it.
`positions` is a `;`-separated list of `name:start-end`, one per sample, in the
sample's own coordinates.

**The first token is the claimant**: the sequence in whose coordinates the
interval is defined and which the scheduler assigned the region to. The
remaining tokens are sorted by sample id. Reading the file "per genome" by its
first token is therefore an order-dependent view — this is the mistake that
invalidated `clade_cons/full_*` in the analysis that consumes this tool.

## The algorithm, in the shape the paper needs

1. **Index.** Every pairwise alignment is indexed; with ≥ 100 files the default
   is one index per file (`--index-mode auto`), which is what 331,575 VGP
   alignments require.
2. **Order.** Sequences are ordered for anchoring: an explicit order first if
   one was given (`--ref-order-file`), then the `--ref` sample's sequences, then
   by alignment degree descending, then by length descending. With `--ref` the
   degree pre-scan is skipped entirely — the reference decides Phase 1.
3. **Phase 1.** Each anchor sequence is walked in 5 MB chunks. For every chunk
   the alignments covering it are collected (CIGAR-precise, hop 0), projected
   into the anchor's coordinates, and swept to produce intervals of constant
   sample set. Everything emitted is marked processed in a shared tracker.
4. **Phase 2.** Regions not claimed in Phase 1 are processed in locality-aware
   waves over the remaining sequences, re-reading the tracker between waves so a
   region discovered by an earlier anchor is never emitted twice.
5. **Fill.** With `--fai-list`, sequence not covered by any alignment is emitted
   at depth 1 — so every base of every input sequence appears exactly once and
   the file accounts for the whole corpus, not only its aligned part.

## What is and is not a function of the order

| quantity | order-dependent? |
|---|---|
| the depth VALUE at a locus | **no** — it is the count of distinct samples whose alignments cover it |
| which sequence claims a region (the partition) | **yes** |
| the coordinate frame a region is written in | **yes** — it is the claimant's |
| which alignment file is consulted for a locus | **yes**, because homology is not transitive |
| per-genome statistics keyed on the first token | **yes** — do not compute them |
| Σ length·depth, the lineage Venn, per-base depth read in a genome's own frame from every record mentioning it | **no** (demonstrated empirically, see below) |
| row order in the file, and the `id` column | **yes, and also thread-dependent** — the same binary run twice does not reproduce them. Never join on `id` |

The last row is not a design choice, it is a measured fact: the tag rebuild of
2026-09-11 produced a binary that is output-identical to the production one once
the `id` column is dropped and the rows sorted, and *the same binary run twice*
agrees only under the same transformation.

## Non-transitivity, concretely

Homology here is pairwise and not transitive: A may align to B and B to C
without A aligning to C. Two consequences the file shows directly:

- a locus of C can appear in a record claimed by A and in another claimed by B,
  so **the all-token view is not a partition** — only the claimant view is;
- which of A or B is consulted for that locus depends on which was scheduled
  first, so the sample list at a locus can legitimately differ between two runs
  that use different orders.

This is why the order-invariance claim has to be *measured* rather than argued,
and `../../VGP/anchor_invariance/` measures it: the same 380 alignments run under
five orders (production anchor, no anchor, the largest genome, an outgroup, and a
seeded random permutation of all 1,033 sequences), compared base by base in each
genome's own frame.

## The partition guarantee, and its measured exception

Claimant intervals are **disjoint**: no base is claimed twice. That is asserted
in code and confirmed on real data by `VGP/verify_depth_tsv.py`.

Completeness is weaker than the design assumed. On a three-genome run, 290 of
304 sequences are tiled exactly by the claimant view; 14 have gaps, totalling
**0.037% of bp**. Reading every token position instead does not fix it — that
view gains overlaps, for the non-transitivity reason above. So the verifier
checks the claimant partition, reports the all-token bp beside it, and names a
sequence that is short in one view and exact in the other as the benign case
rather than counting it as a failure.

## Reproducing the production runs

Two production TSVs underlie the paper (`581_cmaes.depth.tsv`,
`581_11.depth.tsv`), 581 genomes and 331,575 alignments each. The exact command
is from `a.sh`, **not** from `batch_depth.sh`:

```bash
export MALLOC_ARENA_MAX=24
export IMPG_TRANSIENT_HEADER_CACHE_LIMIT=512
export IMPG_SUB_INDEX_CACHE_BYTES=$((16 * 1024 * 1024 * 1024))
export IMPG_CIGAR_QUERY_BATCH=64
export IMPG_FILE_QUERY_CONCURRENCY=16
export IMPG_FILE_QUERY_BYTES=$((8 * 1024 * 1024 * 1024))

impg depth --alignment-list cmaes.txt -t 112 --fai-list fai.list \
  --output-prefix 581_cmaes --ref GCA_009914755.4 --cigar-precise \
  --phase2-max-waves 580 --resume
```

The binary must be built **without jemalloc**: jemalloc's per-arena mmap model
exhausted `vm.max_map_count` (65,530) while holding 331,575 per-file indices,
and `munmap` then returned ENOMEM at ~11.5 GiB RSS. The recipe, including the
three libraries the vendored wfmash C++ needs and the two dangling conda
symlinks that cost four build attempts, is in `VGP/env/BUILD_impg.md`.

Provenance: tag `depth-v1.0.0` is the working tree that produced both production
TSVs (source 2026-07-16 07:52, binary 07:54; runs 2026-07-27 and 2026-08-02).
`rust-toolchain.toml` pins the compiler.

## Tests

`cargo test depth` covers 58 unit tests in `src/commands/depth.rs` and
`depth_checkpoint.rs` plus 17 integration tests across `tests/depth_chain.rs`,
`tests/depth_resume.rs` and `tests/test_depth_transient.rs`. The checkpoint's
invalidation hash includes the requested sequence order, so resuming a run under
a different order is refused rather than silently interleaving two partitions —
with a unit test asserting that a *permutation* of the same names changes the
hash, not merely a different set of names.

## The Methods paragraph

> Pangenome depth was computed with `impg depth` (impg, tag `depth-v1.0.0`), a
> subcommand that counts, for every locus, the number of distinct assemblies
> whose pairwise alignments cover it. No reference genome is required; CHM13v2.0
> was passed to `--ref`, which schedules its sequences as the first anchors. The
> record partition — which assembly claims a region, and therefore the coordinate
> frame it is reported in — depends on that schedule, as does which alignment
> file is consulted at a locus, since pairwise homology is not transitive.
> Statistics reported here are schedule-invariant reads (Σ length·depth, lineage
> co-occurrence, and per-base depth read in each assembly's own frame from every
> record mentioning it); per-assembly statistics keyed on the claiming sequence
> are not used. Two properties were verified empirically over 380 alignments
> among 20 assemblies (see Supplementary). First, holding the scheduling path
> fixed, per-base depth is invariant to the sequence order: a degree-ranked
> order and a seeded random permutation of all 1,033 sequences agree to 7 parts
> in 10^7 of claimed sequence. Second, the depth itself was checked against an
> independent count made from the pairwise alignments alone, with no depth
> output in it, under the identity depth = partners + 1; the totals agree to
> 0.25%, and per base the two agree at a median of 98.8% of positions per
> assembly, within the 5.8% median disagreement between the two directions of a
> single pairwise alignment.
