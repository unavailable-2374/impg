# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`impg` (implicit pangenome graph) is a Rust-based tool for projecting sequence ranges through many-way pairwise alignments (PAF format) without constructing full pangenome graphs. It enables extraction of homologous loci across genomes from specific genomic regions using alignments built by tools like `wfmash` and `minimap2`.

## Common Development Commands

### Build & Test
```bash
# Build the project
cargo build --verbose

# Build release version
cargo build --release

# Run all tests
cargo test --verbose

# Run a specific test
cargo test <test_name>

# Run tests for AGC integration (requires AGC feature)
cargo test --features agc

# Build without AGC support (FASTA-only)
cargo build --no-default-features

# Install from local source
cargo install --force --path .
```

### Running the Tool
```bash
# Basic query command
./target/release/impg query -p alignments.paf -r chr1:1000-2000

# Debug build (for development)
./target/debug/impg <command> [options]

# Partition command example
./target/debug/impg partition -p alignments.paf -w 1000000
```

## Architecture Overview

### Core Data Structures

**Interval Tree Index (`Impg` struct in `src/impg.rs`)**
- Uses `coitrees` (Cache Oblivious Interval Trees) for efficient range lookup over PAF alignments
- Organizes alignments by target sequence ID using a "forest" of interval trees (one tree per target)
- CIGAR strings are converted to compact delta encoding (`CigarOp`) stored as 32-bit values:
  - Upper 3 bits: operation type (=, X, I, D, M)
  - Lower 29 bits: operation length
- Index is serializable with `bincode` for fast loading/saving

**ForestMap (`src/forest_map.rs`)**
- Maps target sequence IDs to byte offsets in the serialized index file
- Enables random access to interval trees without loading entire index
- Used for efficient querying of specific target sequences

**SequenceIndex (`src/seqidx.rs`)**
- Bidirectional mapping between sequence names (strings) and compact IDs (u32)
- Tracks sequence lengths
- Reduces memory usage by avoiding string duplication in interval tree

**PAF Parsing (`src/paf.rs`)**
- `PartialPafRecord`: Lightweight PAF representation using sequence IDs instead of names
- Stores file offset to CIGAR string for lazy loading
- Packs strand information into MSB of 64-bit offset field

**Sequence File Support (`src/sequence_index.rs`, `src/faidx.rs`, `src/agc_index.rs`)**
- `UnifiedSequenceIndex`: Enum wrapping FASTA index or AGC index
- Provides trait-based interface for sequence extraction
- AGC support is optional (controlled by `agc` feature flag)

### Command Architecture

Commands are in `src/commands/`:
- `query`: Extract regions and their transitive alignments
- `partition`: Split alignment into genomic chunks
- `similarity`: Compute pairwise sequence similarity matrices with optional PCA
- `lace`: Combine multiple GFA or VCF files
- `index`: Build IMPG index from PAF files
- `stats`: Print PAF alignment statistics

Each command in `src/main.rs` uses:
1. `PafOpts` for PAF/index file handling
2. `SequenceOpts` for FASTA/AGC file resolution
3. `CommonOpts` for thread count and logging

### Graph Generation (`src/graph.rs`)

- Uses SPOA (partial order alignment) library for MSA-based graph construction
- Converts extracted sequences into GFA/MAF/FASTA formats
- Handles reverse complement sequences by adjusting path orientations in GFA

## Key Implementation Details

### CIGAR String Handling
- PAF files must use `--eqx` CIGAR format (= for match, X for mismatch)
- CIGARs are stored as file offsets and loaded on-demand to reduce memory
- Delta encoding supports forward/reverse strand projection through alignments

### Transitive Closure
- Query command can perform BFS/DFS traversal to find transitively connected alignments
- Configurable depth limit (`-m`) to control expansion
- DFS mode (`--transitive-dfs`) reduces overlapping results but is slower

### Multithreading
- Uses `rayon` for parallel processing of queries and partitions
- Thread count configurable via `-t` flag (default: 4)
- GZI index (bgzip) enables faster multithreaded decompression of `.paf.gz` files

### PanSN Naming Convention
- Supports PanSN format: `sample#haplotype#contig`
- Similarity command can group by delimiter position for sample/haplotype-level analysis
- Path names in GFA output use `NAME:START-END` format where NAME can contain ':'

## Testing

- Integration tests in `tests/test_agc_integration.rs`
- CI pipeline in `.github/workflows/rust_build_test.yml` runs on Ubuntu
- Tests require Linux dependencies: `build-essential cmake zlib1g-dev`

## Output Formats

- **BED/BEDPE**: Simple coordinate ranges
- **PAF**: Alignment records for extracted regions
- **GFA**: Variation graph with POA-based alignment
- **MAF**: Multiple alignment format
- **FASTA**: Extracted sequences (optionally reverse-complemented)
- **FASTA-aln**: POA-aligned FASTA sequences

Scripts in `scripts/`:
- `faln2html.py`: Convert FASTA alignments to interactive HTML visualization
- `partitioning.sh`: Example partitioning workflow
- `plot_partitioning_stats.R`: Visualize partition statistics

## Dependencies

Critical crates:
- `coitrees`: Interval tree data structure
- `noodles`: PAF/BGZF parsing
- `spoa_rs`: Partial order alignment for graph generation
- `handlegraph`: GFA manipulation for lace command
- `rayon`: Parallel iteration
- `clap`: CLI argument parsing
- `agc-rs` (optional): AGC archive support
- 完成修改后可以提交到https://github.com/unavailable-2374/impg.git上