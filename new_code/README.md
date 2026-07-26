# HiCTools C++ — hic_pre & hic_addnorm

C++ implementation of the HiC pre-processing pipeline. Converts aligned read pairs to
V9 `.hic` files (`hic_pre`) and computes normalization vectors in-place (`hic_addnorm`).

## Requirements

| Dependency | Notes |
|------------|-------|
| C++17 compiler | GCC ≥ 7, Clang ≥ 5, or Apple Clang ≥ 10 |
| CMake ≥ 3.13 | |
| zlib | Required for block compression |
| pthreads | Required for parallelism |
| zstd | Optional — reserved for future V10 support |

On macOS with Homebrew:
```
brew install cmake zlib zstd
```

On Ubuntu/Debian:
```
apt install cmake zlib1g-dev libzstd-dev
```

---

## Building

```bash
cd new_code
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Binaries are placed in `build/`:
```
build/hic_pre
build/hic_addnorm
```

**Debug build** (with AddressSanitizer):
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

**Run tests:**
```bash
ctest --test-dir build --output-on-failure
```

**Install** to a prefix (e.g. `~/.local`):
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=~/.local
cmake --build build
cmake --install build
```

---

## hic_pre

Converts a contact pairs file to a HiC V9 `.hic` file.

### Usage

```
hic_pre [options] <input_pairs> <output.hic> <chrom.sizes>
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-r <res,...>` | `2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000` | Comma-separated BP resolutions |
| `-q <mapq>` | `0` | Minimum MAPQ score; reads below threshold are discarded |
| `-t <threads>` | `4` | Worker threads |
| `-T <tmpdir>` | `/tmp` | Directory for temporary files |
| `-g <genome>` | _(empty)_ | Genome ID written to file header (e.g. `hg38`) |
| `-f <format>` | `auto` | Input format: `auto`, `mnd`, `short`, `pairs` |
| `-d <depth>` | `2` | V9 block depth base (1–10) |
| `--intra` | off | Discard inter-chromosomal contacts |
| `--near-diag` | off | Discard contacts > 10 Mb from diagonal |

### chrom.sizes format

Two-column tab-separated file:
```
chr1    248956422
chr2    242193529
chrX    156040895
```
UCSC chrom.sizes files work directly. `chrM` and other contigs can be included or
omitted as needed.

### Input formats

Format is auto-detected from extension (`.pairs`, `.mnd`) or file content. All
formats support `.gz` compression (decompressed via `gzip -dc`).

| Format | Fields | Description |
|--------|--------|-------------|
| Extra short | `chr1 pos1 chr2 pos2` | Minimal 4-column format |
| Extra short w/ score | `chr1 pos1 chr2 pos2 score` | 4-column with contact weight |
| Short (MND) | `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2` | Standard Juicer short format |
| Short w/ score | `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 score` | Short + contact weight |
| Old MND | `…frag2 mapq1 mapq2` | Short + MAPQ columns |
| New MND | `…frag2 mapq1 mapq2 score` | Short + MAPQ + score |
| Medium | `readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2` | Medium format with read name |
| Long (Juicer) | `str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 …` | Full merged_nodups format |
| 4DN pairs | `readID chr1 pos1 chr2 pos2 strand1 strand2 [frag1 frag2 …]` | 4DN DCIC .pairs format |

> **Note:** For formats that include a fragment field (`frag1`/`frag2`), contacts
> where both reads map to the same fragment are discarded. If using dummy fragment
> numbers, use `frag1=0` and `frag2=1` (not the same value) to avoid losing reads.

As with Java `hic_tools pre`, all records for a normalized chromosome pair must
form one contiguous block in the input. Read order within that block does not need
to be coordinate-sorted. If a chromosome pair appears again after a different pair
has begun, `hic_pre` exits with an error instead of producing an incomplete matrix.

### Example

```bash
hic_pre -r 1000000,500000,100000 \
        -q 30 \
        -t 8 \
        -g hg38 \
        -T /scratch/tmp \
        merged_nodups.txt.gz \
        output.hic \
        hg38.chrom.sizes
```

---

## hic_addnorm

Computes normalization vectors and writes them into an existing V9 `.hic` file
**in-place**. Make a backup before running if needed.

Computes:
- **VC** — Vanilla Coverage (row sums of the contact matrix)
- **VC_SQRT** — Square root of VC
- **SCALE** — Matrix balancing via iterative normalization (intra-chromosomal only)

Normalized expected value vectors are also computed and stored for each norm type.

### Usage

```
hic_addnorm [options] <hic_file>
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--no-vc` | off | Skip VC normalization |
| `--no-vc-sqrt` | off | Skip VC_SQRT normalization |
| `--no-scale` | off | Skip SCALE normalization |
| `-t <threads>` | `4` | Threads for SCALE matrix multiply |
| `--tol <val>` | `1e-4` | SCALE convergence tolerance |
| `--iter <n>` | `2000` | SCALE maximum iterations |
| `--min-res <bp>` | `0` | Skip SCALE for resolutions finer than this (e.g. `5000` to skip 1 kb) |

### Example

```bash
# Typical run: all three normalizations, 8 threads
hic_addnorm -t 8 output.hic

# Skip SCALE at fine resolutions where data is sparse
hic_addnorm -t 8 --min-res 25000 output.hic

# VC only
hic_addnorm --no-scale --no-vc-sqrt output.hic
```

---

## Full pipeline example

```bash
# 1. Build contact matrix from aligned pairs
hic_pre -r 2500000,1000000,500000,250000,100000 \
        -q 30 -t 16 -g hg38 \
        merged_nodups.txt.gz \
        sample.hic \
        hg38.chrom.sizes

# 2. Add normalization vectors
hic_addnorm -t 16 --min-res 25000 sample.hic
```

The resulting `sample.hic` can be loaded in Juicebox, queried with straw, or
processed by any tool that supports V9 `.hic` files.

## Java/C++ validation

`validation/run_java_cpp_validation.sh` builds the local Java reference,
`hic_pre`, `hic_addnorm`, and C++ straw, then creates independent Java and C++
V9 files.  It explicitly requests only `VC`, `VC_SQRT`, and `SCALE` (KR is
deprecated), including 1 kb BP resolution by default.

The script does not merely compare the `.hic` file hashes, because valid files
can have different block ordering, compression bytes, metadata positions, and
floating-point roundoff.  It saves:

- file and canonical straw-dump MD5 fingerprints;
- V9 header, footer, matrix-index, expected-vector, and norm-index checks;
- direct Java/C++ normalization-vector correlations and relative errors;
- deterministic random intra- and inter-chromosomal straw comparisons at every
  requested resolution for observed and observed/expected matrices.

Run from the repository root:

```bash
THREADS=8 JAVA_HEAP=16g \
  new_code/validation/run_java_cpp_validation.sh \
  test.txt \
  new_code/validation-results
```

Each run is written to a timestamped directory under the output path.  SCALE
may legitimately be absent for a chromosome/resolution when the matrix is too
sparse to converge; those cases are recorded as missing-vector dump failures
rather than allowed to hang straw.

---

## Project structure

```
new_code/
├── CMakeLists.txt
├── common/
│   ├── hic_file_def.h       # V9 format constants and shared structs
│   ├── little_endian.h      # Portable little-endian I/O helpers
│   ├── thread_pool.h        # Simple work-stealing thread pool
│   ├── genome.{h,cpp}       # Chromosome list and coordinate utilities
│   ├── pair_parser.{h,cpp}  # Input format detection and parsing
│   ├── block_writer.{h,cpp} # Contact binning, block numbering, zlib compression
│   ├── expected_value.{h,cpp} # Distance-decay expected value computation
│   └── hic_writer.{h,cpp}   # HiC V9 binary writer
├── hic_pre/
│   ├── preprocessor.{h,cpp} # Core preprocessing pipeline
│   └── main.cpp
└── hic_addnorm/
    ├── scale_norm.{h,cpp}   # SCALE matrix balancing algorithm
    ├── norm_updater.{h,cpp} # In-place norm vector computation and writing
    └── main.cpp
```

## Design notes

- **RAM- and descriptor-efficient**: contacts are streamed through the input once
  and distributed to per-chromosome-pair temp files (12 bytes/contact). Only the
  current pair's output file is open during distribution; matrix building opens at
  most one contact file per active worker.
- **Parallelism**: `hic_pre` uses a bounded rolling work queue to bin and compress
  multiple chromosome pairs concurrently while ordered output continues on the
  writer thread. `hic_addnorm` uses multi-threaded matrix-vector multiply for SCALE.
- **Single input pass**: the whole-genome matrix and chromosome-pair temp streams
  are accumulated together, avoiding a redundant scan of large MND files.
- **No external reader dependency**: `hic_addnorm` includes its own lightweight
  block decompressor; it does not depend on straw or libcurl.
