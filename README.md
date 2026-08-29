# HiCTools C++ — hic_pre & hic_addnorm

C++ implementation of the HiC pre-processing pipeline. Converts aligned read pairs to
V9 `.hic` files (`hic_pre`) and computes normalization vectors in-place (`hic_addnorm`).

## V10 support

The default build and `hic_pre` output remain **V9**. A separate implementation
under [`v10/`](v10/README.md) provides direct V10 preprocessing, V9 conversion,
and a separate V10 normalization command:

```sh
cmake -S . -B build -DBUILD_V10=ON
cmake --build build -j4
build/hic_v10 pre input.pairs output.v10.hic hg38
build/hic_v10 convert input.v9.hic output.v10.hic
build/hic_v10 addnorm -t 8 output.v10.hic
```

`hic_v10 addnorm` computes VC, VC_SQRT, and SCALE at both materialized and
derived resolutions and writes raw and normalized expected vectors. See the
[V10 documentation](v10/README.md) for details, memory behavior, and legacy
vector length migration. V10 always derives 20, 50, 200, 500, and 2,000 bp from
their fixed finer anchors, materializes 500 kb, and uses adaptive V9-compatible
rotated blocks for cis matrices.

## Requirements

| Dependency | Notes |
|------------|-------|
| C++17 compiler | GCC ≥ 7, Clang ≥ 5, or Apple Clang ≥ 10 |
| CMake ≥ 3.13 | |
| zlib | Required for block compression |
| pthreads | Required for parallelism |
| zstd | Required only for the opt-in `hic_v10` executable (`-DBUILD_V10=ON`) |

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
hic_pre [options] <input_pairs> <output.hic> <chrom.sizes|genome-id>
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-r <res,...>` | `2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000` | Comma-separated BP resolutions |
| `-q <mapq>` | `0` | Minimum MAPQ score; a pair is discarded when either available MAPQ is below the threshold |
| `-t <threads>` | `4` | Worker threads |
| `-T <tmpdir>` | `/tmp` | Directory for temporary files |
| `-g <genome>` | _(empty)_ | Genome ID written to file header (e.g. `hg38`) |
| `-f <format>` | `auto` | Input format: `auto`, `mnd`, `short`, `pairs`, `bin`, `bn`, or `hbs` |
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
omitted as needed. Instead of a file, the final argument may be one of the built-in
cleaned assemblies: `hg19`, `hg38`, `mm9`, or `mm10`. Built-in assemblies contain
only the standard autosomes, sex chromosomes, and mitochondrial chromosome, in
the same order and at the same lengths as the corresponding Java hic-tools
resources. The built-in ID is also written to the `.hic` header unless `-g`
overrides it.

### Input formats

Format is auto-detected from extension (`.pairs`, `.mnd`, `.bin`, `.bn`) or, for
ambiguous text files, from file content. Text formats support `.gz` compression
(decompressed via `gzip -dc`). Binary `.bin` and `.bn` files are uncompressed
little-endian Juicer formats.

| Format | Fields | Description |
|--------|--------|-------------|
| Extra short | `<chr1> <pos1> <chr2> <pos2>` | Minimal 4-column format |
| Extra short with score | `<chr1> <pos1> <chr2> <pos2> <score>` | Extra short plus contact weight |
| Short | `<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>` | Standard 8-column Juicer short format |
| Short with score | `<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>` | Short plus contact weight |
| Old short MND with MAPQ | `<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>` | 10-column short format with MAPQ |
| New short MND with MAPQ and score | `<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2> <score>` | 11-column short format with MAPQ and contact weight |
| Medium | `<readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>` | Medium format with read name |
| Long (Juicer) | `<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>` | Full 16-column merged_nodups format |
| 4DN/DCIC pairs | `<readID> <chr1> <pos1> <chr2> <pos2> <strand1> <strand2> […]` | `.pairs`; optional `frag1`/`frag2` and `mapq1`/`mapq2` pairs are located from `#columns:` and may appear in any order |
| Juicer binary | `.bin` fixed 26-byte records | Strand, chromosome index, position, and fragment for both ends |
| Juicer short binary | `.bn` fixed 20-byte records | Chromosome index and position for both ends, plus a floating-point score |
| Hi-C binary short | `.hbs.gz` compressed binary records | Header chromosome table and resolution; exact integer counts |

### HBS compressed binary input

Both `hic_pre` and `hic_v10 pre` read `.hbs.gz` exports from straw's `dump` or
`subsample` commands. The suffix is auto-detected; `-f hbs` selects it explicitly.

```sh
build/hic_pre -r 1000,5000,10000 sampled.hbs.gz sampled.v9.hic chrom.sizes
build/hic_v10 pre -r 1000,5000,10000 sampled.hbs.gz sampled.v10.hic chrom.sizes
```

HBS embeds chromosome names, lengths, and a BP resolution. Records use uint16
chromosome IDs, uint32 bin indices, and uint16 counts, with a uint64 escape for
counts of 65,535 or more: 14 or 22 bytes before gzip. The format is specified in
[HBS_FORMAT.md](HBS_FORMAT.md). No text conversion or external decompression is
needed. Gzip and record corruption are reported as errors.

The supplied genome may use a different chromosome ordering, but all embedded
names must resolve uniquely and lengths must match. Output resolutions must be
multiples of the source resolution; current builders require reconstructed BP
coordinates to fit signed int32. V9 converts counts to its existing float32
weight path. V10 carries HBS counts as uint64 through parsing and temporary
spooling, with overflow-checked aggregation (unless `--scores` is requested).
HBS has no mapping qualities or fragment metadata; it uses the existing MAPQ
sentinel of 1000 and dummy fragments 0/1.

### Text input details

The accepted text layouts therefore match Java hic-tools: extra-short (4 or 5
fields), short (8 or 9), old short with adjacent MAPQs (10), medium (11 with a
read name), new short with MAPQs and score (11), long (16), and header-described
4DN/DCIC `.pairs`. `.mnd` is an extension hint for the non-DCIC text formats;
plain `.txt` input is detected from its first data line.

`-q` filters formats that carry MAPQ (`.pairs` with both MAPQ columns, old/new
short MND, medium, and long). As in Java, formats without MAPQ—including
extra-short, ordinary short, `.bin`, and `.bn`—use the sentinel value `1000`, so
normal MAPQ thresholds do not discard them. Binary chromosome fields are numeric
indices and must use the same chromosome ordering as the supplied genome.

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
SCALE first attempts to balance all nonempty rows from a `sqrt(VC)` starting
vector. If either vector convergence or the maximum balanced row-sum error
(`0.05`) fails, progressively stricter nonzero-count cutoffs are tried up to the
more permissive of the 20th-percentile coverage bound or a raw row-count
z-score of -1. A converged vector with row-sum outliers gets one targeted retry
with only those offending rows masked before the global cutoff is increased.

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
  most one contact file per active worker. Sorted block spill runs are merged one
  logical block at a time. Compressed blocks are staged on disk, so completed
  workers retain only block offsets while waiting for ordered output.
- **Parallelism**: `hic_pre` uses a bounded rolling work queue to bin and compress
  multiple chromosome pairs concurrently while ordered output continues on the
  writer thread. The writer streams each completed temporary section into the
  final V9 file and backpatches its index. `hic_addnorm` uses multi-threaded
  matrix-vector multiply for SCALE.
- **Single input pass**: the whole-genome matrix and chromosome-pair temp streams
  are accumulated together, avoiding a redundant scan of large MND files.
- **No external reader dependency**: `hic_addnorm` includes its own lightweight
  block decompressor; it does not depend on straw or libcurl.
