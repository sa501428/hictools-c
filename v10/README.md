# V10 tools

This directory implements the consolidated `HiCFormatV10.md` wire format from
`hic-format`: the 88-byte header, numeric matrix directory, Zstandard pages,
checkpoint indexes, sparse/bitmap/dense blocks, explicit integer counts, and
chunked normalization and expected arrays. It is **not** the earlier experimental
format that added delta/all-one flags to V9 blocks.

## Build

The normal project build is unchanged: `hic_pre` produces V9 and `hic_addnorm`
updates V9. Enable the additional executable explicitly:

```sh
cmake -S . -B build -DBUILD_V10=ON
cmake --build build -j4
```

This adds `build/hic_v10`; it requires zstd headers and the zstd library. Neither
existing executable links the new implementation. The reader and writer repositories
can be built independently; hictools-c does not require a straw checkout.

## From pairs, short, or merged-nodups input

```sh
build/hic_v10 pre -r 1000,5000,10000,50000 \
  input.pairs output.v10.hic hg38

build/hic_v10 pre -f mnd -q 30 -r 1000,5000,10000 \
  merged_nodups.txt output.v10.hic chrom.sizes
```

`pre` reuses the existing input parsers, including extra-short, Juicer short,
merged-nodups/medium/long, header-described DCIC pairs, gzip text, `.bin`, and
`.bn`. The default BP resolution set is the same as V9. `-f short` selects the
4/5-column extra-short parser; use auto detection or `-f mnd` for Juicer's
8/9-column short layout. See `hic_v10 pre --help` for filtering and compression
options.

Each normalized chromosome pair must occupy one contiguous input block, as with
`hic_pre`. Positions must be within `[0, chromosomeLength)` after parsing;
coordinates are consumed using the existing parsers' conventions. The new path
uses one temporary pair spool in `-T` (default `/tmp`) and processes one
chromosome-pair resolution at a time. Its cell accumulator must fit in memory;
it does not yet spill binned cells or parallelize compression. The spool is
removed on success and failure.

Positive integral weights use checked `uint64_t` accumulation, so repeated
contacts do not lose precision above the float integer limit. Fractional,
negative, zero, or nonfinite scores select `SCORE_FLOAT32` for that pair. Use
`--scores` to force score storage. Scores are read as `f32` by the shared parser;
single-cell score bits are preserved, and repeated scores are added in input
order in `f64`, then rounded once to `f32`. Nonfinite/negative input scores disable
raw expected generation, rather than advertising an invalid distance curve.

Direct V10 preprocessing writes real chromosomes only, without the legacy
synthetic `ALL` overview. It generates raw expected vectors using the existing
expected-value calculation, extending unavailable terminal distances with NaNs.
It does not compute normalization vectors. To obtain the existing normalization
algorithms, normalize a V9 file and then convert it:

```sh
build/hic_pre -r 1000,5000,10000 input.pairs intermediate.v9.hic hg38
build/hic_addnorm intermediate.v9.hic
build/hic_v10 convert intermediate.v9.hic output.v10.hic
```

## Convert an existing V9 file

```sh
build/hic_v10 convert input.v9.hic output.v10.hic
```

Conversion preserves chromosome order and lengths, genome ID, unknown attributes
(including duplicates and their order), BP/FRAG resolutions, fragment sites,
canonical raw cells, normalization arrays, expected arrays, and chromosome scale
factors. It supports V9 sparse and dense blocks with all short/int coordinate
combinations. V9 float values that are positive integral `u64` values can become
exact integer counts; other values remain scores. `--scores` forces all matrices
to retain float storage. The converter cannot recover precision already lost by
the original V9 producer.

The V9 synthetic `ALL` matrix is preserved. If its resolution was absent from the
V9 header, that resolution is explicitly added to the V10 BP list. Matrices with
no data at an advertised resolution are represented as empty at that resolution.
Statistics are recomputed from decoded cells; unavailable standard deviation and
percentile fields use the required canonical NaN.

### Legacy vector length migration

V9 commonly stores `floor(length / resolution) + 1` normalization bins; V10
requires `ceil(length / resolution)`. V9 expected arrays may also be shorter than
the full distance range required by V10.

The converter preserves all overlapping words exactly, including signed zeros
and NaN payloads. It fills newly addressable entries with canonical NaNs, and
moves any surplus terminal words into ordered header attributes rather than
silently discarding them:

```text
hictools.v9.vector.<kind>.<norm-id>.<chr-id>.<unit>.<resolution-index>
  = <original-value-count>:<surplus-f32-words-as-8-hex-digits-each>
```

Kinds are `0` normalization, `1` raw expected, and `2` normalized expected. IDs
refer to the output header; fields unused by a kind are zero. This migration
metadata preserves the original words, but the advertised V10 vector lengths
follow V10 semantics. Unsupported/malformed input, out-of-range contacts, and
arithmetic overflow produce errors.

## Derived resolutions

By default every resolution is materialized. To omit redundant matrix pages,
explicitly request one or more derivations:

```sh
build/hic_v10 pre -r 100,200,500,1000,2000 \
  --derive 200:100 --derive 500:100 --derive 2000:1000 \
  input.pairs output.v10.hic chrom.sizes

build/hic_v10 convert --derive 500000:100000 input.v9.hic output.v10.hic
```

Targets and sources must already be present, the source must be finer and
materialized, and its bin size must divide the target exactly. Before discarding
a target's pages, the writer compares every target cell against deterministic
source aggregation. A mismatch, including a float rounding difference, fails the
operation. It never silently substitutes an approximate derived matrix. Each
derived resolution retains its own available normalization and expected arrays.

## Storage and output safety

The writer uses rectangular grids, selects sparse, bitmap, or dense blocks and
all-default, default-exception, or direct value streams, and tries RAW,
BYTE_SHUFFLE, and XOR32 vector transforms. Page bytes are contiguous within each
matrix resolution. Defaults are 256-bin blocks, a 128 KiB **uncompressed** page
target, Zstandard level 3, and 65,536 values per vector chunk. The uncompressed
page target bounds working memory; a single large logical block may exceed it.
The V10 reader also supports rotated cis grids written by other producers.

Output is staged beside its destination, backpatched, flushed, and renamed only
after completion. Input/output aliases are rejected. A failed operation leaves
an existing destination untouched. Successful operations replace the named
output, as the existing preprocessing tools do.

## Tests

```sh
ctest --test-dir build --output-on-failure
```

The V10 test uses an independent Python decoder and independent V9 fixtures. It
covers direct input formats, integer precision, score bits, derived resolution
validation, sparse/dense V9 variants, FRAG conversion, duplicate attributes,
normalization bits, multiple page checkpoint groups, and transactional failures.
Python uses the system zstd shared library through `ctypes`.

Optionally test V9/V10 query parity with the updated straw executable:

```sh
cmake -S . -B build -DBUILD_V10=ON \
  -DSTRAW_TEST_EXECUTABLE=/path/to/straw
ctest --test-dir build --output-on-failure
```
