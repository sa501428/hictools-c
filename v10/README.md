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

`pre` reuses the existing input parsers, including HBS (`.hbs.gz`), extra-short, Juicer short,
merged-nodups/medium/long, header-described DCIC pairs, gzip text, `.bin`, and
`.bn`. The default BP resolution set is the same as V9. `-f short` selects the
4/5-column extra-short parser; use auto detection or `-f mnd` for Juicer's
8/9-column short layout. See `hic_v10 pre --help` for filtering and compression
options. Whenever advertised, 20 and 50 bp are derived from 10 bp, 200 and 500
bp from 100 bp, and 2 kb from 1 kb. The 500 kb level is materialized. These are
format requirements and do not need `--derive` flags.

Each normalized chromosome pair must occupy one contiguous input block, as with
`hic_pre`. Positions must be within `[0, chromosomeLength)` after parsing;
coordinates are consumed using the existing parsers' conventions. Each run creates
a private `hic-v10-run-*` workspace beneath `-T` (default `/tmp`). The parser closes
each pair spool there and hands it to a bounded worker pool while it reads ahead
into later chromosome pairs. Workers aggregate only
materialized resolutions into compact, disk-backed matrix sections; mandatory
derived targets are reconstructed from their source without their own input pass
or temporary matrix. The ordered
writer consumes those sections while pair preparation continues. By default up
to `-t` pair jobs are active or waiting; `--read-ahead N` sets a smaller or larger
explicit spool bound. The same `-t` worker pool performs Zstandard page
compression, so the command does not create a second unbounded set of threads.
Each active cell accumulator must fit in memory; it does not yet spill cells.
Logical blocks reuse that accumulator's cell storage rather than copying all
cells into per-block vectors. After all workers have joined, the private workspace
and all pair and matrix spools are removed on success and on ordinary error exits.
As with any process, an uncatchable termination such as `SIGKILL` can prevent
in-process cleanup.

Positive integral weights use checked `uint64_t` accumulation, so repeated
contacts do not lose precision above the float integer limit. Fractional,
negative, zero, or nonfinite scores select `SCORE_FLOAT32` for that pair. Use
`--scores` to force score storage. Non-HBS scores are read as `f32` by the shared parser;
single-cell score bits are preserved, and repeated scores are added in input
order in `f64`, then rounded once to `f32`. Nonfinite/negative input scores disable
raw expected generation, rather than advertising an invalid distance curve.

For HBS binary exports from straw, use:

```sh
build/hic_v10 pre -r 1000,5000,10000 sampled.hbs.gz output.v10.hic chrom.sizes
```

HBS input preserves exact uint64 counts, including single weights above 2^53,
through the parser, spool, and checked cell aggregation. Its chromosome table is
matched by name and length. Requested resolutions must be multiples of the
embedded BP resolution. See [the HBS specification](../HBS_FORMAT.md).

Direct V10 preprocessing writes real chromosomes only, without the legacy
synthetic `ALL` overview. It generates raw expected vectors using the existing
expected-value calculation, extending unavailable terminal distances with NaNs.
It does not compute normalization vectors during preprocessing. Add them to the
completed V10 file with the separate V10 implementation:

```sh
build/hic_v10 addnorm -t 8 output.v10.hic
```

`hic_v10 addnorm` computes VC, VC_SQRT, and SCALE for every advertised BP and
FRAG resolution. For a materialized resolution it reads that resolution's
pages; for a derived resolution it deterministically aggregates its declared
source in memory and normalizes the resulting cells. Raw expected (`EVI0`) is
rebuilt for every resolution, and each enabled normalization gets a normalized
expected (`NEVI`) vector with chromosome scale factors.

The command uses the same normalization controls as the V9 tool:

```sh
build/hic_v10 addnorm --no-scale output.v10.hic
build/hic_v10 addnorm -t 8 --tol 1e-4 --iter 2000 --min-res 25000 output.v10.hic
```

The update is atomic and in place: matrix bytes are copied unchanged, a fresh
vector section is written, and the temporary file replaces the input only after
it is flushed. Repeated runs replace the prior vector section without growing
the file. `pre` leaves the normalization dictionary empty and the fixed NVI
locator `(0, 0)`. When `addnorm` runs, it adds only the requested normalization
names, expands the variable header, relocates the unchanged matrix section, and
fills the NVI locator. This repack is staged in the same atomic temporary file.

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

The five expensive high-resolution intermediates are always derived:

| Target | Required materialized source |
|---:|---:|
| 20 bp | 10 bp |
| 50 bp | 10 bp |
| 200 bp | 100 bp |
| 500 bp | 100 bp |
| 2 kb | 1 kb |

Advertising a target requires advertising its source. The writer rejects a
materialized mandatory target, a different source, and any attempt to derive 500
kb. `--derive T:S` remains available only for additional nonstandard targets.
During V9 conversion, the writer compares every mandatory or requested target
cell against deterministic source aggregation before discarding its pages. A
mismatch, including a float rounding difference, fails conversion. Direct `pre`
defines these targets from exact source aggregation and never constructs or
spools redundant target matrices. Each derived resolution retains its own
available normalization and expected arrays.

## Storage and output safety

The writer requires V9-compatible rotated distance-band grids for cis matrices
and rectangular grids for trans matrices. It selects sparse, bitmap, or dense blocks and
all-default, default-exception, or direct value streams, and tries RAW,
BYTE_SHUFFLE, and XOR32 vector transforms. Page bytes are contiguous within each
matrix resolution. Defaults are a 256-bin minimum block scale, a 512 KiB
**uncompressed** page target (aiming at the specification's 64–256 KiB compressed
range), four workers, Zstandard level 6, and 65,536 values per vector chunk. Use
`-t N` with `pre` to bound pair preparation and page compression together; with
`convert`, it controls page compression. Use `--read-ahead N` when pair
accumulators are large and their concurrency needs a tighter memory bound. Queued
uncompressed pages are bounded to about 64 MiB in addition to active pair
accumulators; a single large logical block may exceed the page target.
The 256-bin block scale is only a lower bound. The V9 adaptive sizing formula
increases it sharply as resolution becomes finer—hg38 chr1 uses roughly 50,000
bins per rotated cis block at 10 bp—and increases it further if a `u32` block
number would overflow. `--block-bins` may request a larger lower bound but cannot
restore tiny fine-resolution blocks. This prevents hundreds of millions of
one- or two-cell logical blocks and their repeated 40-byte headers.

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
