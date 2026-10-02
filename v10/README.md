# V10 tools

This directory implements the consolidated `HiCFormatV10.md` wire format from
`hic-format`: the 88-byte header, numeric matrix directory, independently
compressed Zstandard blocks, exact block indexes, sparse/bitmap/dense payloads, explicit integer counts, and
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
options. By default, whenever advertised, 2 and 5 bp are derived from 1 bp, 20 and
50 bp from 10 bp, 200 and 500 bp from 100 bp, and 2 kb from 1 kb. The 500 kb
level is materialized by default. Use `--derive TARGET:SOURCE` to replace a
default source or add another derived target, and `--materialize RESOLUTION` to
store a target that would otherwise be derived.

Each normalized chromosome pair must occupy one contiguous input block, as with
`hic_pre`. Numeric positions are accepted in `[0, chromosomeLength]` and consumed
as supplied, without requiring, inferring, or converting a coordinate origin. A
position equal to `chromosomeLength` is folded into the final real bin; values
between the endpoints are binned exactly as provided. Each run creates
a private `hic-v10-run-*` workspace beneath `-T` (default `/tmp`). The parser closes
each pair spool there and hands it to a bounded worker pool while it reads ahead
into later chromosome pairs. Workers aggregate only
materialized resolutions into compact, disk-backed matrix sections; derived
targets are reconstructed from their declared source without their own input pass
or temporary matrix. The ordered
writer consumes those sections while pair preparation continues. By default up
to `-t` pair jobs are active or waiting; `--read-ahead N` sets a smaller or larger
explicit spool bound. The same `-t` worker pool performs Zstandard block
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

`hic_v10 addnorm` adds VC, VC_SQRT, and SCALE at advertised BP and FRAG
resolutions, including derived resolutions. Select one or more types explicitly:

```sh
build/hic_v10 addnorm --norm VC output.v10.hic
build/hic_v10 addnorm --norm VC_SQRT,SCALE -t 8 output.v10.hic
build/hic_v10 add-norm --norm SCALE --min-res 25000 output.v10.hic
```

`add-norm` is an alias for `addnorm`. Without `--norm`, all three algorithms are
selected. The existing `--no-vc`, `--no-vc-sqrt`, `--no-scale`, `--tol`, `--iter`,
`--min-res`, and `--level` controls remain available. A `--no-*` option disables
computation for that type; it never removes existing vectors.

**Existing normalizations are retained exactly**, including types the program
does not know how to compute. Dictionary IDs stay stable and new names are
appended. An existing type/unit/resolution bundle (any normalization or normalized
expected entry at that resolution) is retained as a whole, rather than recomputed
or combined with newly computed vectors. This also preserves bundles available
for only some chromosomes. Raw expected (`EVI0`) entries stay exact; missing raw
expected entries are computed. Each new computed normalization gets its own
normalized expected (`NEVI`) entry and chromosome scale factors.

For materialized resolutions, the regular normalizer reads one chromosome matrix
at a time. For derived resolutions, it aggregates the declared source before
normalizing the target. For large BP count files, use the disk-backed
[`hic_v10_large addnorm`](../v10_large/README.md#add-normalizations-to-a-completed-v10-file)
command; it needs only the finished `.hic` file.

The update is atomic and in place. It writes a temporary file beside the input,
copies stored contact blocks and existing vector chunks without recompression,
updates absolute locators in metadata and indexes, and replaces the input after
flushing and syncing the completed file. File permission bits and ordered header
attributes are retained. Growing the dictionary requires expanding the header;
new vector chunks, NVI/EVI/NEVI directories, and the matrix footer are written with
valid relocated positions. Interleaved matrix/vector layouts are supported.
Obsolete directories are omitted, so repeated additions do not accumulate stale
copies of prior vector sections. A failed update leaves the input untouched.
Allow disk space for the replacement file as well as any computation workspace.

### Supply normalization vectors as text

Use `--norm-file` when the divisors have already been calculated elsewhere and
you want the tool to calculate the corresponding normalized expected arrays:

```sh
build/hic_v10 addnorm --norm-file vectors.txt completed.v10.hic
build/hic_v10_large addnorm --norm-file vectors.txt --memory 8GiB \
  --tmp /local/scratch completed.v10.hic
```

The first non-comment line is `HIC_NORM_VECTORS 1`. Each block specifies a type
name, the exact chromosome name from the file header, `BP` or `FRAG`, and an
advertised bin size. Values follow in bin-index order, starting at bin zero,
with no coordinate column. An `end` line terminates each vector. For example,
for a chromosome with three bins at 10 kb and two bins at 20 kb:

```text
HIC_NORM_VECTORS 1
# Values are divisors: normalized(i,j) = raw(i,j) / (N[i] * N[j]).
vector MY_NORM chr1 BP 10000
1.25
0.90
nan
end

vector MY_NORM chr1 BP 20000
1.10
0.95
end
```

Multiple types, chromosomes, and resolutions can share one text file. One value
per line is recommended; multiple whitespace-separated values per line are also
accepted. Blank lines and `#` comments outside quoted names are ignored. Names containing
spaces or `#` must be double-quoted; quotes and backslashes inside quoted names
are backslash-escaped.
Type names are arbitrary except `NONE`. The synthetic `ALL` overview cannot
receive a supplied vector. BP vectors require exactly
`ceil(chromosomeLength / binSize)` values; FRAG vectors require exactly
`ceil((siteCount + 1) / binSize)` values. Every block must target an advertised
resolution, including any derived target. Supply a separate target-resolution
vector; the tool does not derive normalization values from finer vectors.

Values are parsed directly into float32 and stored **without rescaling or
balancing**. Decimal and scientific notation are supported. Use `nan` or `0`
for unavailable bins. Nonpositive or nonfinite divisors are stored but excluded
from expected-value calculations. For exact bits (including signed zero or a
specific NaN payload), use `bits:XXXXXXXX`, for example `bits:7fc01234`.
Unrepresentable overflow is an error; text values follow float32 rounding.

Straw's C++ CLI can export v9 or v10 normalization vectors directly into this
format, including arbitrary stored custom types such as `RU` and `NDSCALE`:

```sh
straw dump-norms source.v9.hic --output-dir norms
# Or select one type and choose its filename:
straw dump-norms source.v9.hic --norm RU --output RU.norm.txt
build/hic_v10 addnorm --norm-file RU.norm.txt completed.v10.hic
# Both importers accept the same files; import every type additively:
for vectors in norms/*.norm.txt; do
  build/hic_v10_large addnorm --norm-file "$vectors" --tmp /local/scratch completed.v10.hic
done
```

Each exported file contains one normalization type across all its indexed
chromosomes, resolutions, and BP/FRAG units, using exact `bits:XXXXXXXX` words.
`NONE` and the synthetic overview chromosome are omitted. The destination must
advertise the corresponding chromosome names, units and resolutions. Importing
vectors recalculates normalized expected arrays from the **destination's raw
contacts**; it does not copy source expected arrays or source scale factors.

V9 often stores `floor(length/binSize)+1` values, whereas v10 requires the ceiling
bin count. For a vector whose source count differs from the v10 geometry, straw
writes `source-length N` immediately after the `vector` line and emits all `N`
original words. This optional directive explicitly enables length adaptation:
keep the overlapping bins exact, pad missing bins with `bits:7fc00000`, and keep
surplus words outside the addressable v10 bins in a header attribute named
`hictools.import.vector.0.<normId>.<chrId>.<unitId>.<resolutionId>`. Its value is
`N:` followed by eight hexadecimal digits per surplus word. No original words
are discarded. The v10 file's header defines the target geometry. Without this
directive, exact target bin counts remain mandatory. Duplicate/late directives,
invalid counts, and mismatches between `N` and supplied values are errors.
Attributes for a skipped existing bundle are not added or overwritten.

For each newly supplied type/unit/resolution bundle, the command calculates a
full-length NEVI array from the corresponding raw cis contacts divided by the
supplied divisors, using the existing expected-value smoothing and chromosome
scale-factor calculation. Only chromosomes supplied for that bundle contribute
to its normalized expected values. A bundle with no usable normalized contacts
gets an all-NaN expected array. Existing EVI entries are preserved; missing EVI
entries at imported resolutions are calculated from all available cis
chromosomes. Raw contact blocks and all previously stored normalization and
expected vectors remain exact.

The existing preservation rule also applies to text imports: if any NVI or NEVI
entry already exists for that type/unit/resolution, the entire supplied bundle
is skipped with a message. Use a different type name to retain both versions of
an externally calculated normalization. This command does not replace or fill
partial existing bundles. Bad lengths, duplicate keys, unknown chromosomes or
resolutions, invalid values, and missing `end` markers fail without changing
the `.hic` file. `--norm-file` cannot be combined with `--norm` or `--vectors`.

Input vectors are staged as temporary sidecars, with one chromosome vector and
one distance accumulator loaded for expected calculation at a time. `--tmp`
sets the scratch directory (default `/tmp`); workspaces are removed on success
or ordinary failure. The regular tool supports BP/FRAG count and score matrices.
The large text importer supports BP/FRAG count matrices, streaming materialized
blocks and using external sort and exact rollups for derived resolutions. It
requires no original build manifests. Per-bin arrays and decoded blocks require
memory beyond the large tool's `--memory` sort-buffer budget. The updated file
is staged beside the original before atomic replacement.

### Import another normalization algorithm

Both V10 executables can import BP vectors from the large builder's sidecar format:

```sh
build/hic_v10 addnorm --vectors /path/to/vectors.manifest output.v10.hic
build/hic_v10_large addnorm --vectors /path/to/vectors.manifest output.v10.hic
```

This imports a complete sidecar set, including any precomputed expected vectors,
and cannot be combined with `--norm` or `--norm-file`. Use the text option above
when you want normalized expected vectors calculated automatically. The manifest may name additional algorithms such as `KR` or a custom
name. Names map to the destination dictionary by string, so their source IDs
need not match the file's IDs. The same preservation policy applies. Imports use
BP bin sizes to resolve destination resolutions, validate vector dimensions and
chromosome IDs, and verify sidecar checksums while streaming exact float32 words.
Imported vectors must have been computed for the destination's contacts and
chromosome order; the sidecar source fingerprint is provenance, not a comparison
against a V10 file. Supply normalized expected entries too when O/E queries are
needed.

The sidecar format is `HIC_V10_LARGE_VECTORS 1`, implemented in
[`vectors.cpp`](../v10_large/vectors.cpp). Its text records are:

```text
source <16-digit-hex-provenance>
norms <count> "<type-name>" ...
vectors <count>
vector <kind> <norm-id> <chr-id> <resolution-index> <bp-bin-size> <word-count> <checksum-hex> "<sidecar-path>" <scale-count> [<chr-id> <f32-bits-hex>] ...
end
```

Kinds are `0` normalization, `1` raw expected, and `2` normalized expected. IDs
for unused fields conventionally use `UINT32_MAX`; normalization IDs index the
manifest's name list. Each H10W sidecar has a 64-byte header: `H10W`, version 1,
kind, normalization ID, chromosome ID, source resolution index, BP bin size,
word width 4 (all seven integers are u32), then u64 word count and u64 FNV-1a
checksum, and 16 reserved zero bytes. It is followed by exactly `word-count`
little-endian float32 words. The checksum covers the payload bytes. Paths are
used as stored; absolute paths avoid working-directory ambiguity. Expected
vectors cover the maximum chromosome bin count, with chromosome scale factors
stored in the manifest. The shared [`VectorOutput`](vector_updater.h) C++ API also
accepts exact words or a chunk loader for future built-in algorithms.

## Convert an existing V9 file

```sh
build/hic_v10 convert input.v9.hic output.v10.hic

# Bound preparation to two chromosome pairs while using eight workers.
build/hic_v10 convert -t 8 --read-ahead 2 -T /local/scratch \
  input.v9.hic output.v10.hic
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

The converter reads these arrays directly from the V9 file in V10-sized chunks.
Their total size is therefore not part of peak memory use; this is important for
files carrying chromosome-length normalization and expected arrays at 1 bp.

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
follow V10 semantics. A contact in V9's redundant terminal bin (present when
chromosome length is exactly divisible by the resolution) is folded into the
final V10 bin; this migrates legacy endpoint coordinates to V10's finite
`ceil(length / resolution)` geometry. Contacts farther outside the chromosome,
unsupported/malformed input, and arithmetic overflow produce errors.

## Derived resolutions

The standard default derives these seven fine-resolution intermediates:

| Target | Default materialized source |
|---:|---:|
| 2 bp | 1 bp |
| 5 bp | 1 bp |
| 20 bp | 10 bp |
| 50 bp | 10 bp |
| 200 bp | 100 bp |
| 500 bp | 100 bp |
| 2 kb | 1 kb |

When the default applies, advertising a target requires advertising its default
source. `--materialize T` forces the target to be stored, while `--derive T:S`
selects any advertised finer materialized divisor as its source. For example,
`-r 125,250,500 --derive 250:125 --derive 500:125` stores 125 bp and derives
both coarser levels from it. Likewise,
`-r 50,100,150,200,300,500 --materialize 50 --derive 100:50 --derive 150:50
--derive 200:50 --derive 300:50 --derive 500:50` stores 50 bp and derives the
other listed levels from it. Derived sources must be direct and materialized;
chains and nonintegral factors are rejected. The 500 kb level may also be
explicitly derived instead of materialized.

During V9 conversion, the writer compares every requested derived target
cell against deterministic source aggregation before discarding its blocks. A
mismatch, including a float rounding difference, fails conversion. Direct `pre`
defines these targets from exact source aggregation and never constructs or
spools redundant target matrices. Each derived resolution retains its own
available normalization and expected arrays.

## Storage and output safety

The writer requires V9-compatible rotated distance-band grids for cis matrices
and rectangular grids for trans matrices. It selects sparse, bitmap, or dense blocks and
all-default, default-exception, or direct value streams, and tries RAW,
BYTE_SHUFFLE, and XOR32 vector transforms. Each logical block is stored as one
`H10B` record with its own Zstandard frame. An `H10I` version-2 index stores the
exact block number, stored length, and absolute position for every block. Defaults
are a 256-bin minimum block scale, four workers, Zstandard level 6, and 65,536
values per vector chunk. With both `pre` and `convert`, `-t N` controls the shared
chromosome-pair preparation and block-compression worker pool. Conversion workers
decode independent V9 chromosome pairs into run-scoped matrix spools beneath
`-T DIR`; the ordered writer reads one resolution at a time while later pairs are
prepared in parallel. `--read-ahead N` bounds the number of prepared or active
pair spools (default `-t`), allowing disk usage and concurrent decoding memory to
be reduced independently of compression parallelism. At most `-t` encoded logical
blocks are queued for ordered output. The private spool workspace is removed on
success and ordinary error exits.
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
normalization bits, exact multi-block indexes, and transactional failures.
Python uses the system zstd shared library through `ctypes`.

Optionally test V9/V10 query parity and normalization export/import through
both v10 addnorm executables with the updated straw executable:

```sh
cmake -S . -B build -DBUILD_V10=ON \
  -DSTRAW_TEST_EXECUTABLE=/path/to/straw
ctest --test-dir build --output-on-failure
```
