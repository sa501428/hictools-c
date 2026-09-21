# `hic_v10_large`: out-of-core V10 construction

`hic_v10_large` is a separate implementation for very large, integer-count HBS
inputs. It does not call or modify `v10/pre.cpp`. Its intermediate formats are
fixed-width, checksummed, restartable, and suitable for jobs sharing a cluster
filesystem.

The design target is a 100-billion-record HBS stream at 1 bp without retaining
all occupied cells in RAM. Input counts stay `uint64_t`; duplicate-cell addition
is checked for overflow. Raw expected histograms and coverage accumulation use
`uint128_t` until their final float vectors are produced.

V10 itself stores each cell count and each matrix count sum as `u64`. A duplicate
cell or complete matrix whose exact total exceeds `UINT64_MAX` therefore fails
instead of wrapping, even though staging can report its larger `uint128` weight.

## Build

```sh
cmake -S . -B build -DBUILD_V10=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
build/hic_v10_large --help
```

This adds `hic_v10_large` alongside, rather than in place of, `hic_v10`.

## Pipeline

Use the same resolution list for every build task. The finest value must equal
the resolution embedded in HBS, and every other value must be an integer
multiple of it.

### 1. Validate, index, and shard the single gzip stream

```sh
build/hic_v10_large stage \
  --chunk-records 250000000 \
  --progress-records 50000000 \
  input.hbs.gz stage-dir

build/hic_v10_large inspect --verify stage-dir/stage.manifest
```

Staging is the only necessarily serial pass over the gzip member. It validates
the complete HBS header, chromosome geometry, records, pair contiguity, escaped
counts, gzip trailer, and CRC. It canonicalizes pair orientation and cis
coordinates and emits record-aligned `H10S` shards. `stage.manifest` records the
exact record interval, pair, checksum, count, and `uint128` weight for every
shard; it is also the requested quick index of contiguous pair boundaries.
`inspect` prints `[firstRecord, endRecord)` for every pair and shard plus each
shard filename.

Each staged record is 16 bytes. Thus 100 billion physical input records require
about 1.6 TB for this durable stage, before filesystem compression.

### 2. Build sorted cell streams

The local orchestrator is:

```sh
R=1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,25000,50000,100000,500000

build/hic_v10_large build-cells -r "$R" \
  --root-only --memory 16GiB --fan-in 128 \
  stage-dir/stage.manifest build-dir
```

For a scheduler, print the deterministic task graph:

```sh
build/hic_v10_large plan -r "$R" --root-only --fan-in 128 \
  stage-dir/stage.manifest build-dir
```

Run its stages in this dependency order:

1. every `map-root` task;
2. every listed `reduce-root` task;
3. every `build-pair` task;
4. one `finalize-build` task.

The corresponding commands are:

```sh
build/hic_v10_large map-root -r "$R" --root-only --memory 16GiB --fan-in 128 \
  stage-dir/stage.manifest build-dir SHARD_ID

build/hic_v10_large reduce-root -r "$R" --root-only --fan-in 128 \
  stage-dir/stage.manifest build-dir PAIR_ID GROUP_ID

build/hic_v10_large build-pair -r "$R" --root-only --fan-in 128 \
  stage-dir/stage.manifest build-dir PAIR_ID

build/hic_v10_large finalize-build -r "$R" --root-only \
  stage-dir/stage.manifest build-dir
```

Map tasks use a stable LSD radix sort on `(y,x)`, aggregate duplicates in each
RAM chunk, spill sorted runs, and perform a map-side combine. If a pair has more
shards than `--fan-in`, `plan` adds parallel reduction groups before its final
merge. Published task manifests make all commands idempotent. A completed pair
reclaims its map and reduction runs only after its durable pair manifest exists.

`--root-only` is the recommended low-scratch mode. It retains only the canonical
finest-resolution stream for each pair. Normalization and assembly materialize
one requested coarser stream from that root, consume it, verify it, and remove
it. Peak retained cell storage is therefore about 16 bytes per occupied finest
cell instead of up to 16 bytes per occupied cell per resolution. The tradeoff is
one sequential root-cell pass per requested coarser resolution in each consumer.

Omit `--root-only` on every build command to select the faster, high-scratch
mode. In that mode coarser resolutions never rescan HBS: the nearest completed
divisor is consumed and retained in a rollup DAG. Both modes aggregate one
coarsened row at a time in canonical order, using a sparse open-addressing row
accumulator rather than a chromosome-width array. Neither mode externally sorts
a coarser resolution.

`--memory` is the combined size of the two radix-sort record arrays, not a
request to allocate that amount twice. Merge readers and library buffers add a
smaller amount controlled by `--fan-in`. On a node with a 150 GB hard limit,
16–32 GiB per concurrent sort job leaves ample headroom. Running several jobs
on one node requires dividing that budget between them.

### 3. Compute normalization sidecars

The local orchestrator is:

```sh
build/hic_v10_large normalize -t 16 --memory 16GiB \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir
```

For a cluster:

```sh
build/hic_v10_large plan-normalize \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir
```

Run every `normalize-chr` task first, every `expected-res` task second, and
`finalize-vectors` last. Pass identical normalization flags to all tasks.

```sh
build/hic_v10_large normalize-chr -t 16 --memory 16GiB \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir CHR_ID

build/hic_v10_large expected-res \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir RESOLUTION_INDEX

build/hic_v10_large finalize-vectors \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir
```

VC and VC_SQRT are computed from exact `uint128` row coverage and converted to
float only for the V10 vector. SCALE first creates a disk-backed directed CSR;
its row offsets and 8-byte `(column,value)` entries are memory mapped. Matrix
multiplication partitions output rows between persistent worker threads, so it
does not allocate one chromosome-sized accumulation array per thread.

Each chromosome task processes resolutions from coarse bins to fine bins. If
SCALE fails its convergence and balanced-row-sum criteria at one resolution,
that chromosome skips SCALE at all finer/smaller-bin resolutions. VC,
VC_SQRT, and every other chromosome continue. The balancing policy mirrors the
existing V10 SCALE cutoff, row-rescue, convergence-rate, and post-scaling logic.

At hg38 chr1 1 bp, normalization vectors have roughly 249 million entries.
SCALE's peak heap is approximately 120 bytes per bin plus mapped CSR pages and
sort memory; expected-value tasks are lower. A 150 GB node is therefore a
realistic ceiling with a conservative `--memory` setting, assuming the CSR
sidecars live on sufficiently fast local or parallel storage.

### 4. Assemble the V10 file

```sh
build/hic_v10_large write --genome hg38 --memory 16GiB -t 16 \
  --tmp /fast/scratch/block-work \
  --vectors vectors-dir/vectors.manifest \
  stage-dir/stage.manifest build-dir/build.manifest output.hic
```

Assembly is transactional: it writes beside the destination, flushes and
`fsync`s it, then renames it. Matrix cells are externally sorted by V10 logical
block number. Block-run inputs are reclaimed after successful merges, and the
potentially very large H10I block index is streamed through a sidecar instead
of being held twice in memory. Zstandard compression uses a bounded, ordered
worker queue with at most `-t` logical blocks in flight. Normalization arrays
are chunked into independent H10V frames.

The writer follows V10's required adaptive block geometry. It does not retain
an adaptive logical block in RAM: workers identify a fixed-record span, scan it
to obtain bounds and encoding sizes, and stream positions and values through
Zstandard into a temporary encoded-block sidecar. Peak encoding memory is thus
independent of the block's occupied-cell count. A block whose raw or compressed
payload exceeds V10's own `u32` length fields fails explicitly. Mandatory V10
derived-resolution declarations are applied, while the independently rolled-up
cells remain available for norms and expected vectors.

## Restart and storage behavior

- Every published stage, run, and vector sidecar has a typed header, record
  count, and checksum. Task manifests are written by atomic rename.
- Re-running a completed task validates its artifacts and returns successfully.
- The final `.hic` destination is unchanged after an ordinary failed assembly.
- Use a shared high-throughput filesystem for `stage-dir` and `build-dir`.
  `--tmp` may point writer/SCALE scratch at node-local NVMe.
- Do not delete staged shards until all map tasks are complete. Do not delete
  canonical `build-dir/cells` or vector sidecars until the final `.hic` has been
  verified.
- In high-scratch mode, sparse 1 bp data can remain nearly unique at many
  resolutions. The worst-case retained cell storage for 17 resolutions is
  about 27.2 TB per 100 billion occupied root cells, in addition to the 1.6 TB
  stage. The one-billion-record validation sample measured 15.7x amplification
  through all 17 resolution streams.
- In `--root-only` mode, the corresponding retained root is at most 1.6 TB per
  100 billion occupied cells. Allow additional scratch for one on-demand
  rollup, the largest active merge level, writer block runs, and SCALE CSR.
  Put `--tmp` on node-local NVMe or a high-throughput parallel filesystem.

## Numerical equivalence

Integer cell aggregation and raw expected histogram numerators are exact and
overflow checked. Stored normalization and expected arrays are float32 because
that is the V10 representation. Parallel/out-of-core SCALE changes the order of
some floating-point additions, so equivalence means the same cells, masks,
convergence policy, and values within normal float rounding—not necessarily
byte-identical float words. The test suite compares the large builder's VC,
VC_SQRT, raw/normalized expected, and SCALE vectors against the existing
in-memory V10 path with explicit tolerances.

## Tests

```sh
ctest --test-dir build -R v10_large --output-on-failure
```

The integration test covers HBS validation, noncontiguous-pair rejection,
terminal-endpoint folding, zero counts, counts above `2^53`, shard map/reduce
resume, duplicate aggregation, streaming rollups, independent V10 decoding,
normalization task resume, and comparison with the in-memory V10 normalizer.
