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

### Slurm: V9 input to a validated V10 file

[`run_v9_to_v10.sbatch`](run_v9_to_v10.sbatch) submits the dump, staging,
build, normalization, write, validation, and comparison jobs in dependency
order. All input, output, work, scratch, source, and binary paths are supplied
at submission time. For example:

```sh
sbatch v10_large/run_v9_to_v10.sbatch \
  --input /path/to/input.v9.hic \
  --work-dir /path/to/shared/work \
  --tmp-dir /path/to/shared/scratch \
  --output /path/to/output.v10.hic \
  --hictools-repo /path/to/hictools-c \
  --straw-repo /path/to/straw \
  --genome hg38
```

Optional time limits apply to each task in the indicated category and are
exported to every child job. Omitted options retain these defaults:

| Option | Phases | Default |
| --- | --- | --- |
| `--time-short` | Dispatchers | `04:00:00` |
| `--time-medium` | Map, reduce, structural validation | `1-00:00:00` |
| `--time-long` | Dump, stage, expected values, pair writing, comparison | `3-00:00:00` |
| `--time-extra-long` | Chromosome normalization, final merge | `5-00:00:00` |
| `--time-build` | Pair cell construction | `2-00:00:00` |
| `--time-finalize` | Build and vector finalization | `08:00:00` |

For example, append `--time-medium 2-00:00:00 --time-long 5-00:00:00
--time-extra-long 7-00:00:00` to the command above. Values use Slurm time syntax
and are validated by `sbatch`. To change the initial setup job's limit, put
`--time=08:00:00` before the script name: `sbatch --time=08:00:00
v10_large/run_v9_to_v10.sbatch ...`. Script arguments can only affect jobs
submitted after setup starts.

Pass `--hic-bin` and `--straw-bin` instead of the corresponding repository
options to use prebuilt executables. With repository options, the setup job
builds both tools inside the work directory. The work and scratch directories
must be accessible from every node, and the work directory needs enough space
for the HBS dump, staged shards, cell streams, vector sidecars, and pair
fragments. The setup job's log uses Slurm's default output location; all later
jobs log to `WORK_DIR/logs`. The script keeps its submitted version in the work
directory so that child jobs run the same code. Published task artifacts can be
reused after a failed array task; an interrupted dump or staging pass starts
again. If a dependency chain stops, inspect the failed job or held array tasks
with `squeue` and `sacct`.

The final validation is structural. The two `straw compare` passes write their
results to the comparison job log; numerical differences can be expected for
normalization vectors, so review the cell comparison there.

### Manual pipeline

Use the same resolution list for every build task. The finest value must equal
the resolution embedded in HBS, and every other value must be an integer
multiple of it. If a resolution from the standard derived set must be retained
as a cell stream for materialized output, pass the same repeatable
`--materialize RESOLUTION` option to every build task and to `write` or each
`write-pair` task. `--derive TARGET:SOURCE` on the writer overrides a default
source or adds another derived target.

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
finest-resolution stream for each pair. Normalization materializes one requested
coarser stream from that root, consumes it, verifies it, and removes it. The
writer rolls up a bounded batch of resolutions during each root-stream pass and
feeds them directly to block builders without writing intermediate cell files.
Peak retained cell storage is therefore about 16 bytes per occupied finest cell
instead of up to 16 bytes per occupied cell per resolution. The tradeoff is a
small number of sequential root-cell passes in the writer and one pass per
requested coarser resolution in normalization.

Omit `--root-only` on every build command to select the faster, high-scratch
mode. In that mode materialized coarser resolutions never rescan HBS: the nearest completed
divisor is consumed and retained in a rollup DAG. Mandatory derived levels,
including 2 and 5 bp from 1 bp, have no retained cell stream. Both modes aggregate one
coarsened row at a time in canonical order, using a sparse open-addressing row
accumulator rather than a chromosome-width array. Neither mode externally sorts
a coarser resolution.

`--memory` is the combined size of the two radix-sort record arrays, not a
request to allocate that amount twice. Merge readers and library buffers add a
smaller amount controlled by `--fan-in`. On a node with a 150 GB hard limit,
16–32 GiB per concurrent sort job leaves ample headroom. Running several jobs
on one node requires dividing that budget between them.

### Add normalizations to a completed V10 file

The finished `.hic` file can be normalized after assembly, even when the original
HBS input, stage/build manifests, and cell runs are no longer available:

```sh
build/hic_v10_large addnorm --norm VC completed.v10.hic
build/hic_v10_large addnorm --norm VC_SQRT,SCALE -t 8 \
  --memory 8GiB --tmp /local/scratch completed.v10.hic
```

`add-norm` is an alias. Omitting `--norm` selects VC, VC_SQRT, and SCALE. Existing
`--no-*`, `--tol`, `--iter`, `--fan-in`, `--no-warm-start`, and
`--warm-vc-exponent` options also work. `--level` controls compression of new
vector chunks (default 3). Existing normalization types, their vectors,
normalized expected values, and raw expected values remain exact. Existing
normalization IDs stay stable. An existing type/resolution bundle is skipped as
a whole, including one present for only a subset of chromosomes; this command
adds new bundles rather than repairing or replacing existing ones.

The command decodes cis matrices one materialized block at a time, externally
sorts cells into disk-backed runs, aggregates derived targets from their declared
sources, then runs the existing large-data VC/VC_SQRT/SCALE pipeline. It does not
retain a full contact matrix in RAM. `--memory` bounds sort buffers, not total
RSS: decoded blocks, block indexes, per-bin vectors, expected accumulators, and
SCALE working arrays require additional memory. The normalization memory estimates
below still apply. Extraction and normalization workspaces are created under
`--tmp` (default `/tmp`) and removed on success or ordinary failure. Scratch must
hold extracted cis cell runs and normalization sidecars. The replacement `.hic`
is staged beside the original, requiring space for a second file on that
filesystem; the original is replaced only after successful completion.

This computation path supports **BP COUNT_UINT files**. For FRAG or score files,
use `hic_v10 addnorm`, which retains its in-memory normalization path. The shared
additive writer preserves all raw blocks and existing vector chunks without
recompression, including custom normalizations.

Externally computed divisors can also be supplied directly as text; the tool
calculates their normalized expected arrays and chromosome scale factors:

```sh
build/hic_v10_large addnorm --norm-file supplied-vectors.txt \
  --memory 8GiB --tmp /local/scratch completed.v10.hic
```

The [text vector format](../v10/README.md#supply-normalization-vectors-as-text)
uses `HIC_NORM_VECTORS 1`, followed by `vector TYPE CHROMOSOME UNIT BIN_SIZE`
blocks, one divisor per bin, and `end` markers. Imported divisors are not rescaled.
The large text path streams materialized count blocks and externally sorts and
rolls up derived targets; it supports BP and FRAG count matrices. Use the regular
text importer for score matrices. Existing bundles are preserved and skipped.

### Distributed import of supplied normalization vectors

For an already built v10 count file, `plan-addnorm` stages the same text dumps
accepted by `addnorm --norm-file` and creates independent expected-value and
compression tasks. It accepts multiple files in one plan, including one straw
dump per custom type. Neither the original build directory nor HBS input is
needed. BP and FRAG, materialized and derived resolutions, custom type names,
exact float32 words, and the `source-length` v9 migration convention all use the
[existing text contract](../v10/README.md#supply-normalization-vectors-as-text).
Score matrices are rejected by this large-data import path.

Build with `-DBUILD_V10=ON`. Prepare once on a shared filesystem:

```sh
# Produce one file per norm type, including custom names such as RU/NDSCALE.
/path/to/straw dump-norms source.v9.hic --output-dir /shared/norms

build/hic_v10_large plan-addnorm \
  --norm-file /shared/norms/RU.norm.txt \
  --norm-file /shared/norms/NDSCALE.norm.txt \
  --rows-per-task 262144 --norms-per-task 1 --chunks-per-task 16 \
  /shared/target.v10.hic /shared/addnorm-work
```

The target's chromosome names, units and advertised resolutions must match the
supplied vectors. Planning validates the text, stages uncompressed vector words,
and records the target's identity. Repeating preparation with identical inputs
and settings reuses the plan; other inputs require a new work directory. Input
text is no longer needed by workers after preparation. Keep the target unchanged
until assembly; mutation invalidates the plan.

Use the included Slurm launcher to create the job arrays and final assembly job:

```sh
python3 scripts/slurm_addnorm_v10_large.py /shared/addnorm-work/addnorm.plan \
  --hic-v10-large /shared/hictools-c/build/hic_v10_large \
  --partition compute --account my_account \
  --job-memory 32G --sort-memory 8GiB --fan-in 64 \
  --array-batch-size 1000 --max-jobs 100 --assemble-threads 16 --dry-run

# Remove --dry-run to submit the displayed dependency graph.
```

Paths to the target, executable, and work directory must be visible from every
node. Python 3 and `sbatch` are needed on the submission host; compute nodes only
need the executable and Bash. `--tmp /local/scratch` can move cell-extraction
scratch to node-local storage; durable outputs always go into the shared work
directory. `--array-batch-size` must fit the cluster's `MaxArraySize`.
`--cell-array-batch-size` also must fit that limit (default 1).
`--max-jobs` throttles each array separately; multiple batches/phases can run
concurrently. The launcher records submitted IDs, dependencies and scripts in
`slurm-submit/jobs.json`, including successful submissions preceding an error.
An append-only, synced `jobs.jsonl` journal in the submission folder preserves
IDs during submission, even if the launcher is killed before its final JSON
snapshot. Inspect the ledger/journal and stop obsolete pending jobs before an
explicit `--resubmit`; retries use new script
paths and retain the previous ledger. `--dry-run` creates scripts but submits no
jobs and creates no job ledger.

The phases are:

| Phase | Independent work | Required results |
| --- | --- | --- |
| `cells` | One chromosome/unit/resolution; extract, externally sort and partition cis contacts by row range | Prepared plan |
| `compress-norm` | Groups of 64K-word chunks of supplied vectors | Prepared plan; runs alongside `cells` |
| `expected` | Chromosome row range × group of normalization types | Its `cells` result |
| `reduce` | One normalization/unit/resolution, or missing raw expected bundle | Its expected partials |
| `compress-expected` | Groups of 64K-word expected-array chunks | Its `reduce` result |
| Assembly | One coordinator with parallel chunk copying | All vector and compression results |

The launcher defaults to `--dependency-mode ready`: it submits cell tasks
individually, releases expected arrays after their particular cell job succeeds,
then releases each reduction/compression group after its own prerequisite arrays
succeed. Arrays can run across phase boundaries. Increasing
`--cell-array-batch-size` batches extraction jobs at the cost of waiting for the
whole batch; `--dependency-mode phase` enables whole-phase barriers. Each array
contains tasks sharing prerequisites, with a generated lookup from its zero-based
Slurm indices to plan task indices. List commands, count tasks, or inspect the
machine-readable DAG for custom scheduling:

```sh
build/hic_v10_large addnorm-tasks --phase expected /shared/addnorm-work/addnorm.plan
build/hic_v10_large addnorm-tasks --phase expected --count /shared/addnorm-work/addnorm.plan
build/hic_v10_large addnorm-tasks --phase expected --dependencies /shared/addnorm-work/addnorm.plan

# In a corresponding Slurm array, after the required phase has completed:
build/hic_v10_large addnorm-run /shared/addnorm-work/addnorm.plan expected "$SLURM_ARRAY_TASK_ID"

# Cell workers also accept the sort budget, merge fan-in and scratch location.
build/hic_v10_large addnorm-run --memory 8GiB --fan-in 64 --tmp /local/scratch \
  /shared/addnorm-work/addnorm.plan cells "$SLURM_ARRAY_TASK_ID"

# Run after all compression jobs complete; this is the only target writer.
build/hic_v10_large addnorm-assemble -t 16 /shared/addnorm-work/addnorm.plan
```

`--dependencies` emits one line per task: its phase-local index followed by
`phase:index` prerequisites. For example, `3 cells:1` means expected task 3
requires cell task 1. `--count` and `--dependencies` are mutually exclusive.

Smaller `--rows-per-task` values create more expected jobs. `--norms-per-task 1`
maximizes independence across types; larger groups reuse each contact read and
reduce task/file counts. `--chunks-per-task 1` exposes every compression chunk
as its own job; larger groups amortize scheduler overhead. The defaults are
1,048,576 rows, one norm per expected task and 64 chunks per compression task.
`--level` sets Zstandard compression level when planning (default 3).

Expected workers hold a full chromosome-distance accumulator plus a full
chromosome normalization vector for each grouped type: budget approximately
`(sizeof(long double) + 4) × chromosome_bins × grouped_types`, plus decoded
records, manifests and allocator overhead. Reductions hold a full longest
chromosome-distance accumulator and expected array. Row sharding distributes
contact processing but does not shrink these arrays. `--sort-memory` only bounds
cell sort buffers; it does not cap process RSS. Smaller shards also create more
partial files and reduction I/O. Whole-chromosome extraction/sorting within a
cell task, global smoothing within a bundle reduction, and final directory/header
publication still each have a coordinator. New vector
compression is distributed. Assembly uses up to `-t` threads to copy raw matrix
ranges, existing vector chunks, and new precompressed chunks to disjoint file
offsets without recompression; pointer relocation is planned by the coordinator.
Its matrix-copy buffers add up to about 8 MiB per thread.

All durable result manifests identify their plan/task and contain checksums.
Successful tasks can be rerun without rewriting their output. Incomplete tasks
can be rerun; missing, corrupt or mismatched completed artifacts cause an error
before the target is replaced. After investigating corruption, remove the
corresponding `.done` marker and rerun that task and any affected downstream
tasks. Keep the work directory until assembly succeeds. Expected partials are
summed in a fixed task order, then the shared smoothing/scaling algorithm runs;
float rounding can differ slightly from the single-process importer.

Assembly stages a complete replacement beside the target and needs space for a
second `.hic`, in addition to the durable runs, vectors and compressed fragments.
It preserves raw block bytes, all existing vector chunks, normalization IDs and
existing expected/scaling values. Existing type/unit/resolution bundles are
skipped as a whole. Only missing supplied bundles and their expected values are
added; missing raw expected arrays are also built. Validated completion markers
make a repeated assembly a no-op, including recovery after replacement. Workers
never write the target, and concurrent assembly attempts for this workflow are
locked. This workflow imports supplied divisors; computed VC/VC_SQRT/SCALE from
raw contacts still use the existing `addnorm`/`normalize` commands.

If sidecars were already computed, attach them without rebuilding raw matrices:

```sh
build/hic_v10_large addnorm --vectors vectors/vectors.manifest completed.v10.hic
```

`--vectors` supports additional normalization names supplied by another
algorithm. See the [sidecar import contract](../v10/README.md#import-another-normalization-algorithm)
for file layout, dimensions, chromosome ordering, and checksum requirements.
The vector-only import can also preserve FRAG sections in the destination; the
imported sidecars describe BP vectors.

### 3. Compute normalization sidecars

The local orchestrator is:

```sh
build/hic_v10_large normalize -t 16 --memory 16GiB \
  stage-dir/stage.manifest build-dir/build.manifest vectors-dir
```

SCALE uses the last successful coarser-resolution vector as the starting point
for the next finer resolution. The stored normalization divisor is lifted by
genomic overlap (including non-integral transitions such as 500 bp to 200 bp),
converted to balancing weights, and adjusted by the finer-resolution VC. Use
`--no-warm-start` to retain the legacy independent `sqrt(VC)` starts, or
`--warm-vc-exponent N` to change the coverage blend from its conservative
default of `0.5`. A warm attempt that stalls automatically retries the legacy
start before any rows are discarded.

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
The cutoff ceiling is assessed from the complete nonzero-count distribution in
two ways: its 20th percentile and `mean - 1 standard deviation`; the more
permissive ceiling is retained. A successful coarse-resolution excluded-row
fraction is also used as a hint after an all-row fine-resolution attempt stalls.
Cutoffs that would produce an unchanged row mask are skipped, and rows with no
remaining active neighbor are removed before another balancing attempt.

With a `--root-only` build, `--cache-rollups` trades temporary disk space for
speed by retaining each normalization rollup until the corresponding expected-
value task consumes it. The default remains low-scratch behavior. Expected
calculation reads each chromosome/resolution cell stream once and accumulates
raw and all normalized expected vectors together.

At hg38 chr1 1 bp, normalization vectors have roughly 249 million entries.
SCALE's peak heap is approximately 120 bytes per bin plus mapped CSR pages and
sort memory; expected-value tasks are lower. A 150 GB node is therefore a
realistic ceiling with a conservative `--memory` setting, assuming the CSR
sidecars live on sufficiently fast local or parallel storage.

### 4. Assemble the V10 file

```sh
build/hic_v10_large write --genome hg38 --memory 16GiB -t 16 \
  --resolution-batch 4 \
  --tmp /fast/scratch/block-work \
  --vectors vectors-dir/vectors.manifest \
  stage-dir/stage.manifest build-dir/build.manifest output.hic
```

For a cluster, make matrix assembly a second map/reduce phase:

```sh
build/hic_v10_large plan-write \
  stage-dir/stage.manifest build-dir/build.manifest pair-parts output.hic

# Run one array job for every write-pair row printed by plan-write.
build/hic_v10_large write-pair --genome hg38 --memory 16GiB -t 16 \
  --resolution-batch 4 --tmp /fast/scratch/block-work \
  stage-dir/stage.manifest build-dir/build.manifest pair-parts PAIR_ID

# Run once after every pair job succeeds. Vectors are written only here.
build/hic_v10_large merge-pairs -t 16 \
  --vectors vectors-dir/vectors.manifest \
  stage-dir/stage.manifest build-dir/build.manifest pair-parts output.hic
```

Pair tasks are transactional and idempotent. Each produces a standalone,
independently valid V10 file with exactly one matrix. `merge-pairs` verifies the
source geometry, exact fragment header, expected pair identity, footer, block
indexes, and relocation fields. It then copies already-compressed matrix bytes
in 8 MiB windows and relocates absolute offsets in flight. Cell blocks are not
decoded or recompressed, so merge memory is bounded and its work is essentially
one sequential read plus one sequential write of the final matrix section.
Fragments must be on storage visible to the merge job. Do not pass `--vectors`
to `write-pair`; the final merge writes each normalization array only once.

Assembly is transactional: it writes beside the destination, flushes and
`fsync`s it, then renames it. Matrix cells are externally sorted by V10 logical
block number. In root-only mode, `--resolution-batch` controls how many
materialized resolutions share one scan of the root stream; their block-builder
memory budgets divide `--memory`. Block-run inputs are reclaimed after
successful merges, and the potentially very large H10I block index is streamed
through a sidecar instead of being held twice in memory. Zstandard compression
uses a bounded, ordered worker queue with at most `-t` logical blocks in flight.
Normalization arrays are chunked into independent H10V frames.

The writer follows V10's required adaptive block geometry. Ordinary logical
blocks with at most 8 MiB of fixed-width input records are cached, encoded, and
compressed directly into their final in-memory bytes, eliminating per-block
temporary-file traffic. Larger blocks use the bounded streaming encoder and a
temporary encoded-block sidecar, so an unusually dense block cannot exhaust
RAM. The cache/raw/compressed buffers are bounded per compression worker. A
block whose raw or compressed payload exceeds V10's own `u32` length fields
fails explicitly. The standard V10 derived-resolution defaults are applied
unless overridden, while derived cells are rolled up from materialized sources when
computing norms and expected vectors.

Use `hic_v10_large validate-v10 output.hic` for a streaming structural check.
Add `--matrix CHR1_ID:CHR2_ID:BIN` one or more times to fully decode selected
matrices and report their occupied-cell count, exact count sum, and checksum.

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

Normalization text exported by `straw dump-norms source.v9.hic --output-dir norms`
can be passed directly to `addnorm --norm-file`, one file per normalization.
Custom names (including `RU` and `NDSCALE`) use the same import path. The optional
`source-length` directive preserves v9 length differences and surplus bits as
metadata while storing vectors with legal v10 bin counts. See the
[shared text format and migration examples](../v10/README.md#supply-normalization-vectors-as-text).
