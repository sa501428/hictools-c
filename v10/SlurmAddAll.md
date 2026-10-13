# Repair stored All on a Slurm cluster

Use the **updated `hic_v10 add-all`** for large files. It now uses a streamed
repair implementation; a separate `hic_v10_large` flavor is unnecessary for
this operation. This does not change the `.hic` format or the memory behavior
of ordinary `hic_v10 pre` or unrelated normalization commands.

For a human V10 map with approximately 25 chromosomes and normal resolution
lists, a **100G allocation should be comfortably sufficient**, including for a
500 GB file with 100B+ original reads. 150G is also supported. This assessment
comes from bounded copy/index buffers, coarse aggregation, metadata inspection,
and regression tests; it is not a measurement on a 500 GB production file.

## What scales and what costs time

The repair aggregates raw canonical cis and trans contacts from the coarsest
stored BP anchor no coarser than the requested source, normally 2,500,000 bp.
It does not read original pairs/HBS or reconstruct fine matrices. Read depth
changes count magnitudes, not the number of iterations through original reads.
If 2.5 Mb is derived, the actual stored source can be finer; preflight reports
that source and the number of occupied source cells.

Original file bytes after the header are copied in parallel 8 MiB ranges to a
new temporary file beside the output. Header growth introduces one constant
offset shift. Matrix indexes are traversed/patched in 1 MiB batches; vector
indexes are processed one entry at a time. The All matrix and a new active footer
are appended. Original compressed block/vector bytes stay exact, and original
chromosome IDs, normalization vectors and expected scales are preserved.

There is no per-block relocation map, no sort of all block records, no
whole-index allocation, and no 512 MiB limit on the **total block index** in
this repair path. Existing limits on individual compressed/decompressed blocks,
header/footer records and individual vector entries remain. Matrix descriptor
caches still scale with chromosome-pair count and resolution count. The
preflight estimate includes these caches plus a deliberately conservative 8 GiB
allowance for block decoding and other allocations. It refuses estimates above
80% of the requested memory; this is a planning estimate, not a formal memory
bound for arbitrary malformed files.

The file is synced and validated before atomic publication. Validation scans
original matrix indexes, checks vector chunk headers, and decodes All at every
BP level; it does not decompress every original fine matrix. The source must
remain unchanged during repair. Its size and modification/change times are
checked before publication. Keep other writers and normalization jobs off the
source and destination throughout this job.

Unreferenced old bytes, including an obsolete footer, are intentionally retained
by this fast repair path. Output is approximately the original size plus All,
header growth and the new footer. This tradeoff avoids sorting, extra vector
copies and billions of tiny I/O operations.

The dominant cost is a full copy: about 500 GB read plus 500 GB written for a
500 GB source. Ideal copy time at **measured end-to-end copy throughput** is:

| Copy throughput | 500 GB copy time |
| --- | --- |
| 0.1 GB/s | 83 minutes |
| 0.25 GB/s | 33 minutes |
| 0.5 GB/s | 17 minutes |
| 1 GB/s | 8.3 minutes |
| 2 GB/s | 4.2 minutes |

These are estimates for copying only. Add aggregation, index updates, sync and
validation time. A rate quoted for aggregate read+write traffic must be divided
by two before using this table. Lustre/GPFS/NFS load, quotas and stripe placement
can dominate; extra RAM cannot solve a storage bottleneck. The command logs
stage transitions and copy progress/rate every approximately 30 seconds.

## Build on Linux

Transfer the updated source repository to the cluster. A macOS executable
cannot run on Linux. Load the cluster's compiler/CMake/zlib/zstd modules as
needed; the build requires C++17, zlib and zstd development headers/libraries.
From the `hictools-c` checkout:

```bash
BUILD_JOBS=8 bash scripts/build_add_all_v10.sh
build/v10-all/hic_v10 --help
```

Use a site-approved build node if compilation is restricted on login nodes.
`build_add_all_v10.sh` accepts an optional build-directory argument.
Both the executable and Python script, as well as input/output/job directories,
must be visible under the same paths on compute nodes. Python 3.8+ is required;
the launcher needs only the standard library.

## Preflight, then submit

Choose an output on a sufficiently fast filesystem. Allow approximately
**525 GB of additional free space** for a 500 GB input (5% headroom; a minimum
2 GiB margin is used for smaller files). Existing input/output files do not count
as free space. This also applies to `--in-place`: atomic replacement still needs
a complete second file. Check your user/project quota separately; `statvfs`
free-space checks cannot establish quota headroom. Set Lustre output-directory
striping according to site policy before launching, if applicable.

Start with 16 CPUs and 100G RAM on one node:

```bash
python3 scripts/slurm_add_all_v10.py \
  /shared/maps/input.v10.hic /shared/maps/input.with-all.v10.hic \
  --hic-v10 "$PWD/build/v10-all/hic_v10" \
  --work-dir /shared/jobs/hic-all \
  --threads 16 --memory 100G --time 12:00:00 \
  --partition YOUR_PARTITION --account YOUR_ACCOUNT \
  --dry-run
```

Omit `--partition`/`--account` if your Slurm defaults are appropriate. The dry run
reads metadata and writes a unique job directory containing `preflight.json`,
`repair.sbatch` and `submission.json`. It does not submit or copy the file.
Inspect the reported source resolution, source cell count, block count, memory
estimate and disk headroom. `--copy-gbps 0.5` changes only the displayed copy-time
estimate, not an I/O speed limit.

Run the same command **without `--dry-run`** to submit. This creates a new job
directory and prints the job ID. To request 150G, use `--memory 150G`.

For approximately native human 2.5 Mb All bins, add `--all-bins 1250`. The default
500-bin budget groups the source to a smaller overview. The data remain stored
at the existing finest advertised synthetic BP resolution in both cases.
`--all-source-resolution 2500000` is the default source preference.

Atomic in-place repair is explicit:

```bash
python3 scripts/slurm_add_all_v10.py /shared/maps/input.v10.hic --in-place \
  --hic-v10 "$PWD/build/v10-all/hic_v10" \
  --work-dir /shared/jobs/hic-all --threads 16 --memory 100G
```

A separate existing output is rejected unless `--overwrite` is supplied.
Existing usable All is preserved; an in-place call on it is a no-op.

## Monitor and assess the result

```bash
squeue -j JOB_ID
sacct -j JOB_ID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

The job directory holds `slurm-JOB_ID.out` and `.err`, a repeated node-side
`preflight-worker.json`, GNU `/usr/bin/time -v` output in `resources.txt` when
available, and `result.json` with exit code, elapsed time and output size. The
launcher checks the input fingerprint again on the compute node and refuses a
file changed since submission. `MaxRSS` may be reported on the `.batch` row.
A successful exit follows sync and structural/All validation before publication.

Start with 8–16 copy workers. Increase to 32 only if measured throughput benefits
on your filesystem. More workers may make a shared storage system slower.
Aggregation and index patching are sequential; parallel workers serve the
bulk-copy phase, which normally dominates. One job per file avoids distributed
coordination and a shared output writer. Multiple nodes/jobs for the *same* file
would still share the storage bottleneck and complicate safe pointer updates, so
this launcher deliberately requests one node. Different files can use separate
jobs when filesystem bandwidth permits.

## Interrupted jobs

Repair is atomic but is not resumable: a retry repeats the copy. A regular error
removes native temporary files and leaves the original/previous output intact.
Slurm termination, node failure or SIGKILL can leave a destination lock
`OUTPUT.add-all.lock`, `OUTPUT.add-all.*` copy temporary, or
`OUTPUT.all-matrix.*` small temporary. The lock records the Slurm job ID/PID.
**Confirm the previous job is no longer running before removing its lock and
its specific temporary files.** Do not use broad wildcards against unrelated
jobs. Then rerun the launcher. Never run two repairs against the same output;
the worker lock prevents concurrent jobs from this launcher.

## Code-only validation

The existing 4.40 GB hg38 sample repaired with 16 workers in approximately
1.7–2.8 seconds in local warm-cache runs; peak resident memory was 166,346,752
bytes (approximately 159 MiB), measured via `getrusage` on macOS. Independent
C++ Straw verification found 82,718 All cells with an exact raw sum of
1,000,062,822. This is evidence for the implementation and memory behavior,
not a prediction of cluster storage speed or a 500 GB benchmark.


CTest `v10_all_cluster` covers 131,075 fine blocks across multiple index batches,
counts above 2^53, byte/data preservation, atomic rejection of a malformed index,
Slurm script quoting, dry-run and mocked submission, worker execution and stale
input rejection. A sparse fixture exercises offsets beyond 4 GiB without copying
its hole, and a sparse index above 512 MiB verifies streamed traversal. Neither
creates a 500 GB file. On Linux, native fixture tests also impose a 256 MiB
address-space limit; on macOS that Linux-only check is skipped.

```bash
cmake --build build/v10-all --target test_v10_streamed_reader
ctest --test-dir build/v10-all -R v10_all_cluster --output-on-failure
```
