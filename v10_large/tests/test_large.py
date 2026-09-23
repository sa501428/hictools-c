#!/usr/bin/env python3
import gzip
import math
import pathlib
import struct
import subprocess
import sys
import tempfile
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "v10" / "tests"))
from inspect_v10 import Hic


def run(args, ok=True):
    p = subprocess.run([str(x) for x in args], text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if (p.returncode == 0) != ok:
        raise AssertionError(f"command returned {p.returncode}: {args}\nstdout={p.stdout}\nstderr={p.stderr}")
    return p


def hbs(path, records):
    names = [(b"chr1", 20), (b"chr2", 12)]
    data = bytearray(b"HICBS\0\r\n")
    data += struct.pack("<HHII", 1, 0, 1, len(names))
    for name, length in names:
        data += struct.pack("<H", len(name)) + name + struct.pack("<Q", length)
    for a, x, b, y, count in records:
        data += struct.pack("<HIHIH", a, x, b, y, min(count, 65535))
        if count >= 65535:
            data += struct.pack("<Q", count)
    path.write_bytes(gzip.compress(bytes(data)))


def read_run(path):
    data = path.read_bytes()
    assert data[:4] == b"H10R"
    assert struct.unpack_from("<I", data, 4)[0] == 1
    resolution = struct.unpack_from("<I", data, 16)[0]
    count = struct.unpack_from("<Q", data, 24)[0]
    records = [struct.unpack_from("<IIQ", data, 64 + 16 * i) for i in range(count)]
    return resolution, records


def float_word(word):
    return struct.unpack("<f", struct.pack("<I", word))[0]


def close_words(actual, expected, tolerance=2e-5):
    assert len(actual) == len(expected)
    for aw, ew in zip(actual, expected):
        a, e = float_word(aw), float_word(ew)
        if math.isnan(a) or math.isnan(e):
            assert math.isnan(a) and math.isnan(e)
        else:
            assert math.isclose(a, e, rel_tol=tolerance, abs_tol=1e-7), (a, e)


def main():
    executable = pathlib.Path(sys.argv[1])
    reference_executable = pathlib.Path(sys.argv[2])
    with tempfile.TemporaryDirectory() as temporary:
        root = pathlib.Path(temporary)
        source = root / "input.hbs.gz"
        # Pair 0x0 contains duplicates and deliberately unordered cells. Pair
        # 0x1 exercises pair transitions and trans coordinate ordering.
        hbs(source, [
            (0, 5, 0, 8, 2),
            (0, 1, 0, 2, 3),
            (0, 5, 0, 8, 4),
            (0, 9, 0, 9, (1 << 53) + 2),
            (0, 10, 0, 11, 0),
            (0, 3, 1, 7, 5),
            (1, 2, 1, 5, 7),
        ])
        stage = root / "stage"
        run([executable, "stage", "--chunk-records", "1", source, stage])
        inspected = run([executable, "inspect", "--verify", stage / "stage.manifest"])
        assert "records\t7" in inspected.stdout
        assert "pairs\t3" in inspected.stdout
        build = root / "build"
        plan = run([executable, "plan", "-r", "1,2,4", "--fan-in", "2",
                    stage / "stage.manifest", build]).stdout.splitlines()
        map_tasks = [int(line.split("\t")[1]) for line in plan if line.startswith("map-root\t")]
        reduce_tasks = [line.split("\t")[1] for line in plan if line.startswith("reduce-root\t")]
        pair_tasks = [int(line.split("\t")[1]) for line in plan if line.startswith("build-pair\t")]
        assert len(map_tasks) == 7 and reduce_tasks and len(pair_tasks) == 3
        for shard in map_tasks:
            run([executable, "map-root", "-r", "1,2,4", "--memory", "1MiB",
                 "--fan-in", "2", stage / "stage.manifest", build, shard])
        for task in reduce_tasks:
            pair, group = task.split(":")
            run([executable, "reduce-root", "-r", "1,2,4", "--fan-in", "2",
                 stage / "stage.manifest", build, pair, group])
        for pair in pair_tasks:
            run([executable, "build-pair", "-r", "1,2,4", "--memory", "1MiB",
                 "--fan-in", "2", stage / "stage.manifest", build, pair])
        run([executable, "finalize-build", "-r", "1,2,4", "--fan-in", "2",
             stage / "stage.manifest", build])
        run([executable, "build-cells", "-r", "1,2,4", "--memory", "1MiB",
             "--fan-in", "2", stage / "stage.manifest", build])
        files = sorted((build / "cells").glob("*-cells.h10r"))
        assert len(files) == 6
        cis = [p for p in files if "pair-00000-00000" in p.name]
        assert len(cis) == 2
        finest = next(p for p in cis if "-r1-" in p.name)
        resolution, records = read_run(finest)
        assert resolution == 1
        assert records == [(1, 2, 3), (5, 8, 6), (9, 9, (1 << 53) + 2)]
        assert not any("-r2-" in p.name for p in files)
        expected_r2 = [(0, 1, 3), (2, 4, 6), (4, 4, (1 << 53) + 2)]
        r4 = next(p for p in cis if "-r4-" in p.name)
        assert read_run(r4) == (4, [(0, 0, 3), (1, 2, 6), (2, 2, (1 << 53) + 2)])
        output = root / "large.hic"
        run([executable, "write", "--genome", "tiny", "--memory", "1MiB",
             stage / "stage.manifest", build / "build.manifest", output])
        hic = Hic(output)
        assert hic.records(0, 0, 1) == records
        assert hic.records(0, 0, 2) == expected_r2
        assert hic.res[0][1][1:] == (1, 1, 0, 0)
        assert hic.records(0, 0, 4) == read_run(r4)[1]
        assert hic.records(0, 1, 1) == [(3, 7, 5)]
        assert hic.records(1, 1, 1) == [(2, 5, 7)]
        assert hic.vector_locs == [(0, 0), (0, 0), (0, 0)]
        assert not list((build / "block-work").glob("*"))

        fine_build = root / "fine-build"
        run([executable, "build-cells", "-r", "1,2,5", "--memory", "1MiB",
             stage / "stage.manifest", fine_build])
        assert all("-r2-" not in cell.name and "-r5-" not in cell.name
                   for cell in (fine_build / "cells").glob("*-cells.h10r"))
        fine_output = root / "fine.hic"
        run([executable, "write", "--genome", "tiny", "--memory", "1MiB",
             stage / "stage.manifest", fine_build / "build.manifest", fine_output])
        fine = Hic(fine_output)
        assert fine.res[0][1][1:] == (1, 1, 0, 0)
        assert fine.res[0][2][1:] == (1, 1, 0, 0)
        assert fine.records(0, 0, 5) == [(0, 0, 3), (1, 1, (1 << 53) + 8)]

        # Matrix assembly can be distributed one chromosome pair per job and
        # merged without decoding or recompressing block payloads.
        pair_parts = root / "pair-parts"
        write_plan = run([executable, "plan-write", stage / "stage.manifest",
                          build / "build.manifest", pair_parts,
                          root / "merged.hic"]).stdout.splitlines()
        assert len([line for line in write_plan if line.startswith("write-pair\t")]) == 3
        assert write_plan[-1].startswith("merge-pairs\t")
        for pair in range(3):
            run([executable, "write-pair", "--genome", "tiny", "--memory", "1MiB",
                 "--resolution-batch", "2", stage / "stage.manifest",
                 build / "build.manifest", pair_parts, pair])
        resumed_pair = run([executable, "write-pair", "--genome", "tiny",
                            "--memory", "1MiB", stage / "stage.manifest",
                            build / "build.manifest", pair_parts, 0])
        assert "already complete" in resumed_pair.stderr
        merged = root / "merged.hic"
        run([executable, "merge-pairs", stage / "stage.manifest",
             build / "build.manifest", pair_parts, merged])
        merged_hic = Hic(merged)
        for pair in ((0, 0), (0, 1), (1, 1)):
            for resolution in (1, 2, 4):
                assert merged_hic.records(*pair, resolution) == \
                       hic.records(*pair, resolution)
        validated = run([executable, "validate-v10", "--matrix", "0:0:1", merged])
        assert "matrices\t3" in validated.stdout
        assert "matrix\t0\t0\t1\t3" in validated.stdout

        # Low-scratch mode retains only the finest cells. Writer and
        # normalization tasks roll up one requested resolution at a time and
        # retire the temporary run immediately.
        root_build = root / "root-build"
        run([executable, "build-cells", "-r", "1,2,4", "--root-only",
             "--memory", "1MiB", stage / "stage.manifest", root_build])
        assert len(list((root_build / "cells").glob("*-cells.h10r"))) == 3
        root_output = root / "root-only.hic"
        run([executable, "write", "--genome", "tiny", "--memory", "1MiB",
             stage / "stage.manifest", root_build / "build.manifest", root_output])
        root_hic = Hic(root_output)
        for pair in ((0, 0), (0, 1), (1, 1)):
            for resolution in (1, 2, 4):
                assert root_hic.records(*pair, resolution) == \
                       hic.records(*pair, resolution)
        assert not list((root_build / "block-work").glob("materialize-*.h10r"))

        vectors = root / "vectors"
        norm_plan = run([executable, "plan-normalize", "--no-scale",
                         stage / "stage.manifest", build / "build.manifest", vectors])
        assert norm_plan.stdout.count("normalize-chr\t") == 2
        assert norm_plan.stdout.count("expected-res\t") == 3
        for chromosome in range(2):
            run([executable, "normalize-chr", "--no-scale", "--memory", "1MiB",
                 stage / "stage.manifest", build / "build.manifest", vectors, chromosome])
        for resolution_index in range(3):
            run([executable, "expected-res", "--no-scale", "--memory", "1MiB",
                 stage / "stage.manifest", build / "build.manifest", vectors,
                 resolution_index])
        run([executable, "finalize-vectors", "--no-scale",
             stage / "stage.manifest", build / "build.manifest", vectors])
        resumed = run([executable, "normalize-chr", "--no-scale",
                       stage / "stage.manifest", build / "build.manifest", vectors, 0])
        assert "already complete" in resumed.stderr
        normalized = root / "large-normalized.hic"
        run([executable, "write", "--genome", "tiny", "--memory", "1MiB",
             "--vectors", vectors / "vectors.manifest",
             stage / "stage.manifest", build / "build.manifest", normalized])
        hic = Hic(normalized)
        assert hic.norms == ["VC", "VC_SQRT"]
        assert all(position and length for position, length in hic.vector_locs)
        for resolution, bins in ((1, 20), (2, 10), (4, 5)):
            assert len(hic.vectors[0, "VC", 0, 0, resolution][0]) == bins
            assert len(hic.vectors[0, "VC_SQRT", 0, 0, resolution][0]) == bins
            assert len(hic.vectors[1, None, None, 0, resolution][0]) == bins
            assert len(hic.vectors[2, "VC", None, 0, resolution][0]) == bins
            assert len(hic.vectors[2, "VC_SQRT", None, 0, resolution][0]) == bins

        merged_normalized = root / "merged-normalized.hic"
        run([executable, "merge-pairs", "--vectors", vectors / "vectors.manifest",
             "-t", "2", stage / "stage.manifest", build / "build.manifest",
             pair_parts, merged_normalized])
        merged_norm_hic = Hic(merged_normalized)
        assert set(merged_norm_hic.vectors) == set(hic.vectors)
        for key in hic.vectors:
            close_words(merged_norm_hic.vectors[key][0], hic.vectors[key][0], 0)
            assert dict(merged_norm_hic.vectors[key][1]) == dict(hic.vectors[key][1])

        root_vectors = root / "root-vectors"
        run([executable, "normalize", "--no-scale", "--cache-rollups", "--memory", "1MiB",
             stage / "stage.manifest", root_build / "build.manifest", root_vectors])
        root_normalized = root / "root-normalized.hic"
        run([executable, "write", "--genome", "tiny", "--memory", "1MiB",
             "--vectors", root_vectors / "vectors.manifest",
             stage / "stage.manifest", root_build / "build.manifest", root_normalized])
        root_norm_hic = Hic(root_normalized)
        assert set(root_norm_hic.vectors) == set(hic.vectors)
        for key in hic.vectors:
            close_words(root_norm_hic.vectors[key][0], hic.vectors[key][0])
            assert dict(root_norm_hic.vectors[key][1]) == dict(hic.vectors[key][1])
        assert not list(root_vectors.rglob("materialize-*.h10r"))
        assert not list(root_vectors.rglob("*.h10r"))

        chromosome_sizes = root / "chrom.sizes"
        chromosome_sizes.write_text("chr1\t20\nchr2\t12\n")
        reference_norm = root / "reference-normalized.hic"
        run([reference_executable, "pre", "-r", "1,2,4", source,
             reference_norm, chromosome_sizes])
        run([reference_executable, "addnorm", "--no-scale", reference_norm])
        reference_no_scale = Hic(reference_norm)
        assert set(hic.vectors) == set(reference_no_scale.vectors)
        for key in hic.vectors:
            close_words(hic.vectors[key][0], reference_no_scale.vectors[key][0])
            actual_scales = dict(hic.vectors[key][1])
            expected_scales = dict(reference_no_scale.vectors[key][1])
            assert set(actual_scales) == set(expected_scales)
            close_words(list(actual_scales.values()), list(expected_scales.values()))

        scale_source = root / "scale.hbs.gz"
        scale_records = []
        for i in range(10):
            scale_records.append((0, i, 0, i, 10))
            if i + 1 < 10:
                scale_records.append((0, i, 0, i + 1, 5))
            if i + 2 < 10:
                scale_records.append((0, i, 0, i + 2, 2))
        hbs(scale_source, scale_records)
        scale_stage = root / "scale-stage"
        scale_build = root / "scale-build"
        scale_vectors = root / "scale-vectors"
        run([executable, "stage", scale_source, scale_stage])
        run([executable, "build-cells", "-r", "1,2", "--memory", "1MiB",
             scale_stage / "stage.manifest", scale_build])
        scale_run = run([executable, "normalize", "--no-vc", "--no-vc-sqrt",
                         "--memory", "1MiB", "-t", "2", "--iter", "500",
                         scale_stage / "stage.manifest", scale_build / "build.manifest",
                         scale_vectors])
        assert "coarse warm start" in scale_run.stderr
        assert "degree cutoffs percentile=" in scale_run.stderr
        scale_hic = root / "scale.hic"
        run([executable, "write", "--vectors", scale_vectors / "vectors.manifest",
             "--memory", "1MiB", scale_stage / "stage.manifest",
             scale_build / "build.manifest", scale_hic])
        scale = Hic(scale_hic)
        assert scale.norms == ["SCALE"]
        assert len(scale.vectors[0, "SCALE", 0, 0, 1][0]) == 20
        assert len(scale.vectors[0, "SCALE", 0, 0, 2][0]) == 10
        assert not [path for path in scale_vectors.rglob("*")
                    if path.is_file() and "manifest" not in path.name and path.suffix not in (".h10w",)]
        reference_hic = root / "scale-reference.hic"
        run([reference_executable, "pre", "-r", "1,2", scale_source,
             reference_hic, chromosome_sizes])
        run([reference_executable, "addnorm", "--no-vc", "--no-vc-sqrt",
             "-t", "2", "--iter", "500", reference_hic])
        reference = Hic(reference_hic)
        for resolution in (1, 2):
            actual = scale.vectors[0, "SCALE", 0, 0, resolution][0]
            expected = reference.vectors[0, "SCALE", 0, 0, resolution][0]
            for a, b in zip(actual, expected):
                if math.isnan(a) or math.isnan(b):
                    assert math.isnan(a) and math.isnan(b)
                else:
                    assert math.isclose(a, b, rel_tol=2e-4, abs_tol=1e-6)

        bad = root / "bad.hbs.gz"
        hbs(bad, [(0, 1, 0, 2, 1), (0, 1, 1, 2, 1), (0, 3, 0, 4, 1)])
        failed = run([executable, "stage", bad, root / "bad-stage"], ok=False)
        assert "not contiguous" in failed.stderr

        endpoint = root / "endpoint.hbs.gz"
        hbs(endpoint, [(0, 20, 0, 20, 9)])
        endpoint_stage = root / "endpoint-stage"
        endpoint_build = root / "endpoint-build"
        run([executable, "stage", endpoint, endpoint_stage])
        run([executable, "build-cells", "-r", "1", "--memory", "1MiB",
             endpoint_stage / "stage.manifest", endpoint_build])
        endpoint_run = next((endpoint_build / "cells").glob("*-cells.h10r"))
        assert read_run(endpoint_run) == (1, [(19, 19, 9)])

        protected = root / "protected.hic"
        protected.write_bytes(b"keep-existing-output")
        damaged = bytearray(endpoint_run.read_bytes())
        damaged[-1] ^= 1
        endpoint_run.write_bytes(damaged)
        corrupt_write = run([executable, "write", "--memory", "1MiB",
                             endpoint_stage / "stage.manifest",
                             endpoint_build / "build.manifest", protected], ok=False)
        assert "checksum mismatch" in corrupt_write.stderr
        assert protected.read_bytes() == b"keep-existing-output"

        overflow = root / "overflow.hbs.gz"
        hbs(overflow, [(0, 1, 0, 2, (1 << 64) - 1), (0, 1, 0, 2, 1)])
        overflow_stage = root / "overflow-stage"
        run([executable, "stage", overflow, overflow_stage])
        overflow_build = run([executable, "build-cells", "-r", "1", "--memory", "1MiB",
                              overflow_stage / "stage.manifest", root / "overflow-build"],
                             ok=False)
        assert "exceeds uint64" in overflow_build.stderr

        truncated = root / "truncated.hbs.gz"
        truncated.write_bytes(source.read_bytes()[:-4])
        gzip_failure = run([executable, "stage", truncated, root / "truncated-stage"], ok=False)
        assert "gzip" in gzip_failure.stderr or "truncated" in gzip_failure.stderr

        # Force several in-shard radix spills and a map-side combine with the
        # test's 1 MiB sort budget (capacity is 32,768 CellRecords).
        spill = root / "spill.hbs.gz"
        hbs(spill, [(0, i % 20, 0, (i * 7) % 20, 1) for i in range(70000)])
        spill_stage = root / "spill-stage"
        spill_build = root / "spill-build"
        run([executable, "stage", "--chunk-records", "100000", spill, spill_stage])
        run([executable, "build-cells", "-r", "1", "--memory", "1MiB",
             "--fan-in", "2", spill_stage / "stage.manifest", spill_build])
        spill_run = next((spill_build / "cells").glob("*-cells.h10r"))
        assert sum(count for _, _, count in read_run(spill_run)[1]) == 70000


if __name__ == "__main__":
    main()
