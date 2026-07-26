#!/usr/bin/env python3
"""Structural, normalization-vector, and straw-dump checks for V9 .hic files."""

import argparse
import json
import math
import os
import statistics
import struct
import subprocess
import sys
from pathlib import Path


class Reader:
    def __init__(self, path):
        self.path = Path(path)
        self.fh = self.path.open("rb")

    def close(self):
        self.fh.close()

    def seek(self, pos):
        self.fh.seek(pos)

    def tell(self):
        return self.fh.tell()

    def take(self, size):
        data = self.fh.read(size)
        if len(data) != size:
            raise EOFError(f"{self.path}: wanted {size} bytes at {self.tell() - len(data)}")
        return data

    def i32(self):
        return struct.unpack("<i", self.take(4))[0]

    def i64(self):
        return struct.unpack("<q", self.take(8))[0]

    def f32(self):
        return struct.unpack("<f", self.take(4))[0]

    def string(self):
        chunks = bytearray()
        while True:
            value = self.take(1)
            if value == b"\0":
                return chunks.decode("utf-8")
            chunks.extend(value)


def skip_expected(reader, normalized):
    count = reader.i32()
    entries = []
    for _ in range(count):
        norm = reader.string() if normalized else "NONE"
        unit = reader.string()
        resolution = reader.i32()
        n_values = reader.i64()
        reader.seek(reader.tell() + 4 * n_values)
        n_factors = reader.i32()
        reader.seek(reader.tell() + 8 * n_factors)
        entries.append((norm, unit, resolution, n_values, n_factors))
    return entries


def parse_hic(path, load_vectors=True):
    reader = Reader(path)
    try:
        if reader.take(4) != b"HIC\0":
            raise ValueError(f"{path}: bad HIC magic")
        version = reader.i32()
        if version != 9:
            raise ValueError(f"{path}: expected V9, found V{version}")
        footer = reader.i64()
        genome = reader.string()
        nvi_position = reader.i64()
        nvi_length = reader.i64()
        attributes = {}
        for _ in range(reader.i32()):
            key = reader.string()
            value = reader.string()
            attributes[key] = value
        chromosomes = []
        for index in range(reader.i32()):
            chromosomes.append((index, reader.string(), reader.i64()))
        resolutions = [reader.i32() for _ in range(reader.i32())]
        frag_resolutions = [reader.i32() for _ in range(reader.i32())]

        file_size = os.path.getsize(path)
        if not (0 < footer < file_size):
            raise ValueError(f"{path}: footer offset {footer} is outside file")
        reader.seek(footer)
        v5_size = reader.i64()
        v5_end = reader.tell() + v5_size
        master = []
        for _ in range(reader.i32()):
            key, position, size = reader.string(), reader.i64(), reader.i32()
            if not (0 <= position < footer):
                raise ValueError(f"{path}: matrix {key} points outside body")
            master.append((key, position, size))
        expected = skip_expected(reader, normalized=False)
        if reader.tell() != v5_end:
            raise ValueError(
                f"{path}: V5 footer length mismatch, expected end {v5_end}, read {reader.tell()}"
            )
        normalized_expected = skip_expected(reader, normalized=True)

        if nvi_position <= 0 or nvi_length <= 0:
            norm_index = []
        else:
            if nvi_position + nvi_length > file_size:
                raise ValueError(f"{path}: normalization index is outside file")
            reader.seek(nvi_position)
            norm_index = []
            seen = set()
            for _ in range(reader.i32()):
                entry = {
                    "type": reader.string(),
                    "chr": reader.i32(),
                    "unit": reader.string(),
                    "resolution": reader.i32(),
                    "position": reader.i64(),
                    "size": reader.i64(),
                }
                key = (entry["type"], entry["chr"], entry["unit"], entry["resolution"])
                if key in seen:
                    raise ValueError(f"{path}: duplicate normalization entry {key}")
                seen.add(key)
                if entry["position"] < 0 or entry["position"] + entry["size"] > file_size:
                    raise ValueError(f"{path}: normalization entry {key} is outside file")
                norm_index.append(entry)

            expected_nvi_end = nvi_position + nvi_length
            if reader.tell() != expected_nvi_end:
                raise ValueError(
                    f"{path}: NVI length mismatch, expected end {expected_nvi_end}, read {reader.tell()}"
                )

        vectors = {}
        if load_vectors:
            chr_lengths = {idx: length for idx, _name, length in chromosomes}
            for entry in norm_index:
                reader.seek(entry["position"])
                n_values = reader.i64()
                if entry["size"] != 8 + 4 * n_values:
                    raise ValueError(f"{path}: bad byte count for normalization entry {entry}")
                values = list(struct.unpack(f"<{n_values}f", reader.take(4 * n_values)))
                if entry["chr"] != 0 and entry["unit"] == "BP":
                    expected_length = chr_lengths[entry["chr"]] // entry["resolution"] + 1
                    if n_values != expected_length:
                        raise ValueError(
                            f"{path}: {entry['type']} chr {entry['chr']} "
                            f"{entry['resolution']} has {n_values} values; expected {expected_length}"
                        )
                vectors[(entry["type"], entry["chr"], entry["unit"], entry["resolution"])] = values

        return {
            "path": str(Path(path).resolve()),
            "size": file_size,
            "version": version,
            "genome": genome,
            "footer": footer,
            "nvi_position": nvi_position,
            "nvi_length": nvi_length,
            "attributes": attributes,
            "chromosomes": chromosomes,
            "resolutions": resolutions,
            "frag_resolutions": frag_resolutions,
            "master_entries": len(master),
            "expected": expected,
            "normalized_expected": normalized_expected,
            "norm_index": norm_index,
            "vectors": vectors,
        }
    finally:
        reader.close()


def correlation(xs, ys):
    if len(xs) < 2:
        return float("nan")
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    xx = sum((x - mx) ** 2 for x in xs)
    yy = sum((y - my) ** 2 for y in ys)
    if xx == 0 or yy == 0:
        return 1.0 if xs == ys else float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(xx * yy)


def vector_metrics(left, right):
    pairs = [
        (a, b)
        for a, b in zip(left, right)
        if math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0
    ]
    xs = [x for x, _ in pairs]
    ys = [y for _, y in pairs]
    rel = [abs(x - y) / max(abs(x), abs(y), 1e-30) for x, y in pairs]
    ratios = [y / x for x, y in pairs]
    return {
        "length_left": len(left),
        "length_right": len(right),
        "common_positive_finite": len(pairs),
        "pearson": correlation(xs, ys),
        "median_ratio_cpp_over_java": statistics.median(ratios) if ratios else float("nan"),
        "median_relative_error": statistics.median(rel) if rel else float("nan"),
        "max_relative_error": max(rel) if rel else float("nan"),
    }


def compare_hic(java_path, cpp_path):
    java = parse_hic(java_path)
    cpp = parse_hic(cpp_path)
    structural = {
        "version_equal": java["version"] == cpp["version"],
        "genome_equal": java["genome"] == cpp["genome"],
        "chromosomes_equal": java["chromosomes"] == cpp["chromosomes"],
        "resolutions_equal": java["resolutions"] == cpp["resolutions"],
        "fragment_resolutions_equal": java["frag_resolutions"] == cpp["frag_resolutions"],
        "master_entry_count_java": java["master_entries"],
        "master_entry_count_cpp": cpp["master_entries"],
    }
    vector_rows = []
    common = sorted(set(java["vectors"]) & set(cpp["vectors"]))
    for key in common:
        metrics = vector_metrics(java["vectors"][key], cpp["vectors"][key])
        vector_rows.append({"key": key, **metrics})
    return {
        "structural": structural,
        "java_norm_count": len(java["vectors"]),
        "cpp_norm_count": len(cpp["vectors"]),
        "missing_from_cpp": sorted(set(java["vectors"]) - set(cpp["vectors"])),
        "extra_in_cpp": sorted(set(cpp["vectors"]) - set(java["vectors"])),
        "vectors": vector_rows,
    }


def canonicalize_dump(text):
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"unexpected straw output: {line}")
        rows.append((int(fields[0]), int(fields[1]), float(fields[2])))
    rows.sort(key=lambda row: (row[0], row[1]))
    return rows


def write_dump(path, rows):
    with Path(path).open("w", encoding="utf-8") as out:
        for x, y, value in rows:
            out.write(f"{x}\t{y}\t{value:.12g}\n")


def run_straw(args):
    if args.norm != "NONE":
        metadata = parse_hic(args.hic, load_vectors=False)
        name_to_index = {name: index for index, name, _length in metadata["chromosomes"]}
        chr1 = args.chr1.split(":", 1)[0]
        chr2 = args.chr2.split(":", 1)[0]
        available = {
            (entry["type"], entry["chr"], entry["unit"], entry["resolution"])
            for entry in metadata["norm_index"]
        }
        for chrom in (chr1, chr2):
            if chrom not in name_to_index:
                raise RuntimeError(f"chromosome {chrom} is not present in {args.hic}")
            key = (args.norm, name_to_index[chrom], "BP", args.resolution)
            if key not in available:
                raise RuntimeError(f"normalization vector is not present: {key}")
    command = [
        args.straw, args.matrix, args.norm, args.hic,
        args.chr1, args.chr2, "BP", str(args.resolution),
    ]
    result = subprocess.run(
        command, text=True, capture_output=True, timeout=args.timeout, check=False
    )
    Path(args.stderr).write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"straw exited {result.returncode}: {' '.join(command)}")
    rows = canonicalize_dump(result.stdout)
    write_dump(args.out, rows)
    print(len(rows))


def compare_dumps(left_path, right_path, norm):
    left = canonicalize_dump(Path(left_path).read_text(encoding="utf-8"))
    right = canonicalize_dump(Path(right_path).read_text(encoding="utf-8"))
    left_map = {(x, y): value for x, y, value in left}
    right_map = {(x, y): value for x, y, value in right}
    common_keys = sorted(set(left_map) & set(right_map))
    metrics = vector_metrics(
        [left_map[key] for key in common_keys],
        [right_map[key] for key in common_keys],
    )
    metrics.update({
        "norm": norm,
        "records_left": len(left),
        "records_right": len(right),
        "coordinate_sets_equal": set(left_map) == set(right_map),
        "exact_values_equal": left == right,
    })
    return metrics


def json_default(value):
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    inspect_parser = sub.add_parser("inspect")
    inspect_parser.add_argument("hic")

    compare_parser = sub.add_parser("compare-hic")
    compare_parser.add_argument("java_hic")
    compare_parser.add_argument("cpp_hic")

    straw_parser = sub.add_parser("straw")
    straw_parser.add_argument("--straw", required=True)
    straw_parser.add_argument("--hic", required=True)
    straw_parser.add_argument("--matrix", default="observed")
    straw_parser.add_argument("--norm", required=True)
    straw_parser.add_argument("--chr1", required=True)
    straw_parser.add_argument("--chr2", required=True)
    straw_parser.add_argument("--resolution", type=int, required=True)
    straw_parser.add_argument("--out", required=True)
    straw_parser.add_argument("--stderr", required=True)
    straw_parser.add_argument("--timeout", type=int, default=300)

    dump_parser = sub.add_parser("compare-dumps")
    dump_parser.add_argument("left")
    dump_parser.add_argument("right")
    dump_parser.add_argument("--norm", required=True)

    args = parser.parse_args()
    if args.command == "inspect":
        result = parse_hic(args.hic, load_vectors=False)
        result.pop("vectors", None)
    elif args.command == "compare-hic":
        result = compare_hic(args.java_hic, args.cpp_hic)
    elif args.command == "straw":
        run_straw(args)
        return
    else:
        result = compare_dumps(args.left, args.right, args.norm)
    print(json.dumps(result, indent=2, sort_keys=True, default=json_default, allow_nan=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
