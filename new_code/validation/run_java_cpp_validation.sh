#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
INPUT="${1:-$REPO_ROOT/test.txt}"
OUTPUT_PARENT="${2:-$REPO_ROOT/new_code/validation-results}"
THREADS="${THREADS:-4}"
JAVA_HEAP="${JAVA_HEAP:-12g}"
JAVA_PRE_THREADS="${JAVA_PRE_THREADS:-1}"
RESOLUTIONS="${RESOLUTIONS:-2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000}"
STRAW_TIMEOUT="${STRAW_TIMEOUT:-300}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="$OUTPUT_PARENT/$RUN_ID"

CPP_BUILD="$REPO_ROOT/new_code/build-validation"
STRAW_BUILD="$REPO_ROOT/Straw_C++/build-validation"
JAVA_BUILD="$REPO_ROOT/Java_Code/build-validation"
JAVA_LIB="$REPO_ROOT/Java_Code/lib/general"
JAVA_CP="$JAVA_BUILD/hic_tools.jar:$JAVA_LIB/*"
CPP_PRE="$CPP_BUILD/hic_pre"
CPP_ADDNORM="$CPP_BUILD/hic_addnorm"
STRAW="$STRAW_BUILD/straw"
TOOLS="$SCRIPT_DIR/compare_hic.py"

mkdir -p "$OUT"/{hic,logs,dumps,reports,tmp,build}

CURRENT_STEP="initialization"
CURRENT_LOG=""
report_failure() {
  status=$?
  trap - ERR
  echo >&2
  echo "Validation failed during: $CURRENT_STEP (exit $status)" >&2
  if [[ -n "$CURRENT_LOG" && -f "$CURRENT_LOG" ]]; then
    echo "Last 40 lines of $CURRENT_LOG:" >&2
    tail -n 40 "$CURRENT_LOG" >&2
  fi
  echo "All logs: $OUT/logs/" >&2
  exit "$status"
}
trap report_failure ERR

echo "Validation output: $OUT"
echo "Input: $INPUT"
echo "Resolutions: $RESOLUTIONS"
echo "C++ threads: $THREADS; Java pre threads: $JAVA_PRE_THREADS; Java heap: $JAVA_HEAP"

if [[ ! -s "$INPUT" ]]; then
  echo "Input is missing or empty: $INPUT" >&2
  exit 2
fi

echo "[1/7] Building C++ pre/addnorm"
CURRENT_STEP="[1/7] Building C++ pre/addnorm"
CURRENT_LOG="$OUT/logs/cmake_cpp.log"
cmake -S "$REPO_ROOT/new_code" -B "$CPP_BUILD" -DCMAKE_BUILD_TYPE=Release \
  >"$OUT/logs/cmake_cpp.log" 2>&1
cmake --build "$CPP_BUILD" -j "$THREADS" >>"$OUT/logs/cmake_cpp.log" 2>&1

echo "[2/7] Building straw"
CURRENT_STEP="[2/7] Building straw"
CURRENT_LOG="$OUT/logs/cmake_straw.log"
cmake -S "$REPO_ROOT/Straw_C++" -B "$STRAW_BUILD" -DCMAKE_BUILD_TYPE=Release \
  >"$OUT/logs/cmake_straw.log" 2>&1
cmake --build "$STRAW_BUILD" -j "$THREADS" >>"$OUT/logs/cmake_straw.log" 2>&1

echo "[3/7] Building Java reference"
CURRENT_STEP="[3/7] Building Java reference"
CURRENT_LOG="$OUT/logs/javac.log"
mkdir -p "$JAVA_BUILD/classes"
find "$JAVA_BUILD/classes" -type f -name '*.class' -delete
find "$REPO_ROOT/Java_Code/src" -name '*.java' -print0 |
  xargs -0 javac -source 8 -target 8 -encoding UTF-8 \
    -cp "$JAVA_LIB/*" -d "$JAVA_BUILD/classes" \
    >"$OUT/logs/javac.log" 2>&1
jar cfe "$JAVA_BUILD/hic_tools.jar" hic.tools.HiCTools -C "$JAVA_BUILD/classes" .
unzip -p "$JAVA_LIB/java-straw.2.21.00.jar" \
  javastraw/reader/basics/chrom/sizes/hg19.chrom.sizes \
  >"$OUT/build/hg19.chrom.sizes"

{
  echo "run_id=$RUN_ID"
  echo "input=$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
  echo "input_bytes=$(wc -c < "$INPUT" | tr -d ' ')"
  echo "input_lines=$(wc -l < "$INPUT" | tr -d ' ')"
  echo "resolutions=$RESOLUTIONS"
  echo "threads=$THREADS"
  echo "java_pre_threads=$JAVA_PRE_THREADS"
  echo "java_heap=$JAVA_HEAP"
  java -version 2>&1 | sed 's/^/java_version=/'
  c++ --version 2>&1 | head -n 1 | sed 's/^/cxx_version=/'
} >"$OUT/reports/manifest.txt"

echo "[4/7] Building Java .hic (pre without norms, then VC/VC_SQRT/SCALE addNorm)"
CURRENT_STEP="[4/7] Building Java .hic (pre)"
CURRENT_LOG="$OUT/logs/java_pre.log"
java "-Xmx$JAVA_HEAP" -cp "$JAVA_CP" hic.tools.HiCTools \
  pre -n -r "$RESOLUTIONS" -j "$JAVA_PRE_THREADS" \
  "$INPUT" "$OUT/hic/java.hic" hg19 \
  >"$OUT/logs/java_pre.log" 2>&1
CURRENT_STEP="[4/7] Building Java .hic (addNorm)"
CURRENT_LOG="$OUT/logs/java_addnorm.log"
java "-Xmx$JAVA_HEAP" -cp "$JAVA_CP" hic.tools.HiCTools \
  addNorm -k VC,VC_SQRT,SCALE -r 1000,1000,1000 -j "$THREADS" \
  "$OUT/hic/java.hic" \
  >"$OUT/logs/java_addnorm.log" 2>&1

echo "[5/7] Building C++ .hic (pre, then VC/VC_SQRT/SCALE addnorm)"
CURRENT_STEP="[5/7] Building C++ .hic (pre)"
CURRENT_LOG="$OUT/logs/cpp_pre.log"
"$CPP_PRE" -r "$RESOLUTIONS" -t "$THREADS" -T "$OUT/tmp" -g hg19 -f mnd \
  "$INPUT" "$OUT/hic/cpp.hic" "$OUT/build/hg19.chrom.sizes" \
  >"$OUT/logs/cpp_pre.log" 2>&1
CURRENT_STEP="[5/7] Building C++ .hic (addnorm)"
CURRENT_LOG="$OUT/logs/cpp_addnorm.log"
"$CPP_ADDNORM" -t "$THREADS" "$OUT/hic/cpp.hic" \
  >"$OUT/logs/cpp_addnorm.log" 2>&1

echo "[6/7] Structural and direct normalization-vector checks"
CURRENT_STEP="[6/7] Structural and direct normalization-vector checks"
CURRENT_LOG=""
python3 "$TOOLS" inspect "$OUT/hic/java.hic" >"$OUT/reports/java_structure.json"
python3 "$TOOLS" inspect "$OUT/hic/cpp.hic" >"$OUT/reports/cpp_structure.json"
python3 "$TOOLS" compare-hic "$OUT/hic/java.hic" "$OUT/hic/cpp.hic" \
  >"$OUT/reports/hic_comparison.json"
md5 "$INPUT" "$OUT/hic/java.hic" "$OUT/hic/cpp.hic" \
  >"$OUT/reports/file_md5.txt"

python3 - "$OUT/build/hg19.chrom.sizes" "$RESOLUTIONS" >"$OUT/reports/queries.tsv" <<'PY'
import random
import sys

sizes_path, resolutions_text = sys.argv[1], sys.argv[2]
sizes = {}
with open(sizes_path, encoding="utf-8") as handle:
    for line in handle:
        chrom, length = line.split()
        sizes[chrom] = int(length)
resolutions = [int(value) for value in resolutions_text.split(",")]
rng = random.Random(8675309)
chroms = [str(i) for i in range(1, 23)] + ["X"]
print("query_id\tresolution\tchr1_region\tchr2_region")
query_id = 0
for resolution in resolutions:
    window = max(25 * resolution, 2_000_000)
    for mode in ("intra", "inter"):
        query_id += 1
        chrom1 = rng.choice(chroms)
        chrom2 = chrom1 if mode == "intra" else rng.choice([c for c in chroms if c != chrom1])
        start1 = rng.randrange(0, max(1, sizes[chrom1] - window))
        start1 = start1 // resolution * resolution
        if mode == "intra":
            start2 = min(sizes[chrom2] - window, max(0, start1 + rng.randint(-5, 5) * window))
        else:
            start2 = rng.randrange(0, max(1, sizes[chrom2] - window))
        start2 = start2 // resolution * resolution
        end1 = min(sizes[chrom1], start1 + window)
        end2 = min(sizes[chrom2], start2 + window)
        print(f"q{query_id:03d}_{mode}\t{resolution}\t{chrom1}:{start1}:{end1}\t{chrom2}:{start2}:{end2}")
PY

echo "[7/7] Deterministic random-region straw comparisons"
CURRENT_STEP="[7/7] Deterministic random-region straw comparisons"
CURRENT_LOG="$OUT/reports/straw_comparisons.tsv"
printf 'query_id\tmatrix\tnorm\tstatus\tcomparison_json\n' >"$OUT/reports/straw_comparisons.tsv"
tail -n +2 "$OUT/reports/queries.tsv" |
while IFS=$'\t' read -r query_id resolution chr1_region chr2_region; do
  mkdir -p "$OUT/dumps/$query_id"
  for matrix in observed oe; do
   for norm in NONE VC VC_SQRT SCALE; do
    java_dump="$OUT/dumps/$query_id/java_${matrix}_${norm}.tsv"
    cpp_dump="$OUT/dumps/$query_id/cpp_${matrix}_${norm}.tsv"
    java_err="$OUT/dumps/$query_id/java_${matrix}_${norm}.stderr"
    cpp_err="$OUT/dumps/$query_id/cpp_${matrix}_${norm}.stderr"
    status=ok
    if ! python3 "$TOOLS" straw --straw "$STRAW" --hic "$OUT/hic/java.hic" \
      --matrix "$matrix" --norm "$norm" --chr1 "$chr1_region" --chr2 "$chr2_region" \
      --resolution "$resolution" --out "$java_dump" --stderr "$java_err" \
      --timeout "$STRAW_TIMEOUT" >/dev/null 2>>"$java_err.runner"; then
      status=java_dump_failed
    fi
    if ! python3 "$TOOLS" straw --straw "$STRAW" --hic "$OUT/hic/cpp.hic" \
      --matrix "$matrix" --norm "$norm" --chr1 "$chr1_region" --chr2 "$chr2_region" \
      --resolution "$resolution" --out "$cpp_dump" --stderr "$cpp_err" \
      --timeout "$STRAW_TIMEOUT" >/dev/null 2>>"$cpp_err.runner"; then
      if [[ "$status" == ok ]]; then
        status=cpp_dump_failed
      else
        status="${status}+cpp_dump_failed"
      fi
    fi
    comparison="$OUT/dumps/$query_id/compare_${matrix}_${norm}.json"
    if [[ "$status" == ok ]]; then
      python3 "$TOOLS" compare-dumps "$java_dump" "$cpp_dump" --norm "$norm" \
        >"$comparison"
    else
      printf '{"status":"%s"}\n' "$status" >"$comparison"
    fi
    compact="$(tr -d '\n' < "$comparison")"
    printf '%s\t%s\t%s\t%s\t%s\n' "$query_id" "$matrix" "$norm" "$status" "$compact" \
      >>"$OUT/reports/straw_comparisons.tsv"
   done
  done
done

find "$OUT/dumps" -type f -name '*.tsv' -print0 |
  xargs -0 md5 >"$OUT/reports/dump_md5.txt"

echo
echo "Validation finished. Review:"
echo "  $OUT/reports/hic_comparison.json"
echo "  $OUT/reports/straw_comparisons.tsv"
echo "  $OUT/reports/file_md5.txt"
echo "  $OUT/logs/"
