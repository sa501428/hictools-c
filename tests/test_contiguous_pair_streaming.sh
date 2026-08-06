#!/usr/bin/env bash
set -euo pipefail

HIC_PRE="$1"
TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/hic_pre_streaming_test.XXXXXX")"
trap 'rm -rf "$TEST_ROOT"' EXIT

CHROM_SIZES="$TEST_ROOT/chrom.sizes"
SORTED_MND="$TEST_ROOT/sorted.mnd"
REPEATED_MND="$TEST_ROOT/repeated.mnd"
SUCCESS_TMP="$TEST_ROOT/success-tmp"
FAILURE_TMP="$TEST_ROOT/failure-tmp"
mkdir -p "$SUCCESS_TMP" "$FAILURE_TMP"

for chromosome in $(seq 1 12); do
  printf '%s\t1000000\n' "$chromosome"
done >"$CHROM_SIZES"

# 78 normalized chromosome pairs: enough to exceed a deliberately low file
# descriptor limit if phase 1 accidentally retains one stream per pair.
for chr1 in $(seq 1 12); do
  for chr2 in $(seq "$chr1" 12); do
    printf '0\t%s\t100\t0\t0\t%s\t200\t1\n' "$chr1" "$chr2"
  done
done >"$SORTED_MND"

(
  ulimit -n 32
  "$HIC_PRE" -r 100000 -t 4 -T "$SUCCESS_TMP" -f mnd \
    "$SORTED_MND" "$TEST_ROOT/sorted.hic" "$CHROM_SIZES"
)

test -s "$TEST_ROOT/sorted.hic"
test -z "$(find "$SUCCESS_TMP" -type f -print -quit)"

{
  printf '0\t1\t100\t0\t0\t1\t200\t1\n'
  printf '0\t1\t100\t0\t0\t2\t200\t1\n'
  printf '0\t1\t300\t0\t0\t1\t400\t1\n'
} >"$REPEATED_MND"

if "$HIC_PRE" -r 100000 -t 2 -T "$FAILURE_TMP" -f mnd \
    "$REPEATED_MND" "$TEST_ROOT/repeated.hic" "$CHROM_SIZES" \
    >"$TEST_ROOT/repeated.stdout" 2>"$TEST_ROOT/repeated.stderr"; then
  echo "Expected repeated chromosome-pair input to fail" >&2
  exit 1
fi

grep -q "appears in multiple blocks" "$TEST_ROOT/repeated.stderr"
test -z "$(find "$FAILURE_TMP" -type f -print -quit)"
