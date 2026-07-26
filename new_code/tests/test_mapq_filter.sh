#!/usr/bin/env bash
set -euo pipefail

hic_pre=$1
input=$2
workdir=$(mktemp -d "${TMPDIR:-/tmp}/hictools-mapq-test.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

"$hic_pre" -q 30 -r 100000 -t 1 \
    "$input" "$workdir/output.hic" hg38 \
    >"$workdir/stdout.log" 2>"$workdir/stderr.log"

test -s "$workdir/output.hic"
grep -q "total pairs read: 2, skipped: 1" "$workdir/stderr.log"
