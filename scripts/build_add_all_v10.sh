#!/usr/bin/env bash
# Build on the cluster (Linux binary); needs CMake, C++17, zlib and zstd headers.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${1:-"$repo/build/v10-all"}
jobs=${BUILD_JOBS:-4}
cmake -S "$repo" -B "$build_dir" -DBUILD_V10=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" --target hic_v10 --parallel "$jobs"
printf 'Built: %s/hic_v10\n' "$build_dir"
