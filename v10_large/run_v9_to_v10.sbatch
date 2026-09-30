#!/usr/bin/env bash
#SBATCH --job-name=v9to10
#SBATCH --partition=largemem,long,commons,weka,weka2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00

# End-to-end, restartable V9 -> HBS -> V10 pipeline. See README.md for usage.
# The first job uses Slurm's default log location. Child logs go to WORK_DIR/logs.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch run_v9_to_v10.sbatch --input INPUT.v9.hic --work-dir WORK_DIR
       (--hictools-repo REPO | --hic-bin EXECUTABLE)
       (--straw-repo REPO | --straw-bin EXECUTABLE)
       [--output OUTPUT.v10.hic] [--tmp-dir SCRATCH_DIR] [--genome NAME]

WORK_DIR and SCRATCH_DIR must be visible to all nodes. The output defaults to
WORK_DIR/INPUT.v10.hic; scratch defaults to WORK_DIR/scratch. Supply source
repositories to build the tools in WORK_DIR/tools, or supply existing binaries.
EOF
}

fail() { echo "ERROR: $*" >&2; exit 1; }
absolute() { realpath "$1"; }

PARTS=largemem,long,commons,weka,weka2
R=1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000
DERIVE=(--derive 2:1 --derive 5:1)
FAN_IN=128
CHUNK_RECORDS=250000000
THROTTLE_MAP=16
THROTTLE_BUILD=12
THROTTLE_NORM=8
THROTTLE_EXPECTED=4
THROTTLE_WRITE=12

PHASE=${PHASE:-setup}
if [[ $PHASE == setup ]]; then
    INPUT= WORK_DIR= OUTPUT= TMP_DIR= GENOME=hg38
    HICTOOLS_REPO= STRAW_REPO= HIC_BIN= STRAW_BIN=
    while (($#)); do
        case $1 in
            --input|--work-dir|--output|--tmp-dir|--genome|--hictools-repo|--straw-repo|--hic-bin|--straw-bin)
                (($# >= 2)) || fail "missing value for $1"
                case $1 in
                    --input) INPUT=$2 ;; --work-dir) WORK_DIR=$2 ;;
                    --output) OUTPUT=$2 ;; --tmp-dir) TMP_DIR=$2 ;;
                    --genome) GENOME=$2 ;; --hictools-repo) HICTOOLS_REPO=$2 ;;
                    --straw-repo) STRAW_REPO=$2 ;; --hic-bin) HIC_BIN=$2 ;;
                    --straw-bin) STRAW_BIN=$2 ;;
                esac
                shift 2 ;;
            -h|--help) usage; exit 0 ;;
            *) fail "unknown argument: $1" ;;
        esac
    done
    [[ -n $INPUT && -n $WORK_DIR ]] || { usage >&2; fail "--input and --work-dir are required"; }
    [[ -n $HICTOOLS_REPO || -n $HIC_BIN ]] || fail "supply --hictools-repo or --hic-bin"
    [[ -n $STRAW_REPO || -n $STRAW_BIN ]] || fail "supply --straw-repo or --straw-bin"
    [[ -z $HICTOOLS_REPO || -z $HIC_BIN ]] || fail "choose one hictools source or binary"
    [[ -z $STRAW_REPO || -z $STRAW_BIN ]] || fail "choose one straw source or binary"
    [[ -r $INPUT ]] || fail "input is not readable: $INPUT"
    INPUT=$(absolute "$INPUT")
    mkdir -p "$WORK_DIR"
    WORK_DIR=$(absolute "$WORK_DIR")
    RUN=$(basename "$INPUT" .hic)
    RUN=${RUN%.v9}
    [[ -n $OUTPUT ]] || OUTPUT="$WORK_DIR/${RUN}.v10.hic"
    mkdir -p "$(dirname "$OUTPUT")"
    OUTPUT=$(absolute "$(dirname "$OUTPUT")")/$(basename "$OUTPUT")
    [[ $OUTPUT != "$INPUT" ]] || fail "output must differ from input"
    [[ -n $TMP_DIR ]] || TMP_DIR="$WORK_DIR/scratch"
    mkdir -p "$TMP_DIR"
    TMP_DIR=$(absolute "$TMP_DIR")
    if [[ -n $HICTOOLS_REPO ]]; then
        HICTOOLS_REPO=$(absolute "$HICTOOLS_REPO")
        [[ -f $HICTOOLS_REPO/CMakeLists.txt ]] || fail "invalid hictools repository"
        HIC_BIN="$WORK_DIR/tools/hictools/hic_v10_large"
    else
        HIC_BIN=$(absolute "$HIC_BIN")
    fi
    if [[ -n $STRAW_REPO ]]; then
        STRAW_REPO=$(absolute "$STRAW_REPO")
        [[ -f $STRAW_REPO/C++/CMakeLists.txt ]] || fail "invalid straw repository"
        STRAW_BIN="$WORK_DIR/tools/straw/straw"
    else
        STRAW_BIN=$(absolute "$STRAW_BIN")
    fi
    export INPUT WORK_DIR OUTPUT TMP_DIR GENOME HICTOOLS_REPO STRAW_REPO HIC_BIN STRAW_BIN RUN
fi

[[ -n ${SLURM_JOB_ID:-} ]] || fail "run this script with sbatch"
: "${INPUT:?}" "${WORK_DIR:?}" "${OUTPUT:?}" "${TMP_DIR:?}" "${HIC_BIN:?}" "${STRAW_BIN:?}"
STAGE="$WORK_DIR/stage"
BUILD="$WORK_DIR/build"
VECTORS="$WORK_DIR/vectors"
PARTS_DIR="$WORK_DIR/pair-parts"
LOGS="$WORK_DIR/logs"
HBS="$WORK_DIR/${RUN}.hbs.gz"
SELF="$WORK_DIR/run_v9_to_v10.sbatch"

submit() {
    local phase=$1 name=$2 dependency=$3
    shift 3
    local args=(--parsable --job-name="$name" --partition="$PARTS"
                --output="$LOGS/%x_%A_%a.out" --error="$LOGS/%x_%A_%a.err"
                --export=ALL,PHASE="$phase")
    [[ -z $dependency ]] || args+=(--dependency="$dependency")
    sbatch "${args[@]}" "$@" "$SELF"
}

task_count() { awk -F '\t' -v type="$2" '$1 == type {n++} END {print n+0}' "$1"; }
task_row() {
    awk -F '\t' -v type="$2" -v idx="$SLURM_ARRAY_TASK_ID" \
        '$1 == type && n++ == idx {print $2; found=1; exit} END {if (!found) exit 1}' "$1"
}
scratch_dir() {
    local path="$TMP_DIR/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID:-0}-$1"
    mkdir -p "$path"
    printf '%s\n' "$path"
}

echo "=== phase=$PHASE job=$SLURM_JOB_ID host=$(hostname) $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
case $PHASE in
    setup)
        mkdir -p "$STAGE" "$BUILD" "$VECTORS" "$PARTS_DIR" "$LOGS" "$WORK_DIR/tools"
        # Slurm runs a spool copy of the script. Save that copy for child jobs.
        cp "$0" "$SELF"
        if [[ -n $STRAW_REPO ]]; then
            cmake -S "$STRAW_REPO/C++" -B "$WORK_DIR/tools/straw" \
                -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DSTRAW_BUILD_C_API=OFF \
                -DCMAKE_DISABLE_FIND_PACKAGE_zstd=ON
            cmake --build "$WORK_DIR/tools/straw" -j "${SLURM_CPUS_PER_TASK:-2}" --target straw_cli
        fi
        if [[ -n $HICTOOLS_REPO ]]; then
            cmake -S "$HICTOOLS_REPO" -B "$WORK_DIR/tools/hictools" \
                -DBUILD_V10=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
            cmake --build "$WORK_DIR/tools/hictools" -j "${SLURM_CPUS_PER_TASK:-2}" --target hic_v10_large
        fi
        [[ -x $HIC_BIN && -x $STRAW_BIN ]] || fail "tool binaries are missing or not executable"
        echo "input: $INPUT"
        echo "output: $OUTPUT"
        echo "work: $WORK_DIR"
        echo "scratch: $TMP_DIR"
        submit dump v9to10_dump '' --cpus-per-task=2 --mem=48G --time=3-00:00:00
        ;;
    dump)
        if [[ ! -s $HBS ]]; then
            partial="$WORK_DIR/.${RUN}.${SLURM_JOB_ID}.hbs.gz"
            trap 'rm -f "$partial"' EXIT
            /usr/bin/time -v "$STRAW_BIN" dump observed NONE "$INPUT" BP 1 "$partial"
            [[ -s $partial ]] || fail "straw produced an empty HBS file"
            mv "$partial" "$HBS"
            trap - EXIT
        fi
        ls -lh "$HBS"
        submit stage v9to10_stage '' --cpus-per-task=4 --mem=32G --time=3-00:00:00
        ;;
    stage)
        /usr/bin/time -v "$HIC_BIN" stage --chunk-records "$CHUNK_RECORDS" \
            --progress-records 500000000 "$HBS" "$STAGE"
        "$HIC_BIN" inspect --verify "$STAGE/stage.manifest"
        "$HIC_BIN" plan -r "$R" --root-only --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD" > "$WORK_DIR/build.tasks.tsv"
        n_map=$(task_count "$WORK_DIR/build.tasks.tsv" map-root)
        n_reduce=$(task_count "$WORK_DIR/build.tasks.tsv" reduce-root)
        n_pair=$(task_count "$WORK_DIR/build.tasks.tsv" build-pair)
        ((n_map > 0 && n_pair > 0)) || fail "empty build plan"
        jid=$(submit map_root v9to10_maproot '' --cpus-per-task=4 --mem=24G --time=1-00:00:00 \
            --array="0-$((n_map - 1))%$THROTTLE_MAP")
        if ((n_reduce > 0)); then
            jid=$(submit reduce_root v9to10_reduceroot "afterok:$jid" \
                --cpus-per-task=4 --mem=24G --time=1-00:00:00 \
                --array="0-$((n_reduce - 1))%$THROTTLE_MAP")
        fi
        jid=$(submit build_pair v9to10_buildpair "afterok:$jid" \
            --cpus-per-task=4 --mem=24G --time=2-00:00:00 \
            --array="0-$((n_pair - 1))%$THROTTLE_BUILD")
        jid=$(submit finalize_build v9to10_finbuild "afterok:$jid" \
            --cpus-per-task=2 --mem=16G --time=08:00:00)
        submit dispatch_normalize v9to10_dispnorm "afterok:$jid" \
            --cpus-per-task=2 --mem=8G --time=04:00:00
        ;;
    map_root)
        id=$(task_row "$WORK_DIR/build.tasks.tsv" map-root)
        "$HIC_BIN" map-root -r "$R" --root-only --memory 16GiB --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD" "$id"
        ;;
    reduce_root)
        entry=$(task_row "$WORK_DIR/build.tasks.tsv" reduce-root)
        "$HIC_BIN" reduce-root -r "$R" --root-only --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD" "${entry%%:*}" "${entry##*:}"
        ;;
    build_pair)
        id=$(task_row "$WORK_DIR/build.tasks.tsv" build-pair)
        "$HIC_BIN" build-pair -r "$R" --root-only --memory 16GiB --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD" "$id"
        ;;
    finalize_build)
        "$HIC_BIN" finalize-build -r "$R" --root-only --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD"
        ;;
    dispatch_normalize)
        "$HIC_BIN" plan-normalize --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD/build.manifest" "$VECTORS" > "$WORK_DIR/vectors.tasks.tsv"
        n_chr=$(task_count "$WORK_DIR/vectors.tasks.tsv" normalize-chr)
        n_res=$(task_count "$WORK_DIR/vectors.tasks.tsv" expected-res)
        ((n_chr > 0 && n_res > 0)) || fail "empty normalization plan"
        jid=$(submit normalize_chr v9to10_normchr '' --cpus-per-task=16 --mem=150G --time=5-00:00:00 \
            --array="0-$((n_chr - 1))%$THROTTLE_NORM")
        jid=$(submit expected_res v9to10_expected "afterok:$jid" \
            --cpus-per-task=16 --mem=120G --time=3-00:00:00 \
            --array="0-$((n_res - 1))%$THROTTLE_EXPECTED")
        jid=$(submit finalize_vectors v9to10_finvec "afterok:$jid" \
            --cpus-per-task=2 --mem=16G --time=08:00:00)
        submit dispatch_write v9to10_dispwrite "afterok:$jid" \
            --cpus-per-task=2 --mem=8G --time=04:00:00
        ;;
    normalize_chr|expected_res)
        if [[ $PHASE == normalize_chr ]]; then
            type=normalize-chr
            key=normchr
        else
            type=expected-res
            key=expected
        fi
        id=$(task_row "$WORK_DIR/vectors.tasks.tsv" "$type")
        tmp=$(scratch_dir "$key-$id")
        trap 'rm -rf "$tmp"' EXIT
        "$HIC_BIN" "$type" --memory 16GiB --fan-in "$FAN_IN" --tmp "$tmp" -t 16 \
            "$STAGE/stage.manifest" "$BUILD/build.manifest" "$VECTORS" "$id"
        ;;
    finalize_vectors)
        "$HIC_BIN" finalize-vectors --fan-in "$FAN_IN" \
            "$STAGE/stage.manifest" "$BUILD/build.manifest" "$VECTORS"
        ;;
    dispatch_write)
        "$HIC_BIN" plan-write "$STAGE/stage.manifest" "$BUILD/build.manifest" \
            "$PARTS_DIR" "$OUTPUT" > "$WORK_DIR/pair-parts.tasks.tsv"
        n_pair=$(task_count "$WORK_DIR/pair-parts.tasks.tsv" write-pair)
        ((n_pair > 0)) || fail "empty write plan"
        jid=$(submit write_pair v9to10_writepair '' --cpus-per-task=16 --mem=32G --time=3-00:00:00 \
            --array="0-$((n_pair - 1))%$THROTTLE_WRITE")
        jid=$(submit merge v9to10_merge "afterok:$jid" --cpus-per-task=16 --mem=64G --time=5-00:00:00)
        jid=$(submit validate v9to10_validate "afterok:$jid" \
            --cpus-per-task=8 --mem=64G --time=1-00:00:00)
        submit compare v9to10_compare "afterok:$jid" \
            --cpus-per-task=8 --mem=150G --time=3-00:00:00
        ;;
    write_pair)
        id=$(task_row "$WORK_DIR/pair-parts.tasks.tsv" write-pair)
        tmp=$(scratch_dir "write-$id")
        trap 'rm -rf "$tmp"' EXIT
        "$HIC_BIN" write-pair --genome "$GENOME" --memory 16GiB --fan-in "$FAN_IN" \
            "${DERIVE[@]}" --resolution-batch 4 --tmp "$tmp" -t 16 \
            "$STAGE/stage.manifest" "$BUILD/build.manifest" "$PARTS_DIR" "$id"
        ;;
    merge)
        /usr/bin/time -v "$HIC_BIN" merge-pairs --vectors "$VECTORS/vectors.manifest" \
            --level 6 -t 16 "$STAGE/stage.manifest" "$BUILD/build.manifest" "$PARTS_DIR" "$OUTPUT"
        ;;
    validate)
        /usr/bin/time -v "$HIC_BIN" validate-v10 "$OUTPUT"
        old_size=$(wc -c < "$INPUT")
        new_size=$(wc -c < "$OUTPUT")
        awk -v a="$old_size" -v b="$new_size" 'BEGIN {
            printf "V9: %.0f bytes (%.1f GiB)\nV10: %.0f bytes (%.1f GiB)\n", a, a/2^30, b, b/2^30
            if (b > 0) printf "V9/V10 size ratio: %.2fx\n", a/b
        }'
        ;;
    compare)
        # V9 and V10 may differ in vector padding and SCALE solver details.
        # Preserve compare output, including non-cancelling cell counts, for review.
        "$STRAW_BIN" compare "$INPUT" "$OUTPUT" \
            --resolution 1000000 --resolution 500000 --resolution 250000 --resolution 100000 \
            --all --norm VC --norm SCALE --no-timing || echo "coarse comparison reported differences or an error"
        "$STRAW_BIN" compare "$INPUT" "$OUTPUT" \
            --resolution 100000 --resolution 10000 --resolution 1000 \
            --resolution 100 --resolution 10 --resolution 1 \
            --norm VC --samples 12 --min-window-bins 256 --max-window-bins 4096 \
            --seed 42 --no-timing || echo "guided comparison reported differences or an error"
        ;;
    *) fail "unknown phase: $PHASE" ;;
esac
echo "=== phase=$PHASE done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
