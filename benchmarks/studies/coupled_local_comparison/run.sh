#!/usr/bin/env bash
# Coupled-case comparison on this host. See README.md for what each stage measures.
#
#   run.sh A <outdir>   32 spacecraft: code before vs after PR #222, every route, thread ladder
#   run.sh B <outdir>   spacecraft-count sweep on the post-#222 code
#
# PRE_TREE and POST_TREE are checkouts of the two commits. The harness files from
# this checkout are copied into both, because the case and the per-call state
# dumps it needs postdate those commits; the harness itself is identical at both
# commits, so the copy is exact.
#
# Every point is one controller invocation with its own state-dump directory,
# capped at 16 GB. A point that finished is skipped on rerun. Two guards keep
# other work off the CPU while a point is timed:
#   1. it holds every lock in LOCKS, the files other sessions serialize their
#      Julia jobs on (the agent lock and the paper repository's heavy-job lock);
#   2. once it holds them, it refuses to start while any Julia process is using
#      CPU, which catches jobs that take no lock at all, or while MemAvailable is
#      under MIN_AVAILABLE_GB. It then releases the locks, waits QUIET_RETRY_S
#      and tries again, for up to QUIET_WAIT_S.
set -euo pipefail

QUIET_SAMPLE_S="${QUIET_SAMPLE_S:-3}"     # CPU-use sampling window
QUIET_MAX_CORES="${QUIET_MAX_CORES:-0.1}" # a Julia process above this many cores counts as busy
MIN_AVAILABLE_GB="${MIN_AVAILABLE_GB:-24}" # below this MemAvailable the point waits (swap distorts timings)

# Prints "pid cores command" for every Julia process using more than
# QUIET_MAX_CORES over QUIET_SAMPLE_S, measured from /proc utime+stime deltas
# (ps %CPU is a lifetime average, so an idle old process can look busy and a
# job that just started can look idle).
busy_julia() {
    local tck pid; tck=$(getconf CLK_TCK)
    declare -A t0
    for pid in $(pgrep -x julia); do
        t0[$pid]=$(awk '{print $14 + $15}' "/proc/$pid/stat" 2>/dev/null) || true
    done
    sleep "$QUIET_SAMPLE_S"
    for pid in "${!t0[@]}"; do
        [[ -n "${t0[$pid]}" ]] || continue
        local t1; t1=$(awk '{print $14 + $15}' "/proc/$pid/stat" 2>/dev/null) || continue
        [[ -n "$t1" ]] || continue
        awk -v d=$((t1 - t0[$pid])) -v tck="$tck" -v s="$QUIET_SAMPLE_S" -v max="$QUIET_MAX_CORES" \
            -v pid="$pid" -v cmd="$(tr '\0' ' ' < "/proc/$pid/cmdline" | cut -c1-150)" \
            'BEGIN { c = d / tck / s; if (c > max) printf "%s %.2f %s\n", pid, c, cmd }'
    done
}

# Runs with every lock held: check the CPU, then run the point. Exit 75 = busy.
if [[ "${1:-}" == __locked_point__ ]]; then
    shift
    tree=$1 dir=$2 case=$3 mode=$4 t=$5
    busy=$(busy_julia)
    avail_gb=$(awk '/^MemAvailable:/ { printf "%d", $2 / 1048576 }' /proc/meminfo)
    if (( avail_gb < MIN_AVAILABLE_GB )); then
        busy+="${busy:+$'\n'}MemAvailable ${avail_gb} GB < MIN_AVAILABLE_GB ${MIN_AVAILABLE_GB} GB"
    fi
    if [[ -n "$busy" ]]; then
        echo "[busy] $(date -Is) Julia processes using CPU:" >&2
        echo "$busy" >&2
        exit 75
    fi
    echo "[quiet] $(date -Is) no busy Julia process; starting" >&2
    cd "$tree"
    exec env SPACEAGORA_PPC_DUMP_STATE_DIR="$dir/states" \
        systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
        julia --startup-file=no --project=. benchmarks/studies/parallelization_performance.jl "$PROFILE" \
            --outdir="$dir" --cases="$case" --modes="$mode" --threads="$t" \
            --repeats="$REPEATS" --warmup="$WARMUP" --parity-cases=none
fi

STAGE="${1:?stage A or B}"
OUT="$(realpath -m "${2:?output directory}")"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
SRC_ROOT="$(cd "$HERE/../../.." && pwd)"
PRE_TREE="${PRE_TREE:?checkout of ea42d2946 (main before PR #222)}"
POST_TREE="${POST_TREE:?checkout of 402018b98 (PR #222 merge)}"
export PROFILE="${PROFILE:-smoke}"   # harness profile: smoke = 120 s mission, test = 10 s (plumbing check)
export REPEATS="${REPEATS:-6}"
export WARMUP="${WARMUP:-3}"
export QUIET_SAMPLE_S QUIET_MAX_CORES MIN_AVAILABLE_GB
THREADS="${THREADS:-2 4 8 12}"   # parallel rungs; 12 = physical cores of the reference workstation
ROUTES="${ROUTES:-inner_only rhs_satellite rhs_per_satellite rhs_flat predictive}"
LOCKS="${LOCKS:-/tmp/claude-1000/spaceagora-agents-julia.lock /tmp/claude-1000/-home-space-falcon-1-Documents-JAIS-2026-SpaceAGORA/julia-local-heavy.lock}"
QUIET_RETRY_S="${QUIET_RETRY_S:-120}"
QUIET_WAIT_S="${QUIET_WAIT_S:-28800}"
export OPENBLAS_NUM_THREADS=1 GKSwstype=100

# Nested `flock <lock>` prefixes, always taken in the order LOCKS lists them.
lock_cmd=()
for lock in $LOCKS; do
    mkdir -p "$(dirname "$lock")"
    lock_cmd+=(flock "$lock")
done

# The harness is copied FROM this checkout, so it must be one that defines the
# case; launching run.sh from a tree with the stock harness would overwrite the
# other tree's harness with it.
grep -q 'e6_actuated_saved' "$SRC_ROOT/benchmarks/studies/parallelization_performance/cases.jl" || {
    echo "run.sh: $SRC_ROOT's harness does not define stack<N>_e6_actuated_saved; run it from a checkout that does" >&2
    exit 2
}
for tree in "$PRE_TREE" "$POST_TREE"; do
    [[ "$(realpath "$tree")" == "$(realpath "$SRC_ROOT")" ]] && continue   # running from that tree
    for f in cases.jl cli.jl execution.jl; do
        cp "$SRC_ROOT/benchmarks/studies/parallelization_performance/$f" \
           "$tree/benchmarks/studies/parallelization_performance/$f"
    done
done

point() {  # tree tag case mode threads
    local tree=$1 tag=$2 case=$3 mode=$4 t=$5
    local dir="$OUT/$tag/${case}_${mode}_t${t}"
    [[ -f "$dir/done" ]] && { echo "[skip] $tag $case $mode t=$t"; return 0; }
    mkdir -p "$dir"
    echo "[run]  $tag $case $mode t=$t  $(date -Is)"
    local start=$SECONDS rc
    while :; do
        rc=0
        "${lock_cmd[@]}" "$SELF" __locked_point__ "$tree" "$dir" "$case" "$mode" "$t" > "$dir/run.log" 2>&1 || rc=$?
        [[ $rc -ne 75 ]] && break
        if (( SECONDS - start >= QUIET_WAIT_S )); then
            echo "[fail] $tag $case $mode t=$t: machine not quiet for ${QUIET_WAIT_S} s (see $dir/run.log)"
            return 0
        fi
        echo "[wait] $tag $case $mode t=$t: $(sed -n 2p "$dir/run.log" | cut -c1-120)"
        sleep "$QUIET_RETRY_S"
    done
    if [[ $rc -eq 0 ]]; then
        touch "$dir/done"
    else
        echo "[fail] $tag $case $mode t=$t exit=$rc (see $dir/run.log)"
    fi
}

ladder() {  # tree tag case
    point "$1" "$2" "$3" serial 1
    for t in $THREADS; do
        for mode in $ROUTES; do point "$1" "$2" "$3" "$mode" "$t"; done
    done
}

case "$STAGE" in
    A)
        ladder "$PRE_TREE" pre stack32_e6_actuated_saved
        ladder "$POST_TREE" post stack32_e6_actuated_saved
        ;;
    B)
        for n in ${SWEEP_N:-32 64 128 256}; do
            point "$POST_TREE" post "stack${n}_e6_actuated_saved" serial 1
            for t in ${SWEEP_THREADS:-12}; do
                for mode in $ROUTES; do point "$POST_TREE" post "stack${n}_e6_actuated_saved" "$mode" "$t"; done
            done
        done
        ;;
    *) echo "unknown stage $STAGE" >&2; exit 2 ;;
esac
echo "[done] $STAGE $(date -Is)"
