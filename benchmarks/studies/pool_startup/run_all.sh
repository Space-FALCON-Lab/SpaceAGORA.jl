#!/usr/bin/env bash
# Process-pool startup, before vs after, every build in one session.
#
#   TRIALS=3 bash run_all.sh <outdir> <label>=<tree> [<label>=<tree> ...]
#
# Every trial is a fresh coordinator (measure.jl). Within each point the build
# order is rotated per trial, so no build always runs first or last. Labels
# starting with "after" also get the PROCESS_WARMUP variants (skip, cheap) in
# the campaign step. Finishes with the pool/process-campaign tests per build.
set -u
OUT=$(realpath -m "$1"); shift
NAMES=(); TREES=()
for a in "$@"; do NAMES+=("${a%%=*}"); TREES+=("$(realpath "${a#*=}")"); done
NB=${#NAMES[@]}
TRIALS=${TRIALS:-3}
SAMPLES=${SAMPLES:-256}
CAMPAIGN_WORKERS=${CAMPAIGN_WORKERS:-32}
M="$(cd "$(dirname "$0")" && pwd)/measure.jl"
mkdir -p "$OUT/state" "$OUT/logs"

# Scratch adaptive-policy state: nothing read from or written to a shared store.
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH="$OUT/state/inner_policy_state.toml"
export SPACEAGORA_COST_CONSTANTS_PATH="$OUT/state/cost_constants.toml"
export SPACEAGORA_OUTER_ROUTE_STATE_PATH="$OUT/state/outer_route_state.toml"
export SPACEAGORA_RHS_CALIBRATION_PATH="$OUT/state/rhs_calibration.toml"
export SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH="$OUT/state/campaign_corrections.toml"
export SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1
unset JULIA_NUM_THREADS SPACEAGORA_PERF_PROCS

log() { echo "[$(date -Iseconds)] $*" | tee -a "$OUT/run.log"; }
# measure <build-index> <logname> <args...>
measure() {
  local i=$1 name=$2; shift 2
  log "start ${NAMES[$i]} $*"
  julia --project="${TREES[$i]}" --threads=1 --startup-file=no "$M" "${NAMES[$i]}" "$@" \
    > "$OUT/logs/$name.log" 2>&1
  log "end   ${NAMES[$i]} exit=$? $(grep '^RESULT' "$OUT/logs/$name.log" | tail -1)"
}
# order <offset>: build indices rotated by offset
order() { local k=$(( $1 % NB )) j; for ((j = 0; j < NB; j++)); do echo $(( (j + k) % NB )); done; }

{
  echo "host=$(hostname -s) nproc=$(nproc) date=$(date -Iseconds)"
  julia --version
  for ((i = 0; i < NB; i++)); do echo "build ${NAMES[$i]} tree=${TREES[$i]}"; done
  free -g | head -2
} > "$OUT/environment.txt"

# Prime: precompile each tree (SpaceAGORA, GRAMSuite and its extension, the
# study files' dependencies) and run one discarded 2-worker bootstrap, so no
# timed trial pays a cache build. Each release gets a private policy-state dir.
for ((i = 0; i < NB; i++)); do
  t=${TREES[$i]}
  [ -L "$t/output/parallel_policy_state" ] && rm "$t/output/parallel_policy_state"
  mkdir -p "$t/output/parallel_policy_state"
  log "precompile ${NAMES[$i]}"
  julia --project="$t" --startup-file=no -e '
    pushfirst!(LOAD_PATH, joinpath(dirname(Base.active_project()), "data", "GRAMSuite.jl"))
    using SpaceAGORA; import GRAMSuite; using CSV, DataFrames
    println("ext=", Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing)' \
    > "$OUT/logs/precompile_${NAMES[$i]}.log" 2>&1
  log "precompile ${NAMES[$i]} exit=$? $(tail -1 "$OUT/logs/precompile_${NAMES[$i]}.log")"
  measure "$i" "prime_${NAMES[$i]}" bootstrap 2 "$OUT/prime"
done

p=0
for n in 8 16 32; do
  for ((t = 0; t < TRIALS; t++)); do
    for i in $(order $((t + p))); do measure "$i" "bootstrap_${NAMES[$i]}_n${n}_t$((t + 1))" bootstrap "$n" "$OUT"; done
  done
  p=$((p + 1))
done

for ((t = 0; t < TRIALS; t++)); do
  for i in $(order "$t"); do measure "$i" "warmup_${NAMES[$i]}_n32_t$((t + 1))" warmup 32 "$OUT"; done
done

CONFIGS=()
for ((i = 0; i < NB; i++)); do
  CONFIGS+=("$i:default")
  [[ ${NAMES[$i]} == after* ]] && CONFIGS+=("$i:skip" "$i:cheap")
done
NC=${#CONFIGS[@]}
for ((t = 0; t < TRIALS; t++)); do
  for ((j = 0; j < NC; j++)); do
    c=${CONFIGS[$(( (j + t * 3) % NC ))]}
    i=${c%%:*}; v=${c#*:}
    measure "$i" "campaign_${NAMES[$i]}_${v}_t$((t + 1))" campaign "$CAMPAIGN_WORKERS" "$OUT" "$v" "$SAMPLES"
  done
done

# Tests that exercise the worker pool or the process campaign route.
TESTS="test/probes/process_pool_probes.jl test/probes/campaign_process_route_probes.jl
test/unit/parallel/mixed_dispatch_tests.jl test/unit/parallel/process_worker_load_path_tests.jl
test/unit/parallel/pool_worker_heap_hint_tests.jl"
echo "build,test,exit_code,seconds" > "$OUT/tests.csv"
for ((i = 0; i < NB; i++)); do
  for f in $TESTS; do
    s=$(date +%s)
    (cd "${TREES[$i]}" && julia --project=. --startup-file=no "$f") > "$OUT/logs/test_${NAMES[$i]}_$(basename "$f" .jl).log" 2>&1
    e=$?
    echo "${NAMES[$i]},$f,$e,$(( $(date +%s) - s ))" >> "$OUT/tests.csv"
    log "test ${NAMES[$i]} $f exit=$e"
  done
done
log "done"
