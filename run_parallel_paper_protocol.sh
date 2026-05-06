#!/usr/bin/env bash
set -euo pipefail

# Always run from the repo root (script location).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ ! -f "Project.toml" || ! -f "benchmarks/studies/performance_smart_parallel_ladder.jl" ]]; then
  echo "[error] run_parallel_paper_protocol.sh must live in the SpaceAGORA repo root."
  exit 1
fi

PROFILE="full"
SEED="20260304"
PASSES="3"
LAYER_ATTR="1"

# Auto-detect logical CPU count; override by setting SPACEAGORA_PERF_NCPU.
NCPU="${SPACEAGORA_PERF_NCPU:-$(nproc)}"
THREADS_MAIN="$NCPU"

# Four evenly-spaced process worker counts from NCPU/4 up to NCPU.
_gen_process_workers() {
  local n=$1 step=$(( $1 / 4 )) arr=()
  [[ $step -lt 1 ]] && step=1
  for i in 1 2 3 4; do
    local w=$(( step * i ))
    [[ $w -gt $n ]] && w=$n
    [[ ${#arr[@]} -eq 0 || ${arr[-1]} -ne $w ]] && arr+=("$w")
  done
  echo "${arr[@]}"
}

# Powers of 2 from 1 up to NCPU, then NCPU itself if not already a power of 2.
_gen_scaling_threads() {
  local n=$1 t=1 arr=()
  while [[ $t -lt $n ]]; do arr+=("$t"); t=$(( t * 2 )); done
  arr+=("$n")
  echo "${arr[@]}"
}

read -ra PROCESS_WORKERS_LIST <<< "$(_gen_process_workers "$NCPU")"
read -ra SCALING_THREADS_LIST <<< "$(_gen_scaling_threads "$NCPU")"

BASE_OUT="output/performance/paper_protocol_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BASE_OUT"

echo "[setup] Instantiating project"
julia --project=. -e '
    using Pkg
    # GRAMSuite is declared as a weakdep + sources path; develop it so
    # "using GRAMSuite" works in the benchmark and the SpaceAGORA extension fires.
    gram_path = joinpath(pwd(), "data", "GRAMSuite.jl")
    if isdir(gram_path) && !haskey(Pkg.project().dependencies, "GRAMSuite")
        @info "Developing GRAMSuite from local source for GRAM benchmark cases"
        Pkg.develop(Pkg.PackageSpec(path=gram_path))
    end
    Pkg.instantiate()
    Pkg.precompile()
'

# Thread/process hygiene.
export OPENBLAS_NUM_THREADS=1
export SPACEAGORA_PERF_MACHINE_LABEL="${SPACEAGORA_PERF_MACHINE_LABEL:-$(hostname -s)}"
export JULIA_NUM_THREADS="$THREADS_MAIN"
# export JULIA_EXCLUSIVE=1   # enable only on a dedicated machine

echo "[config] NCPU=$NCPU  THREADS_MAIN=$THREADS_MAIN  MACHINE=$SPACEAGORA_PERF_MACHINE_LABEL"
echo "[config] PROCESS_WORKERS_LIST=(${PROCESS_WORKERS_LIST[*]})"
echo "[config] SCALING_THREADS_LIST=(${SCALING_THREADS_LIST[*]})"

# Freeze solver axis for parallel-only claims.
export SPACEAGORA_PERF_SOLVER_MODE="auto_stiff"
export SPACEAGORA_PERF_SPLIT_ROLLOUT_GATE=0
export SPACEAGORA_PERF_MULTIRATE_ROLLOUT_GATE=0

run_ladder_once() {
  local outdir="$1"
  local clean="$2"
  local reset="$3"
  local backend="$4"
  local threads="$5"
  local workers="${6:-}"

  local cmd=(
    julia --project=. benchmarks/studies/performance_smart_parallel_ladder.jl
    "$PROFILE"
    --outdir="$outdir"
    --clean="$clean"
    --passes="$PASSES"
    --randomize-rung-order=1
    --seed="$SEED"
    --outer-only-backend="$backend"
    --layer-attribution="$LAYER_ATTR"
    --solver-axis=frozen
    --solver-mode=auto_stiff
  )

  if [[ -n "$workers" ]]; then
    cmd+=(--process-workers="$workers")
  fi

  SPACEAGORA_PERF_OUTER_ROUTE_STATE_RESET="$reset" \
  SPACEAGORA_PARALLEL_POLICY_STATE_RESET="$reset" \
  JULIA_NUM_THREADS="$threads" \
  "${cmd[@]}"
}

run_cold_warm_pair() {
  local outdir="$1"
  local backend="$2"
  local threads="$3"
  local workers="${4:-}"

  echo "=== COLD: $outdir ==="
  run_ladder_once "$outdir" 1 1 "$backend" "$threads" "$workers"

  echo "=== WARM: $outdir ==="
  run_ladder_once "$outdir" 0 0 "$backend" "$threads" "$workers"
}

echo "[run] 1) Main full run (threads outer rung), cold+warm"
run_cold_warm_pair "$BASE_OUT/ladder_full_threads_main" "threads" "$THREADS_MAIN"

echo "[run] 2) Process outer-rung comparison, cold+warm"
for W in "${PROCESS_WORKERS_LIST[@]}"; do
  # Keep worker subprocesses single-threaded at root to avoid oversubscription.
  run_cold_warm_pair "$BASE_OUT/ladder_full_process_w${W}" "process" "1" "$W"
done

echo "[run] 3) Strong-scaling sweep (threads), cold+warm"
for T in "${SCALING_THREADS_LIST[@]}"; do
  run_cold_warm_pair "$BASE_OUT/ladder_scaling_threads_${T}" "threads" "$T"
done

echo "[run] 4) Paper pipeline summary across modes"
JULIA_NUM_THREADS="$THREADS_MAIN" \
julia --project=. benchmarks/scripts/performance_paper_pipeline.jl \
  --profile=full \
  --modes=serial,auto,threads,process \
  --outdir="$BASE_OUT/paper_pipeline_full"

echo "[done] Outputs: $BASE_OUT"
