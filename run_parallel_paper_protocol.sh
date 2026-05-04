#!/usr/bin/env bash
set -euo pipefail

# Always run from the repo root (script location).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ ! -f "Project.toml" || ! -f "benchmarks/studies/performance_smart_parallel_ladder.jl" ]]; then
  echo "[error] run_parallel_paper_protocol.sh must live in the SpaceAGORA repo root."
  exit 1
fi

PROFILE="${SPACEAGORA_PAPER_PROFILE:-full}"
SEED="20260304"
PASSES_MAIN="${SPACEAGORA_PAPER_PASSES_MAIN:-2}"
PASSES_LIGHT="${SPACEAGORA_PAPER_PASSES_LIGHT:-1}"
LAYER_ATTR_MAIN="1"
LAYER_ATTR_LIGHT="0"

THREADS_MAIN="32"
PROCESS_WORKERS_LIST=(32)
SCALING_THREADS_LIST=(1 2 4 8 16 32 64)

PROCESS_RUNGS="r1_b_outer_only_process"
SCALING_RUNGS="r5_full_smart"

BASE_OUT="output/performance/paper_protocol_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BASE_OUT"

echo "[setup] Instantiating root project"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Thread/process hygiene.
export OPENBLAS_NUM_THREADS=1
export SPACEAGORA_PERF_MACHINE_LABEL="${SPACEAGORA_PERF_MACHINE_LABEL:-tower_9985WX}"
export JULIA_NUM_THREADS="$THREADS_MAIN"
export SPACEAGORA_SMART_LADDER_PROJECT="$SCRIPT_DIR"
export SPACEAGORA_PERF_WORKER_PROJECT="$SCRIPT_DIR"
# export JULIA_EXCLUSIVE=1   # enable only on a dedicated machine

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
  local passes="${7:-$PASSES_LIGHT}"
  local layer_attr="${8:-$LAYER_ATTR_LIGHT}"
  local rungs="${9:-}"

  local cmd=(
    julia --project=. benchmarks/studies/performance_smart_parallel_ladder.jl
    "$PROFILE"
    --outdir="$outdir"
    --clean="$clean"
    --passes="$passes"
    --randomize-rung-order=1
    --seed="$SEED"
    --outer-only-backend="$backend"
    --layer-attribution="$layer_attr"
    --solver-axis=frozen
    --solver-mode=auto_stiff
  )

  if [[ -n "$rungs" ]]; then
    cmd+=(--rungs="$rungs")
  fi

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
  local passes="${5:-$PASSES_LIGHT}"
  local layer_attr="${6:-$LAYER_ATTR_LIGHT}"
  local rungs="${7:-}"

  echo "=== COLD: $outdir ==="
  run_ladder_once "$outdir" 1 1 "$backend" "$threads" "$workers" "$passes" "$layer_attr" "$rungs"

  echo "=== WARM: $outdir ==="
  run_ladder_once "$outdir" 0 0 "$backend" "$threads" "$workers" "$passes" "$layer_attr" "$rungs"
}

run_cold_only() {
  local outdir="$1"
  local backend="$2"
  local threads="$3"
  local workers="${4:-}"
  local passes="${5:-$PASSES_LIGHT}"
  local layer_attr="${6:-$LAYER_ATTR_LIGHT}"
  local rungs="${7:-}"

  echo "=== COLD ONLY: $outdir ==="
  run_ladder_once "$outdir" 1 1 "$backend" "$threads" "$workers" "$passes" "$layer_attr" "$rungs"
}

echo "[run] 1) Main full run (threads outer rung), cold+warm"
run_cold_warm_pair "$BASE_OUT/ladder_full_threads_main" "threads" "$THREADS_MAIN" "" "$PASSES_MAIN" "$LAYER_ATTR_MAIN"

echo "[run] 2) Process outer-rung comparison, cold+warm (process rung only)"
for W in "${PROCESS_WORKERS_LIST[@]}"; do
  # Keep worker subprocesses single-threaded at root to avoid oversubscription.
  run_cold_warm_pair "$BASE_OUT/ladder_process_w${W}" "process" "1" "$W" "$PASSES_LIGHT" "$LAYER_ATTR_LIGHT" "$PROCESS_RUNGS"
done

echo "[run] 3) Strong-scaling sweep (threads), cold-only, r5_full_smart only"
for T in "${SCALING_THREADS_LIST[@]}"; do
  run_cold_only "$BASE_OUT/ladder_scaling_threads_${T}" "threads" "$T" "" "$PASSES_LIGHT" "$LAYER_ATTR_LIGHT" "$SCALING_RUNGS"
done

echo "[run] 4) Paper pipeline summary across modes"
JULIA_NUM_THREADS="$THREADS_MAIN" \
julia --project=. benchmarks/scripts/performance_paper_pipeline.jl \
  --profile=full \
  --modes=serial,auto,threads,process \
  --outdir="$BASE_OUT/paper_pipeline_full"

echo "[done] Outputs: $BASE_OUT"
