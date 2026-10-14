#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

PROFILE="${SPACEAGORA_PPC_PROFILE:-full}"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUTDIR="${SPACEAGORA_PPC_OUTDIR:-output/performance/parallelization_performance_${STAMP}}"
REPEATS="${SPACEAGORA_PPC_REPEATS:-3}"
WARMUP="${SPACEAGORA_PPC_WARMUP:-1}"
SEED="${SPACEAGORA_PPC_SEED:-20260615}"
THREAD_ARGS=()
if [[ -n "${SPACEAGORA_PPC_THREADS:-}" ]]; then
  THREAD_ARGS+=(--threads="$SPACEAGORA_PPC_THREADS")
fi

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export GKSwstype="${GKSwstype:-100}"

julia --project=. benchmarks/studies/parallelization_performance.jl "$PROFILE" \
  --outdir="$OUTDIR" \
  "${THREAD_ARGS[@]}" \
  --repeats="$REPEATS" \
  --warmup="$WARMUP" \
  --seed="$SEED" \
  --solver-mode="${SPACEAGORA_PPC_SOLVER_MODE:-auto_stiff}"

echo "[done] Outputs: $OUTDIR"
