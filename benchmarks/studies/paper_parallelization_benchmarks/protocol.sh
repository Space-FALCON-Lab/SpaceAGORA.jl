#!/usr/bin/env bash
# Run the full paper benchmark suite.
#
# Environment variable overrides (all optional):
#   SPACEAGORA_PPB_OUTDIR           — output root directory
#                                     default: <repo>/output/performance/paper_benchmarks
#   SPACEAGORA_PPB_PHASES           — comma-separated phase subset, e.g. B1,B2
#                                     default: all phases (B1–B6)
#   SPACEAGORA_PPB_THREADS          — comma-separated thread ladder override
#                                     default: auto-scaled to Sys.CPU_THREADS
#   SPACEAGORA_PPB_PROCESS_WORKERS  — max process workers for MC phases (B4, B6)
#                                     default: 32
#   SPACEAGORA_PPB_SOLVER_MODE      — ODE solver mode
#                                     default: auto_stiff
#   SPACEAGORA_PPB_SEED             — RNG seed
#                                     default: 20260615
#   SPACEAGORA_PPB_DRY_RUN          — set to 1 to print planned runs without executing
#                                     default: 0
#
# Additional arguments are forwarded to the Julia script.
#
# Example — run all phases with 16 process workers:
#   SPACEAGORA_PPB_PROCESS_WORKERS=16 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh
#
# Example — run only B1 and B5, dry run:
#   SPACEAGORA_PPB_PHASES=B1,B5 SPACEAGORA_PPB_DRY_RUN=1 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SCRIPT="${REPO_ROOT}/benchmarks/studies/paper_parallelization_benchmarks.jl"

OUTDIR="${SPACEAGORA_PPB_OUTDIR:-${REPO_ROOT}/output/performance/paper_benchmarks}"
WORKERS="${SPACEAGORA_PPB_PROCESS_WORKERS:-32}"

export OPENBLAS_NUM_THREADS=1
export GKSwstype=100

cd "${REPO_ROOT}"

ARGS=("--outdir=${OUTDIR}" "--process-workers=${WORKERS}")

if [ -n "${SPACEAGORA_PPB_THREADS:-}" ]; then
    ARGS+=("--threads=${SPACEAGORA_PPB_THREADS}")
fi

if [ -n "${SPACEAGORA_PPB_PHASES:-}" ]; then
    ARGS+=("--phases=${SPACEAGORA_PPB_PHASES}")
fi

if [ -n "${SPACEAGORA_PPB_SOLVER_MODE:-}" ]; then
    ARGS+=("--solver-mode=${SPACEAGORA_PPB_SOLVER_MODE}")
fi

if [ -n "${SPACEAGORA_PPB_SEED:-}" ]; then
    ARGS+=("--seed=${SPACEAGORA_PPB_SEED}")
fi

if [ "${SPACEAGORA_PPB_DRY_RUN:-0}" = "1" ]; then
    ARGS+=("--dry-run")
fi

julia --project=. "${SCRIPT}" "${ARGS[@]}" "$@"
