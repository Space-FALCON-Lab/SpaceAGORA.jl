#!/usr/bin/env bash
# Run the full paper benchmark suite.
#
# Environment variable overrides (all optional):
#   SPACEAGORA_PPB_OUTDIR           — output root directory
#                                     default: <repo>/output/performance/paper_benchmarks
#   SPACEAGORA_PPB_PHASES           — comma-separated phase subset, e.g. B1,B2
#                                     default: all phases (B1–B8)
#   SPACEAGORA_PPB_THREADS          — comma-separated thread ladder override
#                                     default: auto-scaled to the machine's
#                                     PHYSICAL core count (not Sys.CPU_THREADS —
#                                     the RHS is SIMD/FP-bound and regresses hard
#                                     once SMT siblings are oversubscribed, so a
#                                     ladder topped out at logical CPUs makes its
#                                     own last rung the worst point on the curve)
#   SPACEAGORA_PPC_PHYSICAL_CORES   — override physical-core detection, e.g. when
#                                     the run is confined to part of the machine
#                                     via taskset/cgroups
#                                     default: unset (detected from /proc/cpuinfo)
#   SPACEAGORA_PPB_PROCESS_WORKERS  — max process workers for MC phases (B4)
#                                     default: 32
#   SPACEAGORA_PPB_MC_SAMPLES_MAX   — caps each phase's Monte Carlo sample
#                                     ladder (e.g. B4's default
#                                     1,4,16,64,256,1024) at this value; at
#                                     least the smallest sample count always
#                                     runs. Unlike --preview, this only
#                                     affects MC samples — N_sat, workers,
#                                     and repeats are untouched.
#                                     default: unset (full ladder)
#   SPACEAGORA_PPB_SOLVER_MODE      — ODE solver mode
#                                     default: auto_stiff
#   SPACEAGORA_PPB_SEED             — RNG seed
#                                     default: 20260615
#   SPACEAGORA_PPB_CPU_LIST         — taskset-style CPU pool to pin worker
#                                     subprocesses to, e.g. "0-15" or "0-7,16-23".
#                                     Each worker is pinned to the first N cores
#                                     of this pool, where N is that worker's
#                                     thread count. Linux only; unset disables
#                                     pinning (default).
#   SPACEAGORA_PPB_DRY_RUN          — set to 1 to print planned runs without executing
#                                     default: 0
#   SPACEAGORA_PPB_PREVIEW          — set to 1 for a local-PC preview run:
#                                     caps N_sat at 64, MC samples at 16, process
#                                     workers at 4, and repeats at 2.
#                                     Thread ladder still auto-scales to local CPU count.
#                                     default: 0
#   SPACEAGORA_PPB_RESUME           — path to an existing run directory (e.g.
#                                     <outdir>/20260819_195133) to resume into,
#                                     instead of starting a fresh timestamped
#                                     run. Any (case, mode, thread_count,
#                                     mc_samples, repeat) worker CSV already on
#                                     disk under <resume>/<phase>/worker_rows/
#                                     with a successful row is skipped rather
#                                     than re-run — use this to pick a crashed
#                                     or killed run back up without redoing
#                                     everything that already finished.
#                                     default: unset (fresh run)
#
# Additional arguments are forwarded to the Julia script.
#
# Example — full run with 16 process workers:
#   SPACEAGORA_PPB_PROCESS_WORKERS=16 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh
#
# Example — local PC preview run:
#   SPACEAGORA_PPB_PREVIEW=1 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh
#
# Example — dry run to inspect what preview would execute:
#   SPACEAGORA_PPB_PREVIEW=1 SPACEAGORA_PPB_DRY_RUN=1 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh
#
# Example — run only B1 and B5:
#   SPACEAGORA_PPB_PHASES=B1,B5 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh
#
# Example — run a single phase (B1, B2, B3, B4, or B5) by passing its id as a
# bare positional argument; this is forwarded straight to the Julia script
# and is equivalent to SPACEAGORA_PPB_PHASES=B4:
#   bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh B4
#
# Note: an older B6 ("Cross-Machine Validation") was removed from the phase
# catalog; the current B6 ("Expanded Router Evaluation") is a later, unrelated
# addition and is a valid phase id.
#
# Example — cap B4's MC samples at 64 (runs 1,4,16,64 instead of up to 1024),
# leaving everything else (N_sat, workers, repeats) at full scale:
#   SPACEAGORA_PPB_MC_SAMPLES_MAX=64 bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh B4
#
# Example — resume a crashed run instead of starting over:
#   SPACEAGORA_PPB_PHASES=B6 SPACEAGORA_PPB_RESUME=output/performance/paper_benchmarks/20260819_195133 \
#     bash benchmarks/studies/paper_parallelization_benchmarks/protocol.sh

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

if [ -n "${SPACEAGORA_PPB_MC_SAMPLES_MAX:-}" ]; then
    ARGS+=("--mc-samples-max=${SPACEAGORA_PPB_MC_SAMPLES_MAX}")
fi

if [ -n "${SPACEAGORA_PPB_SOLVER_MODE:-}" ]; then
    ARGS+=("--solver-mode=${SPACEAGORA_PPB_SOLVER_MODE}")
fi

if [ -n "${SPACEAGORA_PPB_SEED:-}" ]; then
    ARGS+=("--seed=${SPACEAGORA_PPB_SEED}")
fi

if [ -n "${SPACEAGORA_PPB_CPU_LIST:-}" ]; then
    ARGS+=("--cpu-list=${SPACEAGORA_PPB_CPU_LIST}")
fi

if [ "${SPACEAGORA_PPB_DRY_RUN:-0}" = "1" ]; then
    ARGS+=("--dry-run")
fi

if [ "${SPACEAGORA_PPB_PREVIEW:-0}" = "1" ]; then
    ARGS+=("--preview")
fi

julia --project=. "${SCRIPT}" "${ARGS[@]}" "$@"
