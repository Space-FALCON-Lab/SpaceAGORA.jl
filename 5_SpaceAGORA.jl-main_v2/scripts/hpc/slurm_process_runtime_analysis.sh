#!/usr/bin/env bash
#SBATCH --job-name=spaceagora-runtime
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$REPO_ROOT"

# Use a writable job-local depot overlay instead of mutating a shared depot.
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$PWD/.julia_depot_slurm_${SLURM_JOB_ID:-manual}}"
mkdir -p "$JULIA_DEPOT_PATH"

export SPACEAGORA_PERF_WORKER_PROJECT="${SPACEAGORA_PERF_WORKER_PROJECT:-.}"
export SPACEAGORA_PERF_PROCS="${SPACEAGORA_PERF_PROCS:-${SLURM_CPUS_PER_TASK:-4}}"
export SPACEAGORA_PERF_MACHINE_LABEL="${SPACEAGORA_PERF_MACHINE_LABEL:-slurm-${SLURM_JOB_PARTITION:-default}}"
export SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS="${SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS:-1}"
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH="${SPACEAGORA_PARALLEL_POLICY_STATE_PATH:-output/parallel_policy_state/slurm_hint_state.toml}"

mkdir -p output/parallel_policy_state

./bin/spaceagora benchmark runtime-analysis smoke --output-dir="output/hpc_slurm_runtime_analysis_${SLURM_JOB_ID:-manual}"
