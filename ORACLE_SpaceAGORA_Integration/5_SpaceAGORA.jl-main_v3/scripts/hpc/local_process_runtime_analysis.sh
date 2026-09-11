#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

export SPACEAGORA_PERF_WORKER_PROJECT="${SPACEAGORA_PERF_WORKER_PROJECT:-.}"
export SPACEAGORA_PERF_PROCS="${SPACEAGORA_PERF_PROCS:-4}"
export SPACEAGORA_PERF_MACHINE_LABEL="${SPACEAGORA_PERF_MACHINE_LABEL:-local-process}"
export SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS="${SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS:-1}"
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH="${SPACEAGORA_PARALLEL_POLICY_STATE_PATH:-output/parallel_policy_state/local_process_hint_state.toml}"

mkdir -p output/parallel_policy_state

./bin/spaceagora benchmark runtime-analysis smoke --output-dir=output/hpc_local_runtime_analysis
