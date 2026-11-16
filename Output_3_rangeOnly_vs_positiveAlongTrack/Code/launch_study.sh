#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project="$(cd -- "$code_dir/../.." && pwd)"
output="$project/Output_3_rangeOnly_vs_positiveAlongTrack/Output"
baseline="$project/Output_4_no_laser/Output"
job="$output/_comparison_job"
session="oracle-range-along-no-laser"
mkdir -p -- "$job" "$output/_batch_range_along" "$baseline/_batch_no_laser"

case "${1:-status}" in
  run)
    command -v flock >/dev/null || { echo "flock is required for the job lock." >&2; exit 1; }
    julia_bin="${JULIA_BIN:-$(command -v julia)}"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 GKSwstype=100
    export JULIA_NUM_PRECOMPILE_TASKS=1
    exec < /dev/null
    exec flock --nonblock --no-fork "$job/process.lock" \
      "$julia_bin" --startup-file=no --threads=1 --gcthreads=1 \
      "$code_dir/run_jobs.jl" >>"$job/job.log" 2>&1
    ;;
  start)
    command -v tmux >/dev/null || { echo "tmux is required for a manual SSH-independent launch." >&2; exit 1; }
    tmux new-session -d -s "$session" -- bash "$code_dir/launch_study.sh" run
    echo "Started tmux session $session. Check status to confirm the job is running."
    ;;
  status)
    for path in "$job/status.toml" "$output/_batch_range_along/status.toml" "$baseline/_batch_no_laser/status.toml"; do
      if [[ -f "$path" ]]; then
        printf '\n--- %s ---\n' "$path"
        cat -- "$path"
      fi
    done
    if [[ -f "$job/process.lock" ]]; then
      if flock --nonblock "$job/process.lock" true; then
        echo "Job lock is free: no controller is running."
      else
        echo "Job lock is held: the controller is running."
      fi
    fi
    ;;
  stop-after-case)
    touch -- "$job/STOP_AFTER_CASE" "$output/_batch_range_along/STOP_AFTER_CASE" \
      "$baseline/_batch_no_laser/STOP_AFTER_CASE"
    echo "The job will pause after the current simulation and verification."
    ;;
  resume)
    rm -f -- "$job/STOP_AFTER_CASE" "$output/_batch_range_along/STOP_AFTER_CASE" \
      "$baseline/_batch_no_laser/STOP_AFTER_CASE"
    exec bash "$code_dir/launch_study.sh" start
    ;;
  *)
    echo "Usage: bash launch_study.sh {start|run|status|stop-after-case|resume}" >&2
    exit 2
    ;;
esac
