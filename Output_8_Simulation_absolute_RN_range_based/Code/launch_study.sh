#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
output="$(realpath -m -- "${2:-$code_dir/../Output}")"
campaign="${3:-range_based}"
case "$campaign" in
  range_based)
    batch="$output/_batch_range_based"
    session="oracle-study-range-based"
    ;;
  *)
    echo "Unknown campaign '$campaign'; use range_based." >&2
    exit 2
    ;;
esac

case "${1:-status}" in
  start)
    command -v tmux >/dev/null || { echo "tmux is required for an SSH-independent launch." >&2; exit 1; }
    tmux new-session -d -s "$session" -- bash "$code_dir/launch_study.sh" run "$output" "$campaign"
    echo "Started tmux session $session. Logs: $batch/grid.log"
    ;;
  run)
    mkdir -p -- "$batch"
    julia_bin="${JULIA_BIN:-$(command -v julia)}"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    exec "$julia_bin" --startup-file=no --threads=1 --gcthreads=1 "$code_dir/run_grid.jl" "$output" "$campaign" >>"$batch/grid.log" 2>&1
    ;;
  status)
    if [[ -f "$batch/status.toml" ]]; then
      cat -- "$batch/status.toml"
    else
      echo "No study batch has started in $output."
    fi
    ;;
  stop-after-case)
    mkdir -p -- "$batch"
    touch -- "$batch/STOP_AFTER_CASE"
    echo "The batch will stop after its current case and verification."
    ;;
  resume)
    rm -f -- "$batch/STOP_AFTER_CASE"
    exec bash "$code_dir/launch_study.sh" start "$output" "$campaign"
    ;;
  *)
    echo "Usage: bash launch_study.sh {start|run|status|stop-after-case|resume} [OUTPUT_ROOT] [range_based]" >&2
    exit 2
    ;;
esac
