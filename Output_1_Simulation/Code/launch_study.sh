#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
output="$(realpath -m -- "${2:-$code_dir/../Output}")"
campaign="${3:-main}"
case "$campaign" in
  main)
    batch="$output/_batch"
    session="oracle-study"
    ;;
  helpers400500)
    batch="$output/_batch_helpers400500"
    session="oracle-study-helpers400500"
    ;;
  helpers600800)
    batch="$output/_batch_helpers600800"
    session="oracle-study-helpers600800"
    ;;
  helpers700)
    batch="$output/_batch_helpers700"
    session="oracle-study-helpers700"
    ;;
  helpers9001100km)
    batch="$output/_batch_helpers9001100km"
    session="oracle-study-helpers9001100km"
    ;;
  targets9001100km)
    batch="$output/_batch_targets9001100km"
    session="oracle-study-targets9001100km"
    ;;
  *)
    echo "Unknown campaign '$campaign'; use main, helpers400500, helpers600800, helpers700, targets9001100km or helpers9001100km." >&2
    exit 2
    ;;
esac
mkdir -p -- "$batch"

case "${1:-status}" in
  start)
    command -v tmux >/dev/null || { echo "tmux is required for an SSH-independent launch." >&2; exit 1; }
    tmux new-session -d -s "$session" -- bash "$code_dir/launch_study.sh" run "$output" "$campaign"
    echo "Started tmux session $session. Logs: $batch/grid.log"
    ;;
  run)
    julia_bin="${JULIA_BIN:-$(command -v julia)}"
    exec "$julia_bin" --startup-file=no --threads=1 "$code_dir/run_grid.jl" "$output" "$campaign" >>"$batch/grid.log" 2>&1
    ;;
  status)
    if [[ -f "$batch/status.toml" ]]; then
      cat -- "$batch/status.toml"
    else
      echo "No study batch has started in $output."
    fi
    ;;
  stop-after-case)
    touch -- "$batch/STOP_AFTER_CASE"
    echo "The batch will stop after its current case and verification."
    ;;
  resume)
    rm -f -- "$batch/STOP_AFTER_CASE"
    exec bash "$code_dir/launch_study.sh" start "$output" "$campaign"
    ;;
  *)
    echo "Usage: bash launch_study.sh {start|run|status|stop-after-case|resume} [OUTPUT_ROOT] [main|helpers400500|helpers600800|helpers700|targets9001100km|helpers9001100km]" >&2
    exit 2
    ;;
esac
