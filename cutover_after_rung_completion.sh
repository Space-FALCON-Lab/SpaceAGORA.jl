#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$SCRIPT_DIR}"
SESSION_NAME="${SESSION_NAME:-agora-nightly}"
TARGET_PASS="${TARGET_PASS:-1}"
TARGET_RUNG="${TARGET_RUNG:-la_density_multibody}"
POLL_SECONDS="${POLL_SECONDS:-20}"
OUT_DIR="${OUT_DIR:-$REPO_DIR/output/performance}"
AUTO_RESTART="${AUTO_RESTART:-1}"
NEW_RUN_CMD="${NEW_RUN_CMD:-./run_parallel_paper_protocol.sh}"
RUN_LOG="${RUN_LOG:-}"
SHUTDOWN_TIMEOUT_SECONDS="${SHUTDOWN_TIMEOUT_SECONDS:-600}"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S %Z'
}

log() {
  printf '[%s] %s\n' "$(timestamp)" "$*"
}

pane_command() {
  tmux display-message -p -t "$SESSION_NAME" "#{pane_current_command}" 2>/dev/null || true
}

pane_is_shell_ready() {
  local cmd
  cmd="$(pane_command)"
  case "$cmd" in
    bash|zsh|sh|fish)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

if [[ -z "$RUN_LOG" ]]; then
  RUN_LOG="$(ls -1t "$OUT_DIR"/overnight_*.log 2>/dev/null | head -n1 || true)"
fi

if [[ -z "$RUN_LOG" || ! -f "$RUN_LOG" ]]; then
  log "No overnight run log found in $OUT_DIR (or RUN_LOG is invalid)."
  exit 1
fi

if ! tmux has-session -t "$SESSION_NAME" >/dev/null 2>&1; then
  log "tmux session '$SESSION_NAME' not found."
  exit 1
fi

pattern="[smart-ladder] pass=${TARGET_PASS} rung=${TARGET_RUNG} completed"

mapfile -t tracked_pids < <(pgrep -f 'run_parallel_paper_protocol.sh|performance_smart_parallel_ladder.jl|performance_runtime_analysis.jl' || true)
log "Tracking existing run PIDs: ${tracked_pids[*]:-none}"
log "Monitoring $RUN_LOG for boundary: $pattern"

while true; do
  if grep -Fq "$pattern" "$RUN_LOG"; then
    log "Boundary reached."
    break
  fi
  last_line="$(tail -n 1 "$RUN_LOG" 2>/dev/null || true)"
  if [[ -n "$last_line" ]]; then
    log "Waiting... last line: $last_line"
  else
    log "Waiting... log has no lines yet."
  fi
  sleep "$POLL_SECONDS"
done

log "Sending Ctrl-C to tmux session '$SESSION_NAME'."
tmux send-keys -t "$SESSION_NAME" C-c

# Wait for tracked run processes to exit and pane to return to shell.
shutdown_loops=$((SHUTDOWN_TIMEOUT_SECONDS / 2))
if [[ "$shutdown_loops" -lt 1 ]]; then
  shutdown_loops=1
fi
for _ in $(seq 1 "$shutdown_loops"); do
  any_alive=0
  for pid in "${tracked_pids[@]}"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      any_alive=1
      break
    fi
  done
  if [[ "$any_alive" -eq 0 ]] && pane_is_shell_ready; then
    break
  fi
  sleep 2
done

# One extra interrupt if tracked processes are still alive.
need_second_interrupt=0
for pid in "${tracked_pids[@]}"; do
  if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
    need_second_interrupt=1
    break
  fi
done
if [[ "$need_second_interrupt" -eq 1 ]] || ! pane_is_shell_ready; then
  log "Run still active after first grace period; sending second Ctrl-C."
  tmux send-keys -t "$SESSION_NAME" C-c
fi

for _ in $(seq 1 "$shutdown_loops"); do
  any_alive=0
  for pid in "${tracked_pids[@]}"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      any_alive=1
      break
    fi
  done
  if [[ "$any_alive" -eq 0 ]] && pane_is_shell_ready; then
    break
  fi
  sleep 2
done

if ! pane_is_shell_ready; then
  log "tmux pane is not back at a shell command prompt; refusing auto-restart."
  exit 1
fi

if [[ "$AUTO_RESTART" != "1" ]]; then
  log "AUTO_RESTART=$AUTO_RESTART, stopping after cutover interrupt."
  exit 0
fi

new_log="$OUT_DIR/overnight_$(date +%Y%m%d_%H%M%S)_cutover.log"
launch_cmd="cd \"$REPO_DIR\" && $NEW_RUN_CMD 2>&1 | tee -a \"$new_log\""
log "Starting replacement run in tmux session '$SESSION_NAME'."
log "Launch command: $launch_cmd"
tmux send-keys -t "$SESSION_NAME" "$launch_cmd" C-m
log "Cutover complete. New run log: $new_log"
