#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${REPO_DIR:-/tmp/SpaceAGORA.jl}"
OUT_DIR="${OUT_DIR:-$REPO_DIR/output/performance}"
SESSION_NAME="${SESSION_NAME:-agora-nightly}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-10800}"   # 3 hours
DAY_START_HOUR="${DAY_START_HOUR:-8}"           # inclusive
DAY_END_HOUR="${DAY_END_HOUR:-22}"              # exclusive
TZ_NAME="${TZ_NAME:-America/Detroit}"
EMAIL_TO="${EMAIL_TO:-}"
EMAIL_FROM="${EMAIL_FROM:-spaceagora-monitor@$(hostname -f 2>/dev/null || hostname)}"
EMAIL_SUBJECT_PREFIX="${EMAIL_SUBJECT_PREFIX:-[SpaceAGORA Progress]}"
MAILER_BIN="${MAILER_BIN:-}"
SMTP_HOST="${SMTP_HOST:-}"
SMTP_PORT="${SMTP_PORT:-587}"
SMTP_USER="${SMTP_USER:-}"
SMTP_PASSWORD_FILE="${SMTP_PASSWORD_FILE:-}"
SMTP_PASSWORD_ENV="${SMTP_PASSWORD_ENV:-SPACEAGORA_SMTP_PASSWORD}"
SMTP_STARTTLS="${SMTP_STARTTLS:-1}"
RUN_ONCE="${RUN_ONCE:-0}"
MONITOR_LOG="${MONITOR_LOG:-$OUT_DIR/progress_email_daemon.log}"

mkdir -p "$OUT_DIR"

timestamp() {
  TZ="$TZ_NAME" date '+%Y-%m-%d %H:%M:%S %Z'
}

log() {
  printf '[%s] %s\n' "$(timestamp)" "$*" | tee -a "$MONITOR_LOG"
}

count_matches() {
  local regex="$1"
  local file="$2"
  local count
  count="$(grep -E -c "$regex" "$file" 2>/dev/null || true)"
  if [[ -z "$count" ]]; then
    count=0
  fi
  echo "$count"
}

in_daytime_window() {
  local hour
  hour="$(TZ="$TZ_NAME" date +%H)"
  if ((10#$hour >= DAY_START_HOUR && 10#$hour < DAY_END_HOUR)); then
    return 0
  fi
  return 1
}

detect_mailer() {
  if [[ -n "$MAILER_BIN" ]]; then
    if [[ -x "$MAILER_BIN" || "$(command -v "$MAILER_BIN" 2>/dev/null || true)" != "" ]]; then
      echo "$MAILER_BIN"
      return 0
    fi
  fi
  if command -v sendmail >/dev/null 2>&1; then
    echo "sendmail"
    return 0
  fi
  if command -v msmtp >/dev/null 2>&1; then
    echo "msmtp"
    return 0
  fi
  if command -v mailx >/dev/null 2>&1; then
    echo "mailx"
    return 0
  fi
  if command -v mail >/dev/null 2>&1; then
    echo "mail"
    return 0
  fi
  return 1
}

read_smtp_password() {
  local pass=""
  if [[ -n "$SMTP_PASSWORD_FILE" && -f "$SMTP_PASSWORD_FILE" ]]; then
    pass="$(<"$SMTP_PASSWORD_FILE")"
  elif [[ -n "${!SMTP_PASSWORD_ENV-}" ]]; then
    pass="${!SMTP_PASSWORD_ENV}"
  fi
  pass="${pass%$'\r'}"
  pass="${pass%$'\n'}"
  printf '%s' "$pass"
}

send_email_via_python_smtp() {
  local subject="$1"
  local body="$2"
  local smtp_password
  smtp_password="$(read_smtp_password)"

  if [[ -z "$SMTP_HOST" || -z "$SMTP_USER" || -z "$smtp_password" ]]; then
    return 1
  fi

  SMTP_SUBJECT="$subject" \
  SMTP_BODY="$body" \
  SMTP_TO="$EMAIL_TO" \
  SMTP_FROM="$EMAIL_FROM" \
  SMTP_HOST="$SMTP_HOST" \
  SMTP_PORT="$SMTP_PORT" \
  SMTP_USER="$SMTP_USER" \
  SMTP_PASSWORD="$smtp_password" \
  SMTP_STARTTLS="$SMTP_STARTTLS" \
  python3 - <<'PY'
import os
import smtplib
import ssl
from email.message import EmailMessage

msg = EmailMessage()
msg["To"] = os.environ["SMTP_TO"]
msg["From"] = os.environ["SMTP_FROM"]
msg["Subject"] = os.environ["SMTP_SUBJECT"]
msg.set_content(os.environ["SMTP_BODY"])

host = os.environ["SMTP_HOST"]
port = int(os.environ.get("SMTP_PORT", "587"))
user = os.environ["SMTP_USER"]
password = os.environ["SMTP_PASSWORD"]
starttls = os.environ.get("SMTP_STARTTLS", "1") not in ("0", "false", "False")

with smtplib.SMTP(host, port, timeout=30) as server:
    server.ehlo()
    if starttls:
        server.starttls(context=ssl.create_default_context())
        server.ehlo()
    server.login(user, password)
    server.send_message(msg)
PY
}

build_progress_body() {
  local run_log="$1"
  local last_line rung_order_line rung_order current_started last_completed
  local total_rungs current_pass current_rung completed_current runtime_pid runtime_etime
  local smart_lines total_lines active_cmd

  last_line="$(tail -n 1 "$run_log" 2>/dev/null || true)"
  rung_order_line="$(rg -n "\\[smart-ladder\\] pass=[0-9]+ rung-order=" "$run_log" | tail -n1 || true)"
  rung_order="$(echo "$rung_order_line" | sed -E 's/^[0-9]+:.*rung-order=//')"
  if [[ -n "$rung_order" ]]; then
    total_rungs="$(awk -F'->' '{print NF}' <<< "$rung_order")"
  else
    total_rungs="unknown"
  fi

  current_started="$(rg "\\[smart-ladder\\] pass=[0-9]+ rung=[^ ]+ backend=" "$run_log" | tail -n1 || true)"
  last_completed="$(rg "\\[smart-ladder\\] pass=[0-9]+ rung=[^ ]+ completed" "$run_log" | tail -n1 || true)"

  if [[ -n "$current_started" ]]; then
    current_pass="$(sed -E 's/.*pass=([0-9]+).*/\1/' <<< "$current_started")"
    current_rung="$(sed -E 's/.*rung=([^ ]+).*/\1/' <<< "$current_started")"
  elif [[ -n "$last_completed" ]]; then
    current_pass="$(sed -E 's/.*pass=([0-9]+).*/\1/' <<< "$last_completed")"
    current_rung="$(sed -E 's/.*rung=([^ ]+).*/\1/' <<< "$last_completed")"
  else
    current_pass="unknown"
    current_rung="unknown"
  fi

  if [[ "$current_pass" =~ ^[0-9]+$ ]]; then
    completed_current="$(count_matches "\\[smart-ladder\\] pass=${current_pass} rung=[^ ]+ completed" "$run_log")"
  else
    completed_current="0"
  fi

  runtime_pid="$(pgrep -n -f 'performance_runtime_analysis.jl' || true)"
  if [[ -n "$runtime_pid" ]]; then
    runtime_etime="$(ps -p "$runtime_pid" -o etime= 2>/dev/null | xargs || true)"
  else
    runtime_etime="not_running"
  fi

  smart_lines="$(count_matches "\\[smart-ladder\\]" "$run_log")"
  total_lines="$(wc -l < "$run_log" 2>/dev/null || echo 0)"
  active_cmd="$(tmux display-message -p -t "$SESSION_NAME" "#{pane_current_command}" 2>/dev/null || true)"
  if [[ -z "$active_cmd" ]]; then
    active_cmd="unavailable"
  fi

  cat <<EOF
Timestamp: $(timestamp)
Session: $SESSION_NAME
tmux_active_command: $active_cmd
Run log: $run_log
Log lines: $total_lines
Smart-ladder markers: $smart_lines

Current pass: $current_pass
Current rung: $current_rung
Completed in current pass: $completed_current / $total_rungs
Active runtime PID: ${runtime_pid:-none}
Active runtime elapsed: $runtime_etime

Last smart-ladder start:
${current_started:-<none>}

Last smart-ladder completion:
${last_completed:-<none>}

Latest log line:
${last_line:-<none>}
EOF
}

send_email() {
  local subject="$1"
  local body="$2"
  local mailer

  mailer="$(detect_mailer || true)"
  if [[ -z "$mailer" ]]; then
    if send_email_via_python_smtp "$subject" "$body"; then
      return 0
    fi
    spool_dir="$OUT_DIR/email_spool"
    mkdir -p "$spool_dir"
    spool_file="$spool_dir/$(TZ="$TZ_NAME" date +%Y%m%d_%H%M%S).txt"
    {
      echo "To: $EMAIL_TO"
      echo "From: $EMAIL_FROM"
      echo "Subject: $subject"
      echo
      echo "$body"
    } > "$spool_file"
    log "No mailer found. Spooling email to $spool_file"
    return 1
  fi

  case "$mailer" in
    sendmail|msmtp)
      {
        echo "To: $EMAIL_TO"
        echo "From: $EMAIL_FROM"
        echo "Subject: $subject"
        echo
        echo "$body"
      } | "$mailer" -t
      ;;
    mailx|mail)
      printf '%s\n' "$body" | "$mailer" -s "$subject" "$EMAIL_TO"
      ;;
    *)
      {
        echo "To: $EMAIL_TO"
        echo "From: $EMAIL_FROM"
        echo "Subject: $subject"
        echo
        echo "$body"
      } | "$mailer" -t
      ;;
  esac
}

send_progress_update() {
  local run_log subject body

  if [[ -z "$EMAIL_TO" ]]; then
    log "EMAIL_TO is not set; skipping email send."
    return 1
  fi

  run_log="$(ls -1t "$OUT_DIR"/overnight_*.log 2>/dev/null | head -n1 || true)"
  if [[ -z "$run_log" || ! -f "$run_log" ]]; then
    log "No overnight log found under $OUT_DIR; skipping email send."
    return 1
  fi

  body="$(build_progress_body "$run_log")"
  subject="$EMAIL_SUBJECT_PREFIX $(TZ="$TZ_NAME" date +%Y-%m-%d\ %H:%M\ %Z)"
  if send_email "$subject" "$body"; then
    log "Progress email sent to $EMAIL_TO (log=$(basename "$run_log"))."
    return 0
  fi
  return 1
}

log "Starting progress email daemon: every ${INTERVAL_SECONDS}s, daytime ${DAY_START_HOUR}:00-${DAY_END_HOUR}:00 $TZ_NAME"
if [[ "$RUN_ONCE" == "1" ]]; then
  if in_daytime_window; then
    send_progress_update || true
  else
    log "RUN_ONCE=1 and outside daytime window; skipping."
  fi
  exit 0
fi

while true; do
  if in_daytime_window; then
    send_progress_update || true
  else
    log "Outside daytime window; no email sent."
  fi
  sleep "$INTERVAL_SECONDS"
done
