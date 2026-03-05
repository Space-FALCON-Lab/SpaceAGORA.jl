#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${REPO_DIR:-/tmp/SpaceAGORA.jl}"
OUT_DIR="${OUT_DIR:-$REPO_DIR/output/performance}"
EMAIL_TO="${EMAIL_TO:-falconeg@umich.edu}"
ICLOUD_ACCOUNT="${ICLOUD_ACCOUNT:-spacefalconlab@icloud.com}"
DAY_START_HOUR="${DAY_START_HOUR:-6}"
DAY_END_HOUR="${DAY_END_HOUR:-22}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-10800}"
TZ_NAME="${TZ_NAME:-America/Detroit}"
SMTP_PASSWORD_FILE="${SMTP_PASSWORD_FILE:-$HOME/.spaceagora_icloud_app_password}"

mkdir -p "$OUT_DIR"

if [[ ! -s "$SMTP_PASSWORD_FILE" ]]; then
  if [[ ! -t 0 ]]; then
    echo "Missing SMTP password file: $SMTP_PASSWORD_FILE"
    echo "Run this script in an interactive terminal to enter your iCloud app-specific password."
    exit 1
  fi
  read -r -s -p "Enter iCloud app-specific password for $ICLOUD_ACCOUNT: " icloud_app_password
  echo
  if [[ -z "$icloud_app_password" ]]; then
    echo "Password cannot be empty."
    exit 1
  fi
  umask 177
  printf '%s\n' "$icloud_app_password" > "$SMTP_PASSWORD_FILE"
  chmod 600 "$SMTP_PASSWORD_FILE"
  unset icloud_app_password
  echo "Saved app password to $SMTP_PASSWORD_FILE"
fi

pkill -f progress_email_daemon.sh >/dev/null 2>&1 || true

cd "$REPO_DIR"
nohup env \
  EMAIL_TO="$EMAIL_TO" \
  EMAIL_FROM="$ICLOUD_ACCOUNT" \
  DAY_START_HOUR="$DAY_START_HOUR" \
  DAY_END_HOUR="$DAY_END_HOUR" \
  INTERVAL_SECONDS="$INTERVAL_SECONDS" \
  TZ_NAME="$TZ_NAME" \
  SMTP_HOST="smtp.mail.me.com" \
  SMTP_PORT="587" \
  SMTP_USER="$ICLOUD_ACCOUNT" \
  SMTP_PASSWORD_FILE="$SMTP_PASSWORD_FILE" \
  SMTP_STARTTLS="1" \
  ./progress_email_daemon.sh >> "$OUT_DIR/progress_email_daemon.nohup.log" 2>&1 &

daemon_pid=$!
echo "Started progress daemon pid=$daemon_pid"

env \
  EMAIL_TO="$EMAIL_TO" \
  EMAIL_FROM="$ICLOUD_ACCOUNT" \
  DAY_START_HOUR="$DAY_START_HOUR" \
  DAY_END_HOUR="$DAY_END_HOUR" \
  INTERVAL_SECONDS="$INTERVAL_SECONDS" \
  TZ_NAME="$TZ_NAME" \
  SMTP_HOST="smtp.mail.me.com" \
  SMTP_PORT="587" \
  SMTP_USER="$ICLOUD_ACCOUNT" \
  SMTP_PASSWORD_FILE="$SMTP_PASSWORD_FILE" \
  SMTP_STARTTLS="1" \
  RUN_ONCE="1" \
  ./progress_email_daemon.sh

