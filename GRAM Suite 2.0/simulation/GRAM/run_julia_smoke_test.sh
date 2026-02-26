#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "${SCRIPT_DIR}/gram.env" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/gram.env"
fi

if [[ $# -eq 1 ]]; then
  export GRAM_ROOT="$1"
elif [[ $# -gt 1 ]]; then
  echo "Usage: ./run_julia_smoke_test.sh [GRAM_ROOT]" >&2
  exit 1
fi

if command -v julia >/dev/null 2>&1; then
  JULIA_BIN="$(command -v julia)"
elif [[ -x /opt/homebrew/bin/julia ]]; then
  JULIA_BIN="/opt/homebrew/bin/julia"
else
  echo "julia not found on PATH." >&2
  exit 1
fi

echo "Using Julia: ${JULIA_BIN}"
"${JULIA_BIN}" "${SCRIPT_DIR}/julia_smoke_test.jl"
