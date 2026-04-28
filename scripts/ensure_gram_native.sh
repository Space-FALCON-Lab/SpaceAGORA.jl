#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./scripts/ensure_gram_native.sh [--clean] [GRAM_ROOT]

Ensure a native GRAM shared library exists for the current operating system.

This wrapper checks for the expected platform-native `libGRAM` artifact under
the vendored GRAM tree. If it is missing, it delegates to the vendored
`build_gram.sh` helper to build the correct library for the current host.

Arguments:
  --clean       force a clean rebuild
  GRAM_ROOT     optional GRAM Suite root; defaults to vendored GRAM tree

Outputs:
  - native shared library under `GRAM Suite 2.0/Build/lib/`
  - `simulation/GRAM/gram.env` written by the vendored build helper

Examples:
  ./scripts/ensure_gram_native.sh
  ./scripts/ensure_gram_native.sh --clean
EOF
}

is_windows_like() {
  case "${1}" in
    MINGW*|MSYS*|CYGWIN*|Windows_NT) return 0 ;;
    *) return 1 ;;
  esac
}

expected_lib_ext() {
  local os="$1"
  if [[ "${os}" == "Darwin" ]]; then
    echo "dylib"
  elif is_windows_like "${os}"; then
    echo "dll"
  else
    echo "so"
  fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_GRAM_ROOT="${REPO_ROOT}/data/GRAMSuite.jl/GRAM Suite 2.0"
BUILD_HELPER="${DEFAULT_GRAM_ROOT}/simulation/GRAM/build_gram.sh"

CLEAN_FLAG=0
GRAM_ROOT_ARG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)
      CLEAN_FLAG=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -n "${GRAM_ROOT_ARG}" ]]; then
        echo "Only one GRAM_ROOT argument is allowed." >&2
        usage >&2
        exit 1
      fi
      GRAM_ROOT_ARG="$1"
      ;;
  esac
  shift
done

GRAM_ROOT="${GRAM_ROOT_ARG:-${GRAM_ROOT:-${DEFAULT_GRAM_ROOT}}}"
GRAM_ROOT="$(cd "${GRAM_ROOT}" && pwd)"
[[ -d "${GRAM_ROOT}" ]] || { echo "GRAM root not found: ${GRAM_ROOT}" >&2; exit 1; }
[[ -x "${BUILD_HELPER}" ]] || { echo "GRAM build helper not found: ${BUILD_HELPER}" >&2; exit 1; }

OS="$(uname -s)"
LIB_EXT="$(expected_lib_ext "${OS}")"
EXPECTED_LIB="${GRAM_ROOT}/Build/lib/libGRAM.${LIB_EXT}"

echo "GRAM root:        ${GRAM_ROOT}"
echo "Expected library: ${EXPECTED_LIB}"

if [[ -f "${EXPECTED_LIB}" && "${CLEAN_FLAG}" -eq 0 ]]; then
  echo "Native GRAM library already present for this host."
  exit 0
fi

echo "Native GRAM library missing for this host. Building now..."

BUILD_ARGS=()
if [[ "${CLEAN_FLAG}" -eq 1 ]]; then
  BUILD_ARGS+=(--clean)
fi
BUILD_ARGS+=("${GRAM_ROOT}")

"${BUILD_HELPER}" "${BUILD_ARGS[@]}"

if [[ ! -f "${EXPECTED_LIB}" ]]; then
  echo "Build completed but expected native library is still missing: ${EXPECTED_LIB}" >&2
  exit 1
fi

echo "Native GRAM library ready:"
echo "  ${EXPECTED_LIB}"
