#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./build_gram.sh [--clean] [GRAM_ROOT]

Builds GRAM shared library for simulation integration.

Arguments:
  --clean       run a clean build before shared build
  GRAM_ROOT     optional absolute/relative path to GRAM Suite root

Env:
  GRAM_ROOT     used if argument is not provided
EOF
}

cpu_count() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
    return
  fi
  if [[ -n "${NUMBER_OF_PROCESSORS:-}" ]]; then
    echo "${NUMBER_OF_PROCESSORS}"
    return
  fi
  if [[ "$(uname -s)" == "Darwin" ]]; then
    sysctl -n hw.ncpu
    return
  fi
  echo 4
}

is_windows_like() {
  case "$1" in
    MINGW*|MSYS*|CYGWIN*|Windows_NT) return 0 ;;
    *) return 1 ;;
  esac
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLEAN_BUILD=0
GRAM_ROOT_ARG=""
BUILD_MANIFEST_NAME=".gram-build-manifest"

host_tag() {
  local os="$1"
  local arch
  arch="$(uname -m)"
  case "$os" in
    Darwin) echo "macos-${arch}" ;;
    Linux) echo "linux-${arch}" ;;
    *) echo "windows-${arch}" ;;
  esac
}

expected_lib_ext() {
  local os="$1"
  if [[ "$os" == "Darwin" ]]; then
    echo "dylib"
  elif is_windows_like "$os"; then
    echo "dll"
  else
    echo "so"
  fi
}

write_build_manifest() {
  local manifest_path="$1"
  local root="$2"
  local lib="$3"
  local host="$4"
  cat > "${manifest_path}" <<EOF
GRAM_ROOT="${root}"
GRAM_LIB="${lib}"
GRAM_HOST="${host}"
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)
      CLEAN_BUILD=1
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

GRAM_ROOT="${GRAM_ROOT_ARG:-${GRAM_ROOT:-}}"
if [[ -z "${GRAM_ROOT}" ]]; then
  CANDIDATES=(
    "${SCRIPT_DIR}/../.."
    "${SCRIPT_DIR}/../../GRAM Suite 2.0"
    "${SCRIPT_DIR}/../../GRAM"
    "${SCRIPT_DIR}/../GRAM Suite 2.0"
    "${SCRIPT_DIR}/../GRAM"
  )
  for c in "${CANDIDATES[@]}"; do
    if [[ -d "${c}/Build" && -d "${c}/Julia" ]]; then
      GRAM_ROOT="${c}"
      break
    fi
  done
fi

if [[ -z "${GRAM_ROOT}" ]]; then
  echo "Could not find GRAM root. Pass path as argument or set GRAM_ROOT." >&2
  exit 1
fi

GRAM_ROOT="$(cd "${GRAM_ROOT}" && pwd)"
if [[ ! -d "${GRAM_ROOT}/Build" ]]; then
  echo "Invalid GRAM_ROOT: ${GRAM_ROOT}" >&2
  exit 1
fi

OS="$(uname -s)"
HOST_TAG="$(host_tag "${OS}")"
if [[ "${OS}" == "Darwin" ]]; then
  if command -v gmake >/dev/null 2>&1; then
    MAKE_BIN="$(command -v gmake)"
  elif [[ -x /opt/homebrew/bin/gmake ]]; then
    MAKE_BIN="/opt/homebrew/bin/gmake"
  else
    echo "gmake not found. Install GNU make (brew install make)." >&2
    exit 1
  fi
elif is_windows_like "${OS}"; then
  if command -v mingw32-make >/dev/null 2>&1; then
    MAKE_BIN="$(command -v mingw32-make)"
  elif command -v make >/dev/null 2>&1; then
    MAKE_BIN="$(command -v make)"
  else
    echo "No MinGW make found. Install MSYS2/MinGW and ensure mingw32-make (or make) is on PATH." >&2
    exit 1
  fi
else
  if ! command -v make >/dev/null 2>&1; then
    echo "make not found on PATH." >&2
    exit 1
  fi
  MAKE_BIN="$(command -v make)"
fi

echo "GRAM_ROOT: ${GRAM_ROOT}"
echo "MAKE_BIN:  ${MAKE_BIN}"
echo "OS:        ${OS}"
echo "HOST_TAG:  ${HOST_TAG}"

BUILD_MANIFEST="${SCRIPT_DIR}/${BUILD_MANIFEST_NAME}"
GRAM_LIB="${GRAM_ROOT}/Build/lib/libGRAM.$(expected_lib_ext "${OS}")"
REQUESTED_GRAM_ROOT="${GRAM_ROOT}"
if [[ "${CLEAN_BUILD}" -eq 0 ]]; then
  if [[ -f "${BUILD_MANIFEST}" ]]; then
    # shellcheck disable=SC1090
    source "${BUILD_MANIFEST}"
    if [[ "${GRAM_HOST:-}" != "${HOST_TAG}" || "${GRAM_ROOT:-}" != "${REQUESTED_GRAM_ROOT}" ]]; then
      CLEAN_BUILD=1
      echo "Detected stale build artifacts for host '${GRAM_HOST:-unknown}'. Rebuilding for ${HOST_TAG}."
    fi
    GRAM_ROOT="${REQUESTED_GRAM_ROOT}"
  elif [[ -d "${GRAM_ROOT}/Build/lib" && ! -f "${GRAM_LIB}" ]]; then
    CLEAN_BUILD=1
    echo "Detected existing build outputs without a native ${GRAM_LIB##*/}. Rebuilding clean."
  fi
fi

pushd "${GRAM_ROOT}/Build" >/dev/null
./setup_cspice.sh
if [[ "${CLEAN_BUILD}" -eq 1 ]]; then
  "${MAKE_BIN}" clean
fi
"${MAKE_BIN}" shared -j"$(cpu_count)"
popd >/dev/null

if [[ ! -f "${GRAM_LIB}" ]]; then
  echo "Build finished but shared library not found: ${GRAM_LIB}" >&2
  exit 1
fi

ENV_FILE="${SCRIPT_DIR}/gram.env"
cat > "${ENV_FILE}" <<EOF
export GRAM_ROOT="${GRAM_ROOT}"
export GRAM_LIB="${GRAM_LIB}"
EOF
write_build_manifest "${BUILD_MANIFEST}" "${GRAM_ROOT}" "${GRAM_LIB}" "${HOST_TAG}"

echo "Build complete."
echo "GRAM_LIB: ${GRAM_LIB}"
echo "Wrote:    ${ENV_FILE}"
echo "Wrote:    ${BUILD_MANIFEST}"
