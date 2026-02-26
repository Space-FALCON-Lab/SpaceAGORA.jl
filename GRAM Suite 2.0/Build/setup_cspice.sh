#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_DIR="${REPO_ROOT}/common/cspice/lib"
mkdir -p "${LIB_DIR}"

os="$(uname -s)"
arch="$(uname -m)"

echo "Detected platform: ${os} ${arch}"

copy_cspice() {
  local src="$1"
  local dst="$2"
  cp -f "${src}" "${dst}"
  echo "Installed CSPICE archive:"
  echo "  ${dst}"
}

if [[ "${os}" == "Darwin" ]]; then
  brew_bin=""
  if command -v brew >/dev/null 2>&1; then
    brew_bin="$(command -v brew)"
  elif [[ -x /opt/homebrew/bin/brew ]]; then
    brew_bin="/opt/homebrew/bin/brew"
  elif [[ -x /usr/local/bin/brew ]]; then
    brew_bin="/usr/local/bin/brew"
  fi
  if [[ -z "${brew_bin}" ]]; then
    echo "Homebrew is required on macOS to auto-provision CSPICE."
    echo "Install it, then run: brew install cspice"
    exit 1
  fi
  prefix="$("${brew_bin}" --prefix cspice 2>/dev/null || true)"
  if [[ -z "${prefix}" || ! -f "${prefix}/lib/cspice.a" ]]; then
    echo "Could not find Homebrew cspice. Run: brew install cspice"
    exit 1
  fi
  if [[ "${arch}" == "arm64" ]]; then
    copy_cspice "${prefix}/lib/cspice.a" "${LIB_DIR}/cspice_macos_arm64.a"
  elif [[ "${arch}" == "x86_64" ]]; then
    copy_cspice "${prefix}/lib/cspice.a" "${LIB_DIR}/cspice_macos_x86_64.a"
  else
    echo "Unsupported macOS architecture: ${arch}"
    exit 1
  fi
  exit 0
fi

if [[ "${os}" == "Linux" ]]; then
  case "${arch}" in
    x86_64)
      dst="${LIB_DIR}/cspice_linux_x86_64.a"
      if [[ -f "${LIB_DIR}/cspice_gcc85.a" ]]; then
        copy_cspice "${LIB_DIR}/cspice_gcc85.a" "${dst}"
        exit 0
      fi
      ;;
    aarch64|arm64)
      dst="${LIB_DIR}/cspice_linux_aarch64.a"
      ;;
    *)
      dst="${LIB_DIR}/cspice_linux_${arch}.a"
      ;;
  esac

  candidates=(
    "/usr/lib/libcspice.a"
    "/usr/local/lib/libcspice.a"
    "/usr/lib64/libcspice.a"
    "/usr/lib/x86_64-linux-gnu/libcspice.a"
    "/usr/lib/aarch64-linux-gnu/libcspice.a"
  )
  for c in "${candidates[@]}"; do
    if [[ -f "${c}" ]]; then
      copy_cspice "${c}" "${dst}"
      exit 0
    fi
  done

  echo "Could not auto-find a Linux CSPICE archive."
  echo "Install CSPICE for this machine and then either:"
  echo "  1) Re-run this script, or"
  echo "  2) Build with SPICE_LIB=/absolute/path/to/cspice.a"
  exit 1
fi

if [[ "${os}" == MINGW* || "${os}" == MSYS* || "${os}" == CYGWIN* || "${os}" == "Windows_NT" ]]; then
  mingw_archive="${LIB_DIR}/cspice_mingw64.a"
  if [[ -f "${mingw_archive}" ]]; then
    echo "Using bundled Windows CSPICE archive:"
    echo "  ${mingw_archive}"
    exit 0
  fi

  echo "Could not find MinGW CSPICE archive at:"
  echo "  ${mingw_archive}"
  echo "Provide a MinGW-compatible CSPICE static archive and build with:"
  echo "  SPICE_LIB=/absolute/path/to/cspice_mingw64.a"
  exit 1
fi

echo "Unsupported host OS: ${os}"
exit 1
