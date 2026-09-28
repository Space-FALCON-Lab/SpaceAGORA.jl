#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-regular}"
if [[ "${MODE}" != "regular" && "${MODE}" != "dev" ]]; then
  echo "Usage: $0 [regular|dev]" >&2
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/../.."
REVISION="$(awk -v mode="${MODE}" '$1 == mode && NF == 2 { print $2 }' .github/gramsuite-revisions)"
if [[ ! "${REVISION}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "Missing or invalid ${MODE} GRAMSuite revision in .github/gramsuite-revisions" >&2
  exit 1
fi

SUBMODULE_PATH="data/GRAMSuite.jl"
if [[ "${MODE}" == "dev" ]]; then
  URL="https://github.com/Space-FALCON-Lab/dev-GRAMSuite.jl.git"
else
  URL="https://github.com/Space-FALCON-Lab/GRAMSuite.jl.git"
fi

# Supply credentials only to these Git calls; never put the token in a URL or
# persistent Git configuration. The helper is restricted to GitHub HTTPS.
export GH_TOKEN
with_credentials() {
  if [[ -n "${GH_TOKEN:-}" ]]; then
    git -c credential.helper= \
      -c 'credential.https://github.com.helper=!f() { if [ "$1" = get ]; then printf "%s\n" "username=x-access-token" "password=$GH_TOKEN"; fi; }; f' \
      "$@"
  else
    git "$@"
  fi
}

# Start with an empty repository rather than checking out the public gitlink.
# The private dev history need not contain the public revision.
git submodule init -- "${SUBMODULE_PATH}"
git config submodule.GRAMSuite.jl.url "${URL}"
if [[ ! -e "${SUBMODULE_PATH}/.git" ]]; then
  git init --quiet "${SUBMODULE_PATH}"
  git submodule absorbgitdirs -- "${SUBMODULE_PATH}"
fi
if ! git -C "${SUBMODULE_PATH}" diff --quiet ||
   ! git -C "${SUBMODULE_PATH}" diff --cached --quiet; then
  echo "GRAMSuite has local tracked changes; refusing to replace its revision" >&2
  exit 1
fi
git -C "${SUBMODULE_PATH}" config remote.origin.url "${URL}"

# No branch, remote HEAD, or fallback can replace this exact commit. Avoid
# tags and nested submodules so their moving refs cannot affect retrieval.
with_credentials -C "${SUBMODULE_PATH}" -c protocol.version=2 fetch \
  --depth=1 --no-tags --recurse-submodules=no origin "${REVISION}"
FETCHED="$(git -C "${SUBMODULE_PATH}" rev-parse --verify 'FETCH_HEAD^{commit}')"
if [[ "${FETCHED}" != "${REVISION}" ]]; then
  echo "Fetched GRAMSuite commit does not equal the ${MODE} pin" >&2
  exit 1
fi
with_credentials -C "${SUBMODULE_PATH}" checkout --detach "${REVISION}"
OBSERVED="$(git -C "${SUBMODULE_PATH}" rev-parse --verify HEAD)"
if [[ "${OBSERVED}" != "${REVISION}" ]]; then
  echo "Checked-out GRAMSuite commit does not equal the ${MODE} pin" >&2
  exit 1
fi
printf 'GRAMSuite %s revision: %s\n' "${MODE}" "${OBSERVED}"
