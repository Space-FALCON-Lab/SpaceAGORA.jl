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

# Optional sparse retrieval: GRAMSUITE_SPARSE_PROFILE names a path manifest in
# .github/gramsuite-sparse/ (sparse-checkout patterns, paths only). The fetch
# then downloads the commit's trees but only the blobs those patterns select.
PROFILE="${GRAMSUITE_SPARSE_PROFILE:-}"
PATTERNS=""
if [[ -n "${PROFILE}" ]]; then
  PATTERNS=".github/gramsuite-sparse/${PROFILE}.txt"
  if [[ ! "${PROFILE}" =~ ^[a-z0-9_-]+$ || ! -f "${PATTERNS}" ]]; then
    echo "Unknown GRAMSuite sparse profile: ${PROFILE}" >&2
    exit 1
  fi
fi
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

FILTER=()
if [[ -n "${PROFILE}" ]]; then
  # A partial clone: blobs outside the patterns are never downloaded, and the
  # ones inside are fetched on demand at checkout.
  git -C "${SUBMODULE_PATH}" config core.repositoryformatversion 1
  git -C "${SUBMODULE_PATH}" config extensions.partialClone origin
  git -C "${SUBMODULE_PATH}" config remote.origin.promisor true
  git -C "${SUBMODULE_PATH}" config remote.origin.partialclonefilter blob:none
  git -C "${SUBMODULE_PATH}" config core.sparseCheckout true
  git -C "${SUBMODULE_PATH}" config core.sparseCheckoutCone false
  GIT_DIR_PATH="$(git -C "${SUBMODULE_PATH}" rev-parse --absolute-git-dir)"
  mkdir -p "${GIT_DIR_PATH}/info"
  grep -v '^[[:space:]]*\(#\|$\)' "${PATTERNS}" > "${GIT_DIR_PATH}/info/sparse-checkout"
  FILTER=(--filter=blob:none)
fi

# No branch, remote HEAD, or fallback can replace this exact commit. Avoid
# tags and nested submodules so their moving refs cannot affect retrieval.
# FILTER is empty without a sparse profile; bash before 4.4 (macOS's default
# 3.2) rejects "${FILTER[@]}" of an empty array under set -u, so expand it only
# when set.
with_credentials -C "${SUBMODULE_PATH}" -c protocol.version=2 fetch \
  --depth=1 --no-tags --recurse-submodules=no ${FILTER[@]+"${FILTER[@]}"} origin "${REVISION}"
FETCHED="$(git -C "${SUBMODULE_PATH}" rev-parse --verify 'FETCH_HEAD^{commit}')"
if [[ "${FETCHED}" != "${REVISION}" ]]; then
  echo "Fetched GRAMSuite commit does not equal the ${MODE} pin" >&2
  exit 1
fi
# The full dev tree is ~15 GB packed and 18 GB checked out; writing it out is
# the slow part of this script, so let Git inflate and write files on every
# core.
with_credentials -C "${SUBMODULE_PATH}" -c checkout.workers=0 checkout --detach "${REVISION}"
OBSERVED="$(git -C "${SUBMODULE_PATH}" rev-parse --verify HEAD)"
if [[ "${OBSERVED}" != "${REVISION}" ]]; then
  echo "Checked-out GRAMSuite commit does not equal the ${MODE} pin" >&2
  exit 1
fi
if [[ -n "${PROFILE}" ]]; then
  # A path the profile names but the pinned revision lacks (renamed, moved)
  # must fail here, not make a later test skip. Wildcard and negated patterns
  # are not checked.
  missing=0
  while IFS= read -r pattern; do
    [[ -z "${pattern}" || "${pattern}" == \#* || "${pattern}" == \!* || "${pattern}" == *\** ]] && continue
    rel="${pattern#/}"
    if [[ ! -e "${SUBMODULE_PATH}/${rel%/}" ]]; then
      echo "GRAMSuite sparse profile ${PROFILE}: required path missing: ${rel}" >&2
      missing=1
    fi
  done < "${PATTERNS}"
  [[ "${missing}" -eq 0 ]] || exit 1
  printf 'GRAMSuite sparse profile %s: %s files, %s bytes checked out\n' "${PROFILE}" \
    "$(git -C "${SUBMODULE_PATH}" ls-files -t | grep -c '^H')" \
    "$(cd "${SUBMODULE_PATH}" && git ls-files -t -z | tr '\0' '\n' | sed -n 's/^H //p' | tr '\n' '\0' | xargs -0 stat -c %s | awk '{s+=$1} END {print s+0}')"
fi
printf 'GRAMSuite %s revision: %s\n' "${MODE}" "${OBSERVED}"
