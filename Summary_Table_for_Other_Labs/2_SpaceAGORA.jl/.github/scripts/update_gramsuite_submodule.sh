#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-regular}"

if [[ "${MODE}" != "regular" && "${MODE}" != "dev" ]]; then
  echo "Usage: $0 [regular|dev]" >&2
  exit 1
fi

GH_TOKEN="${GH_TOKEN:-}"
if [[ -n "${GH_TOKEN}" ]]; then
  git config --global url."https://x-access-token:${GH_TOKEN}@github.com/".insteadOf "https://github.com/"
fi

if [[ "${MODE}" == "dev" ]]; then
  git config -f .gitmodules submodule.GRAMSuite.jl.url \
    "https://github.com/Space-FALCON-Lab/dev-GRAMSuite.jl.git"
else
  git config -f .gitmodules submodule.GRAMSuite.jl.url \
    "https://github.com/Space-FALCON-Lab/GRAMSuite.jl.git"
fi

git config -f .gitmodules submodule.GRAMSuite.jl.shallow true
git config -f .gitmodules submodule.GRAMSuite.jl.fetchRecurseSubmodules false

git submodule sync data/GRAMSuite.jl
git -c protocol.version=2 submodule update --init --depth=1 --recommend-shallow --remote data/GRAMSuite.jl
