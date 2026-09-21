#!/usr/bin/env bash
# Sync access-restricted mission telemetry into this checkout.
#
# SpaceAGORA.jl is a public repository; flight telemetry provided under
# data-sharing agreements (CYGNSS, GRIFEX) lives in a separate private repo
# and is gitignored here (see data/telemetry/PRIVATE_TELEMETRY.md). This
# script clones/updates a local cache of that repo and syncs its mission
# directories into data/telemetry/. Requires access to the private repo
# (gh auth or git credentials).
#
# Environment overrides:
#   SPACEAGORA_PRIVATE_TELEMETRY_REPO   owner/name of the private repo
#   SPACEAGORA_PRIVATE_TELEMETRY_CACHE  local cache clone location
set -euo pipefail

REPO_SLUG="${SPACEAGORA_PRIVATE_TELEMETRY_REPO:-gfalcon2/spaceagora-private-telemetry}"
CACHE="${SPACEAGORA_PRIVATE_TELEMETRY_CACHE:-$HOME/.julia/dev/spaceagora-private-telemetry-repo}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

if [ -d "$CACHE/.git" ]; then
    git -C "$CACHE" pull --ff-only
elif command -v gh >/dev/null 2>&1; then
    gh repo clone "$REPO_SLUG" "$CACHE"
else
    git clone "https://github.com/$REPO_SLUG.git" "$CACHE"
fi

synced=0
for mission_dir in "$CACHE"/*/; do
    mission="$(basename "$mission_dir")"
    case "$mission" in
        .git|.github) continue ;;
    esac
    mkdir -p "$ROOT/data/telemetry/$mission"
    rsync -a "$mission_dir" "$ROOT/data/telemetry/$mission/"
    echo "synced $mission -> data/telemetry/$mission/"
    synced=$((synced + 1))
done

if [ "$synced" -eq 0 ]; then
    echo "warning: no mission directories found in $CACHE" >&2
    exit 1
fi
echo "Done. These paths are gitignored - never commit or redistribute them."
