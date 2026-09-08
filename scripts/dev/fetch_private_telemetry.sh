#!/usr/bin/env bash
# Sync verification data that lives outside this public repository into the
# gitignored data/telemetry/ paths SpaceAGORA expects. Two sources:
#
#   references  Simulator-to-simulator reference sets (Basilisk, GMAT parity
#               matrix) from the lab-org repo Space-FALCON-Lab/
#               spaceagora-verification-data. Its manifest.toml names each
#               dataset, its path in that repo, and the SpaceAGORA path it
#               syncs to; every dataset carries SHA256SUMS, verified after sync.
#   telemetry   Flight telemetry provided under data-sharing agreements
#               (CYGNSS, GRIFEX) from the access-restricted repo; per-person
#               access, synced as mission directories.
#
# Usage: scripts/dev/fetch_private_telemetry.sh [references|telemetry|all] [dataset ...]
#   default: all (a source you cannot access is skipped with a warning)
#   dataset: restrict the references source to the named manifest datasets
#
# Environment overrides:
#   SPACEAGORA_VERIFICATION_DATA_REPO    owner/name of the references repo
#   SPACEAGORA_VERIFICATION_DATA_CACHE   local cache clone location
#   SPACEAGORA_PRIVATE_TELEMETRY_REPO    owner/name of the restricted telemetry repo
#   SPACEAGORA_PRIVATE_TELEMETRY_CACHE   local cache clone location
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REF_REPO="${SPACEAGORA_VERIFICATION_DATA_REPO:-Space-FALCON-Lab/spaceagora-verification-data}"
REF_CACHE="${SPACEAGORA_VERIFICATION_DATA_CACHE:-$HOME/.julia/dev/spaceagora-verification-data}"
TEL_REPO="${SPACEAGORA_PRIVATE_TELEMETRY_REPO:-gfalcon2/spaceagora-private-telemetry}"
TEL_CACHE="${SPACEAGORA_PRIVATE_TELEMETRY_CACHE:-$HOME/.julia/dev/spaceagora-private-telemetry-repo}"

SOURCE="${1:-all}"
shift $(( $# > 0 ? 1 : 0 ))
WANTED=("$@")

update_cache() {
    local slug="$1" cache="$2"
    if [ -d "$cache/.git" ]; then
        git -C "$cache" pull -q --ff-only
    elif command -v gh >/dev/null 2>&1; then
        gh repo clone "$slug" "$cache" -- -q
    else
        git clone -q "https://github.com/$slug.git" "$cache"
    fi
}

wanted() {
    [ "${#WANTED[@]}" -eq 0 ] && return 0
    local w
    for w in "${WANTED[@]}"; do [ "$w" = "$1" ] && return 0; done
    return 1
}

sync_references() {
    if ! update_cache "$REF_REPO" "$REF_CACHE"; then
        echo "warning: cannot reach $REF_REPO (no access?); skipping simulator references" >&2
        return 0
    fi
    local manifest="$REF_CACHE/manifest.toml"
    [ -f "$manifest" ] || { echo "error: $manifest not found" >&2; return 1; }
    local synced=0 name path target
    # manifest.toml is a flat list of [[dataset]] tables with quoted string values.
    while IFS=$'\t' read -r name path target; do
        [ -n "$name" ] || continue
        wanted "$name" || continue
        mkdir -p "$ROOT/$target"
        rsync -a --delete "$REF_CACHE/$path/" "$ROOT/$target/"
        if [ -f "$ROOT/$target/SHA256SUMS" ]; then
            (cd "$ROOT/$target" && shasum -a 256 --check --quiet SHA256SUMS) \
                || { echo "error: checksum mismatch after syncing $name" >&2; return 1; }
        fi
        echo "synced $name -> $target/ ($(ls "$ROOT/$target" | grep -c '\.feather$') feather files, checksums ok)"
        synced=$((synced + 1))
    done < <(awk -F'"' '
        /^\[\[dataset\]\]/ { if (name != "") print name "\t" path "\t" target; name=""; path=""; target="" }
        /^name = /            { name=$2 }
        /^path = /            { path=$2 }
        /^spaceagora_path = / { target=$2 }
        END { if (name != "") print name "\t" path "\t" target }' "$manifest")
    [ "$synced" -gt 0 ] || echo "warning: no reference dataset matched" >&2
}

sync_telemetry() {
    if ! update_cache "$TEL_REPO" "$TEL_CACHE"; then
        echo "warning: cannot reach $TEL_REPO (access is granted per person); skipping flight telemetry" >&2
        return 0
    fi
    local synced=0 mission_dir mission
    for mission_dir in "$TEL_CACHE"/*/; do
        mission="$(basename "$mission_dir")"
        case "$mission" in
            .git|.github|paper_drafts) continue ;;
        esac
        mkdir -p "$ROOT/data/telemetry/$mission"
        rsync -a "$mission_dir" "$ROOT/data/telemetry/$mission/"
        echo "synced $mission -> data/telemetry/$mission/"
        synced=$((synced + 1))
    done
    [ "$synced" -gt 0 ] || echo "warning: no mission directories found in $TEL_CACHE" >&2
}

case "$SOURCE" in
    references) sync_references ;;
    telemetry)  sync_telemetry ;;
    all)        sync_references; sync_telemetry ;;
    *) echo "usage: $0 [references|telemetry|all] [dataset ...]" >&2; exit 2 ;;
esac
echo "Done. These paths are gitignored - never commit them; flight telemetry must not be redistributed."
