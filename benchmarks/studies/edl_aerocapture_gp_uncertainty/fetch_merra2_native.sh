#!/usr/bin/env bash
# Fetch MERRA-2 native model-level granules (M2I3NVASM, inst3_3d_asm_Nv) for the
# EDL aerocapture GP study.
#
# These are the day-specific reanalysis granules that `merra2_native.jl` reads.
# They are NOT the climatologies vendored with Earth-GRAM, and they are not in
# this repository: each granule is ~2.1 GB and GES DISC requires a NASA
# Earthdata login.
#
# One-time setup
# --------------
#   1. Register: https://urs.earthdata.nasa.gov/users/new
#   2. Authorize the "NASA GESDISC DATA ARCHIVE" application. One click, while
#      logged in to Earthdata:
#
#        https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ
#
#      Until this is done the server answers 200/206 with an HTML page rather
#      than an error, so a naive download silently writes HTML into a .nc4 file.
#      This script checks the HDF5 magic bytes and refuses such a file.
#   3. Authenticate. A bearer token is preferred over a password: it is
#      revocable on its own, scoped, and expires. Generate one at
#      https://urs.earthdata.nasa.gov/profile -> Generate Token, then:
#
#        umask 077; printf '%s' "PASTE_TOKEN" > ~/.edl_token
#
#      The script also accepts $EARTHDATA_TOKEN, and falls back to ~/.netrc
#      (machine urs.earthdata.nasa.gov login USER password PASS, chmod 600) if
#      no token is present.
#
# Usage
# -----
#   ./fetch_merra2_native.sh 2024-03-20 [2024-03-21 ...]
#   ./fetch_merra2_native.sh --dir /path/to/store 2024-03-20
#   ./fetch_merra2_native.sh --dry-run 2024-03-20
#
# The study's default cases anchor at 2024-03-20, 2024-10-05 and 2025-03-21,
# each 18:00 UTC. A 6 hr aerocapture-to-EDL gap crosses midnight, so fetch the
# following day as well:
#
#   ./fetch_merra2_native.sh 2024-03-20 2024-03-21 2024-10-05 2024-10-06 \
#                            2025-03-21 2025-03-22
#
# That is ~12.6 GB. If bandwidth or disk is tight, GES DISC offers spatial
# subsetting through its Subset/Get Data workflow on the collection landing page
# (https://disc.gsfc.nasa.gov/datasets/M2I3NVASM_5.12.4/summary); a
# corridor-sized lat/lon window is a few tens of MB per day. `merra2_native.jl`
# reads subsetted files unchanged as long as the 72 model levels are kept, but
# note the spec's warning that some subsetters alter metadata or flip the
# vertical grid direction — the reader asserts the level count and the tests
# pin the ordering, so a flipped file will produce an obviously wrong profile
# rather than a silently plausible one. Check with --verify below.

set -euo pipefail

COLLECTION="M2I3NVASM.5.12.4"
BASE="https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/${COLLECTION}"
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)/data/merra2_native"
DRY_RUN=0
DATES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir)     DEST="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,45p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*)        echo "Unknown option: $1" >&2; exit 2 ;;
    *)         DATES+=("$1"); shift ;;
  esac
done

if [[ ${#DATES[@]} -eq 0 ]]; then
  echo "No dates given. Try: $0 2024-03-20 2024-03-21" >&2
  exit 2
fi

# Prefer a bearer token; fall back to ~/.netrc.
#
# --location-trusted is required with a token: GES DISC bounces through
# urs.earthdata.nasa.gov, and plain -L drops the Authorization header on a
# cross-host redirect, which lands you back on the HTML login page. Both hosts
# are NASA EDL endpoints.
TOKEN=""
if [[ -n "${EARTHDATA_TOKEN:-}" ]]; then
  TOKEN="$EARTHDATA_TOKEN"
elif [[ -f "$HOME/.edl_token" ]]; then
  TOKEN="$(tr -d '\r\n' < "$HOME/.edl_token")"
fi

if [[ -n "$TOKEN" ]]; then
  AUTH_ARGS=(--location-trusted -H "Authorization: Bearer ${TOKEN}")
elif [[ -f "$HOME/.netrc" ]]; then
  AUTH_ARGS=(-L -n)
else
  cat >&2 <<'EOF'
No Earthdata credentials found. GES DISC requires authentication.

  1. Register at https://urs.earthdata.nasa.gov/users/new
  2. Authorize: https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ
  3. Generate a token at https://urs.earthdata.nasa.gov/profile, then:
       umask 077; printf '%s' "PASTE_TOKEN" > ~/.edl_token

  (A ~/.netrc with login/password also works, but a token is revocable on its
  own and expires, so prefer it.)
EOF
  exit 1
fi

# MERRA-2 production streams; the number is part of every filename.
stream_for_year() {
  local y="$1"
  if   [[ "$y" -le 1991 ]]; then echo 100
  elif [[ "$y" -le 2000 ]]; then echo 200
  elif [[ "$y" -le 2010 ]]; then echo 300
  else                           echo 400
  fi
}

mkdir -p "$DEST"
COOKIES="$(mktemp)"
trap 'rm -f "$COOKIES"' EXIT

for d in "${DATES[@]}"; do
  if ! date -d "$d" +%Y >/dev/null 2>&1; then
    echo "Not a date: $d" >&2
    exit 2
  fi
  Y=$(date -d "$d" +%Y); M=$(date -d "$d" +%m); STAMP=$(date -d "$d" +%Y%m%d)
  S=$(stream_for_year "$Y")
  FILE="MERRA2_${S}.inst3_3d_asm_Nv.${STAMP}.nc4"
  URL="${BASE}/${Y}/${M}/${FILE}"
  OUT="${DEST}/${FILE}"

  if [[ -s "$OUT" ]]; then
    echo "have  ${FILE}"
    continue
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "would fetch ${URL}"
    continue
  fi

  echo "fetch ${FILE} (~2.1 GB)"
  # The cookie jar carries the URS session through the redirect chain.
  # -C - resumes a partial file.
  if ! curl -fS "${AUTH_ARGS[@]}" -c "$COOKIES" -b "$COOKIES" -C - -o "$OUT.part" "$URL"; then
    rm -f "$OUT.part"
    cat >&2 <<EOF

Download failed for ${FILE}.

HTTP 401 almost always means the "NASA GESDISC DATA ARCHIVE" application has not
been authorized on your Earthdata profile, rather than bad credentials:
  https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ
If you are using a token, check it has not expired.
EOF
    exit 1
  fi
  # An unauthorized fetch does NOT return an HTTP error: GES DISC answers 200/206
  # and serves an HTML page, so curl -f is not enough. netCDF-4 is HDF5, which
  # always begins with the 8-byte signature \x89HDF\r\n\x1a\n.
  MAGIC=$(head -c 4 "$OUT.part" | od -An -tx1 | tr -d ' \n')
  if [[ "$MAGIC" != "89484446" ]]; then
    rm -f "$OUT.part"
    cat >&2 <<EOF

Downloaded ${FILE} is not an HDF5/netCDF-4 file (magic: ${MAGIC:-empty}).

This is almost always the GES DISC application-authorization redirect, which is
served as an HTML page with a 200/206 status rather than an error. Authorize the
archive once, while logged in to Earthdata:

  https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

then re-run this script. To confirm what the server is actually returning:

  curl -sS ${AUTH_ARGS[@]} -r 0-100 -D - -o /dev/null "${URL}" | grep -iE '^HTTP/|^location:'
EOF
    exit 1
  fi
  mv "$OUT.part" "$OUT"
done

cat <<EOF

Granules in ${DEST}:
$(ls -1sh "$DEST" 2>/dev/null | tail -n +2 || echo "  (none)")

Verify the reader against them before trusting any results:

  julia --project=benchmarks/studies/edl_aerocapture_gp_uncertainty \\
        benchmarks/studies/edl_aerocapture_gp_uncertainty/verify_merra2_native.jl

Then run the study with:

  --truth-source merra2_native
EOF
