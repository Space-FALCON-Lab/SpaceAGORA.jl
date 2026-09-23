#!/usr/bin/env bash
# WS11c collector-settings grid driver. Runs collector_grid.jl once per
# (gcthreads, heap-size-hint) combination, at a single thread count, for one
# case, sequentially (one Julia process at a time, each under a hard memory
# cap -- see docs/architecture/heap_contention.md's incident note on why the
# cap is non-negotiable for anything launched from this study).
#
# Usage:
#   benchmarks/studies/heap_contention/collector_grid.sh <case> <threads> <out.csv>
#
# This workstation has 12 physical cores; the common WS11 rule caps any one
# agent at 8 threads, which this driver honors (do not raise --threads past
# 8 here without an explicit release). The 16+-thread regime the contract
# actually asks about is out of reach on this box regardless of the flag --
# see the doc for what TRX50 should run instead.
set -euo pipefail

CASE="${1:?usage: collector_grid.sh <case> <threads> <out.csv>}"
THREADS="${2:?usage: collector_grid.sh <case> <threads> <out.csv>}"
OUT="${3:?usage: collector_grid.sh <case> <threads> <out.csv>}"
REPEATS="${REPEATS:-5}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SCRIPT="$REPO_ROOT/benchmarks/studies/heap_contention/collector_grid.jl"

mkdir -p "$(dirname "$OUT")"
echo "case,threads,gcthreads,heap_size_hint,median_wall_s,min_wall_s,max_wall_s" > "$OUT"

# gcthreads: "" (Julia default) and an explicit split matching this box's
# physical-core count divided roughly in half between mark/sweep.
# heap-size-hint: "" (Julia default, no soft cap) and a hint sized to comfortably
# hold one P4 sample's ~365 MiB allocation several times over without forcing
# early collection.
GCTHREADS_GRID=("" "2")
HEAP_HINT_GRID=("" "4G")

for gct in "${GCTHREADS_GRID[@]}"; do
  for heap in "${HEAP_HINT_GRID[@]}"; do
    uptime
    running=$(ps -eo pcpu,args | grep '[j]ulia' | grep -v 'language-julia\|vscode' | wc -l)
    if [ "$running" -gt 6 ]; then
      echo "[collector-grid] backing off: $running other julia processes already running" >&2
      sleep 20
    fi

    flags=()
    [ -n "$gct" ] && flags+=("--gcthreads=$gct")
    [ -n "$heap" ] && flags+=("--heap-size-hint=$heap")

    echo "[collector-grid] case=$CASE threads=$THREADS gcthreads='${gct}' heap='${heap}'" >&2
    line=$(systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
      julia --project="$REPO_ROOT" --threads="$THREADS" "${flags[@]}" \
      "$SCRIPT" --case="$CASE" --repeats="$REPEATS" | grep '^RESULT ')

    median=$(echo "$line" | grep -oP 'median_wall_s=\K[0-9.]+')
    minv=$(echo "$line" | grep -oP 'min_wall_s=\K[0-9.]+')
    maxv=$(echo "$line" | grep -oP 'max_wall_s=\K[0-9.]+')
    echo "$CASE,$THREADS,\"$gct\",\"$heap\",$median,$minv,$maxv" >> "$OUT"
  done
done

echo "[collector-grid] wrote $OUT" >&2
