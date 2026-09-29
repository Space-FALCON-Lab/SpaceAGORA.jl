#!/usr/bin/env bash
# Process-pool heap-growth verification driver. Runs pool_rss_probe.jl twice,
# back to back, one Julia process at a time, each under a hard memory cap:
# once with the pool-worker heap-size-hint off, once with it on (the
# default). While each run is in flight, samples the RSS of every
# julia-named process every 2 s via `ps`, then reports peak RSS per process
# and peak total RSS for both arms, plus the wall-time ratio, via
# summarize_pool_rss.py.
#
# Usage:
#   benchmarks/studies/heap_contention/pool_rss_probe.sh [case] [workers] [samples]
set -uo pipefail

CASE="${1:-montecarlo_mars_gram_live}"
PROFILE="${4:-smoke}"
WORKERS="${2:-4}"
SAMPLES="${3:-800}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SCRIPT="$REPO_ROOT/benchmarks/studies/heap_contention/pool_rss_probe.jl"
OUTDIR="$REPO_ROOT/benchmarks/studies/heap_contention/results"
mkdir -p "$OUTDIR"

run_arm() {
  local label="$1" hint_off="$2"
  local rss_log="$OUTDIR/pool_rss_${label}.log"
  local run_log="$OUTDIR/pool_rss_${label}_run.log"
  : > "$rss_log"

  ( while true; do
      ts=$(date +%s.%N)
      ps -eo pid,ppid,rss,args | grep '[j]ulia' | while IFS= read -r line; do
        echo "$ts $line" >> "$rss_log"
      done
      sleep 2
    done ) &
  local sampler_pid=$!

  echo "[pool-rss] arm=$label starting (sampler pid=$sampler_pid)" >&2
  if [ "$hint_off" = "1" ]; then
    systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
      env SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT=off \
      julia --project="$REPO_ROOT" --threads=1 "$SCRIPT" \
      --case="$CASE" --workers="$WORKERS" --samples="$SAMPLES" --profile="$PROFILE" > "$run_log" 2>&1
  else
    systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
      julia --project="$REPO_ROOT" --threads=1 "$SCRIPT" \
      --case="$CASE" --workers="$WORKERS" --samples="$SAMPLES" --profile="$PROFILE" > "$run_log" 2>&1
  fi
  local rc=$?

  kill "$sampler_pid" 2>/dev/null
  wait "$sampler_pid" 2>/dev/null
  echo "[pool-rss] arm=$label rc=$rc run_log=$run_log rss_log=$rss_log" >&2
  return "$rc"
}

run_arm "nohint" "1"
run_arm "hint" "0"

python3 "$REPO_ROOT/benchmarks/studies/heap_contention/summarize_pool_rss.py" \
  "$OUTDIR/pool_rss_nohint.log" "$OUTDIR/pool_rss_hint.log" \
  "$OUTDIR/pool_rss_nohint_run.log" "$OUTDIR/pool_rss_hint_run.log" \
  | tee "$OUTDIR/pool_rss_summary.txt"
