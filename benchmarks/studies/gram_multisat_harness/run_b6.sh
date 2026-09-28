#!/bin/bash
# Reproduce the B6 interacting-GRAM case at its real full mission duration,
# with RSS sampled every 5 s so the reported leak is measured, not inferred.
set -u
REPO=/home/space-falcon-1/Documents/SpaceAGORA.jl
OUT="$1"; CASE="$2"; MODE="$3"; THREADS="$4"; LIMIT_S="${5:-900}"
cd "$REPO"
env ${EXTRA_ENV:-FOO=bar} julia --project=. -t "$THREADS" benchmarks/studies/parallelization_performance.jl full \
  --worker --case="$CASE" --mode="$MODE" --thread-count="$THREADS" \
  --repeat=1 --outfile="$OUT.json" > "$OUT.log" 2>&1 &
PID=$!
echo "pid=$PID case=$CASE mode=$MODE threads=$THREADS" > "$OUT.mem"
S=0
while kill -0 $PID 2>/dev/null; do
  RSS=$(awk '/VmRSS/{print $2}' /proc/$PID/status 2>/dev/null)
  echo "t=${S}s rss_mb=$(( ${RSS:-0} / 1024 ))" >> "$OUT.mem"
  sleep 5; S=$((S+5))
  if [ $S -ge $LIMIT_S ]; then echo "TIMEOUT at ${S}s -- killing" >> "$OUT.mem"; kill -9 $PID 2>/dev/null; break; fi
done
wait $PID 2>/dev/null
echo "exit=$? elapsed=${S}s" >> "$OUT.mem"
