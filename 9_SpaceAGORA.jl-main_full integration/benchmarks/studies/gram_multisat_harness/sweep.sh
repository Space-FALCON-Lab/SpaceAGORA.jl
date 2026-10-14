#!/bin/bash
# Loose -> tight, one process per trial, hard cap each.
SCRATCH="$1"; CAP="${2:-400}"
cd /home/space-falcon-1/Documents/SpaceAGORA.jl
: > "$SCRATCH/sweep.out"
for RT in 1e-5 1e-6 1e-7 1e-8 1e-9; do
  R=$(B6_RELTOL=$RT B6_MISSION=60 timeout $CAP julia --project=. -t 1 "$SCRATCH/b6_tolerance.jl" 2>&1 \
      | grep -E "reltol=" | tail -1)
  if [ -z "$R" ]; then R="  reltol=$RT mission=60s     DID NOT COMPLETE within ${CAP}s"; fi
  echo "$R" >> "$SCRATCH/sweep.out"
done
echo "DONE" >> "$SCRATCH/sweep.out"
