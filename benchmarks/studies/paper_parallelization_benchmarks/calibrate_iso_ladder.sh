#!/usr/bin/env bash
# Recalibrate the iso-work L50 ladder (PPC_L50_ISO_MISSION_S in cases.jl) for the
# host it runs on.
#
# P1 holds the total work fixed rather than the mission length, so that every
# rung's serial baseline clears the harness's 3 s measurability floor. Those
# durations are per machine: a faster or slower core moves every rung, and a
# ladder carried over from another host either drops rungs under the floor or
# spends minutes per point measuring nothing extra.
#
# This runs each registered rung once in serial mode and prints the duration
# that would put it at the target, as a table you can paste back into cases.jl.
# It does NOT edit cases.jl -- eyeball the numbers first; a rung whose measured
# and suggested durations differ wildly usually means the machine was busy.
#
#   bash benchmarks/studies/paper_parallelization_benchmarks/calibrate_iso_ladder.sh [target_s]
#
# Takes about 10 minutes: one Julia start-up (~80 s) plus a warm-up and a timed
# solve per rung. Run it on an idle machine -- the harness's own quiet guard
# applies to the controller, not to this.
set -uo pipefail

TARGET_S="${1:-10.0}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT="${TMPDIR:-/tmp}/iso_ladder_calibration_$$"
mkdir -p "$OUT"

# The rungs currently registered, read from cases.jl so the two cannot drift.
mapfile -t RUNGS < <(
  sed -n '/^const PPC_L50_ISO_MISSION_S/,/^]/p' \
    "$REPO_ROOT/benchmarks/studies/parallelization_performance/cases.jl" |
    sed -n 's/^ *(\([0-9]*\), *\([0-9]*\)).*/\1 \2/p'
)

if [ "${#RUNGS[@]}" -eq 0 ]; then
  echo "error: no rungs parsed from PPC_L50_ISO_MISSION_S in cases.jl" >&2
  exit 1
fi

echo "Calibrating ${#RUNGS[@]} rungs against a ${TARGET_S}s serial target on $(hostname)."
echo

printf '%-8s %-14s %-12s %-14s\n' "N" "mission (s)" "serial (s)" "suggested (s)"
for rung in "${RUNGS[@]}"; do
  n="${rung% *}"
  mission="${rung#* }"
  case_name="gravity_${n}sat_l50_vacuum_${mission}s"
  julia --threads=1 --project="$REPO_ROOT" \
    "$REPO_ROOT/benchmarks/studies/parallelization_performance.jl" full --worker \
    --case="$case_name" --mode=serial --thread-count=1 --repeat=1 --worker-repeats=1 \
    --worker-seed=20260616 --worker-mc-samples=1 --warmup=1 --solver-mode=auto_stiff \
    --process-workers=1 --parity-samples=512 --outfile="$OUT/$case_name.csv" --parity=0 \
    > "$OUT/$case_name.log" 2>&1

  wall=$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1{print $h["wall_time_s"]; exit}' \
    "$OUT/$case_name.csv" 2>/dev/null)
  if [ -z "$wall" ]; then
    printf '%-8s %-14s %-12s %-14s  (see %s)\n' "$n" "$mission" "FAILED" "--" "$OUT/$case_name.log"
    continue
  fi
  # Scaling is linear in mission length to within a few percent at these sizes,
  # but only above the floor: a rung measuring well under 3 s is extrapolating
  # from mostly fixed per-solve cost and will land short. Re-run those.
  suggested=$(awk -v w="$wall" -v m="$mission" -v t="$TARGET_S" \
    'BEGIN{printf "%d", (m * t / w) + 0.5}')
  printf '%-8s %-14s %-12.2f %-14s\n' "$n" "$mission" "$wall" "$suggested"
done

echo
echo "Raw rows and logs: $OUT"
echo "Paste the suggested column into PPC_L50_ISO_MISSION_S if it moves a rung"
echo "more than ~20%, then re-run this script once to confirm."
