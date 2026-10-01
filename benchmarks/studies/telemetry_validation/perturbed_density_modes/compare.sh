#!/usr/bin/env bash
# Perturbed-density mode comparison on the shortened Odyssey MarsGRAM arc
# (run from a release root). Every run is single-threaded and they run side by
# side in one job, so wall times are comparable with each other (same machine
# load), not with an idle-machine solve.
#
#   nominal        mean density (the record setting) -- the baseline
#   A_s{11,22,33}  design A (per-accepted-step hold), density scale 1.0
#   B_s{11,22,33}  design B (per-pass look-ahead, 1 s knots), density scale 1.0
#   *_tight        10x tighter reltol/abstol (orbit and atmosphere)
#   *_dt01         in-atmosphere step cap 0.1 s instead of 0.2 s
#   *_rep          identical rerun, for bit-identity
set -u
ORBITS=${ORBITS:-120}
OUT=results/perturbed_density_modes
D=benchmarks/studies/telemetry_validation/perturbed_density_modes/run_mode.jl
mkdir -p "$OUT"
run() { local tag=$1; shift; ( julia --project=. $D --tag=$tag --orbits=$ORBITS --out=$OUT "$@" > "$OUT/$tag.log" 2>&1; echo "$tag EXIT=$? $(date -Iseconds)" >> "$OUT/status" ) & }
echo "START $(date -Iseconds) orbits=$ORBITS" >> "$OUT/status"
run nominal        --mode=off
run nominal_rep    --mode=off
run nominal_tight  --mode=off --tight=true
run nominal_dt01   --mode=off --dt-max-atm=0.1
for s in 11 22 33; do
  run A_s$s --mode=step --seed=$s
  run B_s$s --mode=pass --seed=$s
done
run A_s11_rep   --mode=step --seed=11
run B_s11_rep   --mode=pass --seed=11
run A_s11_tight --mode=step --seed=11 --tight=true
run B_s11_tight --mode=pass --seed=11 --tight=true
run A_s11_dt01  --mode=step --seed=11 --dt-max-atm=0.1
run B_s11_dt01  --mode=pass --seed=11 --dt-max-atm=0.1
wait
echo "DONE $(date -Iseconds)" >> "$OUT/status"
