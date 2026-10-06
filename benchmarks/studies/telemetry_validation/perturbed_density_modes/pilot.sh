#!/usr/bin/env bash
# Pilot for the perturbed-density mode comparison (run from a release root).
#  1. instance-isolation probe (native GRAM);
#  2. default-mode bit-identity: one orbit with mode=off, compared afterwards
#     against the smoke test's pre-change one-orbit record run;
#  3. 5-orbit nominal / A / B / naive runs, concurrently, to check the modes run
#     and to measure per-orbit cost before the arc length is chosen.
set -u
OUT=results/perturbed_density_modes_pilot
D=benchmarks/studies/telemetry_validation/perturbed_density_modes/run_mode.jl
mkdir -p "$OUT"
julia --project=. test/probes/gram_perturbation_walk_probes.jl > "$OUT/probe.log" 2>&1
echo "PROBE_EXIT=$?" >> "$OUT/status"
run() { local tag=$1; shift; ( "$@" > "$OUT/$tag.log" 2>&1; echo "$tag EXIT=$?" >> "$OUT/status" ) & }
run off_1orbit julia --project=. $D --tag=off_1orbit --mode=off --orbits=1 --out=$OUT
run nominal_5 julia --project=. $D --tag=nominal_5 --mode=off --orbits=5 --out=$OUT
run A_s11_5 julia --project=. $D --tag=A_s11_5 --mode=step --seed=11 --orbits=5 --out=$OUT
run B_s11_5 julia --project=. $D --tag=B_s11_5 --mode=pass --seed=11 --orbits=5 --out=$OUT
run naive_s11_5 timeout 1200 julia --project=. $D --tag=naive_s11_5 --mode=naive_rhs --seed=11 --orbits=5 --maxiters=2000000 --out=$OUT
wait
echo "PILOT_DONE $(date -Iseconds)" >> "$OUT/status"
