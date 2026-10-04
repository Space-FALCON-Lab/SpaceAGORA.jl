#!/usr/bin/env bash
# Follow-up to compare.sh: the same arc with the walk reseeded at every
# atmospheric entry (SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_RESEED=1), to test
# whether per-pass seeds stop a numerical change from re-randomizing every later
# pass. The nominal baseline is compare.sh's (unchanged code path).
set -u
ORBITS=${ORBITS:-120}
OUT=results/perturbed_density_modes_reseed
D=benchmarks/studies/telemetry_validation/perturbed_density_modes/run_mode.jl
mkdir -p "$OUT"
run() { local tag=$1; shift; ( julia --project=. $D --tag=$tag --orbits=$ORBITS --out=$OUT --reseed=1 "$@" > "$OUT/$tag.log" 2>&1; echo "$tag EXIT=$? $(date -Iseconds)" >> "$OUT/status" ) & }
echo "START $(date -Iseconds) orbits=$ORBITS" >> "$OUT/status"
run nominal --mode=off
for s in 11 22 33; do run Brs_s$s --mode=pass --seed=$s; done
run Brs_s11_rep   --mode=pass --seed=11
run Brs_s11_tight --mode=pass --seed=11 --tight=true
run Brs_s11_dt01  --mode=pass --seed=11 --dt-max-atm=0.1
run Ars_s11       --mode=step --seed=11
run Ars_s11_tight --mode=step --seed=11 --tight=true
run Ars_s11_dt01  --mode=step --seed=11 --dt-max-atm=0.1
wait
echo "DONE $(date -Iseconds)" >> "$OUT/status"
