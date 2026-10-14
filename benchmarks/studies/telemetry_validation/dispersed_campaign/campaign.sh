#!/usr/bin/env bash
# Dispersed Odyssey MarsGRAM campaign job (run from a release root):
# 1) the 32-member campaign under the adaptive policy, full budget;
# 2) then, alone, the nominal (perturbation off) member.
# ORBITS defaults to 42 = N + 2 for the 40-pass analysis (analyze_campaign.py).
set -u
OUT=results/${CAMPAIGN_TAG:-odyssey_dispersed_campaign}
DRV=benchmarks/studies/telemetry_validation/dispersed_campaign/run_campaign.jl
mkdir -p "$OUT"
echo "CAMPAIGN_START $(date -Iseconds)" >> "$OUT/status"
julia --project=. --threads=${JULIA_NUM_THREADS:-32} $DRV --out=$OUT --orbits=${ORBITS:-42} --seeds=${SEEDS:-101:132} > "$OUT/campaign.log" 2>&1
echo "CAMPAIGN_EXIT=$? $(date -Iseconds)" >> "$OUT/status"
julia --project=. $DRV --out=$OUT --orbits=${ORBITS:-42} --member=nominal > "$OUT/nominal.log" 2>&1
echo "NOMINAL_EXIT=$? $(date -Iseconds)" >> "$OUT/status"
