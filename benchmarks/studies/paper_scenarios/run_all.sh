#!/bin/bash
# Full paper_scenarios protocol. Front-loads the short scenarios so wiring problems
# surface before the long GRAM/Monte Carlo sweeps. Set PS_THREADS to the machine's
# physical core count before a paper-quality run.
set -e
cd "$(dirname "$0")/../../.."

echo "== S1: constellation scaling (vacuum baseline) =="
julia --project=. benchmarks/studies/paper_scenarios/s1_constellation_scaling.jl

echo "== S5: routing profile ladder =="
julia --project=. benchmarks/studies/paper_scenarios/s5_routing_profile_ladder.jl

echo "== S4: hybrid MC x constellation =="
julia --project=. benchmarks/studies/paper_scenarios/s4_hybrid_mc_constellation.jl

echo "== S2: GRAM atmosphere modes =="
julia --project=. benchmarks/studies/paper_scenarios/s2_gram_atmosphere_modes.jl

echo "== S3: Monte Carlo process scaling =="
julia --project=. benchmarks/studies/paper_scenarios/s3_montecarlo_process_scaling.jl

echo "All scenarios complete. Results: benchmarks/studies/paper_scenarios/results/$(hostname)/"
