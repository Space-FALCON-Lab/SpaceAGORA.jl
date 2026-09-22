#!/usr/bin/env bash
# Targeted single-point verification of the policy_v2 consistency fixes: the
# specific points that failed the acceptance criterion, each run as policy_v2
# and as predictive against its best pinned static route, 11 repeats, 32
# threads. Meant to run inside a release on the benchmark box (via
# scripts/remote/spaceagora-remote), one POINT per job, with the calibration
# store set by the launcher beforehand (cold = emptied; converged = the
# archived snapshot restored). Gates a full P1-P5 run: no point spending eleven
# hours to learn what one hour would.
#
# Calibration for the predictive arm. The predictive campaign planner (profile
# R7, SPACEAGORA_CAMPAIGN_PLANNER=predictive) reads the per-machine cost
# constants from
#
#     output/parallel_policy_state/cost_constants_<fingerprint>.toml
#
# written by
#
#     julia --project=. --threads=<T> scripts/calibrate_machine.jl
#
# run at the thread count the jobs will use. The launcher runs that ONCE on the
# box before the first job; this script does not run it, because calibration is
# itself a timed measurement and must not share the machine with a benchmark.
# Without the file the planner still runs, but it models no contention, so the
# predictive rows measure something other than what the launcher intended. The
# check below only warns.
set -euo pipefail
export SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1
if ! ls output/parallel_policy_state/cost_constants_*.toml >/dev/null 2>&1; then
  echo "[targeted] WARNING: no output/parallel_policy_state/cost_constants_*.toml -- the" \
       "predictive rows will run without machine calibration. Run" \
       "'julia --project=. --threads=<T> scripts/calibrate_machine.jl' once on this box first."
fi
OUT=output/performance/targeted_$(date -u +%Y%m%d_%H%M%S); mkdir -p "$OUT"
run() { # case mode workers threads mc label
  local c=$1 m=$2 w=$3 t=$4 mc=$5 lbl=$6
  echo "[targeted] $lbl: case=$c mode=$m workers=$w threads=$t mc=$mc"
  timeout 7200 julia --threads="$t" --project=. benchmarks/studies/parallelization_performance.jl \
    --profile=full --worker --case="$c" --mode="$m" --thread-count="$t" \
    --worker-repeats="${REPEATS:-11}" --worker-mc-samples="$mc" --warmup="${WARMUP:-1}" --process-workers="$w" --parity=0 \
    --outfile="$OUT/${lbl}_${m}.csv" > "$OUT/${lbl}_${m}.log" 2>&1
  tail -3 "$OUT/${lbl}_${m}.log"
}
case "${POINT:?set POINT}" in
  finding8)   # converged store: does R6 re-sweep off cache/heuristic at N=1024?
    run gravity_1024sat_l50_vacuum_24600s policy_v2          32 32 1  f8
    run gravity_1024sat_l50_vacuum_24600s predictive         32 32 1  f8
    run gravity_1024sat_l50_vacuum_24600s outer_inner_static 32 32 1  f8 ;;
  p5_16sat)   # cold store: WS3's inner-split gate at the 1x32 split
    run mcgrid_16sat_8mc policy_v2 1 32 8 p5_16
    run mcgrid_16sat_8mc predictive 1 32 8 p5_16
    run mcgrid_16sat_8mc outer_threads 1 32 8 p5_16 ;;
  finding9)   # cold store: the throw -- error_message + max_sample now recorded
    run mcgrid_8sat_16mc policy_v2 1 32 16 f9
    run mcgrid_8sat_16mc predictive 1 32 16 f9
    run mcgrid_8sat_16mc outer_threads 1 32 16 f9 ;;
  defectA)    # baseline only (no fix yet): mixed dispatch at W=32
    run montecarlo_heavy_aerobraking policy_v2     32 32 32 dA
    run montecarlo_heavy_aerobraking predictive    32 32 32 dA
    run montecarlo_heavy_aerobraking outer_process 32 32 32 dA ;;
  defectA_gc) # cold store: is R7's extra worker-side GC at P4@32 the collection mode?
    SPACEAGORA_POOL_WORKER_GC=full        run montecarlo_heavy_aerobraking predictive    32 32 32 dAgcfull
    SPACEAGORA_POOL_WORKER_GC=full        run montecarlo_heavy_aerobraking outer_process 32 32 32 dAgcfull
    SPACEAGORA_POOL_WORKER_GC=incremental run montecarlo_heavy_aerobraking predictive    32 32 32 dAgcincr
    SPACEAGORA_POOL_WORKER_GC=incremental run montecarlo_heavy_aerobraking outer_process 32 32 32 dAgcincr
    SPACEAGORA_POOL_WORKER_GC=off         run montecarlo_heavy_aerobraking predictive    32 32 32 dAgcoff
    SPACEAGORA_POOL_WORKER_GC=off         run montecarlo_heavy_aerobraking outer_process 32 32 32 dAgcoff ;;
  defectA_33) # cold store, 33 repeats: is R7's higher mid-campaign worker-GC incidence real?
    REPEATS=33 run montecarlo_heavy_aerobraking predictive    32 32 32 dA33
    REPEATS=33 run montecarlo_heavy_aerobraking outer_process 32 32 32 dA33 ;;
  defectA_w3) # cold store, 11 repeats after 3 warm-ups: is the P4@32 gap the heap burn-in?
    WARMUP=3 run montecarlo_heavy_aerobraking predictive    32 32 32 dAw3
    WARMUP=3 run montecarlo_heavy_aerobraking outer_process 32 32 32 dAw3 ;;
  *) echo "unknown POINT=$POINT"; exit 2 ;;
esac
echo "[targeted] done -> $OUT"; ls "$OUT"
