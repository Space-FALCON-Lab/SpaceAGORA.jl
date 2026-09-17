#!/usr/bin/env bash
# Pre-solve sweep vs in-run width trial on the two cold-store single-simulation
# points the TRX50 lost (L10 gram_surrogate 1.28, L14 heavy_1024_6hr 1.25):
# best static route / R6 with the sweep (cold store) / R6 with the sweep off and
# the in-run trial on (cold store). Usage: sweep_vs_trial_probe.sh <threads>
set -u
T=${1:-12}
cd "$(dirname "$0")/../.."
OUT=output/probe_sweep_vs_trial; mkdir -p $OUT
run() { local name=$1 mode=$2 case=$3 extra_env=$4
  local d=$OUT/$name; mkdir -p "$d"
  ( export SPACEAGORA_RHS_CALIBRATION_PATH="$d/store.toml" SPACEAGORA_PARALLEL_POLICY_STATE_PATH="$d/hints.toml" SPACEAGORA_OUTER_ROUTE_STATE_PATH="$d/route.toml"
    [ -n "$extra_env" ] && export $extra_env
    echo "=== $name mode=$mode case=$case env=$extra_env start $(date +%T)"
    julia --threads=$T --project=. benchmarks/studies/parallelization_performance.jl full --worker \
      --case=$case --mode=$mode --thread-count=$T --repeat=1 --worker-repeats=4 --worker-seed=20260615 \
      --worker-mc-samples=1 --warmup=1 --solver-mode=auto_stiff --process-workers=1 --parity-samples=0 \
      --outfile="$d/row.csv" --parity=0 > "$d/stdout.log" 2>&1; echo "  exit=$?  end $(date +%T)" )
}
for case in atmo256_gram_surrogate_10min heavy_1024sat_l50_6hr; do
  run ${case}_static outer_inner_static $case ""
  run ${case}_sweep  policy_v2          $case ""
  run ${case}_trial  policy_v2_nocalib  $case "SPACEAGORA_RHS_IDENTIFY=1"
done
echo "=== ALL DONE $(date +%T)"
