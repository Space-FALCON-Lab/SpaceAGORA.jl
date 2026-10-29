#!/usr/bin/env bash
# Does RHS calibration pay at 32 threads when the region-width tuner is off? Single-simulation vacuum constellation
# cases (mc_samples=1, one 32-thread process each). Modes: inner_only (identical in configuration to the paper's
# outer_threads / outer_inner_static on these cases), predictive, predictive_notuner, predictive_nocalib_notuner.
# Runs INSIDE a Slurm job (julia_job.slurm does the quiet-node gate):
#
#   EXPANSE_ACCOUNT=mia402 ./expanse.sh run --bench --time 01:00:00 -- bash scripts/remote/expanse/cases_notuner_vs_static.sh
#   DRY=1 bash scripts/remote/expanse/cases_notuner_vs_static.sh    # prints the commands, runs no Julia
#   CASES=a,b MODES=m1,m2 REPEATS=6 to narrow (comma lists; spaces would not survive expanse.sh).
#
# Isolation: ONE harness call per (case, mode), SPACEAGORA_PPB_SHARED_PROCESS=0, so every mode runs in its own fresh
# Julia process (harness child) with its own FRESH RHS calibration store, snapshot dir and emptied policy-state dir
# (cost constants put back). No process, store or compiled-code sharing between modes.
# Order: a Latin square. The mode list is rotated by one position per case, so each mode runs in each position
# (1st..4th in the job's case-major order) exactly once across the four cases; position effects are balanced across
# modes but are confounded with case, and are NOT balanced within a case. The cases run in the order given.
# Output: results/expanse_policy/<jobid>_notuner_vs_static/calls/<case>__<mode>{/,.log,_snap/}. Raw CSV per call has
# wall_time_s and final_primary_pos_norm_m / vel_norm_mps / mass_kg / terminal_time_s (compare across modes for
# bit-identical trajectories); _snap/ is the per-solve snapshot (wall time, calibration store). Data only, no plots.
set -u
cd "$(dirname "$0")/../../.." || exit 1
DRY="${DRY:-0}"
JOB="${SLURM_JOB_ID:-local}"
OUT="$PWD/results/expanse_policy/${JOB}_notuner_vs_static"
BUDGET="${BUDGET:-32}"
CASES="${CASES:-gravity_4096sat_l50_vacuum_5800s gravity_1024sat_l50_vacuum_24600s gravity_4096sat_l20_vacuum_19700s gravity_64sat_l50_vacuum_415000s}"
CASES="${CASES//,/ }"
MODES="${MODES:-inner_only,predictive,predictive_notuner,predictive_nocalib_notuner}"
BASE=(${MODES//,/ })
REPEATS="${REPEATS:-6}"; WARMUP="${WARMUP:-1}"
STORE=output/parallel_policy_state
CONSTANTS="$OUT/constants"

SEQ=()   # "<case>:<mode>", case-major, modes rotated by one per case
ci=0
for cs in $CASES; do
  for ((k = 0; k < ${#BASE[@]}; k++)); do SEQ+=("${cs}:${BASE[(k + ci) % ${#BASE[@]}]}"); done
  ci=$((ci + 1))
done
say() { echo "[notuner-vs-static $(date -Iseconds)] $*"; }
run() { say "+ $*"; [ "$DRY" = 1 ] && return 0; "$@"; }

export JULIA_NUM_THREADS="" SPACEAGORA_PERF_PROCS="" SPACEAGORA_PPB_WORKLOAD=0
# Anything that would leak into a mode run (modes set their own SPACEAGORA_RHS_CALIBRATE / REGION_WIDTH via extra_env).
unset JULIA_NUM_GC_THREADS JULIA_THREAD_SLEEP_THRESHOLD SPACEAGORA_RHS_CALIBRATE SPACEAGORA_RHS_CALIBRATE_SEED_TOPK \
      SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER SPACEAGORA_POLICY_REGRET_GUARD \
      SPACEAGORA_POLICY_REGRET_MARGIN SPACEAGORA_POLICY_REGRET_MIN_SOLVES SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR \
      SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN SPACEAGORA_RHS_REGION_WIDTH
export SPACEAGORA_CORE_BUDGET=$BUDGET SPACEAGORA_PPC_PHYSICAL_CORES=$BUDGET SPACEAGORA_PPC_EQUAL_CORE_BUDGET=1
export SPACEAGORA_PPB_SHARED_PROCESS=0 SPACEAGORA_PPB_WARMUP_SPAN=0.1 SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto
export SPACEAGORA_CAMPAIGN_CORRECTIONS=off SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1 SPACEAGORA_PPB_NO_PLOTS=1
export OPENBLAS_NUM_THREADS=1

say "job=$JOB host=$(hostname) commit=${GIT_COMMIT:-unknown} DRY=$DRY repeats=$REPEATS warmup=$WARMUP budget=$BUDGET"
say "sequence: ${SEQ[*]}"

# Not timed: bring both precompile caches up to date for this source tree (the previous job's tree differed), so the
# first mode run does not pay for it. Machine cost constants (usl_alpha_base, usl_beta_alloc): measured once at the
# budget's thread count, put back into the emptied policy store before every call.
run julia --project=. --threads=$BUDGET -e 'using SpaceAGORA'
run julia --project=examples --threads=$BUDGET -e 'using SpaceAGORA'
if [ "${SKIP_CALIBRATE:-0}" != 1 ]; then
  run mkdir -p "$CONSTANTS"
  run julia --project=. --threads=$BUDGET scripts/calibration/calibrate_machine.jl
  if [ "$DRY" != 1 ]; then
    cp -p "$STORE"/cost_constants_*.toml "$CONSTANTS"/ || { say "ERROR: calibration wrote no constants"; exit 5; }
  fi
fi

failed=()
for item in "${SEQ[@]}"; do
  CASE=${item%%:*}; mode=${item#*:}; tag="${CASE}__${mode}"
  store="$OUT/stores/$tag.toml"; snap="$OUT/calls/${tag}_snap"
  say "START $tag load=$(cut -d' ' -f1-3 /proc/loadavg)"
  run mkdir -p "$OUT/stores" "$OUT/calls"
  if [ "$DRY" != 1 ]; then
    rm -rf "$STORE" "$snap"; mkdir -p "$STORE" "$snap"; cp -p "$CONSTANTS"/cost_constants_*.toml "$STORE"/ 2>/dev/null
    rm -f "$store"
  fi
  cmd=(julia --project=examples benchmarks/studies/parallelization_performance.jl full
       --outdir="$OUT/calls/$tag" --cases=$CASE --modes=$mode --threads=$BUDGET --process-workers=$BUDGET
       --repeats=$REPEATS --warmup=$WARMUP --parity-cases=none)
  if [ "$DRY" = 1 ]; then
    echo "  env SPACEAGORA_RHS_CALIBRATION_PATH=$store SPACEAGORA_PPB_SNAPSHOT_DIR=$snap ${cmd[*]}"; continue
  fi
  env SPACEAGORA_RHS_CALIBRATION_PATH="$store" SPACEAGORA_PPB_SNAPSHOT_DIR="$snap" "${cmd[@]}" > "$OUT/calls/$tag.log" 2>&1
  rc=$?
  # the harness can exit 0 with dead workers; a failed child or point is a failed call
  if [ $rc -eq 0 ] && grep -qE '\[child-failed\]|status=error' "$OUT/calls/$tag.log"; then rc=1; fi
  say "END $tag rc=$rc load=$(cut -d' ' -f1-3 /proc/loadavg)"
  [ $rc -eq 0 ] || failed+=("$tag")
done
if [ ${#failed[@]} -gt 0 ]; then say "FAILED: ${failed[*]}"; exit 1; fi
say "done: $OUT"
