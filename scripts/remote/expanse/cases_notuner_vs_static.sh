#!/usr/bin/env bash
# Does the confirmed, tuner-off calibration pay at 32 threads? Single-simulation vacuum constellation cases
# (mc_samples=1, one 32-thread process each). Modes: inner_only (identical in configuration to the paper's
# outer_threads / outer_inner_static on these cases), predictive (default), predictive_notuner_confirm (region-width
# tuner off, SWEEP_NS_FLOOR=1, A/B CONFIRM=1 with K=2 solves per arm; set by modes.jl extra_env).
# Runs INSIDE a Slurm job (julia_job.slurm does the quiet-node gate):
#
#   EXPANSE_ACCOUNT=mia402 ./expanse.sh run --bench --time 01:15:00 -- bash scripts/remote/expanse/cases_notuner_vs_static.sh
#   DRY=1 bash scripts/remote/expanse/cases_notuner_vs_static.sh    # prints the commands, runs no Julia
#   CASES=a,b MODES=m1,m2 REPEATS=6 to narrow (comma lists; spaces would not survive expanse.sh).
#
# Startup cost: all calls run in ONE controller process (--points-file), so the package loads once, not once per call.
# Isolation: ONE harness call per (case, mode), SPACEAGORA_PPB_SHARED_PROCESS=0, so every mode runs in its own fresh
# Julia process (harness child) with its own FRESH RHS calibration store, snapshot dir and emptied policy-state dir
# (cost constants put back). No process, store or compiled-code sharing between modes.
# Order: case-major. Cases 1-3 rotate the mode list by one position each (a full Latin square of 3 modes); case 4,
# which a 3x3 square cannot cover, runs the reversed list. So within-case position counts are 3 modes x 4 cases:
# each mode is 1st/2nd/3rd about 4/3 times; the residual imbalance is one extra count per mode (inner_only 1st 1x,
# 2nd 1x, 3rd 2x; predictive 1/2/1; confirm 2/1/1). Position effects stay confounded with case.
# predictive_notuner_confirm: the sweep runs in the warm-up; timed repeats 1-2 run the pin, 3-4 the heuristic (K=2
# per arm), the verdict follows repeat 4, so repeats 5-8 are the post-verdict state; analyze those separately.
# The last 4 calls are the extra inner_only noise-floor processes (NOISE_FLOOR=0 drops them).
# Output: results/expanse_policy/<jobid>_notuner_vs_static/calls/<case>__<mode>{/,.log,_snap/}. Raw CSV per call has
# wall_time_s and final_primary_pos_norm_m / vel_norm_mps / mass_kg / terminal_time_s (compare across modes for
# bit-identical trajectories); _snap/ is the per-solve snapshot (wall time, calibration store, plan source, confirm state). Data only, no plots.
set -u
cd "$(dirname "$0")/../../.." || exit 1
DRY="${DRY:-0}"
JOB="${SLURM_JOB_ID:-local}"
OUT="$PWD/results/expanse_policy/${JOB}_notuner_vs_static"
BUDGET="${BUDGET:-32}"
CASES="${CASES:-gravity_4096sat_l50_vacuum_5800s gravity_1024sat_l50_vacuum_24600s gravity_4096sat_l20_vacuum_19700s gravity_64sat_l50_vacuum_415000s}"
CASES="${CASES//,/ }"
MODES="${MODES:-inner_only,predictive,predictive_notuner_confirm}"
BASE=(${MODES//,/ })
REPEATS="${REPEATS:-8}"; WARMUP="${WARMUP:-1}"
STORE=output/parallel_policy_state
CONSTANTS="$OUT/constants"

SEQ=()   # "<case>:<mode>:<tag suffix>", case-major, modes rotated by one per case
ci=0
for cs in $CASES; do
  n=${#BASE[@]}
  for ((k = 0; k < n; k++)); do
    if ((ci < n)); then j=$(((k + ci) % n)); else j=$(((n - 1 - k + n) % n)); fi   # ponytail: reversed list for case n+1 only
    SEQ+=("${cs}:${BASE[j]}:")
  done
  ci=$((ci + 1))
done
# Noise floor: a second inner_only process per case, run after all the above, cases rotated by two positions
# (tag suffix _b2, so outputs never collide with the first inner_only call). Each case thus has 2 inner_only processes.
CA=($CASES); nc=${#CA[@]}
if [ "${NOISE_FLOOR:-1}" = 1 ]; then
  for ((k = 0; k < nc; k++)); do SEQ+=("${CA[(k + 2) % nc]}:inner_only:_b2"); done
fi
say() { echo "[notuner-vs-static $(date -Iseconds)] $*"; }
run() { say "+ $*"; [ "$DRY" = 1 ] && return 0; "$@"; }

export JULIA_NUM_THREADS="" SPACEAGORA_PERF_PROCS="" SPACEAGORA_PPB_WORKLOAD=0
# Anything that would leak into a mode run (modes set their own SPACEAGORA_RHS_CALIBRATE / REGION_WIDTH via extra_env).
unset JULIA_NUM_GC_THREADS JULIA_THREAD_SLEEP_THRESHOLD SPACEAGORA_RHS_CALIBRATE SPACEAGORA_RHS_CALIBRATE_SEED_TOPK \
      SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER SPACEAGORA_POLICY_REGRET_GUARD \
      SPACEAGORA_POLICY_REGRET_MARGIN SPACEAGORA_POLICY_REGRET_MIN_SOLVES SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR \
      SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN SPACEAGORA_RHS_REGION_WIDTH SPACEAGORA_RHS_CALIBRATE_CONFIRM \
      SPACEAGORA_RHS_CALIBRATE_CONFIRM_SOLVES SPACEAGORA_RHS_CALIBRATE_CONFIRM_MARGIN
export SPACEAGORA_CORE_BUDGET=$BUDGET SPACEAGORA_PPC_PHYSICAL_CORES=$BUDGET SPACEAGORA_PPC_EQUAL_CORE_BUDGET=1
export SPACEAGORA_PPB_SHARED_PROCESS=0 SPACEAGORA_PPB_WARMUP_SPAN=0.1 SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto
export SPACEAGORA_CAMPAIGN_CORRECTIONS=off SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1 SPACEAGORA_PPB_NO_PLOTS=1
export OPENBLAS_NUM_THREADS=1

say "job=$JOB host=$(hostname) commit=${GIT_COMMIT:-unknown} DRY=$DRY repeats=$REPEATS warmup=$WARMUP budget=$BUDGET"
say "${#SEQ[@]} processes, in order:"
n=0; for item in "${SEQ[@]}"; do n=$((n + 1)); r=${item#*:}; say "  $n. ${item%%:*} ${r%%:*}${r#*:}"; done

# Not timed: bring both precompile caches up to date for this source tree. Normally already done on the `shared`
# partition (precompile_tree.sh, run via `expanse.sh run` without --bench), so each line here is only a cache check
# (~10 s). If it prints `Precompiling`, the exclusive node is paying for it: say so, and the run is still valid
# (every timed child starts after this). Machine cost constants (usl_alpha_base, usl_beta_alloc): measured once at the
# budget's thread count, put back into the emptied policy store before every call.
[ "$DRY" = 1 ] || mkdir -p "$OUT"
for proj in . examples; do
  pl="$OUT/prelude_${proj//./root}.log"
  if [ "$DRY" = 1 ]; then run julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA'; continue; fi
  say "+ julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA'"
  julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA' 2>&1 | tee "$pl"
  [ "${PIPESTATUS[0]}" -eq 0 ] || { say "ERROR: prelude --project=$proj failed"; exit 6; }
  grep -q Precompiling "$pl" && say "WARNING: --project=$proj precompiled on the exclusive node; run precompile_tree.sh on the shared partition first next time"
done
if [ "${SKIP_CALIBRATE:-0}" != 1 ]; then
  run mkdir -p "$CONSTANTS"
  run julia --project=. --threads=$BUDGET scripts/calibration/calibrate_machine.jl
  if [ "$DRY" != 1 ]; then
    cp -p "$STORE"/cost_constants_*.toml "$CONSTANTS"/ || { say "ERROR: calibration wrote no constants"; exit 5; }
  fi
fi


# ONE controller process runs all calls (no package load per call). Points file: tag, case, mode, outdir, store, snapshot dir.
# Per call the controller sets SPACEAGORA_RHS_CALIBRATION_PATH / SPACEAGORA_PPB_SNAPSHOT_DIR, runs the prep command
# (fresh policy-state dir and snapshot dir, cost constants put back, no store), then runs the same harness call as before
# (--cases=<case> --modes=<mode> --outdir=<outdir>) with a fresh child, logging to <outdir>.log. A call fails on a failed
# child or a `Precompiling` line in its log.
POINTS="$OUT/points.tsv"
export SPACEAGORA_PPC_POINT_PREP="rm -rf '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; mkdir -p '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; cp -p '$CONSTANTS'/cost_constants_*.toml '$STORE'/ 2>/dev/null; rm -f \"\$SPACEAGORA_RHS_CALIBRATION_PATH\""
run mkdir -p "$OUT/stores" "$OUT/calls"
pts=""
for item in "${SEQ[@]}"; do
  CASE=${item%%:*}; rest=${item#*:}; mode=${rest%%:*}; tag="${CASE}__${mode}${rest#*:}"
  pts+="$tag\t$CASE\t$mode\t$OUT/calls/$tag\t$OUT/stores/$tag.toml\t$OUT/calls/${tag}_snap\n"
done
if [ "$DRY" = 1 ]; then
  say "points file $POINTS:"; printf "$pts" | sed 's/^/  /'
  say "prep: $SPACEAGORA_PPC_POINT_PREP"
else
  printf "$pts" > "$POINTS"
fi
cmd=(julia --project=examples benchmarks/studies/parallelization_performance.jl full
     --points-file="$POINTS" --threads=$BUDGET --process-workers=$BUDGET
     --repeats=$REPEATS --warmup=$WARMUP --parity-cases=none)
if [ "$DRY" = 1 ]; then echo "  ${cmd[*]}"; say "done (dry run): $OUT"; exit 0; fi
"${cmd[@]}"; rc=$?
[ $rc -eq 0 ] || { say "FAILED (rc=$rc); see [points] lines above and $OUT/calls/*.log"; exit 1; }
say "done: $OUT"
