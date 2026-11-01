#!/usr/bin/env bash
# Final-stack pilot: inner_only (static, the baseline and its own noise floor) vs predictive (the candidate stack: R7,
# region-width tuner off, sweep_ns floor on, plan confirmation on, no multi-sample mixed plans; see modes.jl) at 32 threads
# on four single-simulation vacuum constellation cases (mc_samples=1). 3 fresh processes per (case, mode) = 24 points.
# Runs INSIDE a Slurm job (julia_job.slurm does the quiet-node gate):
#
#   EXPANSE_ACCOUNT=mia402 ./expanse.sh run --bench --time 01:30:00 -- bash scripts/remote/expanse/pilot_final_vs_static.sh
#   DRY=1 bash scripts/remote/expanse/pilot_final_vs_static.sh      # prints the points and commands, runs no Julia
#   CASES=a,b MODES=m1,m2 ROUNDS=3 REPEATS=8 to narrow (comma lists; spaces would not survive expanse.sh).
#
# Startup cost: all calls run in ONE controller process (--points-file), so the package loads once, not once per call.
# Isolation: ONE harness call per point, SPACEAGORA_PPB_SHARED_PROCESS=0, so every point runs in its own fresh Julia
# process (harness child) with its own FRESH RHS calibration store, snapshot dir and emptied policy-state dir
# (cost constants put back). No process, store or compiled-code sharing between points.
# Order: ROUNDS rounds; each round holds every (case, mode) pair once. Pairs are indexed i = 0..P-1 in list order (cases
# in order, modes within case). Round r runs pair i at position (i + 3r) mod P (a cyclic shift by 3 pairs per round), so
# each pair visits three different positions and drift/position effects average out; inner_only is replicated 3x per case,
# so its process-to-process noise is measured in the same job. For the default P=8, R=3 the mean position of a pair is
# 3.3..5.7 (1-based, grand mean 4.5); by mode: inner_only 4.33, predictive 4.67 (DRY=1 prints them). Position balance
# per mode is approximate (3 rounds cannot balance parity of an 8-slot round exactly).
# Output: results/expanse_policy/<jobid>_pilot_final/calls/<case>__<mode>_r<round>{/,.log,_snap/}. Raw CSV per call has
# wall_time_s and final_primary_pos_norm_m / vel_norm_mps / mass_kg / terminal_time_s (compare across modes for
# bit-identical trajectories); _snap/ is the per-solve snapshot (wall time, calibration store, plan source, confirm state).
# Data only, no plots, no in-job summary (analysis happens locally on the pulled CSVs).
set -u
cd "$(dirname "$0")/../../.." || exit 1
DRY="${DRY:-0}"
JOB="${SLURM_JOB_ID:-local}"
OUT="$PWD/results/expanse_policy/${JOB}_pilot_final"
BUDGET="${BUDGET:-32}"
CASES="${CASES:-gravity_4096sat_l50_vacuum_5800s gravity_1024sat_l50_vacuum_24600s gravity_4096sat_l20_vacuum_19700s gravity_64sat_l50_vacuum_415000s}"
CASES="${CASES//,/ }"
MODES="${MODES:-inner_only,predictive}"
BASE=(${MODES//,/ })
ROUNDS="${ROUNDS:-3}"; REPEATS="${REPEATS:-8}"; WARMUP="${WARMUP:-1}"
STORE=output/parallel_policy_state
CONSTANTS="$OUT/constants"

PAIRS=()   # "<case>:<mode>", cases in order, modes within case
for cs in $CASES; do for m in "${BASE[@]}"; do PAIRS+=("$cs:$m"); done; done
P=${#PAIRS[@]}
SEQ=()   # "<case>:<mode>:_r<round>": round r runs pair i at position (i + 3r) mod P
for ((r = 0; r < ROUNDS; r++)); do
  for ((pos = 0; pos < P; pos++)); do
    i=$(( ((pos - 3 * r) % P + P) % P ))
    SEQ+=("${PAIRS[i]}:_r$((r + 1))")
  done
done
say() { echo "[pilot-final $(date -Iseconds)] $*"; }
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

# Position balance per mode (mean 1-based position within a round over all rounds).
for m in "${BASE[@]}"; do
  awk -v m="$m" -v P="$P" -v R="$ROUNDS" -v s="${SEQ[*]}" 'BEGIN { n = split(s, a, " "); for (k = 1; k <= n; k++) { split(a[k], f, ":"); if (f[2] == m) { t += (k - 1) % P + 1; c++ } } printf "[pilot-final] mean position of %s: %.2f over %d processes\n", m, t / c, c }'
done

# Not timed: bring both precompile caches up to date for this source tree. Normally already done on the `shared`
# partition (precompile_tree.sh, run via `expanse.sh run` without --bench), so each line here is only a cache check
# (~10 s). If it prints `Precompiling`, the exclusive node is paying for it: say so, and the run is still valid
# (every timed child starts after this). Machine cost constants (usl_alpha_base, usl_beta_alloc): measured once at the
# budget's thread count, put back into the emptied policy store before every call.
[ "$DRY" = 1 ] || mkdir -p "$OUT"
# The examples manifest is gitignored. Instantiate it first (a no-op when it is complete; if the file is missing, run
# `julia --project=. examples/setup.jl` on the shared partition first, e.g. via precompile_tree.sh).
if [ "$DRY" = 1 ]; then run julia --project=examples -e 'using Pkg; Pkg.instantiate()'
else
  say "+ julia --project=examples -e 'using Pkg; Pkg.instantiate()'"
  julia --project=examples -e 'using Pkg; Pkg.instantiate()' || { say "ERROR: examples instantiate failed (no examples/Manifest.toml? run examples/setup.jl via precompile_tree.sh)"; exit 7; }
fi
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
