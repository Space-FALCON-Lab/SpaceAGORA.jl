#!/usr/bin/env bash
# Second pilot: the two things the single-sample vacuum pilot (pilot_final_vs_static.sh, job 55488643) did not cover.
#   Stage gram: aero_4096sat_l50_gram_lookahead_100s (live native GRAM through the look-ahead cache, single simulation).
#               Modes outer_inner_static, inner_only, predictive. Thread budget GRAM_BUDGETS (default 16, the budget of the
#               paper row whose gap was +9.2%; add 32 to also run the 32-thread budget as a second, separate stage).
#   Stage mc:   mcgrid_16sat_8mc_long and mcgrid_8sat_16mc_long (multi-sample, 10 h missions, both routing axes live).
#               Modes outer_process, outer_threads, predictive at MC_BUDGET (default 32). The process route starts at most
#               min(samples, budget) workers (8 and 16), each with --threads=1.
# 3 fresh processes per (case, mode). Runs INSIDE a Slurm job (julia_job.slurm does the quiet-node gate):
#
#   EXPANSE_ACCOUNT=mia402 ./expanse.sh run --bench --time 01:30:00 -- bash scripts/remote/expanse/pilot2_gram_multisample.sh
#   DRY=1 bash scripts/remote/expanse/pilot2_gram_multisample.sh      # prints the points and commands, runs no Julia
#   STAGES=gram|mc|gram,mc  GRAM_BUDGETS=16|32|"16,32"  MC_BUDGET=32  ROUNDS=3  GRAM_REPEATS=6  MC_REPEATS=8
#   GRAM_CASES= GRAM_MODES= MC_CASES= MC_MODES=   (comma lists; spaces would not survive expanse.sh)
#
# Startup cost: one controller process per stage (--points-file), so the package loads once per stage, not once per call.
# Isolation: ONE harness call per point, SPACEAGORA_PPB_SHARED_PROCESS=0, so every point runs in its own fresh Julia
# process (harness child) with its own FRESH RHS calibration store, snapshot dir and emptied policy-state dir
# (cost constants of that stage's budget put back). No process, store or compiled-code sharing between points.
# Order and position balance: per stage, pairs (case, mode) are indexed i = 0..P-1 (cases in order, modes within case).
# Round r runs pair i at position (i + SHIFT*r) mod P, SHIFT = P/3 when 3 divides P (gram P=3: SHIFT=1; mc P=6: SHIFT=2),
# else 1. For P=3 and P=6 with 3 modes and 3 rounds this is exact: every pair (P=3) or every mode summed over the cases
# (P=6) visits each of the P positions exactly once, so mean position is (P+1)/2 for every mode (DRY=1 prints the per-mode
# position histogram to check it). For other P or ROUNDS the balance is approximate and the printed histogram says how much.
# Output: results/expanse_policy/<jobid>_pilot2/<stage>/calls/<case>__<mode>_r<round>{/,.log,_snap/}. Raw CSV per call has
# wall_time_s and final_primary_pos_norm_m / vel_norm_mps / mass_kg / terminal_time_s (compare across modes for
# bit-identical trajectories); _snap/ is the per-solve snapshot (wall time, calibration store, plan source, confirm state).
# Data only, no plots, no in-job summary (analysis happens locally on the pulled CSVs).
set -u
cd "$(dirname "$0")/../../.." || exit 1
DRY="${DRY:-0}"
JOB="${SLURM_JOB_ID:-local}"
OUT="$PWD/results/expanse_policy/${JOB}_pilot2"
STAGES="${STAGES:-gram,mc}"; STAGES="${STAGES//,/ }"
GRAM_BUDGETS="${GRAM_BUDGETS:-16}"; GRAM_BUDGETS="${GRAM_BUDGETS//,/ }"
GRAM_CASES="${GRAM_CASES:-aero_4096sat_l50_gram_lookahead_100s}"
GRAM_MODES="${GRAM_MODES:-outer_inner_static,inner_only,predictive}"
MC_BUDGET="${MC_BUDGET:-32}"
MC_CASES="${MC_CASES:-mcgrid_16sat_8mc_long,mcgrid_8sat_16mc_long}"
MC_MODES="${MC_MODES:-outer_process,outer_threads,predictive}"
ROUNDS="${ROUNDS:-3}"; GRAM_REPEATS="${GRAM_REPEATS:-6}"; MC_REPEATS="${MC_REPEATS:-8}"; WARMUP="${WARMUP:-1}"
STORE=output/parallel_policy_state

say() { echo "[pilot2 $(date -Iseconds)] $*"; }
run() { say "+ $*"; [ "$DRY" = 1 ] && return 0; "$@"; }

export JULIA_NUM_THREADS="" SPACEAGORA_PERF_PROCS="" SPACEAGORA_PPB_WORKLOAD=0
# Anything that would leak into a mode run (modes set their own SPACEAGORA_RHS_CALIBRATE / REGION_WIDTH via extra_env).
unset JULIA_NUM_GC_THREADS JULIA_THREAD_SLEEP_THRESHOLD SPACEAGORA_RHS_CALIBRATE SPACEAGORA_RHS_CALIBRATE_SEED_TOPK \
      SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER SPACEAGORA_POLICY_REGRET_GUARD \
      SPACEAGORA_POLICY_REGRET_MARGIN SPACEAGORA_POLICY_REGRET_MIN_SOLVES SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR \
      SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN SPACEAGORA_RHS_REGION_WIDTH SPACEAGORA_RHS_CALIBRATE_CONFIRM \
      SPACEAGORA_RHS_CALIBRATE_CONFIRM_SOLVES SPACEAGORA_RHS_CALIBRATE_CONFIRM_MARGIN
export SPACEAGORA_PPB_SHARED_PROCESS=0 SPACEAGORA_PPB_WARMUP_SPAN=0.1 SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto
export SPACEAGORA_CAMPAIGN_CORRECTIONS=off SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1 SPACEAGORA_PPB_NO_PLOTS=1
export SPACEAGORA_PPC_EQUAL_CORE_BUDGET=1 OPENBLAS_NUM_THREADS=1

say "job=$JOB host=$(hostname) commit=${GIT_COMMIT:-unknown} DRY=$DRY stages=[$STAGES] gram_budgets=[$GRAM_BUDGETS] mc_budget=$MC_BUDGET rounds=$ROUNDS warmup=$WARMUP gram_repeats=$GRAM_REPEATS mc_repeats=$MC_REPEATS"

# Prerequisites for the GRAM stage (read-only checks; nothing is fixed here): the data-store link, the native library,
# and the Earth GRAM data and SPICE kernels under it.
GRAMDIR="data/GRAMSuite.jl/GRAM Suite 2.0"
if [[ " $STAGES " == *" gram "* ]]; then
  for p in "$GRAMDIR/Build/lib/libGRAM.so" "$GRAMDIR/SPICE/pck/pck00011.tpc" "$GRAMDIR/Earth/data"; do
    if [ -e "$p" ]; then say "gram prerequisite present: $p"
    elif [ "$DRY" = 1 ]; then say "gram prerequisite (not checked in DRY): $p"
    else say "ERROR: gram prerequisite missing: $p (data store link or native build; see the task notes)"; exit 4; fi
  done
fi
[ "$DRY" = 1 ] || { say "memory: $(grep -E 'MemTotal|MemAvailable' /proc/meminfo | tr -s ' ' | tr '\n' ' ')"; }

# Not timed: bring both precompile caches up to date for this source tree. Normally already done on the `shared`
# partition (precompile_tree.sh, run via `expanse.sh run` without --bench), so each line here is only a cache check
# (~10 s). If it prints `Precompiling`, the exclusive node is paying for it: say so, and the run is still valid
# (every timed child starts after this). The examples manifest is gitignored: instantiate it first (a no-op when
# complete; if missing, run `julia --project=. examples/setup.jl` on the shared partition first, e.g. via precompile_tree.sh).
[ "$DRY" = 1 ] || mkdir -p "$OUT"
PRE_THREADS=$MC_BUDGET
if [ "$DRY" = 1 ]; then run julia --project=examples -e 'using Pkg; Pkg.instantiate()'
else
  say "+ julia --project=examples -e 'using Pkg; Pkg.instantiate()'"
  julia --project=examples -e 'using Pkg; Pkg.instantiate()' || { say "ERROR: examples instantiate failed (no examples/Manifest.toml? run examples/setup.jl via precompile_tree.sh)"; exit 7; }
fi
for proj in . examples; do
  pl="$OUT/prelude_${proj//./root}.log"
  if [ "$DRY" = 1 ]; then run julia --project=$proj --threads=$PRE_THREADS -e 'using SpaceAGORA'; continue; fi
  say "+ julia --project=$proj --threads=$PRE_THREADS -e 'using SpaceAGORA'"
  julia --project=$proj --threads=$PRE_THREADS -e 'using SpaceAGORA' 2>&1 | tee "$pl"
  [ "${PIPESTATUS[0]}" -eq 0 ] || { say "ERROR: prelude --project=$proj failed"; exit 6; }
  grep -q Precompiling "$pl" && say "WARNING: --project=$proj precompiled on the exclusive node; run precompile_tree.sh on the shared partition first next time"
done
if [[ " $STAGES " == *" gram "* ]]; then   # GRAM must load in the native-library path too (a load failure here is cheaper than one per point)
  if [ "$DRY" = 1 ]; then run julia --project=. --threads=$PRE_THREADS -e 'using SpaceAGORA; SpaceAGORA.RuntimeServices._prepare_gramsuite_load_path!(joinpath(pwd(), "data", "GRAMSuite.jl"); prepend=true); @eval import GRAMSuite; ext = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing; println("GRAMSuite loaded, extension attached: ", ext); ext || exit(4)'
  else
    say "+ julia --project=. --threads=$PRE_THREADS -e '<load GRAMSuite via the vendored load path>'"
    julia --project=. --threads=$PRE_THREADS -e 'using SpaceAGORA; SpaceAGORA.RuntimeServices._prepare_gramsuite_load_path!(joinpath(pwd(), "data", "GRAMSuite.jl"); prepend=true); @eval import GRAMSuite; ext = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing; println("GRAMSuite loaded, extension attached: ", ext); ext || exit(4)' 2>&1 | tee "$OUT/prelude_gram.log"
    [ "${PIPESTATUS[0]}" -eq 0 ] || { say "ERROR: GRAMSuite did not load"; exit 3; }
    grep -q Precompiling "$OUT/prelude_gram.log" && say "WARNING: GRAM extension precompiled on the exclusive node; add it to precompile_tree.sh's shared-partition run"
  fi
fi

FAILED=()
# stage <name> <budget> <cases...> -- <modes,comma> <repeats>
stage() {
  local name=$1 budget=$2 modes=$3 repeats=$4; shift 4
  local cases=("$@") sdir="$OUT/${name}_b${budget}"
  local BASE=(${modes//,/ }) PAIRS=() SEQ=() cs m r pos i P SHIFT
  for cs in "${cases[@]}"; do for m in "${BASE[@]}"; do PAIRS+=("$cs:$m"); done; done
  P=${#PAIRS[@]}
  if [ $((P % 3)) -eq 0 ]; then SHIFT=$((P / 3)); else SHIFT=1; fi
  SHIFT="${SHIFT_OVERRIDE:-$SHIFT}"
  # ROUNDS=5 with P=3: rounds 4 and 5 add one forward and one reversed order, which makes the mean position exactly (P+1)/2
  # for every pair (the 3-round cyclic block is exact, forward+reverse sums to P+1). Round order below has no same pair
  # ending one round and starting the next. Any other ROUNDS/P keeps the cyclic rotation.
  local ORDER5=(c0 f c2 c1 r) rk
  for ((r = 0; r < ROUNDS; r++)); do
    rk=cyc; if [ "$ROUNDS" -eq 5 ] && [ "$P" -eq 3 ] && [ -z "${SHIFT_OVERRIDE:-}" ]; then rk=${ORDER5[r]}; fi
    for ((pos = 0; pos < P; pos++)); do
      case $rk in
        cyc) i=$(( ((pos - SHIFT * r) % P + P) % P )) ;;
        c0) i=$pos ;; c1) i=$(( ((pos - 1) % P + P) % P )) ;; c2) i=$(( ((pos - 2) % P + P) % P )) ;;
        f) i=$pos ;; r) i=$((P - 1 - pos)) ;;
      esac
      SEQ+=("${PAIRS[i]}:_r$((r + 1))")
    done
  done
  say "=== stage $name budget=$budget: ${#SEQ[@]} processes, P=$P pairs, shift=$SHIFT, repeats=$repeats, in order:"
  local n=0 item rest; for item in "${SEQ[@]}"; do n=$((n + 1)); rest=${item#*:}; say "  $n. ${item%%:*} ${rest%%:*}${rest#*:}"; done
  for m in "${BASE[@]}"; do   # position histogram per mode (1-based position within a round), summed over cases and rounds
    awk -v m="$m" -v P="$P" -v s="${SEQ[*]}" 'BEGIN { n = split(s, a, " "); for (k = 1; k <= n; k++) { split(a[k], f, ":"); if (f[2] == m) { p = (k - 1) % P + 1; h[p]++; t += p; c++ } }
      printf "[pilot2] position histogram of %s: ", m; for (p = 1; p <= P; p++) printf "%d ", h[p]; printf " mean %.2f over %d processes\n", t / c, c }'
  done
  local constants="$sdir/constants" pts="" points="$sdir/points.tsv" tag
  [ "$DRY" = 1 ] || mkdir -p "$sdir/stores" "$sdir/calls" "$constants"
  # Machine cost constants (usl_alpha_base, usl_beta_alloc), measured at this stage's thread count, put back into the
  # emptied policy store before every call. Not timed.
  if [ "${SKIP_CALIBRATE:-0}" != 1 ]; then
    run env SPACEAGORA_CORE_BUDGET=$budget SPACEAGORA_PPC_PHYSICAL_CORES=$budget julia --project=. --threads=$budget scripts/calibration/calibrate_machine.jl
    if [ "$DRY" != 1 ]; then
      cp -p "$STORE"/cost_constants_*.toml "$constants"/ || { say "ERROR: calibration wrote no constants"; FAILED+=("$name:calibrate"); return 1; }
    fi
  fi
  export SPACEAGORA_CORE_BUDGET=$budget SPACEAGORA_PPC_PHYSICAL_CORES=$budget
  # ONE controller process runs all calls of this stage. Points file: tag, case, mode, outdir, store, snapshot dir. Per call
  # the controller sets SPACEAGORA_RHS_CALIBRATION_PATH / SPACEAGORA_PPB_SNAPSHOT_DIR, runs the prep command (fresh policy-state
  # dir and snapshot dir, cost constants put back, no store), then runs the same harness call as a single call would
  # (--cases=<case> --modes=<mode> --outdir=<outdir>) with a fresh child, logging to <outdir>.log. A call fails on a failed
  # child or a `Precompiling` line in its log.
  export SPACEAGORA_PPC_POINT_PREP="rm -rf '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; mkdir -p '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; cp -p '$constants'/cost_constants_*.toml '$STORE'/ 2>/dev/null; rm -f \"\$SPACEAGORA_RHS_CALIBRATION_PATH\"; bash '$PWD/scripts/remote/expanse/node_snapshot.sh' \"\${SPACEAGORA_PPB_SNAPSHOT_DIR%_snap}_node_before.txt\""
  for item in "${SEQ[@]}"; do
    cs=${item%%:*}; rest=${item#*:}; m=${rest%%:*}; tag="${cs}__${m}${rest#*:}"
    pts+="$tag\t$cs\t$m\t$sdir/calls/$tag\t$sdir/stores/$tag.toml\t$sdir/calls/${tag}_snap\n"
  done
  if [ "$DRY" = 1 ]; then
    say "points file $points:"; printf "$pts" | sed 's/^/  /'
    say "prep: $SPACEAGORA_PPC_POINT_PREP"
  else
    printf "$pts" > "$points"
  fi
  local cmd=(julia --project=examples benchmarks/studies/parallelization_performance.jl full
             --points-file="$points" --threads=$budget --process-workers=$budget
             --repeats=$repeats --warmup=$WARMUP --parity-cases=none)
  if [ "$DRY" = 1 ]; then echo "  ${cmd[*]}"; return 0; fi
  "${cmd[@]}"; local rc=$?
  [ $rc -eq 0 ] || { say "FAILED stage $name (rc=$rc); see [points] lines above and $sdir/calls/*.log"; FAILED+=("$name:$budget"); return 1; }
  say "stage $name budget=$budget done"
}

for s in $STAGES; do
  case $s in
    gram) for b in $GRAM_BUDGETS; do stage gram "$b" "$GRAM_MODES" "$GRAM_REPEATS" ${GRAM_CASES//,/ }; done ;;
    mc)   stage mc "$MC_BUDGET" "$MC_MODES" "$MC_REPEATS" ${MC_CASES//,/ } ;;
    *)    say "ERROR: unknown stage $s"; exit 2 ;;
  esac
done
[ ${#FAILED[@]} -eq 0 ] || { say "FAILED: ${FAILED[*]}"; exit 1; }
say "done: $OUT"
