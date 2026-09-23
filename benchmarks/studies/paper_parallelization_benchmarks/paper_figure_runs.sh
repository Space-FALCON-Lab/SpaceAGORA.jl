#!/usr/bin/env bash
# The paper's figure runs: which phases go on which machine, in what order, from
# what calibration-store state.
#
# protocol.sh runs "the suite"; this runs "the figures". The difference is that a
# figure has a machine, an ordering and a store state attached to it, and getting
# any of those wrong produces a plausible-looking CSV that answers a different
# question. Every step below is therefore named, ordered and refusable rather
# than left to the launcher's memory.
#
#   bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh <target> [--execute]
#   TARGET=<target> bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh [--execute]
#
# Targets:
#   workstation   F4's cross-machine arm: P3 and P5 on this 12-core box, every
#                 mode including predictive (R7) and policy_v2 (R6), 11 repeats,
#                 from a cold calibration store. Runs here, now.
#   trx50         the full benchmark-box sequence (calibration, targeted points,
#                 the cold 11-repeat P1-P5, the converged P1/P5, P6/P6p, and the
#                 archive call after each). Every step is a
#                 scripts/remote/spaceagora-remote push; this target PRINTS them
#                 and, with --execute, offers each one in order for confirmation.
#   calibrate-p6  re-derive P6's per-trace mission lengths on the host it runs
#                 on, the way calibrate_iso_ladder.sh does for PPC_L50_ISO_MISSION_S.
#                 This is itself a timed measurement: benchmark box, idle, alone.
#
# Without --execute every target is a dry description and touches nothing.
#
# Why 11 repeats on the workstation arm. Five repeats put the 90% interval of an
# adaptive point's median at ~7%, which is the same size as the band used to
# decide whether two runs differ at all, so an adaptive point cannot be compared
# against anything; eleven takes it to ~2.2% and it plateaus there. See
# _ppb_min_repeats in main.jl. F4 compares one machine's adaptive arm against
# another's, so it is exactly the comparison that needs the tighter interval.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
STUDY_DIR="${REPO_ROOT}/benchmarks/studies/paper_parallelization_benchmarks"
PPB="${REPO_ROOT}/benchmarks/studies/paper_parallelization_benchmarks.jl"
PPC_LAUNCHER="${REPO_ROOT}/benchmarks/studies/parallelization_performance.jl"
STORE_DIR="${REPO_ROOT}/output/parallel_policy_state"

# Where an archived run goes. A benchmark run lands in a gitignored output/
# directory inside a worktree, so a reboot, a git clean or a deleted worktree
# takes it with it; every figure has to stay regenerable from a directory that
# outlives the worktree. Override with PAPER_DATA_RAW.
PAPER_DATA_RAW="${PAPER_DATA_RAW:-/home/space-falcon-1/Documents/SpaceAGORA-paper-data/data/raw}"
ARCHIVE="${REPO_ROOT}/scripts/archive_paper_run.py"

# The remote, and where its releases live. Both are properties of
# scripts/remote/remotes.conf, repeated here only so the rsync line below can be
# written out in full -- spaceagora-remote's auto-pull takes a release's
# top-level output/ and results/ and nothing else.
REMOTE="${REMOTE:-trx50}"
REMOTE_BASE="${REMOTE_BASE:-~/spaceagora_remote}"

TARGET="${TARGET:-}"
EXECUTE=0
for arg in "$@"; do
  case "$arg" in
    --execute) EXECUTE=1 ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "error: unknown option '$arg'" >&2; exit 2 ;;
    *)  TARGET="$arg" ;;
  esac
done
if [ -z "$TARGET" ]; then
  echo "error: no target. Use 'workstation', 'trx50' or 'calibrate-p6'." >&2
  exit 2
fi

say()  { printf '[figure-runs] %s\n' "$*"; }
step() { printf '\n[figure-runs] ── step %s ── %s\n' "$1" "$2"; }
cmd()  { printf '    %s\n' "$*"; }

physical_cores() {
  if [ -n "${SPACEAGORA_PPC_PHYSICAL_CORES:-}" ]; then
    echo "${SPACEAGORA_PPC_PHYSICAL_CORES}"
    return
  fi
  local n
  n="$(awk -F: '/^core id/{print $2"/"id} /^physical id/{id=$2}' /proc/cpuinfo 2>/dev/null | sort -u | wc -l)"
  [ "${n:-0}" -ge 1 ] || n="$(nproc 2>/dev/null || echo 1)"
  echo "$n"
}

# Refuse to start a timed run on a machine that is already working.
#
# This is the same guard ppc_assert_machine_quiet! applies inside the harness
# (50% load headroom), hoisted to the launcher so the refusal happens before the
# store is rotated and the calibration is spent, not forty minutes into a phase.
# The julia check is the half that matters most here: two benchmark jobs on one
# box do not merely add noise to each other, they make both sets of numbers
# unusable and the machine unusable for its owner.
require_quiet() {
  local cores load busy others
  cores="$(physical_cores)"
  load="$(awk '{print $1}' /proc/loadavg)"
  busy="$(awk -v l="$load" -v c="$cores" 'BEGIN{print (l > 0.5 * c) ? 1 : 0}')"
  # A Julia language server is an editor feature, not a job: it sits near 0% CPU
  # and is always there. Anything else running Julia is treated as a job.
  others="$(ps -eo pcpu,args | grep '[j]ulia' | grep -v 'languageserver' | grep -v 'vsc-jl' | wc -l)"
  say "machine check: load1=${load} physical_cores=${cores} other_julia_processes=${others}"
  if [ "$busy" = "1" ] || [ "${others}" -gt 0 ]; then
    echo "[figure-runs] REFUSING: this machine is busy." >&2
    echo "[figure-runs]   load1=${load} against ${cores} physical cores (needs <= $(awk -v c="$cores" 'BEGIN{print 0.5*c}'))," >&2
    echo "[figure-runs]   and ${others} other Julia process(es) running." >&2
    echo "[figure-runs] Wait for them. Serializing the launches is the only thing that" >&2
    echo "[figure-runs] prevents two jobs contaminating each other's measurements." >&2
    exit 3
  fi
}

# Empty the calibration store, keeping what was there.
#
# A cold arm means the adaptive profiles meet every point without prior hints;
# leaving yesterday's store in place makes the run measure convergence instead.
# The previous contents are MOVED to a dated sibling, never deleted: a converged
# store is hours of measurement and is the input to the converged arm.
cold_store() {
  local stamp backup
  stamp="$(date -u +%Y%m%d_%H%M%S)"
  backup="${STORE_DIR}_backup_${stamp}"
  if [ -d "$STORE_DIR" ] && [ -n "$(ls -A "$STORE_DIR" 2>/dev/null)" ]; then
    say "calibration store -> ${backup}"
    if [ "$EXECUTE" = "1" ]; then
      mv "$STORE_DIR" "$backup" || exit 4
    else
      cmd "mv ${STORE_DIR} ${backup}"
    fi
  else
    say "calibration store is already empty or absent: ${STORE_DIR}"
  fi
  if [ "$EXECUTE" = "1" ]; then
    mkdir -p "$STORE_DIR"
  else
    cmd "mkdir -p ${STORE_DIR}"
  fi
}

run() {
  if [ "$EXECUTE" = "1" ]; then
    say "running: $*"
    "$@"
    local rc=$?
    [ $rc -eq 0 ] || { echo "[figure-runs] step failed with exit ${rc}" >&2; exit "$rc"; }
  else
    cmd "$*"
  fi
}

# Print a step and, under --execute, offer it. Used by the trx50 target, where
# every step is a job on another machine that takes hours and must not be
# launched on top of the previous one.
offer() {
  local label="$1"; shift
  printf '\n'
  cmd "$*"
  [ "$EXECUTE" = "1" ] || return 0
  printf '[figure-runs] run this step now? [y/N] '
  local reply=""
  read -r reply || true
  case "$reply" in
    y|Y|yes|YES) say "launching: ${label}"; bash -c "$*" ;;
    *) say "skipped: ${label}" ;;
  esac
}

# ── workstation: F4's cross-machine arm ──────────────────────────────────────
target_workstation() {
  local cores threads workers
  cores="$(physical_cores)"
  threads="$(seq -s, 1 1 "$cores" >/dev/null 2>&1 && echo "")"
  # The budget ladder the P phases derive from the host (geometric to the
  # physical core count, the count itself always included) -- written out here so
  # the launcher passes exactly what the dry run showed.
  threads="1"
  local b=1
  while [ $((b * 2)) -lt "$cores" ]; do b=$((b * 2)); threads="${threads},${b}"; done
  [ "$b" -eq "$cores" ] || threads="${threads},${cores}"
  workers="$cores"

  say "target       = workstation ($(hostname))"
  say "phases       = P3, P5   (F4: the same two Monte Carlo phases the benchmark box runs)"
  say "threads      = ${threads}"
  say "workers      = ${workers}"
  say "repeats      = 11 (SPACEAGORA_PPB_MIN_REPEATS)"
  say "store        = cold"
  say "modes        = each phase's own full ladder, predictive and policy_v2 included"

  require_quiet

  step 1 "machine calibration for the predictive planner (R7)"
  # The predictive campaign planner reads per-machine cost constants from
  # output/parallel_policy_state/cost_constants_<fingerprint>.toml, written by
  # this script at the thread count the jobs will use. Without the file the
  # planner still runs but models no contention, so the predictive rows measure
  # something other than what was intended -- and on a cross-machine figure that
  # is the difference being reported. It runs BEFORE the store is emptied, then
  # is re-run after, because emptying the store takes the constants with them.
  run julia --project="$REPO_ROOT" --threads="$cores" "${REPO_ROOT}/scripts/calibrate_machine.jl"

  step 2 "cold calibration store"
  cold_store

  step 3 "re-calibrate into the now-empty store"
  run julia --project="$REPO_ROOT" --threads="$cores" "${REPO_ROOT}/scripts/calibrate_machine.jl"

  step 4 "P3 and P5, 11 repeats"
  if [ "$EXECUTE" = "1" ]; then
    SPACEAGORA_PPB_MIN_REPEATS=11 SPACEAGORA_CAMPAIGN_CORRECTIONS=off OPENBLAS_NUM_THREADS=1 GKSwstype=100 \
      julia --project="$REPO_ROOT" "$PPB" \
        --phases=P3,P5 --threads="${threads}" --process-workers="${workers}"
    rc=$?
    [ $rc -eq 0 ] || { echo "[figure-runs] P3/P5 failed with exit ${rc}" >&2; exit "$rc"; }
  else
    cmd "SPACEAGORA_PPB_MIN_REPEATS=11 SPACEAGORA_CAMPAIGN_CORRECTIONS=off OPENBLAS_NUM_THREADS=1 GKSwstype=100 julia --project=. $PPB --phases=P3,P5 --threads=${threads} --process-workers=${workers}"
  fi

  step 5 "archive the run"
  say "run directory: output/performance/paper_benchmarks/<stamp> (the newest one)"
  cmd "python3 ${ARCHIVE} output/performance/paper_benchmarks/<stamp> --archive ${PAPER_DATA_RAW} --machine workstation --store cold --notes 'F4 cross-machine arm, P3+P5, 11 repeats'"

  say ""
  say "Inspect it first: --dry-run on the same command line prints the rungs"
  say "without running anything."
}

# ── trx50: the benchmark-box sequence ────────────────────────────────────────
target_trx50() {
  local remote_sh="${REPO_ROOT}/scripts/remote/spaceagora-remote"
  say "target = ${REMOTE} (the benchmark box), via ${remote_sh}"
  say ""
  say "Every step is one detached job. They are ordered and must not overlap:"
  say "the calibration is itself a timed measurement, the targeted points gate"
  say "the eleven-hour run, and the converged arm needs the store the cold arm"
  say "left behind. Check 'spaceagora-remote status <job-id>' (or list) before"
  say "launching the next one."
  say ""
  say "Flags take a SPACE, not an '=': --threads 1,2,4,8,16,32. The '=' form is"
  say "rejected by the argument parser."
  say ""
  say "The calibration store lives on the remote at ${REMOTE_BASE}/policy_state and"
  say "is symlinked into each release's output/parallel_policy_state, so it"
  say "persists across jobs and is NOT shipped by the push. Its state is set over"
  say "ssh, out of band, before the job that depends on it."

  step 1 "machine calibration at the thread count the jobs use"
  offer "calibrate" \
    "${remote_sh} push --remote ${REMOTE} --threads 32 -- julia --project=. --threads=32 scripts/calibrate_machine.jl"

  step 2 "targeted points — one job per POINT, store set per point"
  say "finding8 wants the converged store; the other three want it empty."
  offer "targeted finding8 (converged store)" \
    "ssh ${REMOTE} 'rm -rf ${REMOTE_BASE}/policy_state && cp -a ${REMOTE_BASE}/policy_state_converged_20260918 ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env POINT=finding8 bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh"
  offer "targeted p5_16sat (cold store)" \
    "ssh ${REMOTE} 'mv ${REMOTE_BASE}/policy_state ${REMOTE_BASE}/policy_state_bak_\$(date -u +%Y%m%d_%H%M%S); mkdir -p ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env POINT=p5_16sat bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh"
  offer "targeted finding9 (cold store)" \
    "ssh ${REMOTE} 'mv ${REMOTE_BASE}/policy_state ${REMOTE_BASE}/policy_state_bak_\$(date -u +%Y%m%d_%H%M%S); mkdir -p ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env POINT=finding9 bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh"
  offer "targeted defectA (cold store)" \
    "ssh ${REMOTE} 'mv ${REMOTE_BASE}/policy_state ${REMOTE_BASE}/policy_state_bak_\$(date -u +%Y%m%d_%H%M%S); mkdir -p ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env POINT=defectA bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh"

  step 3 "the cold 11-repeat P1-P5 run"
  offer "P1-P5 cold, 11 repeats" \
    "ssh ${REMOTE} 'mv ${REMOTE_BASE}/policy_state ${REMOTE_BASE}/policy_state_bak_\$(date -u +%Y%m%d_%H%M%S); mkdir -p ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env SPACEAGORA_PPB_MIN_REPEATS=11 julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl --phases=P1,P2,P3,P4,P5 --threads=1,2,4,8,16,32 --process-workers=32"

  step 4 "the converged P1/P5 arm"
  say "Restore the converged snapshot first; it is the input, not a by-product."
  offer "P1,P5 converged, 11 repeats" \
    "ssh ${REMOTE} 'rm -rf ${REMOTE_BASE}/policy_state && cp -a ${REMOTE_BASE}/policy_state_converged_20260918 ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- env SPACEAGORA_PPB_MIN_REPEATS=11 julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl --phases=P1,P5 --threads=1,2,4,8,16,32 --process-workers=32"

  step 5 "P6 and P6p, 3 repeats (figure F2)"
  say "Recalibrate the mission lengths first if this box has not run P6 before:"
  cmd "${remote_sh} push --remote ${REMOTE} --threads 1 -- bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh calibrate-p6 --execute"
  offer "P6 + P6p cold, 3 repeats" \
    "ssh ${REMOTE} 'mv ${REMOTE_BASE}/policy_state ${REMOTE_BASE}/policy_state_bak_\$(date -u +%Y%m%d_%H%M%S); mkdir -p ${REMOTE_BASE}/policy_state' && ${remote_sh} push --remote ${REMOTE} --threads 1,2,4,8,16,32 --process-workers 32 -- julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl --phases=P6,P6p --threads=1,2,4,8,16,32 --process-workers=32"

  step 6 "pull each run back, then archive it"
  say "spaceagora-remote's auto-pull takes a release's top-level output/ and"
  say "results/ only, and it runs on THIS machine -- if it slept or rebooted"
  say "while the job ran, nothing was pulled. Name the path either way:"
  cmd "rsync -az ${REMOTE}:${REMOTE_BASE}/releases/<job-id>/output/performance/paper_benchmarks/ ${REPO_ROOT}/output/performance/paper_benchmarks/"
  say ""
  say "Then archive, with --store naming the state that job actually ran under:"
  cmd "python3 ${ARCHIVE} output/performance/paper_benchmarks/<stamp> --archive ${PAPER_DATA_RAW} --machine trx50 --store cold      --notes 'P1-P5, 11 repeats, cold store'"
  cmd "python3 ${ARCHIVE} output/performance/paper_benchmarks/<stamp> --archive ${PAPER_DATA_RAW} --machine trx50 --store converged --notes 'P1+P5, 11 repeats, converged store'"
  cmd "python3 ${ARCHIVE} output/performance/paper_benchmarks/<stamp> --archive ${PAPER_DATA_RAW} --machine trx50 --store cold      --notes 'P6+P6p, 3 repeats, figure F2'"
  say ""
  say "And check the archive reads back:"
  cmd "python3 ${ARCHIVE} --verify --archive ${PAPER_DATA_RAW}"
}

# ── calibrate-p6: re-derive the per-trace mission lengths ────────────────────
#
# Same shape as calibrate_iso_ladder.sh, for P6's five new traces instead of the
# L50 iso-work ladder: run each once in serial mode at the full profile, and
# print the mission length that would put its serial baseline at the target.
# Reads the case names out of cli.jl so the two cannot drift.
target_calibrate_p6() {
  local target_s="${P6_TARGET_S:-11.86}"
  local out; out="${TMPDIR:-/tmp}/p6_calibration_$$"

  say "target = calibrate-p6 on $(hostname)"
  say "serial-baseline target = ${target_s} s"
  say "(the default is the measured TRX50 serial baseline of P6's trace 2,"
  say " gravity_4096sat_l50_vacuum_5800s, run 20260918_162845 -- so every trace"
  say " is sized against the rung P1 and P2 already report)"

  local n l20_s vac_s aero_s
  n="$(sed -n 's/^const PPB_P6_N_SAT  *= *\([0-9]*\).*/\1/p' "${STUDY_DIR}/cli.jl" | head -1)"
  l20_s="$(sed -n 's/^const PPB_P6_L20_MISSION_S  *= *\([0-9]*\).*/\1/p' "${STUDY_DIR}/cli.jl" | head -1)"
  vac_s="$(sed -n 's/^const PPB_P6_VACUUM_MISSION_S  *= *\([0-9]*\).*/\1/p' "${STUDY_DIR}/cli.jl" | head -1)"
  aero_s="$(sed -n 's/^const PPB_P6_AERO_MISSION_S  *= *\([0-9]*\).*/\1/p' "${STUDY_DIR}/cli.jl" | head -1)"
  if [ -z "$n" ] || [ -z "$l20_s" ] || [ -z "$vac_s" ] || [ -z "$aero_s" ]; then
    echo "error: could not read the P6 constants out of ${STUDY_DIR}/cli.jl" >&2
    exit 5
  fi

  # (case name, mission seconds, samples). The process trace is a Monte Carlo
  # case and its sample count is the constellation size, so it is calibrated at
  # its own shape rather than as one solve.
  local rungs=(
    "gravity_${n}sat_l20_vacuum_${l20_s}s ${l20_s} 1"
    "gravity_${n}sat_l50_vacuum_${vac_s}s ${vac_s} 1"
    "gravity_${n}sat_l50_srp_nbody_vacuum_${vac_s}s ${vac_s} 1"
    "aero_${n}sat_l50_expatm_${aero_s}s ${aero_s} 1"
    "aero_${n}sat_l50_gram_lookahead_${aero_s}s ${aero_s} 1"
    "aero_${n}sat_l50_gram_process_${aero_s}s ${aero_s} ${n}"
  )

  if [ "$EXECUTE" != "1" ]; then
    say "would time, one at a time (add --execute to run):"
    for rung in "${rungs[@]}"; do cmd "serial  ${rung%% *}"; done
    say ""
    say "This is a timing measurement. Benchmark box, idle, nothing else running."
    return 0
  fi

  require_quiet
  mkdir -p "$out"
  printf '%-46s %-12s %-12s %-14s\n' "case" "mission (s)" "serial (s)" "suggested (s)"
  for rung in "${rungs[@]}"; do
    set -- $rung
    local case_name="$1" mission="$2" samples="$3" wall suggested
    timeout "${P6_CALIBRATION_TIMEOUT_S:-7200}" \
      julia --threads=1 --project="$REPO_ROOT" "$PPC_LAUNCHER" full --worker \
        --case="$case_name" --mode=serial --thread-count=1 --repeat=1 --worker-repeats=1 \
        --worker-seed=20260616 --worker-mc-samples="$samples" --warmup=1 \
        --solver-mode=auto_stiff --process-workers=1 --parity-samples=512 \
        --outfile="${out}/${case_name}.csv" --parity=0 \
        > "${out}/${case_name}.log" 2>&1
    wall="$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1{print $h["wall_time_s"]; exit}' \
      "${out}/${case_name}.csv" 2>/dev/null)"
    if [ -z "$wall" ]; then
      printf '%-46s %-12s %-12s %-14s  (see %s)\n' "$case_name" "$mission" "FAILED" "--" "${out}/${case_name}.log"
      continue
    fi
    # Linear in mission length, and only above the floor: a rung measuring well
    # under 3 s is extrapolating from mostly fixed per-solve cost and will land
    # short. Re-run those rather than trusting the suggestion.
    suggested="$(awk -v w="$wall" -v m="$mission" -v t="$target_s" 'BEGIN{printf "%d", (m * t / w) + 0.5}')"
    printf '%-46s %-12s %-12.2f %-14s\n' "$case_name" "$mission" "$wall" "$suggested"
  done
  echo
  echo "Raw rows and logs: ${out}"
  echo "Move a suggestion into PPB_P6_{L20,VACUUM,AERO}_MISSION_S in"
  echo "${STUDY_DIR}/cli.jl (and the case-catalog names in"
  echo "parallelization_performance/cases.jl, which carry the duration) if it"
  echo "differs by more than ~20%, then re-run this once to confirm."
  echo "The three aero traces share one duration on purpose -- they are"
  echo "single-variable in the density path only -- so size that one from the"
  echo "expatm rung and let the GRAM rungs land where they land."
}

case "$TARGET" in
  workstation)  target_workstation ;;
  trx50)        target_trx50 ;;
  calibrate-p6) target_calibrate_p6 ;;
  *) echo "error: unknown target '$TARGET' (workstation | trx50 | calibrate-p6)" >&2; exit 2 ;;
esac
