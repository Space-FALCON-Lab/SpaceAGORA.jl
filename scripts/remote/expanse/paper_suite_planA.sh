#!/usr/bin/env bash
# Plan A paper suite for the adaptive policy (R7, harness mode `predictive`) on one exclusive Expanse node.
# 19 workloads of the 64-core table, full budget only (32 threads / 32 workers), in THREE jobs (J1, J2, J3).
# Scope, cost model and per-row pilot facts: tasks/2026-10-21_expanse-remote-jobs/paper_suite_scope.md (Plan A),
# paper_suite_incremental_option.md; case names from sf3_plan_audit.md. Notes: planA_scripts_notes.md.
#
# USAGE (runs INSIDE a Slurm job; julia_job.slurm does the quiet-node gate when launched with --bench)
#   precompile first, once per source tree, on the shared partition (not timed):
#     EXPANSE_ACCOUNT=... ./expanse.sh run --cpus 8 --mem 24G --time 00:45:00 -- bash scripts/remote/expanse/precompile_tree.sh
#   then one bench job per part (one timing job at a time, per the standing rule):
#     EXPANSE_ACCOUNT=... ./expanse.sh run --bench --time 04:00:00 -- bash scripts/remote/expanse/paper_suite_planA.sh J1
#     EXPANSE_ACCOUNT=... ./expanse.sh run --bench --time 04:00:00 -- bash scripts/remote/expanse/paper_suite_planA.sh J2
#     EXPANSE_ACCOUNT=... ./expanse.sh run --bench --time 03:00:00 -- bash scripts/remote/expanse/paper_suite_planA.sh J3
#   DRY=1 bash scripts/remote/expanse/paper_suite_planA.sh J1|J2|J3|all    prints every command and the core-hour
#   estimate, runs no Julia (J1|J2|J3 = one job, all = the three jobs plus the grand total). See ESTIMATOR below.
#
# JOBS (ids are the rows of paper_suite_scope.md Section 0; walltime central/high from the DRY estimate):
#   J1  P1 constellation 1/16/64/256/1024/4096 sat, P6s 4096 sat degree-20 / SRP + third body / exponential atm /
#       GRAM look-ahead                                         (10 rows, 50 processes, 6 controller calls)
#   J2  P3 256 independent samples, P4 32 Mars samples, P4g 32 MarsGRAM samples, P6p 4096 GRAM samples,
#       P4gs serial MarsGRAM baseline                           (5 rows, 14 processes, 10 controller calls)
#   J3  P7 one-orbit missions (3 cases), P5f campaigns of constellations 16sat x8 and 8sat x16
#                                                               (5 rows, 34 processes, 7 controller calls)
#
# DEFAULTS FOR THE OPEN SUB-DECISIONS (scope Section 7, the scope doc's recommendation unless marked):
#   1. Plan A: n = 1 fresh Julia process per (row, mode), every mode of a row in the SAME job.
#   2. GRAM look-ahead row (P6s) at the table's 32 threads, all five modes (not the 16-thread pilot data).
#   3. Fresh process per mode (SPACEAGORA_PPB_SHARED_PROCESS=0), cold policy store, cost constants put back.
#   4. The scope's cuts (Section 3): no 16-thread rung of P6s, no P6s degree-50 duplicate (= P1 4096), no P6ps,
#      no ladder rungs below the full budget (serial once), P5f only thread splits with w <= 8 and process w >= 8.
#   5. MarsGRAM serial baseline (P4gs) at 5 repeats.
#   6. The three 12-core rows (P7 on the workstation) are NOT in these jobs. (P7 here means the three one-orbit
#      rows of the 64-core table, run on Expanse.)
#   7. DEVIATION from the scope's suggested split: P7 is in J3 (with P5f), not J1, to balance walltime
#      (J1 with P7 would reach about 3.8 h at the high bound against the 4 h cap).
#   8. DEVIATION on P5f statics: see "P5F" below.
#   9. Tiered repeats: serial 5, static modes 7, predictive and outer_process 9, P7 all modes 11; 3 warm-ups. The
#      harness takes ONE repeat count per controller call (SPACEAGORA_PPB_MODE_REPEATS does not exist in this tree),
#      so the plan makes one controller call per (sample count, repeat tier). Each extra call costs one package
#      load (33.9 s, measured, startup_cost_options.md); the estimate includes it.
#
# ORDER. Within a call the rows run in table order and each row's modes in that tier are rotated left by the row's
# position in the call (mode order rotated across rows). Single-simulation rows (P1, P6s) are split into two halves
# of five rows: half A runs the repeat tiers 5, 7, 9; half B runs 9, 7, 5, so serial and predictive each visit an
# early and a late position in the job and the statics sit in the middle both times. Rows with their own sample
# count run their own tiers, alternating ascending / descending from row to row. DRY=1 prints the exact order.
#
# P5F (campaign of constellations, deviation 8). This tree's harness has no per-width outer_threads / outer_process
# modes (the old phase driver's split grid is not here). Process width is min(samples, 32) = 8 (16sat x8) and
# 16 (8sat x16). The static set built here is serial, predictive, outer_process, outer_threads, outer_inner_static
# and the plan-pinned thread splits force_threads@w8+l0+b4, @w4+l0+b8, @w2+l0+b16, force_none@w1+l0+b32 (and
# force_threads@w16+l0+b2 for the 16-sample case): 9 and 10 processes against the scope's 13 each. Override with
# P5F_MODES=comma,list (modes must exist in modes.jl). The scope's 13-process grid is reproduced by the estimator
# only (SCOPE_P5F=1 in DRY mode) to check that the model matches the scope's totals.
#
# ENV KNOBS (as in the pilots; comma lists, spaces would not survive expanse.sh)
#   ROWS=id,id        only these rows of the job (rerun after a failure; ids are printed by DRY=1)
#   BUDGET=32  WARMUP=3  REP_SERIAL=5  REP_STATIC=7  REP_ADAPT=9  REP_P7=11  P7_MISSION_S=5702
#   P5F_MODES=...     override the P5f mode list
#   MODES=a,b         only these modes in every non-P5f row (P5F_MODES wins for P5f rows); smoke tests
#   PROFILE=smoke     harness profile (default full); smoke shortens every mission, for functional smoke tests only
#   SKIP_CALIBRATE=1  skip the per-job machine-constant calibration
#   P7_MISSION_S is the name suffix the sf3 audit resolved ([A]: not recomputed, check on sf3 first).
#
# OUTPUT: results/expanse_policy/<jobid>_paperA_<J>/
#   plan.tsv                      the call/point plan (C and P lines) this job executed
#   c<NN>/points.tsv              points file of controller call NN (tag, case, mode, outdir, store, snapshot dir)
#   c<NN>/calls/<row>__<mode>{/,.log,_snap/,_node_before.txt}
#                                 per-point raw CSV (all repeats, final-state columns), log, per-solve snapshot
#                                 (calibration store, plan source, confirm state), node state before the point
#   c<NN>/controller.log          the controller's [points ...] START/END lines
#   constants/                    machine cost constants measured at the start of the job
#   summary_points.csv            ANALYSIS-READY: one row per point: job,row,case,mode,tag,repeats,n_timed,
#                                 median_r2plus_s (median of repeats 2..n), min_r2plus_s, max_r2plus_s, repeat1_s,
#                                 pos_norm_m, vel_norm_mps, mass_kg (final state, compare across modes of a row)
#   prelude_*.log                 the untimed prelude (precompile cache checks)
# Data only: no plots, no in-job analysis (SPACEAGORA_PPB_NO_PLOTS=1). Never commit results or logs.
# The controller fails a point whose log contains `Precompiling` (precompile guard): run precompile_tree.sh first.
# Quiet-node gate: julia_job.slurm (BENCH=1), unchanged.
#
# ESTIMATOR (DRY=1). Same model as the scope doc's cost.py/cost3.py, per job:
#   process cost = U(kind) + timed_sum(row, mode) x band;  band = 0.9 / 1.0 / 1.1 (low / central / high) [A]
#   U (untimed fixed cost of a fresh process, s; low / central / high): vac 72.6 / 76.8 / 131, gram 121 / 125 / 144,
#       mc_thr 94 / 94 / 105, mc_proc 191 / 196 / 210 [M for low and central of vac and gram, D for mc_*, A on highs;
#       paper_suite_scope.md Section 2]
#   timed_sum = sum of wall_time_s over the tiered repeats of that mode on the OLD tree 41223f282 (W.json of the
#       scope; P5f from results/expanse_paper/55249215_D), P7 = one third of the three-case total per case [A]
#   controller call = 33.9 / 33.9 / 49.5 s per call (package load before the call [M]; high adds the 15.6 s after-cost)
#   job overhead = 100 / 177 / 355 s (quiet-gate wait + untimed prelude; scope Section 2) [A low, D central, M high]
#   walltime = overhead + sum(process) + sum(call);  core-h = walltime x 128 / 3600 (a full node is billed)
# Calibration (calibrate_machine.jl, once per job) is NOT in the model: its duration is not in the scope doc.
# DRY=1 result (walltime h low/central/high; billed core-h): J1 2.16/2.35/3.24 (277/301/415), J2 2.37/2.59/2.92
# (304/332/374), J3 1.40/1.49/1.90 (179/190/243); three jobs 759/823/1032 core-h for the built plan. With the scope's own
# 13-process P5f grid (SCOPE_P5F=1 DRY=1) 808/874/1087 = the scope's 780/846/1046 plus 28/28/41, which is the 23
# controller calls (package load each) that the scope did not model. Row cells match the scope's to the printed digit.

set -u
cd "$(dirname "$0")/../../.." || exit 1
DRY="${DRY:-0}"
JSEL="${1:-}"
case "$JSEL" in J1|J2|J3) JOBS="$JSEL" ;; all) [ "$DRY" = 1 ] || { echo "all is DRY-only" >&2; exit 2; }; JOBS="J1 J2 J3" ;;
  *) echo "usage: paper_suite_planA.sh J1|J2|J3 (|all with DRY=1)" >&2; exit 2 ;; esac
JOB="${SLURM_JOB_ID:-local}"
BUDGET="${BUDGET:-32}"; WARMUP="${WARMUP:-3}"
REP_SERIAL="${REP_SERIAL:-5}"; REP_STATIC="${REP_STATIC:-7}"; REP_ADAPT="${REP_ADAPT:-9}"; REP_P7="${REP_P7:-11}"
P7S="${P7_MISSION_S:-5702}"
ROWS="${ROWS:-}"
P5F_MODES="${P5F_MODES:-}"
MODES="${MODES:-}"
PROFILE="${PROFILE:-full}"   # harness profile; smoke shortens every mission (functional checks only)
STORE=output/parallel_policy_state
GRAMDIR="data/GRAMSuite.jl/GRAM Suite 2.0"

say() { echo "[planA $(date -Iseconds)] $*"; }
run() { say "+ $*"; [ "$DRY" = 1 ] && return 0; "$@"; }

# ---- the workload table: id|job|case|samples|override tier|mode=timed_sum_s/U-kind ... (see ESTIMATOR) -------------
# samples: 1 = single simulation, N = --mc-samples=N, own = the case carries its own sample count (mcgrid).
rowtab() { cat <<'EOF'
p1_1sat|J1|gravity_1sat_l50_vacuum_4150000s|1|-|serial=88.9/vac outer_threads=124.7/vac inner_only=125.1/vac outer_inner_static=122.2/vac predictive=153.1/vac
p1_16sat|J1|gravity_16sat_l50_vacuum_514000s|1|-|serial=118.4/vac outer_threads=32.3/vac inner_only=31.6/vac outer_inner_static=32.5/vac predictive=42.4/vac
p1_64sat|J1|gravity_64sat_l50_vacuum_415000s|1|-|serial=115.5/vac outer_threads=44.6/vac inner_only=45.2/vac outer_inner_static=46.4/vac predictive=60.2/vac
p1_256sat|J1|gravity_256sat_l50_vacuum_124000s|1|-|serial=115.9/vac outer_threads=37.3/vac inner_only=37.0/vac outer_inner_static=35.8/vac predictive=50.0/vac
p1_1024sat|J1|gravity_1024sat_l50_vacuum_24600s|1|-|serial=101.2/vac outer_threads=24.3/vac inner_only=23.6/vac outer_inner_static=23.5/vac predictive=30.7/vac
p1_4096sat|J1|gravity_4096sat_l50_vacuum_5800s|1|-|serial=106.6/vac outer_threads=17.3/vac inner_only=16.6/vac outer_inner_static=18.4/vac predictive=24.3/vac
p6s_l20|J1|gravity_4096sat_l20_vacuum_19700s|1|-|serial=90.1/vac outer_threads=38.2/vac inner_only=41.3/vac outer_inner_static=47.8/vac predictive=54.2/vac
p6s_srp|J1|gravity_4096sat_l50_srp_nbody_vacuum_5800s|1|-|serial=126.7/vac outer_threads=60.6/vac inner_only=61.0/vac outer_inner_static=56.6/vac predictive=57.8/vac
p6s_expatm|J1|aero_4096sat_l50_expatm_100s|1|-|serial=91.9/vac outer_threads=43.2/vac inner_only=40.7/vac outer_inner_static=40.6/vac predictive=30.3/vac
p6s_gram|J1|aero_4096sat_l50_gram_lookahead_100s|1|-|serial=322.1/gram outer_threads=327.3/gram inner_only=205.2/gram outer_inner_static=205.5/gram predictive=219.3/gram
p7_vacuum|J3|gravity_1sat_l50_vacuum_@P7S@s|1|11|serial=3.5/vac outer_threads=3.5/vac inner_only=3.5/vac outer_inner_static=3.4/vac predictive=3.8/vac
p7_srp|J3|gravity_1sat_l50_srp_nbody_vacuum_@P7S@s|1|11|serial=3.5/vac outer_threads=3.5/vac inner_only=3.5/vac outer_inner_static=3.4/vac predictive=3.8/vac
p7_aero|J3|aero_1sat_l50_expatm_@P7S@s|1|11|serial=3.5/vac outer_threads=3.5/vac inner_only=3.5/vac outer_inner_static=3.4/vac predictive=3.8/vac
p3_indep256|J2|independent_1sat_1hr|256|-|serial=68.4/mc_thr outer_threads=38.4/mc_thr outer_process=6.5/mc_proc predictive=4.3/mc_proc
p4_mars32|J2|montecarlo_heavy_aerobraking|32|-|serial=212.7/mc_thr outer_threads=16.0/mc_thr outer_process=14.7/mc_proc predictive=12.5/mc_proc
p4g_marsgram32|J2|montecarlo_heavy_aerobraking_gram|32|-|outer_process=309.2/mc_proc predictive=306.3/mc_proc
p6p_gram4096|J2|aero_4096sat_l50_gram_process_100s|4096|-|serial=585.8/gram outer_process=37.1/mc_proc predictive=33.9/mc_proc
p4gs_marsgram_serial|J2|montecarlo_heavy_aerobraking_gram|32|-|serial=4971.7/gram
p5f_16sat_x8|J3|mcgrid_16sat_8mc_long|own|-|serial=267.9/mc_thr predictive=57.7/mc_proc outer_process=62.3/mc_proc outer_threads=59.9/mc_thr outer_inner_static=58.9/mc_thr force_threads@w8+l0+b4=96.1/mc_thr force_threads@w4+l0+b8=49.1/mc_thr force_threads@w2+l0+b16=45.8/mc_thr force_none@w1+l0+b32=44.3/mc_thr
p5f_8sat_x16|J3|mcgrid_8sat_16mc_long|own|-|serial=291.8/mc_thr predictive=31.9/mc_proc outer_process=35.1/mc_proc outer_threads=55.0/mc_thr outer_inner_static=55.1/mc_thr force_threads@w8+l0+b4=104.5/mc_thr force_threads@w4+l0+b8=53.8/mc_thr force_threads@w2+l0+b16=29.5/mc_thr force_none@w1+l0+b32=32.5/mc_thr force_threads@w16+l0+b2=104.5/mc_thr
EOF
}
# The scope's 13-process P5f grid (estimator check only; never executed).
rowtab_scope() { rowtab | awk -F'|' '$1 !~ /^p5f_/'; cat <<'EOF'
p5f_16sat_x8|J3|mcgrid_16sat_8mc_long|own|-|serial=267.9/mc_thr predictive=57.7/mc_proc outer_threads@w1=44.3/mc_thr outer_inner_static@w1=44.3/mc_thr outer_threads@w2=45.9/mc_thr outer_inner_static@w2=45.8/mc_thr outer_threads@w4=50.6/mc_thr outer_inner_static@w4=49.1/mc_thr outer_threads@w8=98.7/mc_thr outer_inner_static@w8=96.1/mc_thr outer_process@w8=62.3/mc_proc outer_process@w16=62.6/mc_proc outer_process@w32=63.9/mc_proc
p5f_8sat_x16|J3|mcgrid_8sat_16mc_long|own|-|serial=291.8/mc_thr predictive=31.9/mc_proc outer_threads@w1=32.1/mc_thr outer_inner_static@w1=32.5/mc_thr outer_threads@w2=29.2/mc_thr outer_inner_static@w2=29.5/mc_thr outer_threads@w4=54.1/mc_thr outer_inner_static@w4=53.8/mc_thr outer_threads@w8=104.6/mc_thr outer_inner_static@w8=104.5/mc_thr outer_process@w8=70.2/mc_proc outer_process@w16=35.1/mc_proc outer_process@w32=35.6/mc_proc
EOF
}
# Scope Plan A cells (core-h low|central|high, paper_suite_scope.md Section 4; the scope table's two P5f labels are
# swapped against the data: the 90.1 cell is the 16-sample case). P7 cells are the scope's 40.4/42.9/72.0 split in thirds.
scopecells() { cat <<'EOF'
p1_1sat 32.6 35.5 47.3
p1_16sat 21.1 22.8 33.3
p1_64sat 22.9 24.7 35.5
p1_256sat 21.7 23.5 34.1
p1_1024sat 19.4 20.9 31.2
p1_4096sat 18.8 20.2 30.4
p6s_l20 21.6 23.3 33.9
p6s_srp 24.5 26.5 37.5
p6s_expatm 20.8 22.4 32.9
p6s_gram 62.5 67.7 75.6
p3_indep256 24.0 24.8 27.0
p4_mars32 28.5 29.7 32.4
p4g_marsgram32 33.3 35.8 39.0
p6p_gram4096 38.9 41.7 45.7
p4gs_marsgram_serial 163.4 181.2 199.6
p7_vacuum 13.5 14.3 24.0
p7_srp 13.5 14.3 24.0
p7_aero 13.5 14.3 24.0
p5f_8sat_x16 86.2 90.1 98.9
p5f_16sat_x8 88.9 93.1 102.2
EOF
}

# ---- plan builder: table rows of one job -> C/P lines ---------------------------------------------------------------
#   C <call> <samples> <repeats> <npoints>
#   P <call> <tag> <row> <case> <mode> <timed_s> <ukind>
build_plan() {   # $1 = job; reads the table on stdin
  awk -F'|' -v job="$1" -v SER="$REP_SERIAL" -v STAT="$REP_STATIC" -v ADAPT="$REP_ADAPT" -v P7="$REP_P7" -v p7s="$P7S" \
      -v rows="$ROWS" -v p5f="$P5F_MODES" -v modes="$MODES" '
  function cls(m, ov,   m0) { if (ov != "-") return 4; m0 = m; sub(/@.*/, "", m0)
    if (m0 == "serial") return 1; if (m0 == "predictive" || m0 == "outer_process") return 3; return 2 }
  function rep(c) { return c == 1 ? SER : c == 2 ? STAT : c == 3 ? ADAPT : P7 }
  function emit(list, c,   rl, nr, k, r, j, nm, mm, t, ord, key, buf, npt) {   # one call: rows in list, tier class c
    nr = split(list, rl, " "); npt = 0; buf = ""
    for (k = 1; k <= nr; k++) { r = rl[k]; nm = 0
      for (j = 1; j <= nmod[r]; j++) if (cls(mode[r, j], ov[r]) == c) mm[++nm] = j
      for (j = 1; j <= nm; j++) { t = mm[((j - 1 + k - 1) % nm) + 1]
        npt++; key = mode[r, t]; gsub(/[@+]/, "_", key)
        buf = buf sprintf("P\t%d\t%s__%s\t%s\t%s\t%s\t%s\t%s\n", ncall + 1, id[r], key, id[r], cs[r], mode[r, t], tm[r, t], uk[r, t]) } }
    if (npt == 0) return
    ncall++; printf "C\t%d\t%s\t%d\t%d\n%s", ncall, smp[rl[1]], rep(c), npt, buf }
  function tiers(list, dir,   c) { if (dir == "asc") { for (c = 1; c <= 3; c++) emit(list, c) } else { for (c = 3; c >= 1; c--) emit(list, c) } }
  BEGIN { nsel = split(rows, sel, ",") }
  $2 != job { next }
  { if (nsel > 0 && rows != "") { ok = 0; for (i = 1; i <= nsel; i++) if (sel[i] == $1) ok = 1; if (!ok) next }
    n++; id[n] = $1; cs[n] = $3; gsub(/@P7S@/, p7s, cs[n]); smp[n] = $4; ov[n] = $5
    ml = $6
    flt = ($1 ~ /^p5f_/ && p5f != "") ? p5f : modes
    if (flt != "") { split(flt, want, ","); ml2 = ""; for (w in want) { for (z = 1; z <= split(ml, e, " "); z++) { split(e[z], a, "="); if (a[1] == want[w]) ml2 = ml2 " " e[z] } } ml = ml2 }
    nmod[n] = split(ml, me, " ")
    for (j = 1; j <= nmod[n]; j++) { split(me[j], a, "="); split(a[2], b, "/"); mode[n, j] = a[1]; tm[n, j] = b[1]; uk[n, j] = b[2] } }
  END {
    # pool: single-simulation rows without an override tier -> two halves; p7 group; own-sample rows
    for (r = 1; r <= n; r++) { if (smp[r] == "1" && ov[r] == "-") pool[++np] = r; else if (ov[r] != "-") p7g = p7g " " r; else own[++no] = r }
    half = int((np + 1) / 2); la = ""; lb = ""
    for (i = 1; i <= np; i++) { if (i <= half) la = la " " pool[i]; else lb = lb " " pool[i] }
    if (la != "") tiers(la, "asc"); if (lb != "") tiers(lb, "desc")
    if (p7g != "") emit(p7g, 4)
    for (i = 1; i <= no; i++) tiers(own[i], (i % 2 == 1) ? "asc" : "desc")
  }'
}
# order the output: C line then its P lines
plan_sorted() { awk -F'\t' '{ k = ($1 == "C") ? $2 * 2 : $2 * 2 + 1; printf "%d\t%s\n", k, $0 }' | sort -s -n -k1,1 | cut -f2-; }
job_plan() { # $1 = job -> plan on stdout
  if [ "${SCOPE_P5F:-0}" = 1 ]; then rowtab_scope; else rowtab; fi | build_plan "$1" | plan_sorted; }

# ---- estimator -----------------------------------------------------------------------------------------------------
estimate() {   # $1 = job label; plan on stdin; prints rows + job lines; last line "TOTAL low central high procs calls node_h_c node_h_h"
  local sc; sc=$(scopecells | awk '{printf "%s:%s:%s:%s;", $1, $2, $3, $4}')
  awk -F'\t' -v label="$1" -v sc="$sc" '
  BEGIN { U["vac",1] = 72.6; U["vac",2] = 76.8; U["vac",3] = 131; U["gram",1] = 121; U["gram",2] = 125; U["gram",3] = 144
          U["mc_thr",1] = 94; U["mc_thr",2] = 94; U["mc_thr",3] = 105; U["mc_proc",1] = 191; U["mc_proc",2] = 196; U["mc_proc",3] = 210
          tf[1] = 0.9; tf[2] = 1.0; tf[3] = 1.1; cc[1] = 33.9; cc[2] = 33.9; cc[3] = 49.5; jo[1] = 100; jo[2] = 177; jo[3] = 355
          n = split(sc, a, ";"); for (i = 1; i <= n; i++) { split(a[i], b, ":"); if (b[1] != "") { s1[b[1]] = b[2]; s2[b[1]] = b[3]; s3[b[1]] = b[4] } } }
  $1 == "C" { nc++ }
  $1 == "P" { np++; r = $4; if (!(r in seen)) { seen[r] = 1; ord[++nr] = r; nrp[r] = 0 } nrp[r]++
              for (s = 1; s <= 3; s++) { c = U[$8, s] + $7 * tf[s]; rc[r, s] += c; tot[s] += c } }
  END { printf "  %-22s %5s  %-26s  %-26s\n", "row", "procs", "model core-h low/central/high", "scope cell low/central/high"
        for (i = 1; i <= nr; i++) { r = ord[i]
          printf "  %-22s %5d  %7.1f %7.1f %7.1f      ", r, nrp[r], rc[r,1] * 128 / 3600, rc[r,2] * 128 / 3600, rc[r,3] * 128 / 3600
          if (r in s1) printf "%7.1f %7.1f %7.1f\n", s1[r], s2[r], s3[r]; else printf "     (no scope cell)\n"
          for (s = 1; s <= 3; s++) if (r in s1) { sm[s] += (s == 1 ? s1[r] : s == 2 ? s2[r] : s3[r]) } }
        for (s = 1; s <= 3; s++) { w[s] = jo[s] + tot[s] + nc * cc[s]; ch[s] = w[s] * 128 / 3600 }
        printf "  %s: %d processes, %d controller calls\n", label, np, nc
        printf "  %s model-only sum of rows (no controller calls, no job overhead): %.1f / %.1f / %.1f core-h; scope cells for these rows: %.1f / %.1f / %.1f\n", label, \
          tot[1] * 128 / 3600, tot[2] * 128 / 3600, tot[3] * 128 / 3600, sm[1], sm[2], sm[3]
        printf "  %s walltime %.2f / %.2f / %.2f h (overhead %d/%d/%d s, calls %.0f/%.0f/%.0f s); billed %.1f / %.1f / %.1f core-h\n", label, \
          w[1] / 3600, w[2] / 3600, w[3] / 3600, jo[1], jo[2], jo[3], nc * cc[1], nc * cc[2], nc * cc[3], ch[1], ch[2], ch[3]
        if (w[3] > 3.6 * 3600) printf "  WARNING %s: high-bound walltime above 3.6 h (4 h cap)\n", label
        printf "TOTAL\t%.2f\t%.2f\t%.2f\t%d\t%d\n", ch[1], ch[2], ch[3], np, nc }'
}

# ---- DRY: estimate (and show the plan) for each selected job -----------------------------------------------------------
dry_job() {
  local j=$1 plan; plan=$(job_plan "$j")
  say "=== plan $j: points in execution order (call | repeats | samples | row | mode | case) ==="
  echo "$plan" | awk -F'\t' '$1 == "C" { printf "  call %02d  samples=%s repeats=%d points=%d\n", $2, $3, $4, $5 } $1 == "P" { printf "     %-3d %-22s %-32s %s\n", ++k, $4, $6, $5 }'
  say "=== estimate $j ==="
  local est; est=$(echo "$plan" | estimate "$j")
  echo "$est" | grep -v '^TOTAL'
  echo "$est" | grep '^TOTAL' | sed "s/^/$j\t/" >> "$DRYTOT"
}

# ---- the job -----------------------------------------------------------------------------------------------------------
run_job() {
  local j=$1 OUT="$PWD/results/expanse_policy/${JOB}_paperA_$1" plan
  plan=$(job_plan "$j")
  [ -n "$plan" ] || { say "ERROR: empty plan for $j (ROWS filter?)"; exit 2; }
  export JULIA_NUM_THREADS="" SPACEAGORA_PERF_PROCS="" SPACEAGORA_PPB_WORKLOAD=0
  unset JULIA_NUM_GC_THREADS JULIA_THREAD_SLEEP_THRESHOLD SPACEAGORA_RHS_CALIBRATE SPACEAGORA_RHS_CALIBRATE_SEED_TOPK \
        SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER SPACEAGORA_POLICY_REGRET_GUARD \
        SPACEAGORA_POLICY_REGRET_MARGIN SPACEAGORA_POLICY_REGRET_MIN_SOLVES SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR \
        SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN SPACEAGORA_RHS_REGION_WIDTH SPACEAGORA_RHS_CALIBRATE_CONFIRM \
        SPACEAGORA_RHS_CALIBRATE_CONFIRM_SOLVES SPACEAGORA_RHS_CALIBRATE_CONFIRM_MARGIN
  export SPACEAGORA_PPB_SHARED_PROCESS=0 SPACEAGORA_PPB_WARMUP_SPAN=0.1 SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto
  export SPACEAGORA_CAMPAIGN_CORRECTIONS=off SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1 SPACEAGORA_PPB_NO_PLOTS=1
  export SPACEAGORA_PPC_EQUAL_CORE_BUDGET=1 OPENBLAS_NUM_THREADS=1
  say "job=$JOB part=$j host=$(hostname) commit=${GIT_COMMIT:-unknown} DRY=$DRY budget=$BUDGET warmup=$WARMUP repeats=$REP_SERIAL/$REP_STATIC/$REP_ADAPT/$REP_P7 rows=[${ROWS:-all}]"
  # GRAM prerequisites (read-only checks): J1 has the Earth GRAM look-ahead row, J2 MarsGRAM and Earth GRAM rows.
  if [ "$j" != J3 ]; then
    local need=("$GRAMDIR/Build/lib/libGRAM.so" "$GRAMDIR/SPICE/pck/pck00011.tpc" "$GRAMDIR/Earth/data")
    [ "$j" = J2 ] && need+=("$GRAMDIR/Mars/data")
    for p in "${need[@]}"; do
      if [ -e "$p" ]; then say "gram prerequisite present: $p"
      elif [ "$DRY" = 1 ]; then say "gram prerequisite (not checked in DRY): $p"
      else say "ERROR: gram prerequisite missing: $p (data store link or native build)"; exit 4; fi
    done
  fi
  [ "$DRY" = 1 ] || { mkdir -p "$OUT"; echo "$plan" > "$OUT/plan.tsv"
    say "memory: $(grep -E 'MemTotal|MemAvailable' /proc/meminfo | tr -s ' ' | tr '\n' ' ')"; }
  # Untimed prelude: cache checks only (precompile_tree.sh has done the work on the shared partition).
  if [ "$DRY" = 1 ]; then run julia --project=examples -e 'using Pkg; Pkg.instantiate()'
  else
    say "+ julia --project=examples -e 'using Pkg; Pkg.instantiate()'"
    julia --project=examples -e 'using Pkg; Pkg.instantiate()' || { say "ERROR: examples instantiate failed (run precompile_tree.sh first)"; exit 7; }
  fi
  local proj pl
  for proj in . examples; do
    pl="$OUT/prelude_${proj//./root}.log"
    if [ "$DRY" = 1 ]; then run julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA'; continue; fi
    say "+ julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA'"
    julia --project=$proj --threads=$BUDGET -e 'using SpaceAGORA' 2>&1 | tee "$pl"
    [ "${PIPESTATUS[0]}" -eq 0 ] || { say "ERROR: prelude --project=$proj failed"; exit 6; }
    grep -q Precompiling "$pl" && say "WARNING: --project=$proj precompiled on the exclusive node; run precompile_tree.sh on the shared partition first next time"
  done
  if [ "$j" != J3 ]; then   # GRAM must load through the vendored path; this also settles the first-touch Precompiling line
    local gl='using SpaceAGORA; SpaceAGORA.RuntimeServices._prepare_gramsuite_load_path!(joinpath(pwd(), "data", "GRAMSuite.jl"); prepend=true); @eval import GRAMSuite; ext = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing; println("GRAMSuite loaded, extension attached: ", ext); ext || exit(4)'
    if [ "$DRY" = 1 ]; then run julia --project=. --threads=$BUDGET -e "$gl"
    else
      say "+ julia --project=. --threads=$BUDGET -e '<load GRAMSuite via the vendored load path>'"
      julia --project=. --threads=$BUDGET -e "$gl" 2>&1 | tee "$OUT/prelude_gram.log"
      [ "${PIPESTATUS[0]}" -eq 0 ] || { say "ERROR: GRAMSuite did not load"; exit 3; }
      grep -q Precompiling "$OUT/prelude_gram.log" && say "WARNING: GRAM extension precompiled on the exclusive node; add it to precompile_tree.sh"
    fi
  fi
  # Machine cost constants at the budget, measured once per job, put back into the emptied policy store before every point.
  local constants="$OUT/constants"
  if [ "${SKIP_CALIBRATE:-0}" != 1 ]; then
    run mkdir -p "$constants"
    run env SPACEAGORA_CORE_BUDGET=$BUDGET SPACEAGORA_PPC_PHYSICAL_CORES=$BUDGET julia --project=. --threads=$BUDGET scripts/calibration/calibrate_machine.jl
    if [ "$DRY" != 1 ]; then
      cp -p "$STORE"/cost_constants_*.toml "$constants"/ || { say "ERROR: calibration wrote no constants"; exit 5; }
    fi
  fi
  export SPACEAGORA_CORE_BUDGET=$BUDGET SPACEAGORA_PPC_PHYSICAL_CORES=$BUDGET
  # Per point (controller): fresh policy-state dir and snapshot dir, constants put back, no calibration store, node snapshot.
  export SPACEAGORA_PPC_POINT_PREP="rm -rf '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; mkdir -p '$STORE' \"\$SPACEAGORA_PPB_SNAPSHOT_DIR\"; cp -p '$constants'/cost_constants_*.toml '$STORE'/ 2>/dev/null; rm -f \"\$SPACEAGORA_RHS_CALIBRATION_PATH\"; bash '$PWD/scripts/remote/expanse/node_snapshot.sh' \"\${SPACEAGORA_PPB_SNAPSHOT_DIR%_snap}_node_before.txt\""
  [ "$DRY" = 1 ] && say "prep: $SPACEAGORA_PPC_POINT_PREP"
  local FAILED=() c ncalls cdir pts line smp reps cmd rc t cs m tag
  ncalls=$(echo "$plan" | awk -F'\t' '$1 == "C" { n++ } END { print n + 0 }')
  for ((c = 1; c <= ncalls; c++)); do
    line=$(echo "$plan" | awk -F'\t' -v c="$c" '$1 == "C" && $2 == c')
    smp=$(echo "$line" | cut -f3); reps=$(echo "$line" | cut -f4)
    cdir="$OUT/c$(printf %02d "$c")"
    # P lines: P call tag row case mode timed uk  -> read gives tag=$3, case=$5, mode=$6
    pts=""
    while IFS=$'\t' read -r _ _ tag _ cs m _ _; do
      pts+="${tag}\t${cs}\t${m}\t$cdir/calls/${tag}\t$cdir/stores/${tag}.toml\t$cdir/calls/${tag}_snap\n"
    done < <(echo "$plan" | awk -F'\t' -v c="$c" '$1 == "P" && $2 == c { print }')
    cmd=(julia --project=examples benchmarks/studies/parallelization_performance.jl "$PROFILE"
         --points-file="$cdir/points.tsv" --threads=$BUDGET --process-workers=$BUDGET
         --repeats="$reps" --warmup=$WARMUP --parity-cases=none)
    case "$smp" in 1|own) ;; *) cmd+=(--mc-samples="$smp") ;; esac
    say "=== call $c/$ncalls: samples=$smp repeats=$reps"
    if [ "$DRY" = 1 ]; then
      say "points file $cdir/points.tsv:"; printf "$pts" | sed 's/^/  /'; echo "  ${cmd[*]}"; continue
    fi
    mkdir -p "$cdir/stores" "$cdir/calls"; printf "$pts" > "$cdir/points.tsv"
    "${cmd[@]}" 2>&1 | tee "$cdir/controller.log"; rc=${PIPESTATUS[0]}
    [ $rc -eq 0 ] || { say "FAILED call $c (rc=$rc); see $cdir/calls/*.log"; FAILED+=("c$c"); }
  done
  [ "$DRY" = 1 ] && { say "done (dry run): $OUT"; return 0; }
  summarize "$j" "$OUT"
  [ ${#FAILED[@]} -eq 0 ] || { say "FAILED: ${FAILED[*]}"; exit 1; }
  say "done: $OUT"
}

# summary_points.csv: median of repeats 2..n per point (repeat 1 is kept separately), plus the final-state triple.
summarize() {   # $1 job, $2 out
  local f="$2/summary_points.csv" row tag d
  echo "job,row,case,mode,tag,repeats,n_timed,median_r2plus_s,min_r2plus_s,max_r2plus_s,repeat1_s,pos_norm_m,vel_norm_mps,mass_kg" > "$f"
  while IFS=$'\t' read -r _ _ tag row _ _ _ _; do
    d=$(ls -d "$2"/c*/calls/"$tag" 2>/dev/null | head -1); [ -n "$d" ] || { echo "$1,$row,,,$tag,0,0,,,,,,," >> "$f"; continue; }
    cat "$d"/parallelization_performance_raw_*.csv 2>/dev/null | awk -F, -v job="$1" -v row="$row" -v tag="$tag" '
      $1 == "timestamp_utc" { for (i = 1; i <= NF; i++) h[$i] = i; next }
      $(h["success"]) == "true" { r = $(h["repeat"]) + 0; w = $(h["wall_time_s"]) + 0; cs = $(h["case"]); md = $(h["mode"]); nrep++
        if (r == 1) r1 = w; if (r >= 2) v[++n] = w
        pos = $(h["final_primary_pos_norm_m"]); vel = $(h["final_primary_vel_norm_mps"]); ms = $(h["final_primary_mass_kg"]) }
      END { for (i = 2; i <= n; i++) { x = v[i]; k = i - 1; while (k >= 1 && v[k] > x) { v[k + 1] = v[k]; k-- } v[k + 1] = x }
            med = (n == 0) ? "" : (n % 2 ? v[(n + 1) / 2] : (v[n / 2] + v[n / 2 + 1]) / 2)
            printf "%s,%s,%s,%s,%s,%d,%d,%s,%s,%s,%s,%s,%s,%s\n", job, row, cs, md, tag, nrep, n, med, (n ? v[1] : ""), (n ? v[n] : ""), r1, pos, vel, ms }' >> "$f"
  done < <(echo "$plan" | awk -F'\t' '$1 == "P"')
  say "summary: $f ($(($(wc -l < "$f") - 1)) points)"
}

DRYTOT=$(mktemp)
for j in $JOBS; do
  run_job "$j"
  [ "$DRY" = 1 ] && dry_job "$j"
done
if [ "$DRY" = 1 ] && [ "$JSEL" = all ]; then
  say "=== grand total (three jobs, estimator above) vs the scope's Plan A 780 / 846 / 1046 core-h ==="
  awk -F'\t' '{ l += $3; c += $4; h += $5; p += $6; k += $7 } END { printf "  built plan: %.0f / %.0f / %.0f core-h, %d processes, %d controller calls\n", l, c, h, p, k }' "$DRYTOT"
fi
rm -f "$DRYTOT"
