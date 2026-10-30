"""Stage M analysis: apoapsis-targeting controller. Plan (frozen, with Addenda A and B): tasks/2026-10-15_periapsis-region-border-targeting/apoapsis_controller_plan.md.
Reuses stagei_analysis (tiers, run_metrics) and stagef20_analysis.boot_ci (cell-clustered bootstrap, 10,000 draws, seed 20261020, percentiles 5 and 95, read from stagei_config.toml).
Tiers are the paper's, at the periapsis event (q_eff): low < 0.05, soft 0.25-0.30, medium 0.30-0.45, hard > 0.45 W/cm2. error = r_a,end - r_a,t, signed, before the trim.
  selftest                         synthetic tables and a synthetic file tree; every criterion function is checked against exact cases
  integrity   --results-root R     route/failures per campaign, GRAM epoch, calibration convergence, NaN scan, v5 exceptions (dev: 0 required), horizon hits
  cal         --results-root R     T4 fit-only calibration report: implied first-pass Q, cells below the 0.05 floor, polar cells, flown first-pass Q next to it
  probe       --results-root R     probe runs: passes, wall per pass (contended; not a performance number), v5 exceptions
  held        --results-root R     Mars criteria P1-P3 (S-E), P2/S1/D1/C1 (T4), reported A3-J2, tier counts, error distribution, dV, days, reported-only items
  gen         --results-root R [--planets mars,earth,venus]   G1-G3 at Earth and Venus, and the same rules on Mars broad
  report      --results-root R     integrity + cal + held + gen, one JSON
Layout expected from the Julia side (WP3/WP4), below <R>/stagem/:
  summary_<phase>_<host>.csv, mc_<phase>_<host>.toml (mc_route, mc_failed, mc_samples)
  summary columns: phase, set (se|se_j2|t4|broad|g), planet, arm (M|B1|B2|B3|B4), tag, seed_idx (0 when one seed), split (dev|heldout), run_id, status, retcode, ra_end_m, rp_end_m,
                   err_m, dv_flown, dv_raise, dv_trim, t_end_days, burns, horizon_hit, v5_exceptions, wall_s (optional), cd_factor (T4, optional), mdl_stopped_at_target (optional)
  <set>/<planet>/cells_base.csv (tag, split, inc, raan, omega, epoch_idx, epoch_utc), <set>/<planet>/runs/<run_id>/{passes.csv,plan.csv}
  passes.csv columns: pass, time_s, q_eff_w_cm2, apo_alt_m, phase (optional: main|penultimate|last)
  plan.csv optional columns for the decay-map check: a_m, obs_drop_m, pred_drop_m; Pi_g (last row)
  <set>/<planet>/cal/<tag>/guess.csv columns: status, q_first_w_cm2, q_target_w_cm2, gram_core_epoch, h_ic_m, h_guess_m, hloc_at_hg_m, peri_lat_deg
Usage: python3 -I stagem_analysis.py --phase ... --results-root DIR"""
import argparse, glob, json, pathlib, sys, tempfile, tomllib
from statistics import NormalDist
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagei_analysis as sia
import stagef20_analysis as sf

AN = sia.AN; Q_LOW, Q_HIGH, Q_HARD, Q_TGT = sia.Q_LOW, sia.Q_HIGH, sia.Q_HARD, sia.Q_TGT
SUBDIR = "stagem"
# Frozen thresholds (plan section 5, A.8, D8; ASSUMED unless the plan sources them)
W10, W5 = 10000.0, 5000.0                 # precision windows, m (D7: absolute at every planet)
P1_FRAC10, P1_WILSON5 = 0.90, 0.38        # P1; 0.38 is the paper PR-DRL value (SOURCED)
N_CELLS_10, N_CELLS_5 = 22, 12            # P2 / G1 counts of 24 cells
C1_RATIO, G3_RATIO = 1.10, 1.10
POLAR_LAT_DEG = 75.0                      # ASSUMED: "polar" cells of the T4 cal report (Addendum B: roughly periapsis above about 75 deg N)
ZW = NormalDist().inv_cdf(0.95)           # Wilson 90% (two-sided) bound
ARMS_REPORT = ("M", "B1", "B2", "B3", "B4")
ENDGAME = ("penultimate", "last")


def sdir(root): return pathlib.Path(root) / SUBDIR
def setdir(root, st, pl): return sdir(root) / st / pl


def wilson_lower(k, n, z=ZW):
    if n == 0: return float("nan")
    p = k / n; d = 1 + z * z / n
    return float((p + z * z / (2 * n) - z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d)


def summaries(root):
    fs = sorted(glob.glob(str(sdir(root) / "summary_*.csv"))); return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def walk_in(q):
    """Leading consecutive low-tier passes (ASSUMED reading of 'low-tier passes before the first burn': the paper treats walk-in as out of scope)."""
    q = np.asarray(q, float); nz = np.flatnonzero(q >= Q_LOW); return int(nz[0]) if len(nz) else len(q)


def one_run(root, row):
    """Metrics of one flown run from its summary row and passes.csv (None when missing). total dV = flown + raise + trim (D6)."""
    rd = setdir(root, row["set"], row["planet"]) / "runs" / str(row["run_id"])
    if not (rd / "passes.csv").exists(): return None
    P = pd.read_csv(rd / "passes.csv"); dec = pd.read_csv(rd / "plan.csv") if (rd / "plan.csv").exists() else pd.DataFrame()
    if len(P) == 0: return None
    r = sia.run_metrics(P, dict(row), dec); q = P.q_eff_w_cm2.values
    ph = P["phase"].astype(str).values if "phase" in P else np.array(["main"] * len(P)); lowm = q < Q_LOW
    r.update(low_planned=int((lowm & np.isin(ph, ENDGAME)).sum()), walk_in=walk_in(q), q_first=float(q[0]), q_last=float(q[-1]))
    if len(dec) and {"a_m", "obs_drop_m", "pred_drop_m"} <= set(dec.columns):
        ok = (dec.pred_drop_m > 0) & np.isfinite(dec.obs_drop_m); r["decay_ratios"] = (dec.obs_drop_m[ok] / dec.pred_drop_m[ok]).tolist(); r["decay_a_m"] = dec.a_m[ok].tolist()
    return r


def load_runs(root, sets=None):
    """One row per ok run: summary fields + run metrics + cell geometry. Cells' omega/inc come from cells_base.csv."""
    S = summaries(root)
    if len(S) == 0: return pd.DataFrame()
    S = S[S.status == "ok"] if "status" in S else S
    rows = []; geo = {}
    for _, s in S.iterrows():
        if sets is not None and s["set"] not in sets: continue
        m = one_run(root, s)
        if m is None: continue
        k = (s["set"], s["planet"])
        if k not in geo:
            f = setdir(root, *k) / "cells_base.csv"; geo[k] = pd.read_csv(f).set_index("tag") if f.exists() else pd.DataFrame()
        g = geo[k]; om = float(g.loc[s["tag"], "omega"]) if len(g) and s["tag"] in g.index else np.nan; inc = float(g.loc[s["tag"], "inc"]) if len(g) and s["tag"] in g.index else np.nan
        sd = dict(s); sd.update(m); sd["dv_total"] = float(s["dv_flown"] + s["dv_raise"] + s["dv_trim"]); sd["omega"] = om; sd["inc"] = inc; sd["seed_idx"] = int(s.get("seed_idx", 0) or 0)
        rows.append(sd)
    return pd.DataFrame(rows)


def complete(R, arms):
    """Missing-data rule: a cell missing any arm (or any seed of any arm) is dropped from every arm."""
    R = R[R.arm.isin(arms)]; keep = []
    for t, G in R.groupby("tag"):
        ks = [set(G[G.arm == a].seed_idx) for a in arms]
        if all(ks) and all(k == ks[0] for k in ks): keep.append(t)
    return R[R.tag.isin(keep)]


def sel(R, st, pl, split="heldout"):
    return R[(R["set"] == st) & (R.planet == pl) & (R.split == split)] if len(R) else R


def within(R, w): return (R.err_m.abs() <= w).values & np.isfinite(R.err_m.values)


def cell_vals(R, arm, f):
    """Per-cell mean over seeds of f(run) for one arm, indexed by tag (the cluster is the cell)."""
    A = R[R.arm == arm].copy(); A["_v"] = f(A); return A.groupby("tag")["_v"].mean()


def paired_cell_diff(R, a, b, f):
    d = cell_vals(R, a, f) - cell_vals(R, b, f); return d.dropna().values


def arm_stats(R):
    n = len(R); npass = float(R.n_pass.sum())
    if n == 0: return dict(n=0)
    e = R.err_m.values; w10 = within(R, W10); w5 = within(R, W5)
    return dict(n=n, cells=int(R.tag.nunique()), frac10=float(w10.mean()), frac5=float(w5.mean()), n10=int(w10.sum()), n5=int(w5.sum()), mean_abs_err_m=float(np.nanmean(np.abs(e))), signed_mean_m=float(np.nanmean(e)),
                signed_median_m=float(np.nanmedian(e)), n_over=int((e > 0).sum()), n_under=int((e < 0).sum()), passes=int(npass), low=int(R.low.sum()), low_planned_endgame=int(R.low_planned.sum()),
                low_other=int(R.low.sum() - R.low_planned.sum()), walk_in=int(R.walk_in.sum()), soft=int(R.soft.sum()), medium=int(R.medium.sum()), hard=int(R.hard.sum()), over=int(R.over.sum()),
                frac_over=float(R.over.sum() / npass), dv_flown_mean=float(R.dv_flown.mean()), dv_raise_mean=float(R.dv_raise.mean()), dv_trim_mean=float(R.dv_trim.mean()), dv_total_mean=float(R.dv_total.mean()),
                days_mean=float(R.t_end_days.mean()), days_median=float(R.t_end_days.median()), q_last_median=float(R.q_last.median()), horizon_hits=int(R.get("horizon_hit", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()))


def breakdowns(R):
    """Reported only: dV and precision per AOP quadrant and per inclination (input to the initial-condition stage)."""
    out = {}
    if len(R) == 0 or R.omega.isna().all(): return out
    R = R.copy(); R["aop_quadrant"] = (np.mod(R.omega, 360.0) // 90).astype(int); R["inc_deg"] = R.inc.round(0)
    for k in ("aop_quadrant", "inc_deg"):
        out[k] = {str(g): dict(n=len(G), frac5=float(within(G, W5).mean()), frac10=float(within(G, W10).mean()), dv_total_mean=float(G.dv_total.mean())) for g, G in R.groupby(k)}
    return out


def decay_map(R):
    r = [x for L in R.get("decay_ratios", pd.Series(dtype=object)).dropna() for x in L]; a = [x for L in R.get("decay_a_m", pd.Series(dtype=object)).dropna() for x in L]
    if not r: return None
    r = np.array(r); a = np.array(a); out = dict(n=len(r), median_obs_over_pred=float(np.median(r)), frac_within_5pct=float((np.abs(r - 1) <= 0.05).mean())); bins = [0, 5e6, 8e6, 12e6, 1e12]
    out["by_a"] = {f"{bins[i] / 1e6:g}-{bins[i + 1] / 1e6:g} Mm": dict(n=int(m.sum()), median=float(np.median(r[m])), frac_within_5pct=float((np.abs(r[m] - 1) <= 0.05).mean())) for i in range(4) for m in [(a >= bins[i]) & (a < bins[i + 1])] if m.any()}
    return out


def vs_other(R, arm, rng):
    """Reported: M against another arm, paired by cell."""
    d5 = paired_cell_diff(R, "M", arm, lambda A: within(A, W5).astype(float)); dv = paired_cell_diff(R, "M", arm, lambda A: A.dv_total)
    dd = paired_cell_diff(R, "M", arm, lambda A: A.over / A.n_pass); return dict(dfrac5=dict(mean=float(np.mean(d5)), ci90=sf.boot_ci(d5, rng, AN["bootstrap_draws"])), ddv_total=dict(mean=float(np.mean(dv)), ci90=sf.boot_ci(dv, rng, AN["bootstrap_draws"])),
                                                                                  dover_frac=dict(mean=float(np.mean(dd)), ci90=sf.boot_ci(dd, rng, AN["bootstrap_draws"])))


def crit_P1(RM):
    n = len(RM); k10 = int(within(RM, W10).sum()); k5 = int(within(RM, W5).sum()); wl = wilson_lower(k5, n)
    return dict(n=n, frac10=k10 / n if n else np.nan, k5=k5, wilson90_lower_frac5=wl, passed=bool(n and k10 / n >= P1_FRAC10 and wl > P1_WILSON5))


def crit_P2(RM):
    k10 = int(within(RM, W10).sum()); k5 = int(within(RM, W5).sum()); return dict(n=len(RM), k10=k10, k5=k5, passed=bool(k10 >= N_CELLS_10 and k5 >= N_CELLS_5))


def crit_P3(R, rng):
    d = paired_cell_diff(R, "M", "B1", lambda A: within(A, W5).astype(float)); lo, hi = sf.boot_ci(d, rng, AN["bootstrap_draws"])
    return dict(n_cells=len(d), mean=float(np.mean(d)) if len(d) else np.nan, ci90=[lo, hi], passed=bool(np.isfinite(lo) and lo > 0))


def crit_S1(R):
    M, B = R[R.arm == "M"], R[R.arm == "B1"]; hm, hb = int(M.hard.sum()), int(B.hard.sum()); fm, fb = M.over.sum() / M.n_pass.sum(), B.over.sum() / B.n_pass.sum()
    return dict(hard_M=hm, hard_B1=hb, frac_over_M=float(fm), frac_over_B1=float(fb), passed=bool(hm <= hb and fm <= fb))


def crit_D1(R, rng):
    d = paired_cell_diff(R, "M", "B1", lambda A: A.dv_total); lo, hi = sf.boot_ci(d, rng, AN["bootstrap_draws"])
    return dict(n_cells=len(d), mean_mps=float(np.mean(d)) if len(d) else np.nan, ci90=[lo, hi], passed=bool(np.isfinite(hi) and hi < 0))


def crit_C1(R):
    m = R[R.arm == "M"].groupby("tag").t_end_days.mean(); b = R[R.arm == "B1"].groupby("tag").t_end_days.mean(); r = (m / b).dropna()
    return dict(n_cells=len(r), median_ratio=float(r.median()) if len(r) else np.nan, passed=bool(len(r) and r.median() <= C1_RATIO))


def crit_G(R, rng):
    """G1 (as P2), G2 (Stage I J1: paired over fraction lower bound not > 0, hard <= B1), G3 (mean total dV M <= 1.10 x B1)."""
    G1 = crit_P2(R[R.arm == "M"]); d = paired_cell_diff(R, "M", "B1", lambda A: A.over / A.n_pass); lo, hi = sf.boot_ci(d, rng, AN["bootstrap_draws"])
    hm, hb = int(R[R.arm == "M"].hard.sum()), int(R[R.arm == "B1"].hard.sum()); G2 = dict(over_diff_mean=float(np.mean(d)) if len(d) else np.nan, ci90=[lo, hi], hard_M=hm, hard_B1=hb, passed=bool(not (lo > 0 or hm > hb)))
    mM, mB = R[R.arm == "M"].dv_total.mean(), R[R.arm == "B1"].dv_total.mean(); G3 = dict(dv_M=float(mM), dv_B1=float(mB), passed=bool(mM <= G3_RATIO * mB))
    return dict(G1=G1, G2=G2, G3=G3, all=bool(G1["passed"] and G2["passed"] and G3["passed"]))


def set_report(R, arms):
    out = dict(arms={a: arm_stats(R[R.arm == a]) for a in arms if (R.arm == a).any()}, decay_map=decay_map(R[R.arm == "M"]) if (R.arm == "M").any() else None, breakdowns_M=breakdowns(R[R.arm == "M"]))
    return out


def phase_held_tables(R, rng):
    """Mars criteria; R is the full run table. Returns the result dict."""
    res = {}
    se = complete(sel(R, "se", "mars"), ["M", "B1"])
    if len(se):
        res["S-E"] = dict(P1=crit_P1(se[se.arm == "M"]), P3=crit_P3(se, rng), **set_report(se, ("M", "B1")))
    t4 = complete(sel(R, "t4", "mars"), ["M", "B1"])
    if len(t4):
        res["T4"] = dict(P2=crit_P2(t4[t4.arm == "M"]), S1=crit_S1(t4), D1=crit_D1(t4, rng), C1=crit_C1(t4), **set_report(sel(R, "t4", "mars"), ARMS_REPORT))
        for a in ("B3",):
            t = complete(sel(R, "t4", "mars"), ["M", a])
            if len(t): res["T4"]["vs_" + a] = vs_other(t, a, rng)
    br = sel(R, "broad", "mars")
    if len(br):
        c = complete(br, ["M", "B1"]); res["Mars_broad"] = dict(**crit_G(c, rng), **set_report(br, ARMS_REPORT))
        for a in ("B3", "B4"):
            t = complete(br, ["M", a])
            if len(t): res["Mars_broad"]["vs_" + a] = vs_other(t, a, rng)
    j2 = complete(sel(R, "se_j2", "mars"), ["M"])
    if len(j2) and "S-E" in res:                       # A3-J2, reported only: J2-only M minus 20x20 M, paired by (cell, seed)
        M20 = se[se.arm == "M"]; common = sorted(set(j2.tag) & set(M20.tag))
        if common:
            f5 = lambda A: within(A, W5).astype(float); f10 = lambda A: within(A, W10).astype(float)
            d5 = (cell_vals(j2[j2.tag.isin(common)], "M", f5) - cell_vals(M20[M20.tag.isin(common)], "M", f5)).values; d10 = (cell_vals(j2[j2.tag.isin(common)], "M", f10) - cell_vals(M20[M20.tag.isin(common)], "M", f10)).values
            res["A3_J2"] = dict(n_cells=len(common), dfrac5=dict(mean=float(d5.mean()), ci90=sf.boot_ci(d5, rng, AN["bootstrap_draws"])), dfrac10=dict(mean=float(d10.mean()), ci90=sf.boot_ci(d10, rng, AN["bootstrap_draws"])),
                                abs_le_0p05=bool(abs(d5.mean()) <= 0.05), reported_only=True, stats_j2=arm_stats(j2[j2.arm == "M"]))
    k = ("P1", "S-E"), ("P3", "S-E"), ("P2", "T4"), ("S1", "T4"), ("D1", "T4"), ("C1", "T4")
    ids = {f"{c}": res.get(s, {}).get(c, {}).get("passed") for c, s in k}
    res["verdict"] = "TARGETING_PASS" if all(v is True for v in ids.values()) else "FAILED(" + ",".join(c for c, v in ids.items() if v is not True) + ")"; res["criteria_flags"] = ids
    return res


def phase_gen_tables(R, rng, planets):
    res = {}
    for pl in planets:
        if pl == "mars": continue
        g = sel(R, "g", pl)
        if len(g) == 0: continue
        res[pl] = dict(**crit_G(complete(g, ["M", "B1"]), rng), **set_report(g, ("M", "B1")))
    done = [res[p]["all"] for p in ("earth", "venus") if p in res]
    res["verdict"] = "GENERALIZES" if len(done) == 2 and all(done) else ("INCOMPLETE" if len(done) < 2 else "DOES_NOT_GENERALIZE")
    return res


def cal_report(G, cols=None):
    """T4 fit-only calibration. Implied first-pass Q = Q_tgt exp((h_guess - h_ic)/H_local) (Addendum B); counts below the 0.05 floor, above 0.25, and polar-latitude cells."""
    q = Q_TGT * np.exp((G.h_guess_m - G.h_ic_m) / G.hloc_at_hg_m); lat = G.peri_lat_deg; pol = lat >= POLAR_LAT_DEG; below = q < Q_LOW
    out = dict(n=len(G), q_implied_min=float(q.min()), q_implied_median=float(q.median()), q_implied_max=float(q.max()), n_below_floor=int(below.sum()), frac_below_floor=float(below.mean()), n_above_0p25=int((q > Q_HIGH).sum()),
               n_polar=int(pol.sum()), polar_lat_deg=POLAR_LAT_DEG, n_polar_below_floor=int((pol & below).sum()), n_nonpolar_below_floor=int((~pol & below).sum()), lat_min=float(lat.min()), lat_max=float(lat.max()))
    out["by_lat_bin"] = {lab: dict(n=int(m.sum()), n_below_floor=int((m & below).sum())) for lab, m in (("<70", lat < 70), ("70-80", (lat >= 70) & (lat < 80)), (">=80", lat >= 80)) if m.any()}
    return out, pd.Series(q.values, index=G.tag.values)


def phase_cal(root):
    fs = sorted(glob.glob(str(setdir(root, "t4", "mars") / "cal" / "*" / "guess.csv")))
    if not fs: print("no T4 cal tables"); return {}
    G = pd.concat([pd.read_csv(f).assign(tag=pathlib.Path(f).parent.name) for f in fs], ignore_index=True); G = G[G.status == "ok"]
    out, qi = cal_report(G); R = load_runs(root, {"t4"})
    if len(R):
        M = R[(R.arm == "M")].drop_duplicates("tag").set_index("tag"); c = [t for t in M.index if t in qi.index]
        if c:
            fl = M.loc[c, "q_first"]; out["flown_first_pass_q"] = dict(n=len(c), median=float(fl.median()), tiers=sia.tiers(fl.values), median_flown_over_implied=float((fl / qi.loc[c]).median()),
                                                                          n_walk_in_cells=int((M.loc[c, "walk_in"] > 0).sum()))
    if "cd_factor" in summaries(root): out["cd_factor_logged"] = True
    return out


def phase_integrity(root):
    out = dict(campaigns={}, ok=True); bad = []
    for f in sorted(glob.glob(str(sdir(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); c = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"]); out["campaigns"][pathlib.Path(f).name] = c
        if c["route"] != "process" or c["failed"] != 0: bad.append(f"{pathlib.Path(f).name}: route {c['route']}, failed {c['failed']}")
    for d in sorted(glob.glob(str(sdir(root) / "*" / "*" / "cal"))):
        cd = pathlib.Path(d); key = f"{cd.parent.parent.name}/{cd.parent.name}"; fs = sorted(cd.glob("*/guess.csv"))
        if not fs: continue
        G = pd.concat([pd.read_csv(g).assign(tag=g.parent.name) for g in fs], ignore_index=True); fit_only = cd.parent.parent.name == "t4"
        e = dict(n_cal=len(G), all_ok=bool((G.status == "ok").all()))
        if not fit_only: e["max_rel_q_err"] = float((G.q_first_w_cm2 / Q_TGT - 1).abs().max()); e["converged_0p1pct"] = bool(e["max_rel_q_err"] <= 1e-3)
        cb = cd.parent / "cells_base.csv"
        if cb.exists() and "epoch_utc" in pd.read_csv(cb, nrows=0).columns and "gram_core_epoch" in G:
            C = pd.read_csv(cb).set_index("tag"); dt = [abs((pd.Timestamp(r.gram_core_epoch) - pd.Timestamp(C.loc[r.tag, "epoch_utc"])).total_seconds()) for r in G.itertuples() if r.tag in C.index]
            e["epochs_match"] = bool(dt and max(dt) < 1.0)
        else: e["epochs_match"] = None                         # cannot be verified: counts as a failure
        out[key] = e
        if not e["all_ok"] or e.get("converged_0p1pct") is False or e["epochs_match"] is not True: bad.append(f"{key}: cal {e}")
    nan = 0; nruns = 0
    for p in glob.glob(str(sdir(root) / "*" / "*" / "runs" / "*" / "passes.csv")):
        T = pd.read_csv(p); nruns += 1; nan += int(T.drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum())
    out["runs"] = nruns; out["nan_cells_excl_h_qhigh"] = nan
    if nan: bad.append(f"{nan} NaN cells")
    S = summaries(root)
    if len(S):
        if "v5_exceptions" in S:
            M = S[S.arm == "M"]; dv = M[(M.split == "dev") & (M.v5_exceptions > 0)]; hl = M[(M.split != "dev") & (M.v5_exceptions > 0)]
            out["v5_exceptions_dev_M"] = dv[["set", "planet", "tag"]].to_dict("records"); out["v5_exceptions_held_M_listed"] = hl[["set", "planet", "tag", "v5_exceptions"]].to_dict("records")
            if len(dv): bad.append(f"v5 exceptions in M on {len(dv)} dev runs")
        if "horizon_hit" in S: out["horizon_hits_listed"] = S[S.horizon_hit.fillna(False).astype(bool)][["set", "planet", "arm", "tag", "seed_idx"]].to_dict("records")
        if "status" in S: out["rows_not_ok"] = int((S.status != "ok").sum())
        if (S.get("status", pd.Series(["ok"])) != "ok").any(): bad.append("rows with status != ok")
    out["failures"] = bad; out["ok"] = not bad
    return out


def phase_probe(root):
    S = summaries(root); S = S[S.phase == "probe"] if len(S) and "phase" in S else S
    if len(S) == 0: return {}
    R = load_runs(root); R = R[R.phase == "probe"] if len(R) else R
    out = [dict(set=r["set"], arm=r.arm, tag=r.tag, passes=int(r.n_pass), days=float(r.t_end_days), v5_exceptions=int(r.get("v5_exceptions", 0)), wall_s_per_pass_contended=float(r.wall_s / r.n_pass) if "wall_s" in r and r.n_pass else None,
                err_m=float(r.err_m)) for _, r in R.iterrows()]
    return dict(runs=out, note="wall times are contended and are not performance numbers")


def write(root, name, obj):
    d = sdir(root) / "analysis"; d.mkdir(parents=True, exist_ok=True); s = json.dumps(obj, indent=1, default=str); (d / name).write_text(s); print(s)


def selftest():
    rng = np.random.default_rng(1)
    # tiers (exact; strict inequalities as Stage I) and Wilson bound
    assert sia.tiers([0.04, 0.05, 0.25, 0.26, 0.31, 0.46]) == dict(low=1, soft=1, medium=1, hard=1)
    assert abs(wilson_lower(10, 10) - 1 / (1 + ZW ** 2 / 10)) < 1e-12 and abs(ZW - 1.6448536269514722) < 1e-12
    assert abs(wilson_lower(95, 100) - 0.90084) < 1e-4 and np.isnan(wilson_lower(0, 0))
    assert walk_in([0.02, 0.03, 0.1, 0.02]) == 2 and walk_in([0.1]) == 0 and walk_in([0.01, 0.02]) == 2

    def mk(st, pl, arm, tag, err, dv=(1.0, 20.0, 1.0), days=20.0, n=10, over=0, hard=0, seed=0, split="heldout", low=0, low_pl=0):
        return dict(**{"set": st, "planet": pl}, arm=arm, tag=tag, seed_idx=seed, split=split, err_m=err, dv_flown=dv[0], dv_raise=dv[1], dv_trim=dv[2], dv_total=sum(dv), t_end_days=days, n_pass=n, low=low, low_planned=low_pl, walk_in=0,
                    soft=over, medium=0, hard=hard, over=over, q_last=0.03, omega=70.0, inc=93.5)
    # S-E: 24 cells x 4 seeds; M within 5 km everywhere, B1 at 8 km everywhere -> P1 pass, P3 point mass at +1 (lower bound > 0)
    rows = [mk("se", "mars", a, f"c{c}", e, seed=s) for c in range(24) for s in range(4) for a, e in (("M", 2000.0), ("B1", 8000.0))]
    R = pd.DataFrame(rows); se = complete(sel(R, "se", "mars"), ["M", "B1"])
    p1 = crit_P1(se[se.arm == "M"]); assert p1["passed"] and p1["n"] == 96 and p1["frac10"] == 1.0 and p1["wilson90_lower_frac5"] > 0.96
    p3 = crit_P3(se, rng); assert p3["passed"]
    assert p3["mean"] == 1.0 and p3["ci90"] == [1.0, 1.0]
    # P1 fails on the Wilson bound: 36 of 96 within 5 km (0.375) with all within 10 km; and on the 10 km fraction
    e5 = pd.DataFrame([mk("se", "mars", "M", f"c{i}", 2000.0 if i < 36 else 8000.0) for i in range(96)]); assert not crit_P1(e5)["passed"]
    e10 = pd.DataFrame([mk("se", "mars", "M", f"c{i}", 2000.0 if i < 85 else 12000.0) for i in range(96)]); assert crit_P1(e10)["frac10"] < 0.90 and not crit_P1(e10)["passed"]
    # P3 fails when B1 is as precise as M
    eq = pd.DataFrame([mk("se", "mars", a, f"c{c}", 2000.0) for c in range(24) for a in ("M", "B1")]); assert not crit_P3(eq, rng)["passed"]
    # P2 counts: exactly 22 within 10 km and 12 within 5 km passes; one fewer of either fails
    def p2(n10, n5): return crit_P2(pd.DataFrame([mk("t4", "mars", "M", f"c{i}", 2000.0 if i < n5 else (8000.0 if i < n10 else 20000.0)) for i in range(24)]))["passed"]
    assert p2(22, 12) and not p2(21, 12) and not p2(22, 11)
    # T4: S1 (hard <= B1 and over fraction <= B1), D1 (paired mean upper bound < 0), C1 (median ratio <= 1.10)
    t4 = pd.DataFrame([mk("t4", "mars", "M", f"c{c}", 1000.0, dv=(1.0, 20.0, 0.4), days=21.0, over=1) for c in range(24)] + [mk("t4", "mars", "B1", f"c{c}", 6000.0, dv=(1.0, 20.0, 1.5), days=20.0, over=1) for c in range(24)])
    assert crit_S1(t4)["passed"] and crit_D1(t4, rng)["passed"] and abs(crit_D1(t4, rng)["mean_mps"] + 1.1) < 1e-12 and crit_C1(t4)["passed"] and abs(crit_C1(t4)["median_ratio"] - 1.05) < 1e-12
    bad = t4.copy(); bad.loc[(bad.arm == "M") & (bad.tag == "c0"), "hard"] = 1; assert not crit_S1(bad)["passed"]
    bad = t4.copy(); bad.loc[bad.arm == "M", "over"] = 2; assert not crit_S1(bad)["passed"]
    bad = t4.copy(); bad.loc[bad.arm == "M", "dv_total"] += 1.2; assert not crit_D1(bad, rng)["passed"]
    bad = t4.copy(); bad.loc[bad.arm == "M", "t_end_days"] = 23.0; assert not crit_C1(bad)["passed"]
    # G rules: M no worse; G2 fails when M is worse in every cell (lower bound > 0); G3 at the 10% edge
    g = pd.DataFrame([mk("g", "earth", a, f"c{c}", 2000.0, over=0 if a == "M" else 1, dv=(1.0, 20.0, 1.0)) for c in range(24) for a in ("M", "B1")]); cg = crit_G(g, rng); assert cg["all"] and cg["G2"]["passed"]
    g2 = pd.DataFrame([mk("g", "earth", a, f"c{c}", 2000.0, over=1 if a == "M" else 0) for c in range(24) for a in ("M", "B1")]); assert not crit_G(g2, rng)["G2"]["passed"]
    g3 = pd.DataFrame([mk("g", "earth", a, f"c{c}", 2000.0, dv=(0.0, 0.0, 22.0) if a == "M" else (0.0, 0.0, 20.0)) for c in range(24) for a in ("M", "B1")]); assert crit_G(g3, rng)["G3"]["passed"]
    g3b = pd.DataFrame([mk("g", "earth", a, f"c{c}", 2000.0, dv=(0.0, 0.0, 22.1) if a == "M" else (0.0, 0.0, 20.0)) for c in range(24) for a in ("M", "B1")]); assert not crit_G(g3b, rng)["G3"]["passed"]
    # missing-data rule and total dV
    m = pd.concat([t4, pd.DataFrame([mk("t4", "mars", "M", "c99", 1000.0)])]); assert "c99" not in set(complete(m, ["M", "B1"]).tag) and len(complete(m, ["M", "B1"])) == 48
    assert mk("t4", "mars", "M", "c", 0.0, dv=(1.5, 19.82, 0.44))["dv_total"] == 1.5 + 19.82 + 0.44
    # full Mars table and verdict, then the same with D1 broken
    full = pd.concat([R, t4]); h = phase_held_tables(full, rng); assert h["verdict"] == "TARGETING_PASS", h["verdict"]
    full.loc[(full["set"] == "t4") & (full.arm == "M"), "dv_total"] += 5.0; assert phase_held_tables(full, rng)["verdict"] == "FAILED(D1)"
    # A3-J2: J2-only M 0.02 worse than 20x20 M at 5 km (one of 4 seeds in 24 cells is out -> -0.25 per cell: flagged), reported only
    j = pd.DataFrame([mk("se_j2", "mars", "M", f"c{c}", 2000.0 if s else 6000.0, seed=s) for c in range(24) for s in range(4)]); a3 = phase_held_tables(pd.concat([R, t4, j]), rng)["A3_J2"]
    assert abs(a3["dfrac5"]["mean"] + 0.25) < 1e-12 and not a3["abs_le_0p05"] and a3["reported_only"]
    # error distribution counts, signed errors, low-tier split
    s = arm_stats(pd.DataFrame([mk("se", "mars", "M", "a", -9000.0, low=3, low_pl=2), mk("se", "mars", "M", "b", 4000.0, low=1, low_pl=0), mk("se", "mars", "M", "c", 12000.0)]))
    assert s["n10"] == 2 and s["n5"] == 1 and s["n_under"] == 1 and s["n_over"] == 2 and s["low_planned_endgame"] == 2 and s["low_other"] == 2 and abs(s["signed_mean_m"] - 7000 / 3) < 1e-9
    # T4 cal report: implied Q = 0.15 exp((h_guess - h_ic)/H); one cell 1.8 H low (below floor, polar), one at target
    G = pd.DataFrame(dict(tag=["a", "b", "c"], h_ic_m=[100000.0 + 1.8 * 6000, 100000.0, 100000.0 - 6000], h_guess_m=[100000.0, 100000.0, 100000.0], hloc_at_hg_m=[6000.0] * 3, peri_lat_deg=[88.0, 64.0, 80.0]))
    c, qi = cal_report(G); assert c["n_below_floor"] == 1 and c["n_polar"] == 2 and c["n_polar_below_floor"] == 1 and abs(qi["a"] - 0.15 * np.exp(-1.8)) < 1e-12 and abs(qi["b"] - 0.15) < 1e-12 and c["n_above_0p25"] == 1
    # file layout round trip: a synthetic tree read by load_runs, integrity and cal
    with tempfile.TemporaryDirectory() as td:
        sd = sdir(td); rows = []
        for st, pl, arm in (("t4", "mars", "M"), ("t4", "mars", "B1")):
            d = setdir(td, st, pl); (d / "runs" / f"{arm}_a").mkdir(parents=True, exist_ok=True); (d / "cal" / "a").mkdir(parents=True, exist_ok=True)
            pd.DataFrame(dict(**{"pass": range(5)}, time_s=np.arange(5) * 86400.0, q_eff_w_cm2=[0.02, 0.1, 0.3, 0.5, 0.04], apo_alt_m=[5e6, 4.9e6, 4.8e6, 4.7e6, 4.6e6], phase=["main"] * 3 + ["penultimate", "last"])).to_csv(d / "runs" / f"{arm}_a" / "passes.csv", index=False)
            pd.DataFrame(dict(status=["ok"], q_first_w_cm2=[0.15], q_target_w_cm2=[0.15], gram_core_epoch=["2001-12-05 00:00:00"], h_ic_m=[110000.0], h_guess_m=[108000.0], hloc_at_hg_m=[6000.0], peri_lat_deg=[80.0])).to_csv(d / "cal" / "a" / "guess.csv", index=False)
            pd.DataFrame(dict(tag=["a"], split=["heldout"], inc=[93.5], raan=[115.0], omega=[75.0], epoch_idx=[0], epoch_utc=["2001-12-05 00:00:00"])).to_csv(d / "cells_base.csv", index=False)
            rows.append(dict(phase="theld", **{"set": st, "planet": pl}, arm=arm, tag="a", seed_idx=0, split="heldout", run_id=f"{arm}_a", status="ok", retcode="Success", err_m=3000.0, dv_flown=2.0, dv_raise=20.0, dv_trim=0.5,
                             t_end_days=4.0, burns=2, horizon_hit=False, v5_exceptions=0))
        pd.DataFrame(rows).to_csv(sd / "summary_theld_ws.csv", index=False); (sd / "mc_theld_ws.toml").write_text('mc_route = "process"\nmc_workers = 4\nmc_samples = 2\nmc_failed = 0\n')
        L = load_runs(td); assert len(L) == 2 and (L.dv_total == 22.5).all()
        assert (L.low == 2).all() and (L.low_planned == 1).all() and (L.walk_in == 1).all() and (L.hard == 1).all() and (L.medium == 0).all() and (L.soft == 1).all(), L[["low", "low_planned", "walk_in", "hard", "medium", "soft"]]
        I = phase_integrity(td); assert I["ok"], I["failures"]
        cr = phase_cal(td); assert cr["n"] == 1 and cr["flown_first_pass_q"]["tiers"]["low"] == 1 and cr["n_below_floor"] == 0
        pd.DataFrame(rows).assign(**{"mc": 0}).iloc[:1].to_csv(sd / "summary_theld_ws.csv", index=False); (sd / "mc_theld_ws.toml").write_text('mc_route = "threads"\nmc_workers = 4\nmc_samples = 2\nmc_failed = 1\n')
        assert not phase_integrity(td)["ok"]
    # optional cross-check that the WP1 library carries the fixed names (not required for the analysis, which reads dv_raise/dv_trim from the run summary)
    try:
        import stagem_lib as sl; miss = [n for n in ("raise_dv", "trim_dv", "endgame_decision", "q_end") if not hasattr(sl, n)]; assert not miss, miss; print("stagem_lib names present")
    except ImportError: print("stagem_lib not importable here: library cross-check skipped")
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--planets", default="mars,earth,venus"); a = ap.parse_args()
    if a.phase == "selftest": return selftest()
    root = a.results_root; rng = np.random.default_rng(AN["bootstrap_seed"]); pls = a.planets.split(",")
    if a.phase == "integrity": return write(root, "m_integrity.json", phase_integrity(root))
    if a.phase == "cal": return write(root, "m_cal_t4.json", phase_cal(root))
    if a.phase == "probe": return write(root, "m_probe.json", phase_probe(root))
    R = load_runs(root)
    if a.phase == "held": return write(root, "m_held.json", phase_held_tables(R, rng))
    if a.phase == "gen": return write(root, "m_gen.json", phase_gen_tables(R, rng, pls))
    if a.phase == "report":
        return write(root, "m_report.json", dict(integrity=phase_integrity(root), cal_t4=phase_cal(root), held=phase_held_tables(R, rng), gen=phase_gen_tables(R, rng, pls)))
    raise SystemExit(f"unknown phase {a.phase}")


if __name__ == "__main__": main()
