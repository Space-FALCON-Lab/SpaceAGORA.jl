"""WP-P3 analysis of the planner arms (criteria D1a, D1b, D2, D3, D5, D6, D7 and C1, C2, C4 to C7; the attribution algebra; the OL / OL-R paired metrics). Python only; reads flight output, runs no Julia.
Plan: tasks/2026-10-15_periapsis-region-border-targeting/analytic_framework_and_transfer_plan.md R1.4 (C1, C2, C4 to C7), R2.2 (ladder, factorial), R4.6 (D1a to D5), R7.4 (D6), R8.2 (OL against OL-R, D7).

  selftest                                   synthetic rows covering every criterion, the ladder telescoping, the factorial estimator, the loader and the per-run metrics (no flights, no Julia)
  metrics  --results-root R [--stageo-root S]   one row per run (dV, burns, days, passes, Q_s violations, tracking against the tube centerline, final apoapsis error)
  report   --results-root R [--stageo-root S] [--task-dir T]   all criteria as JSON (held seeds), paired by seed, 90% percentile-bootstrap intervals over seeds

Definitions fixed here (each a reading of the plan text, flagged for the coordinator in planner_build.md):
 * Q_s is the pass maximum on the atmosphere-relative speed (U15): 0.5 rho_peak V_rel^3 / 1e4 W/cm2; the inertial form is reported next to it. Violations are counts of passes above 0.32 and above 0.65 (R8.2).
 * Tube. q_mean (the passes.csv q_w_cm2, mean-field heat rate at the periapsis event on V_rel) against the planner's tube [q_cap exp(-W), q_cap]; the centerline is ln q_cap - W/2. Periapsis tracking error in meters is
   -H_local x (ln q - ln q_center) (one scale height per e-fold of the heat rate); rms and max of |error| per run; "in tube" is the fraction of passes with q inside the tube.
 * Final apoapsis error = r_a of the last flown pass minus the target 7,408 km (runs stop at the first apoapsis at or below it, so a reached run has an error <= 0; a run that never reaches it has an error > 0).
 * C1 uses D_ref per seed = min(days B-ODY, days F-O) at the Odyssey geometry, for every arm including the geometry arms (R4.5). C2: the pooled count of passes with Q_s > 0.54 over the arm's n seeds has an exact Poisson
   90% interval containing E x n, and no pass above 0.65.
 * Attribution uses per-seed paired MEAN differences (so the ladder telescopes exactly per seed); the pass criteria of D1 use paired MEDIANS of dV ratios (R4.6). P_g is the geometry effect at the full planner,
   dV(L3@ODY) - dV(L3@G1) = main effect + interaction / 2 with main = mean of the geometry effects at L1 and L3 and interaction = their difference (R2.2: the interaction is counted in P).
 * D6 verdict (R7.4): supported if L3 has lower dV at equal-or-lower violations, or fewer violations at equal dV; "lower" = the 90% interval of (comparator - L3) lies above 0, "equal-or-lower" = its upper bound is >= 0.
Usage: python3 -I planner_analysis.py selftest | metrics ... | report ...
"""
import argparse, glob, json, math, pathlib, sys, tempfile
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB

CFG = PLB.load_cfg(); AN = CFG["planner"]["analysis"]
def ra_target_m(override=None):
    """Closing apoapsis RADIUS of the campaign (m). Default: planner_config.toml [planner.analysis] ra_target_m = 7,408,000 (SOURCED stageo_config.toml, Table III end of main phase 2). Revision 25 parameter: the environment variable
    SPACEAGORA_RA_TARGET_M (read once at import, so every importing module sees the same value) or the argument override; unset = the default and every output is unchanged."""
    import os
    v = override if override is not None else os.environ.get("SPACEAGORA_RA_TARGET_M")
    return float(AN["ra_target_m"]) if v in (None, "") else float(v)


MU = AN["mu_mars_m3s2"]; RA_TARGET = ra_target_m(); LEVELS = AN["violation_levels_qs"]
RATIO = CFG["planner"]["burn_response_ratio"]


# ------------------------------------------------------------------------------------------------ statistics

def boot(vals, stat, rng, n=None, ci=None):
    """Percentile bootstrap over seeds of `stat` (np.median or np.mean): dict(point, lo, hi, n)."""
    v = np.asarray(vals, float); v = v[np.isfinite(v)]; n = n or AN["bootstrap_draws"]; ci = ci or AN["ci_percent"]
    if len(v) == 0:
        return dict(point=float("nan"), lo=float("nan"), hi=float("nan"), n=0)
    if len(v) == 1:
        return dict(point=float(v[0]), lo=float("nan"), hi=float("nan"), n=1)
    idx = rng.integers(0, len(v), size=(n, len(v))); s = stat(v[idx], axis=1); a = (100 - ci) / 2
    return dict(point=float(stat(v)), lo=float(np.percentile(s, a)), hi=float(np.percentile(s, 100 - a)), n=int(len(v)))


def poisson_ci(k, ci=None):
    """Exact (Garwood) two-sided interval for a Poisson mean given k events."""
    from scipy.stats import chi2
    ci = ci or AN["ci_percent"]; a = (100 - ci) / 200
    return (0.0 if k == 0 else float(chi2.ppf(a, 2 * k) / 2)), float(chi2.ppf(1 - a, 2 * k + 2) / 2)


def rng_of():
    return np.random.default_rng(AN["bootstrap_seed"])


def pair(A, B, f):
    """Per-seed f(a, b) over the seeds both dicts have (missing-seed rule: a seed absent from either is dropped). Returns (seeds, array)."""
    s = sorted(set(A) & set(B)); return s, np.array([f(A[k], B[k]) for k in s], float)


def get(key):
    return lambda r: r[key]


# ------------------------------------------------------------------------------------------------ per-run metrics

def qs_rel(rho_peak, vrel):
    return 0.5 * np.asarray(rho_peak, float) * np.asarray(vrel, float) ** 3 / 1e4


def qs_inertial(rho_peak, rp, ra):
    rp = np.asarray(rp, float); ra = np.asarray(ra, float)
    return 0.5 * np.asarray(rho_peak, float) * (MU * (2 / rp - 2 / (rp + ra))) ** 1.5 / 1e4


def tube_of(pred_arm, cfg=CFG):
    """Tube [ln q_lo, ln q_hi] and centerline of an arm from its prediction row (tube_W_ln)."""
    hi = math.log(cfg["planner"]["q_cap_mean_w_cm2"]); lo = hi - float(pred_arm["tube_W_ln"]); return dict(lo=lo, hi=hi, center=0.5 * (lo + hi))


def run_metrics(P, D, srow, tube=None, tab=None, ra_target=RA_TARGET, ratio=RATIO):
    """Metrics of one flight. P: passes.csv (pass, rho_peak_kg_m3, vrel_mps, q_w_cm2, rp_m, ra_m, h_local_m, [rho]); D: decisions.csv (target_pass, delta_v_mps) or empty; srow: summary row (total_dv_mps, burns, t_end_days,
    stop_reason). tube: dict(lo, hi, center) in ln q for the tracking metrics; tab: ground table dict (path) for the predicted 10-pass drift of C5."""
    P = P.sort_values("pass").reset_index(drop=True); P = P[np.isfinite(P.rho_peak_kg_m3) & np.isfinite(P.vrel_mps)].reset_index(drop=True)
    ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= ra_target))[0]
    reached = bool(len(hit)) or str(srow.get("stop_reason", "")) == "end_radius"; n = int(hit[0] + 1) if len(hit) else len(P); Pn = P.iloc[:n]
    qs = qs_rel(Pn.rho_peak_kg_m3, Pn.vrel_mps); lnq = np.log(Pn.q_w_cm2.values.astype(float)); H = Pn.h_local_m.values.astype(float)
    out = dict(n_passes=n, reached=reached, days=float(srow.get("t_end_days", float("nan"))), stop_reason=str(srow.get("stop_reason", "")),
               qs_max=float(qs.max()) if n else float("nan"), v054=int((qs > 0.54).sum()), qs_inertial_max=float(np.nanmax(qs_inertial(Pn.rho_peak_kg_m3, Pn.rp_m, Pn.ra_m))) if n and np.isfinite(Pn.ra_m.values.astype(float)).any() else float("nan"))
    for lv in LEVELS:
        out[f"v{int(round(lv * 100)):03d}"] = int((qs > lv).sum())
    out["v032"] = int((qs > 0.32).sum()); out["v065"] = int((qs > 0.65).sum())
    out["lnratio"] = np.log(qs / Pn.q_w_cm2.values.astype(float))                         # ln(Q_s / q_mean) per pass (C7)
    fin = ra[:n][np.isfinite(ra[:n])]; ra_final = float(fin[-1]) if len(fin) else float("nan")
    out["ra_final_m"] = ra_final; out["final_apo_err_m"] = ra_final - ra_target; out["final_apo_abs_err_m"] = abs(ra_final - ra_target)
    Dd = D if D is not None and len(D) else pd.DataFrame({"target_pass": [], "delta_v_mps": []})
    dvp = np.zeros(n + 2)
    for r in Dd.itertuples():
        if 1 <= int(r.target_pass) <= n:
            dvp[int(r.target_pass)] += float(r.delta_v_mps)
    dv_sum = float(np.abs(dvp).sum()); dv = float(srow["total_dv_mps"]) if np.isfinite(srow.get("total_dv_mps", float("nan"))) else dv_sum
    out["dv"] = dv; out["burns"] = int(srow["burns"]) if "burns" in srow and np.isfinite(srow["burns"]) else int((np.abs(dvp) > 0).sum()); out["dv_decisions_sum"] = dv_sum
    if tube is not None and n:
        e_ln = lnq - tube["center"]; e_m = -H * e_ln
        out.update(track_rms_m=float(np.sqrt(np.mean(e_m ** 2))), track_max_m=float(np.abs(e_m).max()), track_rms_ln=float(np.sqrt(np.mean(e_ln ** 2))), in_tube=float(np.mean((lnq >= tube["lo"]) & (lnq <= tube["hi"]))))
    else:
        out.update(track_rms_m=float("nan"), track_max_m=float("nan"), track_rms_ln=float("nan"), in_tube=float("nan"))
    # natural path and C5 quantities: burn between passes i-1 and i acts through the state of pass i-1; cancel_i = dv_i / w_i ln units removed
    rp = Pn.rp_m.values.astype(float); w = np.zeros(n); cancel = np.zeros(n)
    for i in range(1, n):
        if np.isfinite(ra[i - 1]) and np.isfinite(rp[i - 1]):
            a = 0.5 * (ra[i - 1] + rp[i - 1]); e = (ra[i - 1] - rp[i - 1]) / (ra[i - 1] + rp[i - 1]); w[i] = PLB.burn_weight(H[i], ratio, MU, a, e)
            cancel[i] = dvp[i + 1] / w[i]                                                  # decisions.csv target_pass is 1-based: the burn before 0-based pass i has target_pass i + 1
    nat = lnq + np.cumsum(cancel)                                                                   # remove the burns: a burn that cancels heat (dv > 0) lowered later ln q
    out["net_nat_ln"] = float(nat[-1] - nat[0]) if n else float("nan"); out["cancel_sum_ln"] = float(cancel.sum())
    opp = tot = 0.0
    if tab is not None:
        path = np.asarray(tab["path"], float)
        for i in range(1, n):
            dvi = dvp[i + 1]
            if dvi == 0.0:
                continue
            j0 = min(i, len(path) - 1); j1 = min(i + AN["c5_lookahead_passes"], len(path) - 1); drift = path[j1] - path[j0]
            tot += abs(dvi); opp += abs(dvi) if (dvi * drift < 0.0) else 0.0
    out["oppose_dv_frac"] = (opp / tot) if tot > 0 else (0.0 if tab is not None else float("nan"))
    return out


# ------------------------------------------------------------------------------------------------ criteria

def c1(runs, body, fo, rng):
    """C1 time no worse (R1.4): upper 90% bound of the paired median ratio days(P) / D_ref and passes(P) / min passes, D_ref = per seed min(B-ODY, F-O) at Odyssey geometry, <= tolerance."""
    seeds = sorted(set(runs) & set(body) & set(fo)); tol = AN["c1_tolerance"]
    dr = {k: min(body[k]["days"], fo[k]["days"]) for k in seeds}; pr = {k: min(body[k]["n_passes"], fo[k]["n_passes"]) for k in seeds}
    rd = np.array([runs[k]["days"] / dr[k] for k in seeds]); rpn = np.array([runs[k]["n_passes"] / pr[k] for k in seeds])
    bd, bp = boot(rd, np.median, rng), boot(rpn, np.median, rng)
    return dict(days=bd, passes=bp, tolerance=tol, n_seeds=len(seeds), not_reached=sum(not runs[k]["reached"] for k in seeds), **{"pass": bool(bd["hi"] <= tol and bp["hi"] <= tol)})


def c2(runs, E=None):
    """C2 risk (R1.4): pooled count of passes with Q_s > 0.54 has an exact Poisson 90% interval containing E x n_seeds, and no pass above 0.65."""
    E = CFG["planner"]["E_budget"] if E is None else E; n = len(runs); k = int(sum(r["v054"] for r in runs.values())); lo, hi = poisson_ci(k); k65 = int(sum(r["v065"] for r in runs.values()))
    return dict(n_seeds=n, count_gt054=k, poisson_ci=[lo, hi], expected=E * n, within_budget=bool(lo <= E * n <= hi), count_gt065=k65, **{"pass": bool(lo <= E * n <= hi and k65 == 0)})


def c4(runs, pred_dv, floors=None):
    """C4 floor and prediction (R1.4): median dV in [median own net-drift floor, 2 x the prediction]. floors: {seed: floor} from the arm's own flown paths (None -> NEEDS_FLOOR; the upper bound is still tested)."""
    med = float(np.median([r["dv"] for r in runs.values()])) if runs else float("nan"); hi = 2.0 * pred_dv
    lo = float(np.median([floors[k] for k in runs if k in floors])) if floors else None
    ok_hi = bool(med <= hi); ok_lo = None if lo is None else bool(med >= lo)
    return dict(median_dv=med, lower=lo, upper=hi, upper_ok=ok_hi, lower_ok=ok_lo, **{"pass": (ok_hi and ok_lo) if ok_lo is not None else None})


def c4abs(runs, floor_pred):
    """C4abs (Revision 12, REPORTED ONLY, does not replace C4): the median dV lies within +-c4abs_tol_mps of the predicted floor (the ranking's ensemble floor of the arm's geometry)."""
    med = float(np.median([r["dv"] for r in runs.values()])) if runs else float("nan"); tol = AN["c4abs_tol_mps"]
    return dict(median_dv=med, predicted_floor=floor_pred, tolerance=tol, difference=med - floor_pred, **{"pass": bool(abs(med - floor_pred) <= tol)})


def c5(runs):
    """C5 burning with the drift (R1.4): per run the burn-sum (ln units cancelled) equals the net natural drift within +-30%, and fewer than 20% of burn dV opposes the predicted 10-pass drift. Returns per-run flags and the share of runs passing."""
    flags = {}
    for k, r in runs.items():
        net = r["net_nat_ln"]; ok_net = bool(abs(r["cancel_sum_ln"] - net) <= AN["c5_net_tolerance"] * abs(net)) if np.isfinite(net) and net != 0 else None
        ok_opp = bool(r["oppose_dv_frac"] < AN["c5_oppose_max"]) if np.isfinite(r["oppose_dv_frac"]) else None
        flags[k] = dict(net_ok=ok_net, oppose_ok=ok_opp, **{"pass": (ok_net and ok_opp) if ok_net is not None and ok_opp is not None else None})
    sc = [v["pass"] for v in flags.values() if v["pass"] is not None]
    return dict(per_run=flags, share_pass=(float(np.mean(sc)) if sc else None), n_scored=len(sc))


def vsum_ev(obs, pred, burn_mask=None):
    """Explained variance 1 - sum (obs - pred)^2 / sum obs^2 of per-pass changes on non-burn passes (R1.4 P7 / C6; the V-SUM statistic of analytic_results.md)."""
    obs = np.asarray(obs, float); pred = np.asarray(pred, float); m = np.isfinite(obs) & np.isfinite(pred)
    if burn_mask is not None:
        m &= ~np.asarray(burn_mask, bool)
    return float(1.0 - np.sum((obs[m] - pred[m]) ** 2) / np.sum(obs[m] ** 2)) if m.any() and np.sum(obs[m] ** 2) > 0 else float("nan")


def c6(runs):
    """C6 Part 1 maps in closed loop (R1.4 / P7): median over runs of the V-SUM explained variance >= 0.85. Reads r['vsum_ev'] (attached from the analytic predictions of the flown pass pairs, WP-A tools); NOT_COMPUTED if absent."""
    v = [r.get("vsum_ev") for r in runs.values() if r.get("vsum_ev") is not None and np.isfinite(r.get("vsum_ev"))]
    if not v:
        return dict(status="NOT_COMPUTED", **{"pass": None})
    m = float(np.median(v)); return dict(status="OK", median_ev=m, n=len(v), threshold=AN["c6_min_ev"], **{"pass": bool(m >= AN["c6_min_ev"])})


def c7(runs):
    """C7 U1 statistics transfer (R1.4): q99 of ln(Q_s / q_mean) pooled over the arm's passes within +-0.05 of the dev 0.707; Poisson check of the tail above 0.65 against the dev per-pass rate (reported, not scored, while the rate is unset)."""
    x = np.concatenate([r["lnratio"] for r in runs.values()]) if runs else np.array([]); q99 = float(np.percentile(x, 99)) if len(x) else float("nan")
    npass = int(sum(r["n_passes"] for r in runs.values())); k65 = int(sum(r["v065"] for r in runs.values())); rate = AN["c7_dev_rate_gt065"]
    out = dict(q99=q99, dev_q99=AN["c7_q99_dev"], tolerance=AN["c7_tolerance"], q99_ok=bool(abs(q99 - AN["c7_q99_dev"]) <= AN["c7_tolerance"]) if np.isfinite(q99) else None, n_passes=npass, count_gt065=k65)
    if rate >= 0:
        from scipy.stats import poisson
        lam = rate * npass; out.update(dev_rate=rate, expected=lam, tail_ok=bool(k65 <= poisson.ppf(0.95, lam)))
    else:
        out.update(dev_rate=None, tail_ok=None, note="dev rate of passes above 0.65 not measured: the Poisson tail check is reported, not scored")
    out["pass"] = out["q99_ok"] if out["tail_ok"] is None else bool(out["q99_ok"] and out["tail_ok"]); return out


def d1(runs, bases, body, fo, rng, c1r, c2r):
    """D1a / D1b: upper 90% bound of the paired median dV ratio arm / X < 1 for X in {B-ODY, F-O, N2-I} at Odyssey geometry, with C1 and C2 holding for the arm. bases: {name: runs}."""
    ratios = {}
    for name, X in bases.items():
        s, r = pair(runs, X, lambda a, b: a["dv"] / b["dv"]); ratios[name] = boot(r, np.median, rng)
    ok = all(v["hi"] < 1.0 for v in ratios.values())
    return dict(ratios=ratios, c1=c1r["pass"], c2=c2r["pass"], **{"pass": bool(ok and c1r["pass"] and c2r["pass"])})


# ------------------------------------------------------------------------------------------------ attribution (R2.2)

def attribution(dv, rng, firsts=("N2I", "BODY", "FO")):
    """Paired per-seed mean differences of dV. dv: {'arm@geom': {seed: dV}}. Components against each first rung X: C_X = dV(X@ODY) - dV(L0), N = L0 - L1, P_ff = L1 - L2, P_t = L2 - L3 (all at ODY), P_g = dV(L3@ODY) - dV(L3@G1)
    (geometry effect at the full planner = main effect + interaction / 2); the main effect is the mean of the geometry effects at L1 and L3, the interaction their difference. total_X = dV(X@ODY) - dV(L3@G1) = C_X + N + P_ff + P_t + P_g
    holds exactly per seed. P = P_ff + P_t + P_g."""
    need = ["L0@ODY", "L1@ODY", "L2@ODY", "L3@ODY", "L1@G1", "L3@G1"]; out = {"components": {}}
    for X in firsts:
        names = need + [f"{X}@ODY"]
        if any(n not in dv for n in names):
            out["components"][X] = dict(status="MISSING_ARMS", missing=[n for n in names if n not in dv]); continue
        seeds = sorted(set.intersection(*[set(dv[n]) for n in names])); g = lambda n: np.array([dv[n][k] for k in seeds], float)
        C = g(f"{X}@ODY") - g("L0@ODY"); N = g("L0@ODY") - g("L1@ODY"); Pff = g("L1@ODY") - g("L2@ODY"); Pt = g("L2@ODY") - g("L3@ODY")
        gL1 = g("L1@ODY") - g("L1@G1"); gL3 = g("L3@ODY") - g("L3@G1"); Pg = gL3; main = 0.5 * (gL1 + gL3); inter = gL3 - gL1; P = Pff + Pt + Pg; total = g(f"{X}@ODY") - g("L3@G1")
        assert np.allclose(total, C + N + Pff + Pt + Pg, atol=1e-9) and np.allclose(Pg, main + 0.5 * inter, atol=1e-9)
        comp = {k: boot(v, np.mean, rng) for k, v in dict(C=C, N=N, P_ff=Pff, P_t=Pt, P_g=Pg, P_g_main=main, P_g_interaction=inter, P=P, total=total).items()}
        share = comp["P"]["point"] / comp["total"]["point"] if comp["total"]["point"] not in (0.0,) and np.isfinite(comp["total"]["point"]) else float("nan")
        out["components"][X] = dict(n_seeds=len(seeds), **comp, P_share_of_total=float(share), telescopes=True)
    return out


def d2(attr, X="N2I"):
    """D2 perturbation attribution (R4.6): P = P_ff + P_t + P_g has a 90% interval above 0 AND P >= 25% of the total saving L3@G1 against N2-I (threshold ASSUMED)."""
    c = attr["components"].get(X)
    if c is None or "P" not in c:
        return dict(status="MISSING_ARMS", **{"pass": None})
    ok = bool(c["P"]["lo"] > 0 and c["P_share_of_total"] >= AN["d2_min_share"])
    return dict(P=c["P"], share=c["P_share_of_total"], min_share=AN["d2_min_share"], **{"pass": ok}, fallback_framing=None if ok else "R2.1 fallback: noise tolerance makes aerobraking cheaper (not a perturbation result)")


def d3(attr, X="N2I"):
    c = attr["components"].get(X)
    return {k: c[k] for k in ("C", "N", "P_ff", "P_t", "P_g", "P_g_main", "P_g_interaction") if k in c} if c and "C" in c else dict(status="MISSING_ARMS")


def d5(body_g1, body_ody, rng):
    """D5 (reported): dV(B-ODY@G1) against dV(B-ODY@Odyssey), paired: median difference and ratio."""
    s, d = pair(body_g1, body_ody, lambda a, b: a["dv"] - b["dv"]); _, r = pair(body_g1, body_ody, lambda a, b: a["dv"] / b["dv"])
    return dict(n_seeds=len(s), diff_median=boot(d, np.median, rng), ratio_median=boot(r, np.median, rng))


def d6(comp, l3, rng):
    """D6 (R7.4): paired median dV and violation (count over 0.32; over 0.65 reported) differences comparator - L3. Verdict per the module docstring."""
    out = {}
    for key, getter in (("dv", get("dv")), ("v032", get("v032")), ("v065", get("v065"))):
        s, d = pair(comp, l3, lambda a, b: getter(a) - getter(b)); out[key] = boot(d, np.median, rng)
    dvb = out["dv"]; vb = out["v032"]
    better = lambda b: bool(b["lo"] > 0); not_worse = lambda b: bool(b["hi"] >= 0) if np.isfinite(b["hi"]) else None
    sup = bool((better(dvb) and not_worse(vb)) or (better(vb) and not_worse(dvb)))
    lower = lambda b: bool(b["hi"] < 0); not_higher = lambda b: bool(b["lo"] <= 0)             # comparator significantly lower than L3 / not significantly higher
    against = bool((lower(dvb) and not_higher(vb)) or (lower(vb) and not_higher(dvb)))
    out["verdict"] = "SUPPORTED" if sup else ("AGAINST" if against else "NULL")
    out["reading"] = "SUPPORTED: L3 has lower dV at equal-or-lower violations, or fewer violations at equal dV; NULL: reported as such, the paper then claims the planner and not the autonomy"
    return out


D7_METRICS = [("track_rms_m", "periapsis tracking error vs tube centerline, rms (m)"), ("track_max_m", "tracking error, max (m)"), ("in_tube", "fraction of passes in the tube"), ("final_apo_abs_err_m", "final apoapsis |error| (m)"),
              ("final_apo_err_m", "final apoapsis error, signed (m)"), ("v032", "passes with Q_s > 0.32"), ("v065", "passes with Q_s > 0.65"), ("dv", "total dV (m/s)"), ("burns", "burn count")]


def d7(ol, olr, rng):
    """D7 (R8.2): paired median differences OL minus OL-R with 90% intervals; reported, no threshold."""
    out = {}
    for key, label in D7_METRICS:
        s, d = pair(ol, olr, lambda a, b: a[key] - b[key]); out[key] = dict(label=label, n_seeds=len(s), **boot(d, np.median, rng))
    return out


# ------------------------------------------------------------------------------------------------ D8 replan deviation (Revision 11)

MIN_BURN = CFG["planner"]["dp"]["min_burn_mps"]


def replan_metrics(R, run_wall_s=None, min_burn=MIN_BURN):
    """Per-flight D8 metrics from replan.csv (one row per apoapsis: pass, apriori_dv_mps, executed_dv_mps, dp_solve_s, apriori_build_s). The a priori decision is the table schedule's burn for that orbit, the replan decision
    the executed burn; a decision is 'burn' when |dv| >= min_burn. Returns the number of orbits with a differing burn / no-burn decision, rms and absolute sum of the per-orbit dv difference (executed - a priori), the burn-count
    and total-dV differences, the pass of the first divergence (NaN if none), and the planner cost (DP wall per orbit, total, fraction of the flight wall)."""
    ap = R.apriori_dv_mps.values.astype(float); ex = R.executed_dv_mps.values.astype(float); ps = R["pass"].values
    ba = np.abs(ap) >= min_burn; be = np.abs(ex) >= min_burn; div = ba != be; d = ex - ap
    dp = R.dp_solve_s.values.astype(float); tot = float(np.nansum(dp) + (np.nansum(R.apriori_build_s.values.astype(float)) if "apriori_build_s" in R else 0.0))
    return dict(n_orbits=int(len(R)), n_diverge=int(div.sum()), dv_diff_rms=float(np.sqrt(np.mean(d ** 2))) if len(d) else float("nan"), dv_diff_abs_sum=float(np.abs(d).sum()), dv_diff_signed_sum=float(d.sum()),
                burn_count_diff=int(be.sum() - ba.sum()), total_dv_diff=float(np.abs(ex).sum() - np.abs(ap).sum()), first_diverge_pass=float(ps[np.argmax(div)]) if div.any() else float("nan"),
                n_burn_replan=int(be.sum()), n_burn_apriori=int(ba.sum()), dp_median_s=float(np.nanmedian(dp)) if len(dp) else float("nan"), dp_total_s=tot,
                dp_fraction=(tot / run_wall_s) if run_wall_s else float("nan"))


D8_METRICS = ["n_diverge", "dv_diff_rms", "dv_diff_abs_sum", "burn_count_diff", "total_dv_diff", "first_diverge_pass", "dp_median_s", "dp_total_s", "dp_fraction"]


def d8(runs, rng):
    """D8 (Revision 11, reported, no threshold): medians over paired seeds of the per-flight replan deviation metrics with 90% intervals. NaN first-divergence entries (no divergence) are dropped from that median
    and counted in n_no_divergence."""
    rs = [r["replan"] for r in runs.values() if r.get("replan")]
    out = {k: boot([r[k] for r in rs], np.median, rng) for k in D8_METRICS}; out["n_flights"] = len(rs); out["n_no_divergence"] = int(sum(not np.isfinite(r["first_diverge_pass"]) for r in rs))
    return out


def d8_ol(runs):
    """OL logs the a priori schedule only: burn count and total dV of the schedule over the orbits flown (median)."""
    rs = [r["replan"] for r in runs.values() if r.get("replan")]
    return dict(n_flights=len(rs), n_burn_apriori_median=float(np.median([r["n_burn_apriori"] for r in rs])) if rs else float("nan"), n_diverge_median=float(np.median([r["n_diverge"] for r in rs])) if rs else float("nan"))


# ------------------------------------------------------------------------------------------------ loading

def task_dir(arg=None):
    return pathlib.Path(arg) if arg else PLB.task_dir()


def gate_roles(td):
    """gid -> token (ODY, G1, G2, Gworst) from the hash-pinned geometry table."""
    T = json.loads((td / "geometry_gate_table.json").read_text()); m = {"ODY": "ODY"}
    for g in T["geometries"]:
        m[g["id"]] = g["role"]
    return m


def _run_dir(root_sub, run_id):
    for c in (root_sub / "runs" / run_id, root_sub / "mars" / "runs" / run_id, root_sub / "dev" / run_id, root_sub / "smoke" / run_id):
        if (c / "passes.csv").exists():
            return c
    return None


def collect(roots, roles, preds=None, tabs=None, phases=("held",)):
    """{'ARM@TOK': {seed: metrics}} from the summaries below each root (planner/ and, for F-O and N2-I, stageo/). Arm keys: BODY, L0..L3, OL, OLR from the planner (summary column arm), FO and N2I from the Stage O summaries
    (arm F-O -> FO, N2-I -> N2I). Geometry token from the summary's geometry column (planner) or ODY (Stage O). Tube and ground table are attached per planner arm and geometry when preds / tabs are given."""
    out = {}
    for root in roots:
        root = pathlib.Path(root)
        for ph in phases:
            for f in sorted(glob.glob(str(root / f"summary_{ph}_*.csv"))):
                S = pd.read_csv(f); S = S[S.status == "ok"] if "status" in S else S
                for r in S.itertuples():
                    arm = str(r.arm).replace("-", ""); arm = {"FO": "FO", "N2I": "N2I"}.get(arm, arm); tok = roles.get(getattr(r, "geometry", "ODY"), "ODY")
                    d = _run_dir(root, r.run_id)
                    if d is None:
                        continue
                    seed = int(str(r.run_id).split("_pert_", 1)[1].split("_")[1].lstrip("s")) if "_pert_" in str(r.run_id) else -1
                    P = pd.read_csv(d / "passes.csv"); D = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() and (d / "decisions.csv").stat().st_size > 2 else None
                    key = f"{arm}@{tok}"; pa = (preds or {}).get("arms", {}).get(key); tube = tube_of(pa) if pa and "tube_W_ln" in pa and arm not in ("FO", "N2I") else None
                    tab = ((tabs or {}).get("geometries", {}) or {}).get(getattr(r, "geometry", "ODY"))
                    m = run_metrics(P, D, r._asdict(), tube=tube, tab=tab)
                    if arm in ("L1", "L2", "L3", "OL", "OLR") and (d / "replan.csv").exists() and (d / "replan.csv").stat().st_size > 2:
                        m["replan"] = replan_metrics(pd.read_csv(d / "replan.csv"), float(getattr(r, "run_wall_s", float("nan"))) if np.isfinite(getattr(r, "run_wall_s", float("nan"))) else None)
                    out.setdefault(key, {})[seed] = m
    return out


def report(runs, preds=None, floors=None):
    rng = rng_of(); R = {}; body, fo, n2 = runs.get("BODY@ODY", {}), runs.get("FO@ODY", {}), runs.get("N2I@ODY", {})
    bases = {k: v for k, v in (("BODY", body), ("FO", fo), ("N2I", n2)) if v}
    R["integrity"] = {k: dict(n=len(v), not_reached=sum(not r["reached"] for r in v.values()), stop_reasons=sorted({r["stop_reason"] for r in v.values()})) for k, v in runs.items()}
    for key, rr in runs.items():
        arm, tok = key.split("@")
        if arm in ("FO", "N2I") or not rr:
            continue
        R[f"C1:{key}"] = c1(rr, body, fo, rng) if body and fo else dict(status="NEEDS_BODY_AND_FO"); R[f"C2:{key}"] = c2(rr); R[f"C5:{key}"] = c5(rr); R[f"C6:{key}"] = c6(rr); R[f"C7:{key}"] = c7(rr)
        pa = (preds or {}).get("arms", {}).get(key)
        if pa and arm not in ("BODY",):
            R[f"C4:{key}"] = c4(rr, pa["dv_model_world_mps"], (floors or {}).get(key))
            if pa.get("floor_ens_mps") is not None:
                R[f"C4abs:{key}"] = c4abs(rr, pa["floor_ens_mps"])
    for tok, name in (("G1", "D1a"), ("ODY", "D1b")):
        l3 = runs.get(f"L3@{tok}")
        if l3 and bases and body and fo:
            R[name] = d1(l3, bases, body, fo, rng, R[f"C1:L3@{tok}"], R[f"C2:L3@{tok}"])
    dv = {k: {s: r["dv"] for s, r in v.items()} for k, v in runs.items()}
    R["attribution"] = attribution(dv, rng); R["D2"] = d2(R["attribution"]); R["D3"] = d3(R["attribution"])
    if runs.get("BODY@G1") and body:
        R["D5"] = d5(runs["BODY@G1"], body, rng)
    for tok in ("ODY", "G1"):
        l3 = runs.get(f"L3@{tok}")
        for comp in ("OL", "OLR"):
            if l3 and runs.get(f"{comp}@{tok}"):
                R[f"D6:{comp}@{tok}"] = d6(runs[f"{comp}@{tok}"], l3, rng)
        for comp in ("L1", "L2", "L3", "OLR"):
            if runs.get(f"{comp}@{tok}"):
                R[f"D8:{comp}@{tok}"] = d8(runs[f"{comp}@{tok}"], rng)
        if runs.get(f"OL@{tok}"):
            R[f"D8:OL@{tok}"] = d8_ol(runs[f"OL@{tok}"])
        if runs.get(f"OL@{tok}") and runs.get(f"OLR@{tok}"):
            R[f"D7:{tok}"] = d7(runs[f"OL@{tok}"], runs[f"OLR@{tok}"], rng)
    return R


# ------------------------------------------------------------------------------------------------ selftest

def _synth_passes(n=40, ra0=2.85e7, ra_end=7.3e6, q=0.2, rho_ratio=1.5, vrel=4700.0, H=6450.0, h_burn=None):
    ra = np.linspace(ra0, ra_end, n); rp = np.full(n, 3.4838e6)
    P = pd.DataFrame(dict(**{"pass": np.arange(1, n + 1)}, time_s=np.arange(1, n + 1) * 8.0e4, q_w_cm2=np.full(n, q), rho=2 * q * 1e4 / vrel ** 3, rho_peak_kg_m3=rho_ratio * 2 * q * 1e4 / vrel ** 3,
                          vrel_mps=vrel, rp_m=rp, ra_m=ra, h_local_m=H, alt_m=1.0e5)); return P


def _run(dv, v054=0, v032=0, v065=0, days=64.0, passes=200, **kw):
    r = dict(dv=dv, v054=v054, v032=v032, v065=v065, days=days, n_passes=passes, reached=True, stop_reason="end_radius", lnratio=np.array([0.3, 0.4]), net_nat_ln=1.0, cancel_sum_ln=1.0, oppose_dv_frac=0.0, burns=5,
             track_rms_m=1000.0, track_max_m=3000.0, in_tube=0.9, final_apo_abs_err_m=1e4, final_apo_err_m=-1e4); r.update(kw); return r


def selftest():
    rng = rng_of(); seeds = list(range(13, 29))
    # statistics
    b = boot(np.arange(100.0), np.median, rng, 2000); assert abs(b["point"] - 49.5) < 1e-9 and b["lo"] < b["point"] < b["hi"]
    lo, hi = poisson_ci(8); assert lo < 8 < hi and abs(lo - 3.98) < 0.05 and abs(hi - 14.44) < 0.05 and poisson_ci(0)[0] == 0.0
    # per-run metrics on a synthetic flight: Q_s = 1.5 x 0.2 on V_rel; two burns; tube tracking
    P = _synth_passes(); P.loc[10:, "ra_m"] = np.linspace(2.0e7, 7.3e6, 30); D = pd.DataFrame(dict(target_pass=[5, 12], delta_v_mps=[0.1, -0.05]))
    tube = dict(lo=math.log(0.234) - 0.72, hi=math.log(0.234)); tube["center"] = 0.5 * (tube["lo"] + tube["hi"])
    m = run_metrics(P, D, dict(total_dv_mps=0.15, burns=2, t_end_days=64.0, stop_reason="end_radius"), tube=tube)
    qs = 0.3; assert abs(m["qs_max"] - qs) < 1e-9 and m["v032"] == 0 and m["v054"] == 0 and m["reached"] and m["burns"] == 2 and abs(m["dv"] - 0.15) < 1e-12
    assert abs(m["track_rms_ln"] - abs(math.log(0.2) - tube["center"])) < 1e-9 and abs(m["track_rms_m"] - 6450.0 * abs(math.log(0.2) - tube["center"])) < 1e-6 and m["in_tube"] == 1.0
    assert m["final_apo_err_m"] <= 0 and abs(m["final_apo_abs_err_m"] - abs(m["final_apo_err_m"])) < 1e-9 and abs(m["dv_decisions_sum"] - 0.15) < 1e-12
    P2 = P.copy(); P2.loc[3, "rho_peak_kg_m3"] *= 3.0; P2.loc[4, "rho_peak_kg_m3"] *= 5.0                       # Q_s 0.9 and 1.5: over 0.32 and 0.65
    m2 = run_metrics(P2, D, dict(total_dv_mps=0.15, burns=2, t_end_days=64.0, stop_reason="end_radius")); assert m2["v032"] == 2 and m2["v065"] == 2 and m2["v054"] == 2 and np.isnan(m2["track_rms_m"])
    short = P.iloc[:10].copy(); short["ra_m"] = 2.5e7; ms = run_metrics(short, None, dict(total_dv_mps=0.0, burns=0, t_end_days=10.0, stop_reason="safety_q")); assert not ms["reached"] and ms["final_apo_err_m"] > 0
    # C5 quantities: a flat flown path with a burn up (dv > 0 removes heat): net natural drift = +cancel (the burn hid a rise)
    Pf = _synth_passes(); Pf["q_w_cm2"] = 0.2 * np.exp(np.where(Pf["pass"] >= 6, -0.1, 0.0)); Df = pd.DataFrame(dict(target_pass=[6], delta_v_mps=[1.0]))
    mf = run_metrics(Pf, Df, dict(total_dv_mps=1.0, burns=1, t_end_days=64.0, stop_reason="end_radius")); assert abs(mf["net_nat_ln"] - (mf["cancel_sum_ln"] - 0.1)) < 1e-9 and mf["cancel_sum_ln"] > 0     # net natural drift = flown change (-0.1) + the heat the burn removed
    tab_up = dict(path=list(np.linspace(0, 2.0, 40)))                                                            # predicted drift +: heat rises, burn up (dv > 0) is WITH the plan
    tab_dn = dict(path=list(np.linspace(0, -2.0, 40)))
    assert run_metrics(Pf, Df, dict(total_dv_mps=1.0, burns=1, t_end_days=64.0, stop_reason="end_radius"), tab=tab_up)["oppose_dv_frac"] == 0.0
    assert run_metrics(Pf, Df, dict(total_dv_mps=1.0, burns=1, t_end_days=64.0, stop_reason="end_radius"), tab=tab_dn)["oppose_dv_frac"] == 1.0
    # C1
    body = {k: _run(3.0, days=64.0, passes=200) for k in seeds}; fo = {k: _run(3.9, days=70.0, passes=240) for k in seeds}
    ok = {k: _run(2.0, days=64.5, passes=201) for k in seeds}; bad = {k: _run(2.0, days=70.0, passes=220) for k in seeds}
    assert c1(ok, body, fo, rng)["pass"] and not c1(bad, body, fo, rng)["pass"]
    # C2: Poisson interval around E x n (E 0.5, n 16 -> 8), and the 0.65 rule
    rr = {k: _run(2.0, v054=(1 if i < 8 else 0)) for i, k in enumerate(seeds)}; assert c2(rr)["pass"] and c2(rr)["expected"] == 8.0
    assert not c2({k: _run(2.0, v054=0) for k in seeds})["pass"] and not c2({k: _run(2.0, v054=(1 if i < 8 else 0), v065=(1 if i == 0 else 0)) for i, k in enumerate(seeds)})["pass"]
    assert not c2({k: _run(2.0, v054=3) for k in seeds})["pass"]
    # C4
    assert c4({k: _run(2.0) for k in seeds}, 1.5, {k: 1.2 for k in seeds})["pass"] and not c4({k: _run(4.0) for k in seeds}, 1.5, {k: 1.2 for k in seeds})["pass"] and not c4({k: _run(1.0) for k in seeds}, 1.5, {k: 1.2 for k in seeds})["pass"]
    assert c4({k: _run(2.0) for k in seeds}, 1.5)["pass"] is None and c4({k: _run(2.0) for k in seeds}, 1.5)["upper_ok"]
    # C4abs (reported only): within +-0.3 m/s of the predicted floor; C4 itself unchanged
    assert c4abs({k: _run(1.5) for k in seeds}, 1.3)["pass"] and not c4abs({k: _run(1.7) for k in seeds}, 1.3)["pass"] and c4abs({k: _run(0.2) for k in seeds}, 0.0)["pass"] and not c4({k: _run(0.2) for k in seeds}, 0.0)["upper_ok"]
    assert abs(c4abs({k: _run(1.5) for k in seeds}, 1.3)["difference"] - 0.2) < 1e-12
    # C5
    good5 = {k: _run(2.0, net_nat_ln=1.0, cancel_sum_ln=1.1, oppose_dv_frac=0.05) for k in seeds}; bad5 = {k: _run(2.0, net_nat_ln=1.0, cancel_sum_ln=1.5, oppose_dv_frac=0.05) for k in seeds}
    bad5b = {k: _run(2.0, net_nat_ln=1.0, cancel_sum_ln=1.0, oppose_dv_frac=0.3) for k in seeds}
    assert c5(good5)["share_pass"] == 1.0 and c5(bad5)["share_pass"] == 0.0 and c5(bad5b)["share_pass"] == 0.0
    # C6
    assert c6({k: _run(2.0, vsum_ev=0.9) for k in seeds})["pass"] and not c6({k: _run(2.0, vsum_ev=0.8) for k in seeds})["pass"] and c6({k: _run(2.0) for k in seeds})["status"] == "NOT_COMPUTED"
    assert abs(vsum_ev([1, 2, 3], [1, 2, 3]) - 1.0) < 1e-12 and abs(vsum_ev([1, 2, 3, 9], [1, 2, 3, 0], [False, False, False, True]) - 1.0) < 1e-12
    # C7
    rs = np.random.default_rng(3); x = rs.normal(0.30, 0.188, 20000); assert abs(np.percentile(x, 99) - 0.737) < 0.03
    assert c7({1: _run(2.0, lnratio=np.full(500, 0.707))})["q99_ok"] and not c7({1: _run(2.0, lnratio=np.full(500, 0.80))})["q99_ok"] and c7({1: _run(2.0, lnratio=np.full(500, 0.707))})["tail_ok"] is None
    # ladder telescoping and the factorial estimator on a synthetic additive world with a geometry x feed-forward interaction
    rs = np.random.default_rng(11); base = {k: 4.0 + rs.normal(0, 0.3) for k in seeds}; eff = dict(C=0.5, N=0.4, Pff=0.3, Pt=0.2, gL1=0.1, gL3=0.9)
    dvs = {}
    for k in seeds:
        n2i = base[k]; l0 = n2i - eff["C"]; l1 = l0 - eff["N"]; l2 = l1 - eff["Pff"]; l3 = l2 - eff["Pt"]
        for name, val in (("N2I@ODY", n2i), ("L0@ODY", l0), ("L1@ODY", l1), ("L2@ODY", l2), ("L3@ODY", l3), ("L1@G1", l1 - eff["gL1"]), ("L3@G1", l3 - eff["gL3"]), ("BODY@ODY", n2i - 0.2), ("FO@ODY", n2i - 0.1)):
            dvs.setdefault(name, {})[k] = val + (0.0)
    at = attribution(dvs, rng); cN = at["components"]["N2I"]
    assert cN["telescopes"] and abs(cN["C"]["point"] - 0.5) < 1e-9 and abs(cN["N"]["point"] - 0.4) < 1e-9 and abs(cN["P_ff"]["point"] - 0.3) < 1e-9 and abs(cN["P_t"]["point"] - 0.2) < 1e-9
    assert abs(cN["P_g"]["point"] - 0.9) < 1e-9 and abs(cN["P_g_main"]["point"] - 0.5) < 1e-9 and abs(cN["P_g_interaction"]["point"] - 0.8) < 1e-9
    assert abs(cN["total"]["point"] - (0.5 + 0.4 + 0.3 + 0.2 + 0.9)) < 1e-9 and abs(cN["P"]["point"] - 1.4) < 1e-9 and abs(cN["P_share_of_total"] - 1.4 / 2.3) < 1e-9
    assert abs(at["components"]["BODY"]["C"]["point"] - 0.3) < 1e-9 and abs(at["components"]["FO"]["total"]["point"] - (2.3 - 0.1)) < 1e-9       # first rung replaced: BODY is 0.2 below N2-I, FO 0.1 below
    assert d2(at)["pass"] and d2(at)["share"] > 0.25 and cN["P"]["lo"] > 0
    # D2 fallbacks: P share below 25%, and P interval touching zero
    eff2 = dict(eff, Pff=0.0, Pt=0.0, gL1=0.0, gL3=0.05); dv2 = {n: {k: v for k, v in d.items()} for n, d in dvs.items()}
    for k in seeds:
        l1 = dv2["L1@ODY"][k]; dv2["L2@ODY"][k] = l1; dv2["L3@ODY"][k] = l1; dv2["L1@G1"][k] = l1; dv2["L3@G1"][k] = l1 - 0.05
    a2 = attribution(dv2, rng); assert not d2(a2)["pass"] and d2(a2)["fallback_framing"]
    assert attribution({k: v for k, v in dvs.items() if k != "L2@ODY"}, rng)["components"]["N2I"]["status"] == "MISSING_ARMS"
    # missing-seed rule: a seed absent from one arm is dropped from the decomposition
    dv3 = {n: dict(d) for n, d in dvs.items()}; del dv3["L2@ODY"][13]; assert attribution(dv3, rng)["components"]["N2I"]["n_seeds"] == 15
    # D1a / D1b
    l3 = {k: _run(1.5) for k in seeds}; bases = dict(BODY={k: _run(3.0) for k in seeds}, FO={k: _run(3.9) for k in seeds}, N2I={k: _run(4.2) for k in seeds})
    c1r = c1(l3, {k: _run(2.0, days=64.0, passes=200) for k in seeds}, {k: _run(2.0, days=70.0, passes=240) for k in seeds}, rng); c1r["pass"] = True
    c2r = dict(c2({k: _run(1.5, v054=(1 if i < 8 else 0)) for i, k in enumerate(seeds)}))
    r1 = d1(l3, bases, None, None, rng, c1r, c2r); assert r1["pass"] and all(v["hi"] < 1 for v in r1["ratios"].values())
    assert not d1({k: _run(3.5) for k in seeds}, bases, None, None, rng, c1r, c2r)["pass"] and not d1(l3, bases, None, None, rng, dict(c1r, **{"pass": False}), c2r)["pass"]
    # D5
    r5 = d5({k: _run(1.0) for k in seeds}, {k: _run(3.0) for k in seeds}, rng); assert abs(r5["diff_median"]["point"] + 2.0) < 1e-12 and abs(r5["ratio_median"]["point"] - 1 / 3) < 1e-12
    # D6: OL worse than L3 in dV and violations -> SUPPORTED; identical -> NULL; OL better -> AGAINST
    rs = np.random.default_rng(5); l3d = {k: _run(1.5 + rs.normal(0, .05), v032=1) for k in seeds}
    ol_worse = {k: _run(l3d[k]["dv"] + 0.8 + rs.normal(0, .05), v032=3) for k in seeds}; assert d6(ol_worse, l3d, rng)["verdict"] == "SUPPORTED"
    ol_same = {k: _run(l3d[k]["dv"] + rs.normal(0, .05), v032=1) for k in seeds}; assert d6(ol_same, l3d, rng)["verdict"] == "NULL"
    ol_better = {k: _run(l3d[k]["dv"] - 0.8 + rs.normal(0, .05), v032=0) for k in seeds}; assert d6(ol_better, l3d, rng)["verdict"] == "AGAINST"
    ol_viol_only = {k: _run(l3d[k]["dv"] + rs.normal(0, .05), v032=4) for k in seeds}; assert d6(ol_viol_only, l3d, rng)["verdict"] == "SUPPORTED"      # fewer violations at equal dV
    # D7: OL minus OL-R, every metric of R8.2 present and signed correctly
    olr = {k: _run(1.6, track_rms_m=800.0, track_max_m=2000.0, in_tube=0.95, final_apo_abs_err_m=5e3, final_apo_err_m=-5e3, v032=1, v065=0, burns=8) for k in seeds}
    ol = {k: _run(1.4, track_rms_m=1500.0, track_max_m=5000.0, in_tube=0.7, final_apo_abs_err_m=2e4, final_apo_err_m=2e4, v032=5, v065=1, burns=4) for k in seeds}
    r7 = d7(ol, olr, rng); assert set(r7) == {k for k, _ in D7_METRICS} and abs(r7["track_rms_m"]["point"] - 700.0) < 1e-9 and abs(r7["dv"]["point"] + 0.2) < 1e-9 and r7["in_tube"]["point"] < 0 and r7["burns"]["point"] == -4
    assert abs(r7["final_apo_abs_err_m"]["point"] - 1.5e4) < 1e-6 and r7["v032"]["point"] == 4 and r7["v065"]["point"] == 1
    # D8: replan deviation on synthetic replan.csv rows: 10 orbits, a priori burns at passes 3 and 7 (0.2 each); the replan moves the first to pass 4 (0.25) and skips the second
    ap_ = np.zeros(10); ex_ = np.zeros(10); ap_[2] = 0.2; ap_[6] = 0.2; ex_[3] = 0.25
    Rp = pd.DataFrame(dict(**{"pass": np.arange(2, 12)}, apriori_dv_mps=ap_, executed_dv_mps=ex_, dp_solve_s=0.01, apriori_build_s=[0.5] + [0.0] * 9))
    mr = replan_metrics(Rp, run_wall_s=100.0)
    assert mr["n_diverge"] == 3 and mr["first_diverge_pass"] == 4 and mr["n_burn_apriori"] == 2 and mr["n_burn_replan"] == 1 and mr["burn_count_diff"] == -1 and abs(mr["total_dv_diff"] - (0.25 - 0.4)) < 1e-12
    assert abs(mr["dv_diff_abs_sum"] - (0.2 + 0.25 + 0.2)) < 1e-12 and abs(mr["dv_diff_signed_sum"] - (0.25 - 0.4)) < 1e-12 and abs(mr["dv_diff_rms"] - math.sqrt((0.04 + 0.0625 + 0.04) / 10)) < 1e-12
    assert abs(mr["dp_total_s"] - (0.1 + 0.5)) < 1e-12 and abs(mr["dp_fraction"] - 0.006) < 1e-12 and abs(mr["dp_median_s"] - 0.01) < 1e-12
    same = replan_metrics(Rp.assign(executed_dv_mps=ap_)); assert same["n_diverge"] == 0 and np.isnan(same["first_diverge_pass"]) and same["dv_diff_rms"] == 0.0 and same["burn_count_diff"] == 0
    tiny = replan_metrics(Rp.assign(executed_dv_mps=ap_ + 1e-4)); assert tiny["n_diverge"] == 0                       # below min_burn is no burn
    runs8 = {k: dict(_run(2.0), replan=replan_metrics(Rp, 100.0 + k)) for k in seeds[:4]}; runs8[seeds[4]] = dict(_run(2.0), replan=same)
    r8 = d8(runs8, rng); assert r8["n_flights"] == 5 and r8["n_no_divergence"] == 1 and r8["n_diverge"]["point"] == 3 and r8["first_diverge_pass"]["point"] == 4 and set(D8_METRICS) <= set(r8)
    assert d8_ol(runs8)["n_burn_apriori_median"] == 2
    # loader on a synthetic results tree (planner summary + Stage O summary), roles from the gate table
    td = task_dir(); roles = gate_roles(td) if td else {"ODY": "ODY", "G5863": "G1"}; assert roles["G5863"] == "G1" and roles["ODY"] == "ODY"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp); pl = tmp / "planner"; so = tmp / "stageo"
        rows = []
        for arm, gid, root, sub in (("L3", "G5863", pl, "runs"), ("OL", "ODY", pl, "runs"), ("F-O", None, so, "mars/runs")):
            for s in (13, 14):
                rid = f"{gid or 'ODY'}_X_pert_ODY_s{s}" if gid else f"FO_fixed1000_pert_ODY_s{s}"; d = root / sub / rid; d.mkdir(parents=True)
                _synth_passes().to_csv(d / "passes.csv", index=False)
                if root is pl:
                    Rp.to_csv(d / "replan.csv", index=False)
                pd.DataFrame(dict(target_pass=[5], delta_v_mps=[0.1])).to_csv(d / "decisions.csv", index=False)
                row = dict(run_id=rid, arm=arm, status="ok", total_dv_mps=0.1, burns=1, t_end_days=64.0, stop_reason="end_radius"); row.update(geometry=gid) if gid else None
                rows.append((root, row))
        for root in (pl, so):
            rr_ = [r for rt, r in rows if rt == root]; pd.DataFrame(rr_).to_csv(root / "summary_held_x.csv", index=False)
        R = collect([pl, so], roles, phases=("held",))
        assert set(R) == {"L3@G1", "OL@ODY", "FO@ODY"} and sorted(R["L3@G1"]) == [13, 14] and R["FO@ODY"][13]["dv"] == 0.1 and R["OL@ODY"][14]["burns"] == 1 and R["L3@G1"][13]["replan"]["n_diverge"] == 3 and "replan" not in R["FO@ODY"][13]
    print("planner_analysis selftest ok: statistics, per-run metrics (Q_s on V_rel, tube tracking, final apoapsis error, natural path and C5 quantities), C1, C2, C4, C5, C6, C7, ladder telescoping, factorial estimator "
          "(main effect, interaction), D1a / D1b, D2 (incl. fallback), D3, D5, D6 (supported / null / against), D7 (all R8.2 metrics), missing-seed rule, loader")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "metrics", "report"]); ap.add_argument("--results-root"); ap.add_argument("--stageo-root"); ap.add_argument("--task-dir"); ap.add_argument("--phases", default="held")
    a = ap.parse_args()
    if a.cmd == "selftest":
        return selftest()
    td = task_dir(a.task_dir); roles = gate_roles(td)
    preds = json.loads((td / CFG["planner"]["predictions"]["path"]).read_text()) if (td / CFG["planner"]["predictions"]["path"]).exists() else None
    tabs = json.loads((td / CFG["planner"]["ground_table"]["path"]).read_text()) if (td / CFG["planner"]["ground_table"]["path"]).exists() else None
    a.results_root or ap.error("--results-root required")
    roots = [pathlib.Path(a.results_root) / CFG["stagef"]["results_subdir"]] + ([pathlib.Path(a.stageo_root)] if a.stageo_root else [pathlib.Path(a.results_root) / "stageo"])
    runs = collect(roots, roles, preds, tabs, tuple(a.phases.split(",")))
    if a.cmd == "metrics":
        for k, v in sorted(runs.items()):
            for s, r in sorted(v.items()):
                print(k, s, {m: (round(x, 4) if isinstance(x, float) else x) for m, x in r.items() if m != "lnratio"})
        return
    print(json.dumps(report(runs, preds), indent=1, default=float))


if __name__ == "__main__":
    main()
