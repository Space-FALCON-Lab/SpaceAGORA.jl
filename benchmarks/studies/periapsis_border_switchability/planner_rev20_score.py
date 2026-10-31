"""Revision 20 flight scoring (analysis only; reads flight directories, no Julia, writes no other file).
Block rev20olt (Revision 20b, arm OLT at Odyssey, 13 flights): per-flight success/hard/error/dV like results_rev17.md, the identity gate (OLT L0 s16 vs Rev 17 OL L0 s16), the 20-flight OLT table (13 flown + the 8
determined OL outcomes, L0 s16 counted once), P20-7 and the interpretation-rule split.
Block rev20val (Revision 20c, OLRT3-W and OLRT3 at Odyssey/G1, 29 flights): success (Rev 14/13 rule), arrival time against t*, window (W_max = 3 final-orbit periods), dV with the pre-burn and FB(h_des), AC2/AC3/AC4 per
cell, pairing against Rev 17 OLRT3 (Odyssey) and the new OLRT3 (G1), P20-9, and AC7 (regression flight OLRT3@ODY@L4 s14 against Rev 17, bit for bit).
Revision 20d (no final periapsis maneuver): the main total dV is the aerobraking dV to the last counted pass plus the pre-burn; FB(h_des) moved to a side table. The window is scored first-passage (as registered, P20-9
unchanged) and with the wait rule of planner_rev20d.wait_rule (extra final orbits while the apoapsis error stays within 25 km, next drops = the last flown drop, ASSUMED). AC5 corridor columns: last-pass mean q in
[0.05, 0.25] W/cm2 (SOURCED), last-pass Q_s <= 0.65, final periapsis altitude.
Missing or unfinished flights (a flight is finished when passes/decisions/measured/replan.csv all exist) are listed and everything else is scored.
Revision 21 (WP21-B): rendezvous position error (planner_rev21_position: great-circle angle between the flight apoapsis direction at the arrival and the J2-propagated nominal apoapsis direction at t*; decomposition (a) timing, (b) predictor residual,
(c) plan-versus-flight history; tol_pos_deg = 1.0 ASSUMED, 0.5 and 2.0 shown) is added as the last columns of the rev20val per-flight and per-cell tables and as Part C (position for rev20val, rev20olt, rev17 and every directory given by --blocks).
--blocks B1 B2 ...: any block directory under the results root (e.g. rev21c rev21v): scored with the same flight metrics (success, window against the nominal t*, dV) and the position columns, without the config flight list.
Usage: python3 -I planner_rev20_score.py --results-root R/planner --out T/results_rev20_flights.md [--json T/results_rev20_flights.json] [--blocks rev21c rev21v]   (json defaults to the --out path with .json)"""
import argparse, glob, json, math, pathlib, re, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_rev20_window as W
import planner_rev20d as D20
import planner_rev21_position as PS

MU = A.MU; TARGET = 7.408e6; VE = 220.0 * 9.80665; H_DES_M = 250e3; W_MAX_PERIODS = 3.0; TOL_KM = 25.0
LEVELS = ["L0", "L1", "L2", "L3", "L4"]; SEEDS = [13, 14, 15, 16]
DETERMINED = {("L0", 16), ("L1", 13), ("L2", 13), ("L3", 14), ("L4", 13), ("L4", 14), ("L4", 15), ("L4", 16)}      # Revision 20b list (OL stopped on safety_q > 107.5 km above the target)
GEOM = {"ODY": "ODY", "G5863": "G1", "G0975": "G0975", "G0063": "G0063", "G0390": "G0390", "G1216": "G1216", "G0528": "G0528", "G0890": "G0890", "G0975M": "G0975M", "G0063M": "G0063M", "G0486": "G0486", "G1230": "G1230", "G1038": "G1038", "G0140_AODY": "G0140_AODY", "G0042_A2000": "G0042_A2000", "G3257_AODY": "G3257_AODY", "G3161_A2000": "G3161_A2000"}; GID = {"ODY": "ODY", "G1": "G5863", "G0975": "G0975", "G0063": "G0063", "G0390": "G0390", "G1216": "G1216", "G0528": "G0528", "G0890": "G0890", "G0975M": "G0975M", "G0063M": "G0063M", "G0486": "G0486", "G1230": "G1230", "G1038": "G1038", "G0140_AODY": "G0140_AODY", "G0042_A2000": "G0042_A2000", "G3257_AODY": "G3257_AODY", "G3161_A2000": "G3161_A2000"}
FILES = ("passes.csv", "decisions.csv", "measured.csv", "replan.csv")
CFG = tomllib.load(open(HERE / "planner_config.toml", "rb"))["planner"]
OFFSET_D = {k: float(v) for k, v in CFG["guard"]["t_offset_d"].items()}       # ASSUMED aim offsets (days), Revision 20c


def nan_none(o):
    if isinstance(o, dict): return {k: nan_none(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [nan_none(v) for v in o]
    if isinstance(o, (np.floating, float)): return None if not math.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    return o


def f2(x, p=2):
    return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{p}f}"


def med(x):
    x = [v for v in x if v is not None and np.isfinite(v)]; return float(np.median(x)) if x else float("nan")


def table(head, rows):
    return ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]


def parse_id(name):
    m = re.match(r"^(\w+?)_PL_(.+)_pert_[A-Za-z0-9]+_s(\d+)(?:_(L\d\w*))?$", name)
    if not m: return None
    return dict(gid=m.group(1), geom=GEOM.get(m.group(1), m.group(1)), arm=re.sub(r"[^A-Za-z0-9]", "", m.group(2)).upper(), seed=int(m.group(3)), level=m.group(4) or "")


def arm_key(arm): return re.sub(r"[^A-Za-z0-9]", "", arm).upper()


def load_summary(root, pat):
    S = {}
    for f in sorted(glob.glob(str(root / pat))):
        for r in pd.read_csv(f).itertuples(): S[r.run_id] = r._asdict()
    return S


def finished(d): return d is not None and all((d / f).exists() for f in FILES)


def infer_stop(reached, qmax):
    return "end_radius" if reached else ("safety_q" if qmax > 0.65 else "horizon_or_other")


def first_burn_dv(ra, rp, rp_des):
    va = lambda r_p: math.sqrt(MU * (2 / ra - 2 / (ra + r_p)))      # planner_rev18_flights.first_burn_dv
    return va(rp_des) - va(rp)


def flight(d, srow, dv_pre_fallback=float("nan"), target=TARGET, onesided=False):
    """Metrics of one finished flight directory. srow = summary row dict or None (stop inferred then). target = the target apoapsis radius (m, default 7,408,000, SOURCED); onesided (Revision 25, ASSUMED sign): the closing pass is the first with
    r_a <= target + TOL_KM and the apoapsis error must lie in [0, +TOL_KM] (the circle at the target radius is reachable only when r_a >= target); default: first pass with r_a <= target, error in [-TOL_KM, 0] (registered)."""
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= target + (TOL_KM * 1e3 if onesided else 0.0)))[0]
    n = int(hit[0] + 1) if len(hit) else len(P); i = n - 1; Pn = P.iloc[:n]; fin = ra[:n][np.isfinite(ra[:n])]; m = P.mass_kg.values.astype(float); reached = bool(len(hit))
    qs = 0.5 * Pn.rho_peak_kg_m3.values * Pn.vrel_mps.values ** 3 / 1e4; qmax = float(np.nanmax(qs))
    stop = str(srow["stop_reason"]) if srow is not None and isinstance(srow.get("stop_reason"), str) else infer_stop(reached, qmax); stop_src = "summary" if srow is not None and isinstance(srow.get("stop_reason"), str) else "inferred"
    hard = bool(qmax > 0.65 or "safety_q" in stop); err = (fin[-1] - target) / 1e3 if len(fin) else float("nan"); rlo, rhi = (0.0, TOL_KM) if onesided else (-TOL_KM, 0.0)
    D = pd.read_csv(d / "decisions.csv"); dv = D[D.target_pass <= n].delta_v_mps.fillna(0).values; nz = dv[np.abs(dv) > 0]
    a_m = float(P.a_m.iloc[i]); e = float(P.e.iloc[i]); Tn = 2 * math.pi * math.sqrt(a_m ** 3 / MU); t_arr = float(P.time_s.iloc[i]) + 0.5 * Tn
    rp_f = float(P.rp_m.iloc[i]); alt_f = float(P.alt_m.iloc[i]); ra_f = float(P.ra_m.iloc[i]); om = np.unwrap(P.omega_rad.values.astype(float)); inc = float(P.i_rad.iloc[i])
    dv_pre = float("nan"); lvl_pre = "n/a"
    if srow is not None and "ic_dv_pre_mps" in srow and np.isfinite(srow["ic_dv_pre_mps"]): dv_pre = float(srow["ic_dv_pre_mps"]); lvl_pre = "summary"
    elif srow is None and np.isfinite(dv_pre_fallback): dv_pre = dv_pre_fallback; lvl_pre = "IC file"
    elif srow is not None and "ic_dv_pre_mps" not in srow: dv_pre = 0.0; lvl_pre = "no IC column (0)"
    fb = first_burn_dv(ra_f, rp_f, rp_f - alt_f + H_DES_M) if reached else float("nan"); dv_sum = VE * math.log(m[0] / m[i])
    last_drop = float((ra[i - 1] - ra[i]) / 1e3) if i >= 1 and np.isfinite(ra[i - 1]) and np.isfinite(ra[i]) else float("nan")      # km, last flown apoapsis drop (Revision 20d wait rule)
    q_last = float(P.q_w_cm2.iloc[i]); qs_last = float(qs[-1])
    rate = float(W.apsidal_per_orbit(ra_f, rp_f, inc) / Tn) if reached else float("nan")           # rad/s, secular J2 rate of the final orbit
    R = pd.read_csv(d / "replan.csv"); tt = R["ttg_tt_s"].dropna() if "ttg_tt_s" in R else pd.Series(dtype=float)
    return dict(n=n, reached=reached, hard=hard, err=err, ok=bool(reached and not hard and rlo <= err <= rhi), apo_ok=bool(reached and not hard and rlo <= err <= rhi), qmax=qmax, stop=stop, stop_src=stop_src, passes=len(P),
                dv_sum=dv_sum, dv_net=abs(float(dv.sum())), dv_pre=dv_pre, dv_pre_src=lvl_pre, fb=fb, total=dv_sum + (dv_pre if np.isfinite(dv_pre) else float("nan")), total_fb=dv_sum + (dv_pre if np.isfinite(dv_pre) else float("nan")) + fb, last_drop_km=last_drop, q_last=q_last, qs_last=qs_last,
                alt_f_km=alt_f / 1e3, corridor=bool(0.05 <= q_last <= 0.25), qs_last_ok=bool(qs_last <= 0.65), burns=int(len(nz)),
                reversals=int(np.sum(np.sign(nz[1:]) != np.sign(nz[:-1]))) if len(nz) > 1 else 0, soft=int((qs > 0.54).sum()), days=float(P.time_s.iloc[i] / 86400), t_arr_d=t_arr / 86400, T_final_s=Tn, wmax_d=W_MAX_PERIODS * Tn / 86400,
                omega_f_deg=math.degrees(float(om[i])), omega_rot_deg=math.degrees(float(om[i] - om[0])), apse_rate_deg_per_d=math.degrees(rate) * 86400 if reached else float("nan"),
                tt_aim_d=float(tt.iloc[0]) / 86400 if len(tt) else float("nan"), tt_spread_s=float(tt.max() - tt.min()) if len(tt) else float("nan"))


def compare_passes(da, db):
    """Bit-for-bit comparison of the shared passes.csv columns (read as text). Returns dict(equal, ncols, rows_a, rows_b, first=(row, col, a, b))."""
    a = pd.read_csv(da / "passes.csv", dtype=str, keep_default_na=False); b = pd.read_csv(db / "passes.csv", dtype=str, keep_default_na=False); cols = [c for c in a.columns if c in b.columns]
    out = dict(ncols=len(cols), only_a=[c for c in a.columns if c not in b.columns], only_b=[c for c in b.columns if c not in a.columns], rows_a=len(a), rows_b=len(b), equal=False, first=None)
    k = min(len(a), len(b)); X = a[cols].iloc[:k].values != b[cols].iloc[:k].values
    if X.any():
        r, c = np.argwhere(X)[0]; out["first"] = (int(r), cols[c], a[cols].iloc[r, c], b[cols].iloc[r, c])
    out["equal"] = bool(not X.any() and len(a) == len(b)); return out


def expected(block):
    """Run-id keys (arm_key, geom, level, seed) of a block from planner_config.toml cell_seeds."""
    for blk in CFG["terminal_flight"]:
        if blk["block"] == block:
            out = []
            for cell, seeds in blk["cell_seeds"].items():
                arm, tok, lv = cell.split("@"); out += [(arm, GEOM.get(GID[tok], tok), lv, s) for s in seeds]
            return out
    raise KeyError(block)


def scan(root, block):
    """Directories of a block keyed by (arm, geom, level, seed)."""
    D = {}
    for p in sorted((root / block).glob("*")):
        k = parse_id(p.name)
        if k: D[(k["arm"], k["geom"], k["level"], k["seed"])] = p
    return D


def ic_dv_pre(geom):
    try:
        f = next(q for q in (HERE.parents[i] / "tasks" / "2026-10-15_periapsis-region-border-targeting" / CFG["ic"]["path"] for i in (2, 5)) if q.exists()); J = json.load(open(f))
        c = [p for p in J["plans"] if p["geom"] == geom and p["scenario"] == "S-fixed" and p["fb_term"] is False and p["horizon"] == "arrival" and abs(p["delta_ln"]) < 1e-12 and abs(p["h_des_m"] - 250000.0) < 1e-6 and p["t_star_d"] is None]
        return float(c[0]["dv_pre_mps"]) if len(c) == 1 else float("nan")
    except Exception:
        return float("nan")


def cellname(k): return f"{k[0]}@{k[1]}@{k[2]}"


POSKEYS = ("err_deg", "arc_km", "a_timing_deg", "b_pred_deg", "c_hist_deg", "closure_deg", "a_analytic_deg")      # planner_rev21_position.flight_position fields kept on a flight record


def posrow(f): return [f2(f.get("err_deg"), 2), f2(f.get("arc_km"), 0), f2(f.get("a_timing_deg"), 2), f2(f.get("b_pred_deg"), 2), f2(f.get("c_hist_deg"), 2)]


POSHEAD = ["position error (deg)", "arc at r_a (km)", "(a) timing (deg)", "(b) predictor (deg)", "(c) history (deg)"]


def posstats(fs):
    e = np.array([f["err_deg"] for f in fs if np.isfinite(f.get("err_deg", float("nan")))])
    if not len(e): return dict(n=0, med=float("nan"), p90=float("nan"), **{f"within_{t:g}": float("nan") for t in PS.TOL_DEG})
    return dict(n=len(e), med=float(np.median(e)), p90=float(np.percentile(e, 90)), **{f"within_{t:g}": float((e <= t).mean()) for t in PS.TOL_DEG})


def posstat_row(fs):
    p = posstats(fs); return [f"{p['n']}", f2(p["med"], 2), f2(p["p90"], 2), f"{p['within_0.5']:.2f}", f"{p['within_1']:.2f}", f"{p['within_2']:.2f}"] if p["n"] else ["0", "n/a", "n/a", "n/a", "n/a", "n/a"]


POSSTATHEAD = ["position n", "median pos err (deg)", "p90 (deg)", "within 0.5 deg", "within 1.0 deg (ASSUMED tol)", "within 2.0 deg"]


def selftest():
    """Scorer selftest (no flight data needed): run-id parsing (L2p, W3/WF3, Revision 23 and 25 tokens, rejects), target radii, default and one-sided tolerance rules and closing-pass detection on a synthetic flight, and the circular-target wrapper
    arguments. Also runs planner_rev21_position and planner_rev22_rendezvous selftests."""
    import tempfile
    for nm, arm, lv, sd, gm in [("G0975_PL_OLRT3_pert_ODY_s13_L2p", "OLRT3", "L2p", 13, "G0975"), ("ODY_PL_OLRT3-WF3_pert_ODY_s14_L0", "OLRT3WF3", "L0", 14, "ODY"), ("G0486_PL_OLRT3-W3_pert_ODY_s16_L0", "OLRT3W3", "L0", 16, "G0486"),
                                ("G0063M_PL_BODYT_pert_ODY_s15_L0", "BODYT", "L0", 15, "G0063M"), ("G1038_PL_OLRT3-W3_pert_ODY_s13_L0", "OLRT3W3", "L0", 13, "G1038"), ("G5863_PL_OLRT3_pert_ODY_s1", "OLRT3", "", 1, "G1")]:
        k = parse_id(nm); assert k and (arm_key(k["arm"]), k["geom"], k["level"], k["seed"]) == (arm, gm, lv, sd), (nm, k)
    assert parse_id("junk_dir") is None and parse_id("G0975_PL_OLRT3_pert_ODY_s13_Lx") is None and GEOM.get("ZZZ") is None
    assert PS.target_radius("ODY") == 7.408e6 and abs(PS.target_radius("G0486") - 5396190.0) < 1 and abs(PS.target_radius("G1230") - 4396190.0) < 1 and abs(PS.target_radius("G1038") - 3896190.0) < 1
    PS.TARGET_OVERRIDE["G0063"] = 7.0e6; assert PS.target_radius("G0063") == 7.0e6; del PS.TARGET_OVERRIDE["G0063"]
    with tempfile.TemporaryDirectory() as td:
        d = pathlib.Path(td); n = 6; ra = np.array([7.6e6, 7.5e6, 7.43e6, 7.41e6, 7.39e6, 7.30e6])      # default target 7.408e6: closes at pass 5 (7.390e6, 18 km below)
        P = pd.DataFrame({"pass": np.arange(1, n + 1), "time_s": np.arange(n) * 12000.0 + 6000.0, "ra_m": ra, "a_m": (ra + 3.5e6) / 2, "e": (ra - 3.5e6) / (ra + 3.5e6), "i_rad": 1.6, "omega_rad": 1.9, "raan_rad": 0.5, "rp_m": 3.5e6, "alt_m": 1.0e5, "mass_kg": 400.0,
                          "rho_peak_kg_m3": 1e-8, "vrel_mps": 4000.0, "q_w_cm2": 0.1}); P.to_csv(d / "passes.csv", index=False); pd.DataFrame({"target_pass": [2], "delta_v_mps": [0.1]}).to_csv(d / "decisions.csv", index=False)
        pd.DataFrame({"t_s": [0.0]}).to_csv(d / "measured.csv", index=False); pd.DataFrame({"x": [0.0]}).to_csv(d / "replan.csv", index=False)
        f = flight(d, None); assert f["reached"] and f["n"] == 5 and abs(f["err"] + 18.0) < 1e-9 and f["ok"], f
        f = flight(d, None, target=7.38e6); assert f["reached"] and f["n"] == 6 and abs(f["err"] + 80.0) < 1e-9 and not f["ok"], f           # default rule: first pass below 7.38e6 is pass 6, 80 km below: fails
        f = flight(d, None, target=7.38e6, onesided=True); assert f["reached"] and f["n"] == 5 and abs(f["err"] - 10.0) < 1e-9 and f["ok"], f     # one-sided: first pass <= target + 25 km is pass 5 (7.39e6, +10 km): ok
        f = flight(d, None, target=7.35e6, onesided=True); assert f["n"] == 6 and not f["ok"], f                                                      # 7.30e6 is 50 km below the target and 7.39e6 is above target + 25 km: not in [0, 25]
    PS.selftest()
    import planner_rev22_rendezvous as R22; assert R22.selftest()
    print("planner_rev20_score selftest ok")
    return True


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root"); ap.add_argument("--out"); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--json"); ap.add_argument("--blocks", nargs="*", default=[]); ap.add_argument("--one-sided", nargs="*", default=[], help="blocks scored with the one-sided apoapsis tolerance (r_a in [target, target + 25 km], ASSUMED sign)")
    ap.add_argument("--target-radius", nargs="*", default=[], help="TOKEN=metres overrides of the target apoapsis radius (default: table/target_orbits.json, else 7,408,000)"); a = ap.parse_args()
    for kv in a.target_radius: PS.TARGET_OVERRIDE[kv.split("=")[0]] = float(kv.split("=")[1])
    if a.selftest: return selftest()
    if not a.results_root or not a.out: ap.error("--results-root and --out are required (unless --selftest)")
    root = pathlib.Path(a.results_root); out_md = pathlib.Path(a.out); out_js = pathlib.Path(a.json) if a.json else out_md.with_suffix(".json"); J = {}; L = []
    S17, S20o, S20v = load_summary(root, "summary_rev17_*.csv"), load_summary(root, "summary_rev20olt*.csv"), load_summary(root, "summary_rev20val*.csv")

    # ---------------- Rev 17 flights (OL, OLRT3 are the comparators) ----------------
    D17 = scan(root, "rev17"); F17 = {}
    for k, d in D17.items():
        if finished(d): F17[k] = flight(d, S17.get(d.name))
    ol = lambda lv, s: F17.get(("OL", "ODY", lv, s)); olrt3 = lambda lv, s: F17.get(("OLRT3", "ODY", lv, s))

    # ---------------- rev20olt ----------------
    exp_o = expected("rev20olt"); Do = scan(root, "rev20olt"); Fo = {}; err_o = []; miss_o = []
    for k in exp_o:
        d = Do.get((arm_key(k[0]),) + k[1:]); s = S20o.get(d.name) if d else None
        if s is not None and s.get("status") != "ok": err_o.append((k, s.get("status"), s.get("error")))
        elif finished(d): Fo[k] = flight(d, s)
        else: miss_o.append(k)
    L += ["# Results: Revision 20 flight blocks (rev20olt: OLT; rev20val: OLRT3-W), Odyssey and G1, held seeds 13-16", "",
          "Scored by planner_rev20_score.py from the flight directories under `rev20olt/`, `rev20val/`, `rev17/` and the summary CSVs. Success = reached, final apoapsis error in [-25, 0] km (ASSUMED tolerance), no hard violation (pass-max Q_s on V_rel > 0.65 W/cm2 or a safety_q stop; Revision 13); dV = dV to the last counted pass = ve ln(m0/m_last) (Revision 17a). "
          "A flight is scored once passes/decisions/measured/replan.csv all exist; the stop reason is read from the summary CSV when it exists and otherwise inferred from the data (end_radius if the target apoapsis was reached, else safety_q if pass-max Q_s > 0.65). 4 seeds per cell: preliminary.", "",
          "Labels. MEASURED: read directly from the flight logs (final apoapsis error, Q_s, burns, arrival time, dV, stop). DERIVED: computed from measured values with a stated formula (t_arr = time at the closing pass + half the final-orbit period, W_max, FB, J2 apse rotation, nominal t* = logged aim t* minus the aim offset, counts and differences). "
          "ASSUMED: the 25 km tolerance, W_max = 3 final-orbit periods, h_des = 250 km, the 0.15 d Odyssey aim offset, and the AC thresholds scaled to 4 flights per cell (below).", "", "## Part A. rev20olt: OLT at Odyssey", ""]
    mlist = lambda ks: ", ".join(f"{cellname(k)} s{k[3]}" for k in ks) or "none"
    L += [f"Finished flights: {len(Fo)} of {len(exp_o)} (MEASURED). Missing or running: {mlist(miss_o)}. Run errors: {'; '.join(f'{cellname(k)} s{k[3]}: {s} {e}' for k, s, e in err_o) or 'none'}.", ""]
    rows = []
    for k in exp_o:
        f = Fo.get(k)
        if f: rows.append([cellname(k).replace("@", " "), k[3], "OK" if f["ok"] else "FAIL" + (" HARD" if f["hard"] else ""), f2(f["err"], 1), f2(f["dv_sum"]), f2(f["dv_net"]), f["burns"], f["soft"], f["stop"] + ("" if f["stop_src"] == "summary" else " (inferred)"), f2(f["days"], 1)])
        else: rows.append([cellname(k).replace("@", " "), k[3], "missing", "", "", "", "", "", "", ""])
    L += table(["arm level", "seed", "result", "final error (km)", "sum dV (m/s)", "net dV (m/s)", "burns", "passes Q_s>0.54", "stop", "days"], [[r[0].replace("ODY ", "")] + r[1:] for r in rows]) + [""]
    # identity gate
    gk = ("OLT", "ODY", "L0", 16); gd = Do.get(("OLT", "ODY", "L0", 16)); od = D17.get(("OL", "ODY", "L0", 16)); gate = dict(status="NOT RUN (flight missing)")
    if finished(gd) and finished(od):
        c = compare_passes(gd, od); fg, fo = Fo.get(gk), F17.get(("OL", "ODY", "L0", 16)); same = bool(fg and fo and fg["stop"] == fo["stop"])
        gate = dict(status="PASS" if c["equal"] and same else "FAIL", cmp=c, stop_olt=fg["stop"] if fg else None, stop_ol=fo["stop"] if fo else None, same_stop=same)
        L += ["### Identity gate (Revision 20b): OLT L0 s16 against Rev 17 OL L0 s16", "",
              f"**{gate['status']}**: passes.csv shared columns {c['ncols']}, rows OLT {c['rows_a']} / OL {c['rows_b']}, all equal as text: {c['equal']}; first difference: {c['first'] or 'none'}; columns only in OLT: {c['only_a'] or 'none'}; only in OL: {c['only_b'] or 'none'}; stop OLT {gate['stop_olt']} ({fg['stop_src']}), OL {gate['stop_ol']} ({fo['stop_src']}), same stop: {same} (MEASURED).", ""]
    else:
        L += ["### Identity gate (Revision 20b): OLT L0 s16 against Rev 17 OL L0 s16", "", "**NOT RUN**: OLT L0 s16 or the Rev 17 OL L0 s16 flight is not finished.", ""]
    J["rev20olt"] = dict(finished=len(Fo), expected=len(exp_o), missing=[list(k) for k in miss_o], errors=[[list(k), s, e] for k, s, e in err_o], gate=gate, flights={cellname(k) + f"@s{k[3]}": v for k, v in Fo.items()})
    # 20-flight table
    L += ["### 20-flight OLT table (13 OLT flights plus the 8 determined OL Revision 17 outcomes; L0 s16 is both and counted once) ", ""]
    T20 = {}; rows = []
    for lv in LEVELS:
        for s in SEEDS:
            det = (lv, s) in DETERMINED; fo_ = Fo.get(("OLT", "ODY", lv, s)); fl = ol(lv, s)
            if fo_: T20[(lv, s)] = ("OLT flown" + (" (gate)" if det else ""), fo_)
            elif det and fl: T20[(lv, s)] = ("determined (OL Rev 17)", fl)
            else: T20[(lv, s)] = ("missing", None)
    sm = lambda f: "missing" if f is None else ("OK" if f["ok"] else "FAIL" + (" HARD" if f["hard"] else ""))
    per = {}
    for lv in LEVELS:
        ks = [(lv, s) for s in SEEDS]; sc = [T20[k][1] for k in ks if T20[k][1]]; so = [ol(lv, s) for s in SEEDS]; s3 = [olrt3(lv, s) for s in SEEDS]
        cnt = lambda fs: sum(1 for f in fs if f and f["ok"])
        per[lv] = dict(OL=cnt(so), OLT=cnt(sc), OLT_scored=len(sc), OLRT3=cnt(s3), OL_hard=sum(1 for f in so if f and f["hard"]), OLT_hard=sum(1 for f in sc if f["hard"]), OLRT3_hard=sum(1 for f in s3 if f and f["hard"]))
        rows.append([lv] + [f"{s}: {sm(T20[(lv, s)][1])} [{T20[(lv, s)][0]}]" for s in SEEDS])
    L += table(["level"] + [f"seed {s}" for s in SEEDS], rows) + [""]
    tot = {k: sum(per[lv][k] for lv in LEVELS) for k in ("OL", "OLT", "OLT_scored", "OLRT3", "OL_hard", "OLT_hard", "OLRT3_hard")}
    up = [(lv, s) for lv in LEVELS for s in SEEDS if (lv, s) not in DETERMINED]; upF = [Fo[("OLT", "ODY", lv, s)] for lv, s in up if ("OLT", "ODY", lv, s) in Fo]; nsu = sum(f["ok"] for f in upF); nhu = sum(f["hard"] for f in upF); rem = 12 - len(upF)
    if len(upF) == 12: p7 = "HELD" if nsu >= 9 and nhu <= 1 else ("FALSIFIED" if nsu < 7 or nhu >= 3 else "NOT MET, NOT FALSIFIED")
    else: p7 = "FALSIFIED (already decided)" if nhu >= 3 or nsu + rem < 7 else ("PENDING: " + ("held if the rest do not add a hard violation and 9+ succeed" if nhu <= 1 else "undecided"))
    pr = []
    for lv in LEVELS:
        p = per[lv]; pr.append([lv, f"{p['OL']}/4", f"{p['OLT']}/4 ({p['OLT_scored']} scored)", f"{p['OLRT3']}/4", f"{p['OLT'] - p['OL']:+d}", f"{p['OLRT3'] - p['OLT']:+d}", f"{p['OL_hard']} / {p['OLT_hard']} / {p['OLRT3_hard']}"])
    pr.append(["total", f"{tot['OL']}/20", f"{tot['OLT']}/20 ({tot['OLT_scored']} scored)", f"{tot['OLRT3']}/20", f"{tot['OLT'] - tot['OL']:+d}", f"{tot['OLRT3'] - tot['OLT']:+d}", f"{tot['OL_hard']} / {tot['OLT_hard']} / {tot['OLRT3_hard']}"])
    L += ["### Interpretation-rule split per level (successes out of 4; DERIVED)", ""] + table(["level", "OL (Rev 17)", "OLT", "OLRT3 (Rev 17)", "terminal-aim part: OLT - OL", "replanning part: OLRT3 - OLT", "hard violations OL / OLT / OLRT3"], pr) + [""]
    L += [f"Upper bound if the identity gate passes (Revision 20b, DERIVED): OLT at most 3, 3, 3, 3, 0 of 4 at L0-L4 (12 of 20); scored OLT total now {tot['OLT']} with {tot['OLT_scored']} of 20 cells scored.", "",
          f"### P20-7 (disposition exactly as registered): of the 12 undetermined OLT flights at least 9 succeed and at most 1 has a hard violation; falsified if fewer than 7 succeed or 3 or more have a hard violation.", "",
          f"**{p7}**: {len(upF)} of 12 undetermined flights scored; {nsu} succeed, {nhu} hard violations (MEASURED)." + ("" if len(upF) == 12 else f" {rem} not yet scored."), ""]
    if p7 == "HELD": L += [f"Interpretation rule: P20-7 holds, so the OL-to-OLRT3 success difference ({tot['OLRT3'] - tot['OL']:+d} of 20) splits into a terminal-aim part (OLT minus OL, {tot['OLT'] - tot['OL']:+d}) and a replanning part (OLRT3 minus OLT, {tot['OLRT3'] - tot['OLT']:+d}); the replanning part is attributed to heating protection under table error.", ""]
    elif p7.startswith("PENDING"): L += ["Interpretation rule: provisional until all 12 undetermined flights are scored; the split above is the current state, and it is attributed only if P20-7 holds.", ""]
    else: L += ["Interpretation rule: P20-7 does not hold, so the attribution of the OLRT3-over-OL difference is reported as UNRESOLVED; the split above is descriptive only.", ""]
    J["rev20olt"].update(table20={f"{lv}_s{s}": dict(source=T20[(lv, s)][0], ok=(T20[(lv, s)][1] or {}).get("ok"), hard=(T20[(lv, s)][1] or {}).get("hard")) for lv in LEVELS for s in SEEDS}, per_level=per, totals=tot, p20_7=dict(disposition=p7, scored=len(upF), success=nsu, hard=nhu))

    # ---------------- rev20val ----------------
    exp_v = expected("rev20val"); Dv = scan(root, "rev20val"); Fv = {}; err_v = []; miss_v = []; dvpre = {g: ic_dv_pre(g) for g in ("ODY", "G1")}
    for k in exp_v:
        d = Dv.get((arm_key(k[0]),) + k[1:]); s = S20v.get(d.name) if d else None
        if s is not None and s.get("status") != "ok": err_v.append((k, s.get("status"), s.get("error")))
        elif finished(d): Fv[k] = flight(d, s, dv_pre_fallback=dvpre[k[1]] if k[0] == "OLRT3-W" else 0.0)
        else: miss_v.append(k)
    # t* per geometry from the guard arm's own logging (replan.csv ttg_tt_s); the unguarded OLRT3 uses the same value at that geometry (the guard's t* depends on geometry and the pinned IC only, planner_ttg_setup)
    tstar = {}
    for g in ("ODY", "G1"):
        v = [f["tt_aim_d"] for k, f in Fv.items() if k[0] == "OLRT3-W" and k[1] == g and np.isfinite(f["tt_aim_d"])]
        if v: tstar[g] = dict(aim_d=float(np.median(v)), nominal_d=float(np.median(v)) - OFFSET_D.get(GID[g], 0.0), spread_d=float(max(v) - min(v)), n=len(v), offset_d=OFFSET_D.get(GID[g], 0.0))
    def win(f, g, ts=None):
        tstar_ = tstar if ts is None else {g: ts}
        if g not in tstar_ or not f["reached"]: return dict(dt_d=float("nan"), dt_over_w=float("nan"), window=None, window_aim=None, window_wait=None, k_wait=-1, k_max=0)
        dt = f["t_arr_d"] - tstar_[g]["nominal_d"]; dta = f["t_arr_d"] - tstar_[g]["aim_d"]
        w = D20.wait_rule(f["err"], f["t_arr_d"], f["last_drop_km"], f["T_final_s"] / 86400.0, tstar_[g]["nominal_d"], W_MAX_PERIODS)      # Revision 20d: wait rule (extra final orbits within the 25 km tolerance)
        return dict(dt_d=dt, dt_over_w=dt / f["wmax_d"], window=bool(abs(dt) <= f["wmax_d"]), window_aim=bool(abs(dta) <= f["wmax_d"]), window_wait=bool(w["met_wait"]), k_wait=int(w["k_used"]) if w["met_wait"] else -1, k_max=int(w["k_max"]))
    for k, f in list(Fv.items()) + [(k, f) for k, f in F17.items() if k[0] == "OLRT3"]:
        f.update(win(f, k[1])); f["mission"] = bool(f["ok"] and f["window"]) if f["window"] is not None else None; f["mission_wait"] = bool(f["ok"] and f["window_wait"]) if f["window_wait"] is not None else None
        f["apse_timing_deg"] = f["apse_rate_deg_per_d"] * f["dt_d"] if f["window"] is not None else float("nan")
    def add_pos(F_, dirs_, keyf):
        """Position fields (planner_rev21_position) of each flight with a directory; t* = the logged nominal of the geometry (the same t* as the window), else the module's nominal. Failure leaves the fields at NaN."""
        for k_, f_ in F_.items():
            f_.update({q: float("nan") for q in POSKEYS}); d_ = dirs_.get(keyf(k_))
            if d_ is None or not f_["reached"]: continue
            try:
                tk_ = PS.GEOM_ID.get(GID.get(k_[1], k_[1]), k_[1]); ts_ = f_.get("tstar_nom_d"); ts_ = tstar.get(k_[1], {}).get("nominal_d") if ts_ is None or not np.isfinite(ts_) else ts_
                xx_ = f_.get("ic_x0_nominal"); r_ = PS.flight_position(d_, tk_, ts_, PS.nominal(tk_, xx_) if xx_ is not None and np.isfinite(xx_) else None, ra_target=f_.get("target_m"), stop_extra_m=f_.get("stop_extra_m", 0.0))
                if r_: f_.update({q: r_[q] for q in POSKEYS})
            except Exception as ex_:
                f_["pos_error"] = f"{type(ex_).__name__}: {ex_}"
    add_pos(Fv, Dv, lambda k_: (arm_key(k_[0]),) + k_[1:]); add_pos({k_: f_ for k_, f_ in F17.items() if k_[0] == "OLRT3"}, D17, lambda k_: k_)
    L += ["## Part B. rev20val: OLRT3-W (selected IC S-fixed + guard v2) against OLRT3", "",
          f"Finished flights: {len(Fv)} of {len(exp_v)} (MEASURED). Missing or running: {mlist(miss_v)}. Run errors: {'; '.join(f'{cellname(k)} s{k[3]}: {s} {e}' for k, s, e in err_o + err_v if k in exp_v) or 'none'}.", ""]
    if tstar:
        L += ["t* (MEASURED from the guard flights' replan.csv column ttg_tt_s, the aim the guard flew; DERIVED nominal = aim minus the configured aim offset, which is the ic_plan nominal arrival of AC4):", ""] + table(["geometry", "OLRT3-W flights used", "logged aim t* (d)", "aim offset (d, ASSUMED)", "nominal t* (d)", "spread across flights (d)"], [[g, v["n"], f2(v["aim_d"], 4), f2(v["offset_d"], 2), f2(v["nominal_d"], 4), f2(v["spread_d"], 6)] for g, v in tstar.items()]) + [""]
    L += ["Decisions: (1) the window is scored against the NOMINAL t* (AC4 text: t* = the ic_plan nominal arrival); the count against the aimed t* (nominal plus the offset) is given beside it. (2) The unguarded OLRT3 has no logged t*; it uses the OLRT3-W value at the same geometry, which is exact because the guard's t* is computed in planner_ttg_setup from the geometry and the pinned IC alone (no seed or level dependence; the spread column checks this). If no OLRT3-W flight of a geometry is finished, its windows are n/a. "
          "(3) t_arr = flown time at the closing pass + half the Kepler period of the closing orbit; W_max = 3 closing-orbit periods of the flight's own final orbit (as planner_rev20_replay). (4) Pre-burn dV: summary ic_dv_pre_mps when the summary exists, otherwise the pinned ic_plan value (S-fixed, w_pre |nu_pre|) for OLRT3-W and 0 for OLRT3; it is not debited from mass, so total = dV(last counted pass) + pre-burn (Revision 20d: no final periapsis maneuver; FB(h_des 250 km) is in a side table only). "
          "(6) Revision 20d: the window is also scored with the wait rule (planner_rev20d.wait_rule: further final orbits k while the apoapsis error minus k times the last flown drop stays >= -25 km, ASSUMED), labeled 'wait'; the first-passage column is the registered one and P20-9 uses it. AC5: last-pass mean q in [0.05, 0.25] W/cm2 (SOURCED) and last-pass Q_s <= 0.65 are reported, outside the corridor is not a failure. "
          "(5) Apse position: the target apse position is defined in no file, so the table reports the argument of periapsis at the closing pass and the DERIVED apse rotation implied by the timing error alone (J2 rate of the final orbit times t_arr - t*); neither is controlled.", ""]
    keys = [k for k in exp_v if k in Fv]; rows = []
    for k in keys:
        f = Fv[k]; rows.append([k[0], k[1], k[2], k[3], "OK" if f["ok"] else "FAIL", "HARD" if f["hard"] else "", f2(f["err"], 1), f2(f["t_arr_d"], 3), f2(tstar.get(k[1], {}).get("nominal_d"), 3), f2(f["dt_over_w"], 2), {True: "yes", False: "no", None: "n/a"}[f["window"]],
                                {True: "yes", False: "no", None: "n/a"}[f["mission"]], {True: "yes", False: "no", None: "n/a"}[f["window_wait"]], f["k_wait"] if f["window_wait"] else "-", f2(f["dv_sum"]), f2(f["dv_pre"], 3), f2(f["total"]), f2(f["q_last"], 3), "yes" if f["corridor"] else "NO", f2(f["qs_last"], 3), f2(f["alt_f_km"], 1), f["burns"], f["reversals"], f2(f["omega_f_deg"], 1), f2(f["apse_timing_deg"], 2), f["stop"] + ("" if f["stop_src"] == "summary" else "*")] + posrow(f))
    L += ["### Per flight (OLRT3-W and OLRT3 flights of rev20val)", ""] + table(["arm", "geom", "level", "seed", "success", "hard", "final error (km)", "t_arr (d)", "t* (d)", "dt / W_max", "window met", "mission (success and window)", "window (wait, Rev 20d)", "extra orbits k", "campaign dV (m/s)", "pre-burn (m/s)", "aero dV + pre-burn (m/s)", "last-pass q (W/cm2)", "q in corridor", "last-pass Q_s", "final periapsis alt (km)", "burns", "reversals", "omega_f (deg)", "apse shift from dt (deg)", "stop"] + POSHEAD, rows) + ["", "A star on the stop reason marks a stop inferred from the data (no summary row yet)."] + [""]
    # cells
    cells = {}
    for k in exp_v:
        cells.setdefault(k[:3], []).append(k)
    def cmp_flights(cell):
        arm, g, lv = cell
        if arm != "OLRT3-W": return None
        return ("OLRT3 Rev 17" if g == "ODY" else "OLRT3 new", [(F17.get(("OLRT3", g, lv, s)) if g == "ODY" else Fv.get(("OLRT3", g, lv, s))) for s in SEEDS])
    cres = {}; rows = []; Lc = []
    for cell, ks in cells.items():
        fs = [Fv[k] for k in ks if k in Fv]; n_exp = len(ks); n = len(fs); c = dict(expected=n_exp, scored=n, success=sum(f["ok"] for f in fs), hard=sum(f["hard"] for f in fs), apo_ok=sum(f["apo_ok"] for f in fs), window=sum(1 for f in fs if f["window"]), window_aim=sum(1 for f in fs if f["window_aim"]),
                                                       window_wait=sum(1 for f in fs if f["window_wait"]), mission=sum(1 for f in fs if f["mission"]), mission_wait=sum(1 for f in fs if f["mission_wait"]), med_total_succ=med([f["total"] for f in fs if f["ok"]]), med_total_fb_succ=med([f["total_fb"] for f in fs if f["ok"]]),
                                                       med_campaign=med([f["dv_sum"] for f in fs]), corridor=sum(1 for f in fs if f["reached"] and f["corridor"]), q_below=sum(1 for f in fs if f["reached"] and f["q_last"] < 0.05), q_above=sum(1 for f in fs if f["reached"] and f["q_last"] > 0.25), qs_over=sum(1 for f in fs if f["reached"] and not f["qs_last_ok"]),
                                                       med_alt=med([f["alt_f_km"] for f in fs if f["reached"]]), reached=sum(1 for f in fs if f["reached"]))
        cres[cell] = c; rows.append([cellname(cell).replace("@", " "), f"{n}/{n_exp}", c["success"], c["hard"], c["window"], c["window_wait"], c["window_aim"], c["mission"], c["mission_wait"], f2(c["med_total_succ"]), f2(c["med_campaign"])] + posstat_row(fs))
    # Rev 17 OLRT3 comparators at Odyssey (cells not in rev20val except L4 s14)
    for lv in LEVELS:
        fs = [olrt3(lv, s) for s in SEEDS]; fs = [f for f in fs if f]
        if fs: rows.append([f"OLRT3 Rev 17 ODY {lv} (comparator)", f"{len(fs)}/4", sum(f["ok"] for f in fs), sum(f["hard"] for f in fs), sum(1 for f in fs if f["window"]), sum(1 for f in fs if f["window_wait"]), sum(1 for f in fs if f["window_aim"]), sum(1 for f in fs if f["mission"]), sum(1 for f in fs if f["mission_wait"]), f2(med([f["total"] for f in fs if f["ok"]])), f2(med([f["dv_sum"] for f in fs]))] + posstat_row(fs))
    L += ["### Per cell (counts out of the flights finished; the Rev 17 OLRT3 comparator rows use t* of the Odyssey OLRT3-W flights and a pre-burn of 0)", ""] + table(["cell", "finished", "success", "hard", "window met (first passage, registered)", "window met (wait rule, Rev 20d)", "window met (vs aimed t*, first passage)", "mission (success and window, first passage)", "mission (wait rule)", "median aero dV + pre-burn of successes (m/s)", "median campaign dV (m/s)"] + POSSTATHEAD, rows) + [""]
    # side table: FB (labeled side analysis, Revision 20d) and AC5 corridor
    rows = []; rows5 = []
    for cell, c in cres.items():
        rows.append([cellname(cell).replace("@", " "), f2(c["med_total_succ"]), f2(c["med_total_fb_succ"]), f2(c["med_total_fb_succ"] - c["med_total_succ"]) if np.isfinite(c["med_total_fb_succ"]) and np.isfinite(c["med_total_succ"]) else "n/a"])
        rows5.append([cellname(cell).replace("@", " "), c["reached"], c["corridor"], c["q_below"], c["q_above"], c["qs_over"], f2(c["med_alt"], 1)])
    L += ["### Side table (Revision 20d): first burn FB(h_des = 250 km, ASSUMED) is out of scope for the main results; shown only as a side analysis", ""] + table(["cell", "median aero dV + pre-burn of successes (m/s)", "same + FB 250 km (m/s)", "FB share (m/s)"], rows) + [""]
    L += ["### AC5 (Revision 20d): final orbit corridor per cell: last-pass mean q in [0.05, 0.25] W/cm2 (SOURCED), last-pass Q_s <= 0.65 (SOURCED), final periapsis altitude (MEASURED). Outside the corridor is reported, not a failure", ""] + table(["cell", "reached", "q in corridor", "q below 0.05", "q above 0.25", "last-pass Q_s > 0.65", "median final periapsis alt (km)"], rows5) + [""]
    # AC2/3/4
    rows = []; ac = {}
    for cell, c in cres.items():
        if cell[0] != "OLRT3-W": continue
        cname, cf = cmp_flights(cell); cfin = [f for f in cf if f]; cp_ok = len(cfin) == 4; n = c["scored"]; complete = n == c["expected"]
        # AC2
        if c["hard"] > 0 and cell[2] == "L0": v2 = "FAIL"
        elif cp_ok and c["hard"] > sum(f["hard"] for f in cfin): v2 = "FAIL"
        elif complete and cp_ok: v2 = "PASS"
        else: v2 = "PENDING"
        # AC3: apoapsis success not below the comparator by more than n/16 flights
        cs = sum(f["apo_ok"] for f in cfin); v3 = ("PASS" if c["apo_ok"] >= cs - c["expected"] / 16 else "FAIL") if complete and cp_ok else "PENDING"
        # AC4: window count threshold scaled from 14/16 (L0) and 12/16 (L1-L4)
        need = math.ceil((14 if cell[2] == "L0" else 12) / 16 * c["expected"]); v4 = "PASS" if c["window"] >= need else ("FAIL" if c["window"] + (c["expected"] - n) < need else "PENDING")
        needw = c["window_wait"] >= need; v4w = "PASS" if needw else ("FAIL" if c["window_wait"] + (c["expected"] - n) < need else "PENDING")
        ac[cellname(cell)] = dict(AC4_wait=v4w, window_wait=c["window_wait"], AC2=v2, AC3=v3, AC4=v4, need_window=need, comparator=cname, comparator_hard=sum(f["hard"] for f in cfin), comparator_apo_ok=cs, comparator_finished=len(cfin))
        rows.append([cellname(cell).replace("@", " "), f"{c['hard']} vs {sum(f['hard'] for f in cfin)}", v2, f"{c['apo_ok']} vs {cs}", v3, f"{c['window']} (need {need}); mission {c['mission']}", v4, f"{c['window_wait']} (need {need}); mission {c['mission_wait']}", v4w, f"{cname} ({len(cfin)}/4)"])
    L += ["### AC2 / AC3 / AC4 per cell (OLRT3-W against its paired comparator)", "",
          "ASSUMED scaling to 4 flights per cell: AC2 = no more hard violations than the paired OLRT3 and zero at L0; AC3 = apoapsis success (reached, error in [-25, 0] km, no hard) not below OLRT3's by more than 1/16 of the flights, i.e. not below at all for 4 flights; "
          "AC4 = window met in at least ceil(14/16 n) flights at L0 and ceil(12/16 n) at L1-L4 (4 of 4 and 3 of 4); AC4 mission success (success and window) is reported beside the window count. PENDING = cell or comparator not complete and not yet decided.", ""] + table(["cell", "hard W vs comparator", "AC2", "apoapsis success W vs comparator", "AC3", "window met (first passage)", "AC4", "window met (wait rule, Rev 20d)", "AC4 (wait, Rev 20d)", "comparator"], rows) + [""]
    # pairing
    rows = []; pr_ = []
    for k in keys:
        if k[0] != "OLRT3-W": continue
        c_ = F17.get(("OLRT3", k[1], k[2], k[3])) if k[1] == "ODY" else Fv.get(("OLRT3", k[1], k[2], k[3])); f = Fv[k]
        if c_ is None: rows.append([k[1], k[2], k[3], "comparator missing"] + [""] * 6); continue
        d = lambda x, y: f2(x - y) if np.isfinite(x) and np.isfinite(y) else "n/a"
        rows.append([k[1], k[2], k[3], "OK" if f["ok"] else "FAIL" + (" HARD" if f["hard"] else ""), "OK" if c_["ok"] else "FAIL" + (" HARD" if c_["hard"] else ""), d(f["dv_sum"], c_["dv_sum"]), d(f["total"], c_["total"]), f["burns"] - c_["burns"], f2(f["dt_over_w"], 2), f2(c_["dt_over_w"], 2)])
        pr_.append(dict(geom=k[1], level=k[2], seed=k[3], d_campaign=f["dv_sum"] - c_["dv_sum"], d_total=f["total"] - c_["total"], d_burns=f["burns"] - c_["burns"]))
    L += ["### Pairing by seed and level: OLRT3-W minus OLRT3 (Odyssey: Rev 17 flights; G1 L0: new flights). Total dV of OLRT3-W is aerobraking dV plus the pre-burn (no FB, Revision 20d); the comparator has none", ""] + table(["geom", "level", "seed", "OLRT3-W", "OLRT3", "d campaign dV (m/s)", "d aero dV + pre-burn (m/s)", "d burns", "dt/W OLRT3-W", "dt/W OLRT3"], rows) + [""]
    # P20-9
    wc = {lv: cres.get(("OLRT3-W", "ODY", lv)) for lv in LEVELS}; g1 = cres.get(("OLRT3-W", "G1", "L0")); allfin = all(c and c["scored"] == c["expected"] for c in list(wc.values()) + [g1])
    h = lambda lv: wc[lv]["hard"] if wc[lv] else 0
    fals = any(wc[lv] and wc[lv]["window"] <= 2 and wc[lv]["scored"] == wc[lv]["expected"] for lv in LEVELS) or h("L0") > 0
    held = all(wc[lv] and wc[lv]["window"] >= 3 for lv in LEVELS) and g1 and g1["window"] == 4 and h("L0") + h("L1") + h("L3") == 0 and h("L2") + h("L4") <= 1
    p9 = "FALSIFIED" if fals else ("HELD" if held and allfin else ("PENDING" if not allfin else "NOT MET, NOT FALSIFIED"))
    L += ["### P20-9 (disposition exactly as registered): OLRT3-W meets the window in at least 3 of 4 flights per Odyssey level and 4 of 4 at G1 L0, with no hard violation at L0, L1, L3 and at most 1 at L2 and L4 together; falsified if any Odyssey level has 2 or fewer or a hard violation occurs at L0.", "",
          f"**{p9}**: window met by Odyssey level " + ", ".join(f"{lv} {wc[lv]['window']}/{wc[lv]['scored']}" for lv in LEVELS if wc[lv]) + f"; G1 L0 {g1['window'] if g1 else 'n/a'}/{g1['scored'] if g1 else 0}; hard L0 {h('L0')}, L1 {h('L1')}, L3 {h('L3')}, L2+L4 {h('L2') + h('L4')} (MEASURED)." + ("" if allfin else " Block not complete."), "",
          "Revision 20d (labeled; P20-9 above is the registered first-passage count and is unchanged): window met with the wait rule by Odyssey level " + ", ".join(f"{lv} {wc[lv]['window_wait']}/{wc[lv]['scored']}" for lv in LEVELS if wc[lv]) + f"; G1 L0 {g1['window_wait'] if g1 else 'n/a'}/{g1['scored'] if g1 else 0} (MEASURED/DERIVED, ASSUMED 25 km tolerance and last-drop extension).", ""]
    # AC7
    k7 = ("OLRT3", "ODY", "L4", 14); d7 = Dv.get(k7); r7 = D17.get(k7); ac7 = dict(status="NOT RUN (flight missing)")
    if finished(d7) and finished(r7):
        c = compare_passes(d7, r7); f7 = Fv.get(k7) or flight(d7, S20v.get(d7.name)); fr = F17[k7]; same = f7["stop"] == fr["stop"]
        ac7 = dict(status="PASS" if c["equal"] and same else "FAIL", cmp=c, stop_new=f7["stop"], stop_rev17=fr["stop"], same_stop=same)
        L += ["### AC7: regression flight OLRT3 @ ODY @ L4 s14 (all Revision 20 switches off) against the Rev 17 flight", "",
              f"**{ac7['status']}**: shared passes.csv columns {c['ncols']}, rows new {c['rows_a']} / Rev 17 {c['rows_b']}, all equal as text (bit for bit): {c['equal']}; first difference (row index from 0, column, new, Rev 17): {c['first'] or 'none'}; columns only in new: {c['only_a'] or 'none'}; only in Rev 17: {c['only_b'] or 'none'}; stop new {f7['stop']} ({f7['stop_src']}), Rev 17 {fr['stop']} ({fr['stop_src']}), same stop: {same} (MEASURED).", ""]
    else: L += ["### AC7: regression flight OLRT3 @ ODY @ L4 s14 against the Rev 17 flight", "", "**NOT RUN**: the rev20val regression flight is not finished.", ""]
    # ---------------- Part C: rendezvous position (Revision 21, WP21-B) and any extra block directory
    XB = {}; XI = {}; XD = {}
    for blk in a.blocks:
        bdir = root / blk; allD = sorted(q for q in bdir.glob("*") if q.is_dir()) if bdir.exists() else []; Db = scan(root, blk); Sb = load_summary(root, f"summary_{blk}_*.csv"); Fb = {}
        unparsed = [q.name for q in allD if not parse_id(q.name)]; unfin = [d.name for d in Db.values() if not finished(d)]; errs = [(r["run_id"], r.get("status"), r.get("error")) for r in Sb.values() if r.get("status") != "ok" and r["run_id"] in {d.name for d in Db.values()}]
        raw = {}
        for k, d in Db.items():
            if not finished(d): continue
            os_ = blk in a.one_sided; tk0 = PS.GEOM_ID.get(GID.get(k[1], k[1]), k[1]); f = flight(d, Sb.get(d.name), dv_pre_fallback=0.0, target=PS.target_radius(tk0), onesided=os_); f["target_m"] = PS.target_radius(tk0); f["stop_extra_m"] = TOL_KM * 1e3 if os_ else 0.0; sr = Sb.get(d.name) or {}; f["ic_scenario"] = str(sr.get("ic_scenario", "n/a")); icx = sr.get("ic_x0", float("nan")); f["ic_x0"] = float(icx) if icx is not None and np.isfinite(icx) else float("nan")
            off = OFFSET_D.get(GID.get(k[1], k[1]), 0.0); f["tstar_nom_d"] = f["tt_aim_d"] - off if np.isfinite(f["tt_aim_d"]) else float("nan"); f["ic_x0_nominal"] = f["ic_x0"] if f["ic_scenario"] == "S-free" else float("nan"); raw[k] = f
        # geometry t*: the guard's own logged nominal (S-fixed or reference IC flights), else the module's nominal (guard recipe for the table version), else the earlier global value (ODY, G1 from rev20val)
        gts = {}
        for tok in sorted({k[1] for k in raw}):
            lg = [f["tstar_nom_d"] for k, f in raw.items() if k[1] == tok and f["ic_scenario"] != "S-free" and np.isfinite(f["tstar_nom_d"])]; off = OFFSET_D.get(GID.get(tok, tok), 0.0)
            try: nd = float(PS.nominal(PS.GEOM_ID.get(GID.get(tok, tok), tok))["days"])
            except Exception: nd = float("nan")
            if lg: gts[tok] = dict(aim_d=float(np.median(lg)) + off, nominal_d=float(np.median(lg)), source="logged ttg (replan.csv)", module_nominal_d=nd, n=len(lg))
            elif np.isfinite(nd): gts[tok] = dict(aim_d=nd + off, nominal_d=nd, source="planner_rev21_position.nominal (guard recipe, DERIVED)", module_nominal_d=nd, n=0)
            elif tok in tstar: gts[tok] = dict(tstar[tok], source="rev20val logged", module_nominal_d=nd, n=tstar[tok].get("n", 0))
        for k, f in raw.items():
            ts = gts.get(k[1]); own = f["tstar_nom_d"]
            if np.isfinite(own): ts = dict(aim_d=f["tt_aim_d"], nominal_d=own)
            f.update(win(f, k[1], ts) if ts else win(f, "_none")); f["mission"] = bool(f["ok"] and f["window"]) if f["window"] is not None else None; f["mission_wait"] = bool(f["ok"] and f["window_wait"]) if f["window_wait"] is not None else None
            f["apse_timing_deg"] = f["apse_rate_deg_per_d"] * f["dt_d"] if f["window"] is not None else float("nan")
            if f["stop_extra_m"] > 0: f["window_wait"] = None; f["k_wait"] = -1; f["mission_wait"] = None      # the wait rule is defined for the registered two-sided tolerance only
            if ts: f["tstar_nom_d"] = ts["nominal_d"]
            Fb[k] = f
        add_pos(Fb, Db, lambda k_: k_); XB[blk] = Fb; XD[blk] = Db; XI[blk] = dict(onesided=blk in a.one_sided, targets={t_: PS.target_radius(PS.GEOM_ID.get(GID.get(t_, t_), t_)) for t_ in sorted({k_[1] for k_ in Fb})}, dirs=len(allD), parsed=len(Db), finished=len(Fb), unparsed=unparsed, unfinished=unfin, errors=errs, geometry_tstar=gts, reached=sum(1 for f in Fb.values() if f["reached"]), missing_dir=not bdir.exists())
    L += ["## Part C. Rendezvous position error (Revision 21, WP21-B)", "",
          "Method (planner_rev21_position.py): target = apoapsis direction at t* of the model-world nominal plan (initial elements of the geometry, SOURCED, advanced by the secular J2 rates along the nominal plan's (a, e) history, i = i0; t* = the logged guard t* of the geometry); flight = apoapsis direction at the arrival "
          "from the flight's closing-pass osculating (omega, Omega, i) advanced half a J2 orbit. Error = great-circle angle (DERIVED), arc = angle x 7,408 km. Decomposition: (a) timing = J2 nominal at t_arr minus at t*; (b) predictor residual = actual minus the J2 prediction from the flight's own (a, e, i) history; "
          "(c) history = that prediction minus the nominal at t_arr; a + b + c = total to second order in the angles (closure residual in the JSON). tol_pos_deg = 1.0 deg is ASSUMED (0.5 and 2.0 shown); the position is reported, not controlled (AC9 is a reporting criterion).", ""]
    def pos_blocks():
        yield "rev20val", {k: f for k, f in Fv.items()}
        yield "rev20olt", {k: f for k, f in Fo.items()}
        yield "rev17 (all arms)", {k: f for k, f in F17.items()}
        for b_, Fb_ in XB.items(): yield b_, Fb_
    for name, Fx in pos_blocks():
        if name in ("rev20olt", "rev17 (all arms)"):
            add_pos(Fx, Do if name == "rev20olt" else D17, (lambda k_: (arm_key(k_[0]),) + k_[1:]) if name == "rev20olt" else (lambda k_: k_))
        cl = {}
        for k, f in Fx.items(): cl.setdefault(cellname(k[:3]), []).append(f)
        rows = [[c] + posstat_row(fs) + [f2(med([f["a_timing_deg"] for f in fs if np.isfinite(f.get("a_timing_deg", float("nan")))]), 2), f2(med([f["b_pred_deg"] for f in fs if np.isfinite(f.get("b_pred_deg", float("nan")))]), 2),
                                           f2(med([f["c_hist_deg"] for f in fs if np.isfinite(f.get("c_hist_deg", float("nan")))]), 2), f2(max([f["closure_deg"] for f in fs if np.isfinite(f.get("closure_deg", float("nan")))] or [float("nan")]), 4)] for c, fs in cl.items()]
        L += [f"### Position per cell: {name}", ""] + table(["cell"] + POSSTATHEAD + ["median (a) timing", "median (b) predictor", "median (c) history", "max closure residual (deg)"], rows) + [""]
        if name not in ("rev20val", "rev17 (all arms)", "rev20olt"):
            rows = [[k[0], k[1], k[2], k[3], "OK" if f["ok"] else "FAIL" + (" HARD" if f["hard"] else ""), f2(f["err"], 1), f2(f["t_arr_d"], 3), f2(f["dt_over_w"], 2), {True: "yes", False: "no", None: "n/a"}[f["window"]], {True: "yes", False: "no", None: "n/a"}[f["window_wait"]], f2(f["dv_sum"]), f2(f["dv_pre"], 3)] + posrow(f)
                    for k, f in sorted(Fx.items())]
            L += [f"### Per flight: {name} (success rule, window against the nominal t* and the wait rule as in Part B; dV = aerobraking dV to the last counted pass, pre-burn from the summary when present)", ""] + table(["arm", "geom", "level", "seed", "success", "final error (km)", "t_arr (d)", "dt / W_max", "window met", "window (wait)", "campaign dV (m/s)", "pre-burn (m/s)"] + POSHEAD, rows) + [""]
    # ---------------- Part D: extra blocks, integrity, per-cell summary and paired comparisons
    if XB:
        L += ["## Part D. Extra blocks: integrity, per-cell summary and paired comparisons (Revision 21)", "",
              "Integrity: every directory of the block is accounted for. unparsed = name does not match <gid>_PL_<arm>_pert_<field>_s<seed>[_<level>] with a known geometry prefix; unfinished = one of passes/decisions/measured/replan.csv is missing (running or failed); run errors = summary status not ok. Nothing is dropped silently.", ""]
        for blk, I in XI.items():
            L += [f"### Integrity: {blk}", "", f"directories {I['dirs']}{' (BLOCK DIRECTORY MISSING)' if I['missing_dir'] else ''}; target radius per geometry (m): {', '.join(f'{g} {r:.0f}' for g, r in I['targets'].items()) or 'n/a'}; tolerance rule: {'one-sided, r_a in [target, target + 25 km] (ASSUMED sign), wait rule not scored' if I['onesided'] else 'registered, r_a in [target - 25 km, target]'}; parsed {I['parsed']}; finished and scored {I['finished']} (reached {I['reached']}); unparsed {len(I['unparsed'])}: {', '.join(I['unparsed']) or 'none'}; unfinished {len(I['unfinished'])}: {', '.join(I['unfinished']) or 'none'}; run errors: {'; '.join(f'{r} {st} {er}' for r, st, er in I['errors']) or 'none'}.", ""]
            if I["geometry_tstar"]: L += table(["geometry", "t* nominal used (d)", "source", "module nominal (guard recipe, d)", "guard flights logged"], [[g, f2(v["nominal_d"], 3), v["source"], f2(v.get("module_nominal_d"), 3), v.get("n", 0)] for g, v in I["geometry_tstar"].items()]) + [""]
        def cellsum(Fx):
            cl = {}
            for k, f in Fx.items(): cl.setdefault(cellname(k[:3]), []).append(f)
            rows = []
            for c, fs in cl.items():
                rr = [f for f in fs if f["reached"]]
                rows.append([c.replace("@", " "), len(fs), sum(f["ok"] for f in fs), sum(f["hard"] for f in fs), sum(1 for f in fs if f["window"]), sum(1 for f in fs if f["window_wait"]), f2(med([f["dv_sum"] for f in fs]), 2), f2(med([f["dv_pre"] for f in fs if np.isfinite(f["dv_pre"])]), 3),
                             sum(1 for f in rr if f["corridor"]), sum(1 for f in rr if not f["qs_last_ok"]), f2(med([f["alt_f_km"] for f in rr]), 1)] + posstat_row(fs))
            return rows
        CSH = ["cell", "flights", "success", "hard", "window (first)", "window (wait)", "median aerobraking dV (m/s)", "median pre-burn (m/s)", "AC5 q in corridor", "last-pass Q_s > 0.65", "median final periapsis alt (km)"] + POSSTATHEAD
        for blk, Fx in XB.items(): L += [f"### Per-cell summary: {blk}", ""] + table(CSH, cellsum(Fx)) + [""]
        pool = {}
        for src, FF in (("rev21c", XB.get("rev21c", {})), ("rev21v", XB.get("rev21v", {})), ("rev20val", {(arm_key(k[0]),) + k[1:]: f for k, f in Fv.items()}), ("rev17", {(arm_key(k[0]),) + k[1:]: f for k, f in F17.items()})):
            for k, f in FF.items(): pool.setdefault(k, []).append((src, f))
        for blk, Fx in XB.items():
            if blk not in ("rev21c", "rev21v"):
                for k, f in Fx.items(): pool.setdefault(k, []).append((blk, f))
        def pair_rows(Xa, comp, name_a, name_b):
            rows = []; agg = {}
            for k, f in sorted(Xa.items()):
                cm = comp(k)
                if cm is None: continue
                src, g = cm; dd = f["dv_sum"] - g["dv_sum"] if np.isfinite(f["dv_sum"]) and np.isfinite(g["dv_sum"]) else float("nan")
                rows.append([k[1], k[2], k[3], "OK" if f["ok"] else "FAIL" + (" HARD" if f["hard"] else ""), f"{src}: " + ("OK" if g["ok"] else "FAIL" + (" HARD" if g["hard"] else "")), f2(dd), f["burns"] - g["burns"], f2(f["dt_over_w"], 2), f2(g["dt_over_w"], 2), f2(f.get("err_deg"), 2), f2(g.get("err_deg"), 2)])
                agg.setdefault((k[1], k[2]), []).append((f, g, dd))
            ar = [[gk[0], gk[1], len(v), f"{sum(x[0]['ok'] for x in v)} vs {sum(x[1]['ok'] for x in v)}", f"{sum(x[0]['hard'] for x in v)} vs {sum(x[1]['hard'] for x in v)}", f"{sum(1 for x in v if x[0]['window'])} vs {sum(1 for x in v if x[1]['window'])}", f2(med([x[2] for x in v]), 2), f2(med([x[0].get('err_deg', float('nan')) for x in v]), 2), f2(med([x[1].get('err_deg', float('nan')) for x in v]), 2)] for gk, v in sorted(agg.items())]
            return rows, ar
        PH = ["geom", "level", "seed", "success", "comparator", "d aerobraking dV (m/s)", "d burns", "dt/W_max", "dt/W_max comparator", "pos err (deg)", "pos err comparator (deg)"]; AH = ["geom", "level", "pairs", "success (A vs B)", "hard (A vs B)", "window first (A vs B)", "median d dV (A minus B, m/s)", "median pos err A (deg)", "median pos err B (deg)"]
        allx = {k: f for FF in XB.values() for k, f in FF.items()}
        w3 = {k: f for k, f in allx.items() if k[0] == "OLRT3W3"}; wf3 = {k: f for k, f in allx.items() if k[0] == "OLRT3WF3"}
        def comp_olrt3(k):
            for src, g in pool.get(("OLRT3",) + k[1:], []):
                if g is not w3.get(k): return src, g
            return None
        def comp_w3(k):
            g = w3.get(("OLRT3W3",) + k[1:]); return ("OLRT3-W3 (S-fixed)", g) if g is not None else None
        for title, Xa, comp, A_, B_ in (("OLRT3-W3 (guard v3) against OLRT3 by geometry, level and seed (A = W3, B = OLRT3 from rev21c, then rev21v, rev20val, rev17)", w3, comp_olrt3, "W3", "OLRT3"),
                                        ("OLRT3-WF3 (S-free) against OLRT3-W3 (S-fixed) (A = WF3, B = W3)", wf3, comp_w3, "WF3", "W3")):
            rows, ar = pair_rows(Xa, comp, A_, B_)
            L += [f"### Paired: {title}", ""] + (table(AH, ar) + ["", "Per pair:", ""] + table(PH, rows) if rows else ["No pairs yet (no flights of arm A with a comparator in the finished data)."]) + [""]
        bl = {}
        for k, f in allx.items():
            if (k[1], k[2]) in (("G0975", "L2p"), ("G0063", "L0")): bl.setdefault((k[1], k[2], k[0]), []).append(f)
        rows = []
        for (g, lv, arm), fs in sorted(bl.items()):
            rr = [f for f in fs if f["reached"]]; rows.append([g, lv, arm, len(fs), sum(f["ok"] for f in fs), sum(f["hard"] for f in fs), sum(1 for f in fs if f["window"]), sum(1 for f in fs if f["window_wait"]), f2(med([f["dv_sum"] for f in fs]), 2), f2(med([f["burns"] for f in fs]), 0), f2(med([f.get("err_deg", float("nan")) for f in fs]), 2)])
        L += ["### Baselines (BODYT, N2IT) against planner arms at G0975 L2p and G0063 L0 (all finished flights of the given blocks)", ""] + (table(["geom", "level", "arm", "flights", "success", "hard", "window (first)", "window (wait)", "median aerobraking dV (m/s)", "median burns", "median pos err (deg)"], rows) if rows else ["No flights at these cells in the given blocks."]) + [""]
    # ---------------- Part E: circular target (Revision 22 item 1): R1 = circle at the target radius, R3 = circle at the mean radius of the nominal final ellipse, both in the nominal plane at t*
    CIRC = {}
    if a.blocks:
        import planner_rev22_rendezvous as R22
        def circ(d_, tok_, rt_, extra_):
            """planner_rev22_rendezvous.score_flight (read only) with its target radius and closing-pass radius set for this flight (module globals restored afterwards)."""
            def fe(dd):
                P_ = pd.read_csv(dd / "passes.csv").sort_values("pass").reset_index(drop=True); ra_ = P_.ra_m.values.astype(float); h_ = np.where(np.isfinite(ra_) & (ra_ <= rt_ + extra_))[0]
                if not len(h_): return None
                n_ = int(h_[0]); row = P_.iloc[n_]; rp_ = float(row.rp_m); rar = float(row.ra_m); a_, e_ = R22.ell(rp_, rar); a0, e0, i0 = float(row.a_m), float(row.e), float(row.i_rad); dom, dOm = PS.j2_orbit(a0, e0, i0)
                return dict(a=a_, e=e_, i=i0, Om=float(row.raan_rad) + 0.5 * dOm, om=float(row.omega_rad) + 0.5 * dom, rp=rp_, ra=rar, t_arr=float(row.time_s) + 0.5 * R22.period(a0), n_pass=n_ + 1)
            old = (R22.RA_T, R22.final_ellipse); R22.RA_T = rt_; R22.final_ellipse = fe
            try: return R22.score_flight(d_, tok_)
            finally: R22.RA_T, R22.final_ellipse = old
        CK = ("R1_miss_best_km", "R1_nu_r_deg", "R1_cost_best_mps", "R1_cost_apo_mps", "R1_miss_apo_km", "R1_timing_part_km", "R3_miss_best_km", "R3_cost_best_mps")
        for name_, Fx_, Dx_ in [("rev20val", Fv, {k_: Dv.get((arm_key(k_[0]),) + k_[1:]) for k_ in Fv})] + [(b_, Fb_, XD[b_]) for b_, Fb_ in XB.items()]:
            CIRC[name_] = {}
            for k_, f_ in Fx_.items():
                d_ = Dx_.get(k_) if name_ == "rev20val" else Dx_.get(k_)
                if d_ is None or not f_["reached"]: continue
                tk_ = PS.GEOM_ID.get(GID.get(k_[1], k_[1]), k_[1])
                try:
                    r_ = circ(d_, tk_, PS.target_radius(tk_), f_.get("stop_extra_m", 0.0))
                    if r_: CIRC[name_][k_] = {q: r_.get(q) for q in CK}
                except Exception as ex_: CIRC[name_][k_] = dict(error=f"{type(ex_).__name__}: {ex_}")
        L += ["## Part E. Circular target (Revision 22 item 1)", "",
              "From planner_rev22_rendezvous.score_flight (read only): the target is the circle of radius r_t (the target apoapsis radius of the geometry: 7,408,000 m by default, SOURCED; target_orbits.json for the low-altitude tokens) in the NOMINAL plan's plane at t* (R1), and the circle at the mean radius of the nominal final ellipse (R3, ASSUMED definition of Revision 22a). "
              "The rendezvous point nu_r is on the flown final ellipse after the closing pass; miss = nearest point of the target circle (phase free), best nu_r minimizes it, cost = |v_target - v_final| there (DERIVED, outside the planner's control); timing part = arc between the nominal plane normals at t* and at the flight's arrival. "
              "The closing pass is the first with r_a <= r_t (plus 25 km for one-sided blocks).", ""]
        for name_, dd_ in CIRC.items():
            cl_ = {}
            for k_, v_ in dd_.items():
                if "error" not in v_ and v_.get("R1_miss_best_km") is not None: cl_.setdefault(cellname(k_[:3]), []).append(v_)
            rows_ = [[c_.replace("@", " "), len(vs_), f2(med([v["R1_miss_best_km"] for v in vs_]), 1), f2(float(np.percentile([v["R1_miss_best_km"] for v in vs_], 90)), 1), f2(med([v["R1_nu_r_deg"] for v in vs_]), 1), f2(med([v["R1_cost_best_mps"] for v in vs_]), 1), f2(med([v["R1_cost_apo_mps"] for v in vs_]), 1),
                      f2(med([v["R1_timing_part_km"] for v in vs_]), 1), f2(med([v["R3_miss_best_km"] for v in vs_]), 1), f2(med([v["R3_cost_best_mps"] for v in vs_]), 1)] for c_, vs_ in cl_.items()]
            errs_ = [f"{cellname(k_)} s{k_[3]}: {v_['error']}" for k_, v_ in dd_.items() if "error" in v_]
            L += [f"### Circular target per cell: {name_}", ""] + (table(["cell", "n", "R1 median miss at best nu_r (km)", "R1 p90 miss (km)", "R1 median best nu_r (deg from the closing apoapsis = 180)", "R1 median cost at best nu_r (m/s)", "R1 median cost at apoapsis (m/s)", "R1 median timing part (km)", "R3 median miss (km)", "R3 median cost (m/s)"], rows_) if rows_ else ["No scored flights."]) + [""] + ([f"Errors: {'; '.join(errs_)}", ""] if errs_ else [])
            if name_ not in ("rev20val",):
                rows_ = [[k_[0], k_[1], k_[2], k_[3], f2(PS.target_radius(PS.GEOM_ID.get(GID.get(k_[1], k_[1]), k_[1])) / 1e3, 1), f2(v_.get("R1_miss_best_km"), 1), f2(v_.get("R1_nu_r_deg"), 1), f2(v_.get("R1_cost_best_mps"), 1), f2(v_.get("R1_cost_apo_mps"), 1), f2(v_.get("R1_miss_apo_km"), 1), f2(v_.get("R1_timing_part_km"), 1)] for k_, v_ in sorted(dd_.items()) if "error" not in v_]
                L += [f"### Circular target per flight: {name_}", ""] + table(["arm", "geom", "level", "seed", "r_t (km)", "R1 miss at best nu_r (km)", "best nu_r (deg)", "cost at best nu_r (m/s)", "cost at apoapsis (m/s)", "miss at apoapsis (km)", "timing part (km)"], rows_) + [""]
    J["position"] = dict(tol_deg_assumed=PS.TOL_BASE, tols=list(PS.TOL_DEG), cells={name: {c: posstats(fs) for c, fs in {cellname(k[:3]): [f2_ for k2_, f2_ in Fx.items() if cellname(k2_[:3]) == cellname(k[:3])] for k in Fx}.items()} for name, Fx in pos_blocks()}, extra_blocks={b_: len(F_) for b_, F_ in XB.items()}, integrity=XI,
                         flights={name: {cellname(k) + f"@s{k[3]}": {q: f.get(q) for q in POSKEYS} for k, f in Fx.items()} for name, Fx in pos_blocks()})
    J["rev20val"] = dict(finished=len(Fv), expected=len(exp_v), missing=[list(k) for k in miss_v], errors=[[list(k), s, e] for k, s, e in err_v], tstar=tstar, flights={cellname(k) + f"@s{k[3]}": v for k, v in Fv.items()},
                         rev17_olrt3={cellname(k) + f"@s{k[3]}": v for k, v in F17.items() if k[0] == "OLRT3"}, cells={cellname(c): v for c, v in cres.items()}, ac=ac, pairs=pr_, p20_9=dict(disposition=p9), ac7=ac7)
    if CIRC: J["position"]["circular"] = {n_: {cellname(k_) + f"@s{k_[3]}": v_ for k_, v_ in d_.items()} for n_, d_ in CIRC.items()}
    out_md.write_text("\n".join(L) + "\n"); out_js.write_text(json.dumps(nan_none(J), indent=1))
    print(f"rev17 {len(F17)} finished | rev20olt {len(Fo)}/{len(exp_o)} | rev20val {len(Fv)}/{len(exp_v)} | gate {gate['status']} | P20-7 {p7} | P20-9 {p9} | AC7 {ac7['status']} -> {out_md}, {out_js}")


if __name__ == "__main__":
    main()
