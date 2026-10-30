"""Revision 20 WP3: baselines and planner arms (flown, WITHOUT an in-campaign time guard) under the rendezvous window and the end-of-campaign toolbox of planner_rev20_window.py. Python only, no Julia, no flights.

For each flown flight (dv_correction/flights.json, passes.csv) and each t* of a grid: window met without toolbox? cheapest toolbox plan, total dV = corrected campaign dV (Revision 17a) + toolbox + first burn FB(h_des),
apse mismatch. Mission success (Revision 13 rule, ASSUMED tolerance): reached, final apoapsis error in [-25, 0] km, no hard violation (Q_s on V_rel > 0.65 W/cm2 at any pass or a safety_q stop). The window statistics are over
SUCCESSFUL flights (survivor bias: failures are counted per cell and are not scored for dV). Row schema (also the interface for a later planner+guard table):
  rows: DataFrame with columns geom, arm, run_id, phase, success, t_star_d, off_d (t* minus the cell's median arrival of successes), met_none, feasible, tb_cost, fb, dv, total, kind, k, k_ph, pos_err_km, t_c_d, dt_d.
  summarize(rows, by=("geom", "arm", "bin")) -> table; cell_table(rows, flights) -> per cell flight / failure counts.
Usage: python3 -I planner_rev20_baselines.py --out DIR [--selftest] [--dh-table]"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev18_flights as F
import planner_rev20_window as W

import planner_lib as _PLB
TASK = pathlib.Path(_PLB.task_dir())
RA_T = W.A.RA_TARGET
TOK = F.TOK
ERR_LO_KM, ERR_HI_KM = -25.0, 0.0            # ASSUMED (Revision 16b): final apoapsis error tolerance
Q_HARD = 0.65                                # SOURCED: planner_config q_never_exceed_qs_w_cm2 (Revision 13)
NOMINAL_D = {"ODY": 61.6, "G1": 75.2}        # DERIVED (Revision 19a, model-world lambda = 0 plan); GC: none (the generic table has no calibrated kappa; table campaign length given for reference only)
ARMS = {"ODY": ["BODY", "N2I", "OLRT3", "OLR", "OL", "L3"], "G1": ["BODY", "OLRT3", "OLR", "OL", "L3"], "GC": ["BODYT", "N2IT", "OLRT3", "OLRT5", "OLRT35"]}   # as requested; BODY@G1 added for comparison
LABEL = {"BODY": "B-ODY", "N2I": "N2-I", "BODYT": "B-ODY-T", "N2IT": "N2-I-T"}
BASELINES = {"BODY", "N2I", "BODYT", "N2IT"}
STEP_D = 1.0                                  # numerical (ASSUMED): t* grid step in days
PAD_D = 2.0                                   # ASSUMED: grid padding beyond the flown range
BIN_EDGES = [-1e9, -20, -10, -5, -2, -0.5, 0.5, 2, 5, 10, 20, 1e9]
BIN_LABEL = ["<-20", "-20..-10", "-10..-5", "-5..-2", "-2..-0.5", "-0.5..0.5", "0.5..2", "2..5", "5..10", "10..20", ">20"]
H_DES_SENS_KM = [120.0, 200.0, 250.0, 300.0, 400.0]       # ASSUMED sensitivity set (Revision 18)
DH_SENS_KM = [40.0, 70.0, 100.0]                           # ASSUMED sensitivity set for the pop-up height


def load_flights(flights_json=TASK / "dv_correction" / "flights.json"):
    """Read the arms of ARMS from the 366 indexed flights; returns list of dicts with pass arrays and success flags (no toolbox yet)."""
    R = json.load(open(flights_json)); out = []; sums = {}
    for r in R:
        if r.get("status") == "not_corrected": continue
        g = TOK[str(r["geometry"])]
        if r["arm"] not in ARMS.get(g, []): continue
        d = F.run_dir(r["root"], r["phase"], r["run_id"])
        if d is None: continue
        P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True)
        key = (r["root"], r["summary"])
        if key not in sums: sums[key] = pd.read_csv(pathlib.Path(r["root"]) / r["summary"]).set_index("run_id")
        S = sums[key].loc[r["run_id"]]; ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= RA_T))[0]
        reached = len(hit) > 0; i = int(hit[0]) if reached else len(P) - 1
        qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; hard = bool(np.nanmax(qs) > Q_HARD or "safety_q" in str(S.stop_reason))
        a = P.a_m.values.astype(float); Tn = 2 * math.pi * math.sqrt(a[i] ** 3 / W.MU); err = (ra[i] - RA_T) / 1e3
        out.append(dict(geom=g, arm=r["arm"], run_id=r["run_id"], phase=r["phase"], dv=float(r["dv_corrected"]), reached=reached, hard=hard, err_km=err,
                        success=bool(reached and not hard and ERR_LO_KM <= err <= ERR_HI_KM), relaxed=bool(reached and not hard), t_c=float(P.time_s.iloc[i]) + 0.5 * Tn, T_f=Tn,
                        ra=ra[:i], rp=P.rp_m.values[:i].astype(float), inc=P.i_rad.values[:i].astype(float), ra_f=float(ra[i]), rp_f=float(P.rp_m.iloc[i]), alt_f=float(P.alt_m.iloc[i]), inc_f=float(P.i_rad.iloc[i]),
                        omega_f=float(np.unwrap(P.omega_rad.values.astype(float))[i])))
    return out


def prep(f, h_des_m):
    ok = np.isfinite(f["ra"]) & np.isfinite(f["rp"])
    return W.prepare(f["ra"][ok], f["rp"][ok], f["inc"][ok], f["ra_f"], f["rp_f"], f["alt_f"], f["inc_f"], h_des_m)


def tstar_grid(flights, g, key="success"):
    ts = np.array([f["t_c"] / 86400 for f in flights if f["geom"] == g and f[key]])
    base = np.arange(math.floor(ts.min()) - PAD_D, ts.max() + PAD_D + STEP_D, STEP_D)
    extra = [NOMINAL_D[g]] if g in NOMINAL_D else []
    cells = [float(np.median([f["t_c"] / 86400 for f in flights if f["geom"] == g and f["arm"] == a and f["success"]])) for a in ARMS[g] if any(f["geom"] == g and f["arm"] == a and f[key] for f in flights)]
    return np.unique(np.round(np.concatenate([base, extra, cells, [float(np.median(ts))]]), 6))


def evaluate(flights, h_des_m=250e3, dh=None, only_median=False, geoms=None, key="success"):
    """All (successful flight x t*) pairs; only_median -> t* = the cell median arrival only. Returns the row DataFrame."""
    rows = []; geoms = geoms or list(ARMS)
    omega = {}
    for g in geoms:
        S = [f for f in flights if f["geom"] == g and f[key]]
        for a in ARMS[g]:
            c = [f for f in S if f["arm"] == a]
            if c: omega[(g, a)] = (float(np.median([f["t_c"] / 86400 for f in c])), float(np.median([f["omega_f"] for f in c])))
    for g in geoms:
        grid = tstar_grid(flights, g, key)
        for f in flights:
            if f["geom"] != g or not f[key]: continue
            tmed, wmed = omega[(g, f["arm"])]; st = prep(f, h_des_m); fb = float(W.burn(f["ra_f"], f["rp_f"], st["rp_des"])); Wd = W.W_MAX_PERIODS * f["T_f"]
            disp = (f["omega_f"] - wmed + math.pi) % (2 * math.pi) - math.pi
            for ts in ([tmed] if only_median else grid):
                plan = W.cheapest(f["t_c"], ts * 86400.0, Wd, st, dh)
                pos = float("nan") if plan is None else f["ra_f"] * abs(disp + plan["rot"]) / 1e3
                rows.append(dict(geom=g, arm=f["arm"], run_id=f["run_id"], phase=f["phase"], success=f["success"], rule=key, t_star_d=ts, off_d=ts - tmed, met_none=abs(ts * 86400 - f["t_c"]) <= Wd, feasible=plan is not None,
                                 tb_cost=float("nan") if plan is None else plan["cost"], fb=fb, dv=f["dv"], total=float("nan") if plan is None else f["dv"] + plan["cost"] + fb,
                                 kind="fail" if plan is None else plan["kind"], k=-1 if plan is None else plan["k"], k_ph=-1 if plan is None else plan["k_ph"], pos_err_km=pos,
                                 t_c_d=f["t_c"] / 86400, dt_d=ts - f["t_c"] / 86400, W_d=Wd / 86400, h_des_km=h_des_m / 1e3, dh_km=(W.DH_POPUP_M if dh is None else dh) / 1e3))
    return pd.DataFrame(rows)


def summarize(rows, by=("geom", "arm", "bin")):
    """Table of window success and total dV per group. rows: schema in the module docstring (successful flights only; an off_d column is binned when 'bin' is in by).
    Columns: frac_early / frac_late (t* later / earlier than the closing apoapsis by more than W), median toolbox cost and total dV among pairs that needed the toolbox, n_pairs, n_flights, met_none (rate), feasible (rate with toolbox), median toolbox cost (all feasible pairs), median total dV among feasible, median pos error."""
    R = rows.copy()
    if "bin" in by: R["bin"] = pd.cut(R.off_d, BIN_EDGES, labels=BIN_LABEL, right=False)
    R["early"] = R.dt_d > R.W_d; R["late"] = R.dt_d < -R.W_d; R["need_tb"] = R.feasible & ~R.met_none
    R["tb_needed"] = R.tb_cost.where(R.need_tb); R["total_needed"] = R.total.where(R.need_tb)
    g = R.groupby(list(by), observed=True)
    T = g.agg(frac_early=("early", "mean"), frac_late=("late", "mean"), tb_median_needed=("tb_needed", "median"), total_median_needed=("total_needed", "median"), n_pairs=("run_id", "size"), n_flights=("run_id", "nunique"), met_none=("met_none", "mean"), feasible=("feasible", "mean"), tb_median=("tb_cost", "median"),
              total_median=("total", "median"), dv_median=("dv", "median"), fb_median=("fb", "median"), pos_median_km=("pos_err_km", "median")).reset_index()
    return T


def cell_table(flights):
    """Per (geometry, arm): flights, successes, hard violations, not reached, out-of-tolerance, median arrival (successes), median campaign dV (successes)."""
    rows = []
    for g in ARMS:
        for a in ARMS[g]:
            c = [f for f in flights if f["geom"] == g and f["arm"] == a]; s = [f for f in c if f["success"]]
            rows.append(dict(geom=g, arm=a, flights=len(c), success=len(s), hard=sum(f["hard"] for f in c), not_reached=sum(not f["reached"] for f in c),
                             out_tol=sum(f["reached"] and not f["hard"] and not (ERR_LO_KM <= f["err_km"] <= ERR_HI_KM) for f in c),
                             t_med_d=float(np.median([f["t_c"] / 86400 for f in s])) if s else float("nan"), dv_med=float(np.median([f["dv"] for f in s])) if s else float("nan")))
    return pd.DataFrame(rows)


def dh_table(flights_json=TASK / "dv_correction" / "flights.json"):
    """Inputs of the default DH_POPUP_M: per-pass apoapsis drop distribution and local scale height over all indexed flights."""
    R = json.load(open(flights_json)); drops = []; Hs = []
    for r in R:
        d = F.run_dir(r["root"], r["phase"], r["run_id"])
        if d is None: continue
        P = pd.read_csv(d / "passes.csv").sort_values("pass"); x = -np.diff(P.ra_m.values.astype(float)); drops.append(x[np.isfinite(x)]); Hs.append(float(P.h_local_m.median()))
    a = np.concatenate(drops); q = np.percentile(a, [50, 95, 99, 100]); H = float(np.median(Hs))
    return dict(n_passes=int(len(a)), drop_km_p50_p95_p99_max=[float(v / 1e3) for v in q], H_km=H / 1e3, drop_tol_m=W.DROP_TOL_M, dh_p99_km=H * math.log(q[2] / W.DROP_TOL_M) / 1e3, dh_max_km=H * math.log(q[3] / W.DROP_TOL_M) / 1e3, dh_default_km=W.DH_POPUP_M / 1e3)


def selftest(flights):
    cells = cell_table(flights); assert cells.flights.sum() == len(flights) and (cells.success <= cells.flights).all()
    # cross-check against the Revision 18 table: same arrival time and first burn for the same run
    J = {x["run_id"]: x for x in json.load(open(TASK / "rev18_data" / "rev18_flights.json"))["flights"]}; n = 0
    for f in flights:
        if f["run_id"] in J:
            x = J[f["run_id"]]; assert abs(f["t_c"] / 86400 - x["t_arr_d"]) < 1e-6, (f["run_id"], f["t_c"] / 86400, x["t_arr_d"]); n += 1
            st = prep(f, 250e3); assert abs(float(W.burn(f["ra_f"], f["rp_f"], st["rp_des"])) - x["fb_250"]) < 1e-9; assert abs(f["dv"] - x["dv"]) < 1e-9
    assert n >= 0.9 * len(flights), n
    lev = np.array([float(W.lever_T(f["ra_f"], f["rp_f"])) for f in flights if f["success"]]); assert 8.5 < np.median(lev) < 9.3, np.median(lev)       # Revision 18: 8.896 s per m/s
    # t* = own arrival -> window met without toolbox, zero toolbox cost, total = dV + FB
    S = [f for f in flights if f["success"]][:6]; R = evaluate(S, only_median=False, geoms=sorted({f["geom"] for f in S}))
    f0 = S[0]; st = prep(f0, 250e3); p = W.cheapest(f0["t_c"], f0["t_c"], 3 * f0["T_f"], st); assert p["cost"] == 0 and p["kind"] == "none"
    assert (R[R.met_none].tb_cost == 0).all() and np.allclose(R[R.feasible].total, R[R.feasible].dv + R[R.feasible].tb_cost + R[R.feasible].fb)
    # an early arrival costs more than zero, a late arrival beyond the phasing authority fails
    e = W.cheapest(f0["t_c"], f0["t_c"] + 10 * 86400, 3 * f0["T_f"], st); assert e is not None and e["cost"] > 0 and abs(e["t_arr"] - f0["t_c"] - 10 * 86400) <= 3 * f0["T_f"] + 1
    assert W.cheapest(f0["t_c"], f0["t_c"] - 10 * 86400, 3 * f0["T_f"], st) is None
    # summarize schema
    T = summarize(R); assert {"geom", "arm", "bin", "met_none", "feasible", "total_median"} <= set(T.columns)
    print("selftest OK:", len(flights), "flights,", int(cells.success.sum()), "successes; median final-orbit lever %.3f s per m/s; cross-checked %d runs against rev18_flights.json" % (np.median(lev), n))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(TASK / "rev20_data" / "baselines")); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--dh-table", action="store_true"); a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True); flights = load_flights()
    if a.selftest: selftest(flights); return
    if a.dh_table:
        t = dh_table(); json.dump(t, open(out / "dh_table.json", "w"), indent=1); print(t); return
    cells = cell_table(flights); cells.to_csv(out / "cells.csv", index=False); print(cells.to_string())
    R = evaluate(flights); R.to_csv(out / "rows.csv", index=False)
    T = summarize(R); T.to_csv(out / "summary_by_offset.csv", index=False)
    M = evaluate(flights, only_median=True); M.to_csv(out / "rows_median_tstar.csv", index=False)
    Tm = summarize(M, by=("geom", "arm")); Tm.to_csv(out / "summary_median_tstar.csv", index=False); print(Tm.to_string())
    RR = evaluate(flights, key="relaxed"); RR.to_csv(out / "rows_relaxed.csv", index=False); summarize(RR).to_csv(out / "summary_by_offset_relaxed.csv", index=False)
    MR = evaluate(flights, only_median=True, key="relaxed"); summarize(MR, by=("geom", "arm")).to_csv(out / "summary_median_tstar_relaxed.csv", index=False); print(summarize(MR, by=("geom", "arm")).to_string())
    sens = []
    for h in H_DES_SENS_KM:
        S = summarize(evaluate(flights, h_des_m=h * 1e3, only_median=True), by=("geom", "arm")); S["h_des_km"] = h; S["dh_km"] = W.DH_POPUP_M / 1e3; sens.append(S)
    for dh in DH_SENS_KM:
        S = summarize(evaluate(flights, dh=dh * 1e3, only_median=True), by=("geom", "arm")); S["h_des_km"] = 250.0; S["dh_km"] = dh; sens.append(S)
    pd.concat(sens).to_csv(out / "sensitivity_median_tstar.csv", index=False)
    json.dump(dict(W_MAX_PERIODS=W.W_MAX_PERIODS, H_DES_M=W.H_DES_M, DH_POPUP_M=W.DH_POPUP_M, DROP_TOL_M=W.DROP_TOL_M, H_PH_MAX_M=W.H_PH_MAX_M, K_PH_MAX=W.K_PH_MAX, ERR_TOL_KM=[ERR_LO_KM, ERR_HI_KM], Q_HARD=Q_HARD, RELAXED_RULE='reached and no hard violation (final-apoapsis tolerance ignored)', NOMINAL_D=NOMINAL_D,
                   STEP_D=STEP_D, PAD_D=PAD_D, BIN_EDGES=BIN_EDGES), open(out / "parameters.json", "w"), indent=1)


if __name__ == "__main__":
    main()
