"""Revision 19 offline evaluation (analysis only, Python): a terminal-periapsis term in the DP. The first burn (periapsis raise at the closing apoapsis to h_des) is exact from vis-viva; its local slope per unit periapsis radius is checked against the existing burn lever.
(1) slope check, (2) model-world open-loop plans on the ground tables with and without the terminal term (linear cost c_T nu_N on the level at the arrival pass), (3) an upper bound on what the 366 flights could save if their final level were moved to the lower tube edge by an early burn.
Usage: python3 -I planner_rev19_offline.py --flights rev18_flights.json --out F.json [--h-des-km 250]"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_analysis as _AN
import planner_rev18_offline as O
import planner_rev18_flights as FL
MU = O.MU; P = O.P; H_REF = 105319.6      # SOURCED geometry_rank_config.toml h_ref_m: the altitude at nu = 0 (the held tube level)
RES = P["fast_reserve_ln"]


def fb(ra, rp, h_des_m, alt):                 # first burn (m/s): raise from rp (altitude alt) to altitude h_des at the same apse point
    return FL.first_burn_dv(ra, rp, rp - alt + h_des_m)


def slope_check(h_des_m=250e3):
    """dFB/drp (m/s per km, finite difference, exact vis-viva) against 1/lever (existing burn lever at the closing-apoapsis orbit), for the final states of the flights."""
    rows = []
    for ra, rp in ((_AN.RA_TARGET, 3.480e6), (_AN.RA_TARGET, 3.484e6), (_AN.RA_TARGET, 3.49e6)):
        a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp); lev = PLB.lever_m_per_mps(MU, a, e)
        d = (FL.first_burn_dv(ra, rp + 500.0, rp + 500.0 + 150e3) - FL.first_burn_dv(ra, rp - 500.0, rp - 500.0 + 150e3)) / 1000.0   # m/s per m, both ends moved with the same altitude gap: NOT the quantity; see below
        d2 = (FL.first_burn_dv(ra, rp + 500.0, rp + 150e3) - FL.first_burn_dv(ra, rp - 500.0, rp + 150e3)) / 1000.0                   # rp varies, target rp_des fixed: dFB/drp at fixed desired periapsis
        rows.append(dict(ra=ra, rp=rp, lever_m_per_mps=lev, inv_lever_mps_per_km=1000.0 / lev, dFB_drp_fixed_target_mps_per_km=-d2 * 1000.0, dFB_dalt_mps_per_km=-d * 1000.0))
    return rows


def model_plan(tab, with_term, h_des_m=250e3, iters=3):
    """OL-style plan on the table; with_term adds c_T nu_N at the arrival pass N (linear). Returns dict with campaign dV up to the arrival pass, final level, rp, FB, totals, days."""
    inp = PLB.ol_inputs(dict(), tab, P, math.log(O.Q_FIRST)); lvl = np.array(inp["level"]); lo = np.array(inp["lo"]); hi = np.array(inp["hi"]); w = np.array(inp["w"]); tw = PLB.tie_weights(w, 1e-6)
    lin = np.zeros(len(lvl)); N = None; nu = None
    for it in range(iters + 1):
        nu = O.dp_lin(lo - lvl, hi - lvl, tw, lin); x = np.concatenate(([math.log(O.Q_FIRST)], lvl + nu)); r = O.arrival(tab, x); Nn = r["N"]
        if not with_term or Nn == N: N = Nn; break
        N = Nn; kk = min(N, tab["n"]); ra_f = _AN.RA_TARGET; rp_t = float(tab["rp"][kk - 1]); H = float(tab["H"][kk - 1])
        d = (FL.first_burn_dv(ra_f, rp_t - 250.0, rp_t - 250.0 + 150e3) - FL.first_burn_dv(ra_f, rp_t + 250.0, rp_t + 250.0 + 150e3)) / 500.0     # m/s per m, alt fixed gap: dFB/drp at fixed gap
        dFB_drp = (FL.first_burn_dv(ra_f, rp_t - 250.0, rp_t + 150e3) - FL.first_burn_dv(ra_f, rp_t + 250.0, rp_t + 150e3)) / 500.0       # m/s per m of periapsis radius lost, target fixed
        lin = np.zeros(len(lvl)); lin[max(0, N - 2)] = H * dFB_drp                                                                      # nu up -> rp down by H nu -> first burn up by H dFB_drp nu
    du = np.diff(np.concatenate(([0.0], nu))); steps = np.arange(2, 2 + len(nu)); mk = steps <= N; dv_c = float(np.sum((w * np.abs(du))[mk]))
    kk = min(N, tab["n"]); nuN = float(nu[N - 2]) if N >= 2 else 0.0; rp_f = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nuN; alt_f = H_REF - float(tab["H"][kk - 1]) * nuN
    f = fb(_AN.RA_TARGET, rp_f, h_des_m, alt_f); x = np.concatenate(([math.log(O.Q_FIRST)], lvl + nu)); xf = float(x[min(N, len(x)) - 1])
    return dict(N=N, days=O.arrival(tab, x)["days"], dv_campaign=dv_c, nuN=nuN, ln_q_final=xf, rp_final=rp_f, alt_final_km=alt_f / 1e3, first_burn=f, total=dv_c + f, burns=int(((np.abs(du) * w >= P["dp_min_burn_mps"]) & mk).sum()))


def flights_bound(Fl, h_des_m):
    rows = []; tabs = {g: O.load_tab(g) for g in ("ODY", "G1")}
    for f in Fl["flights"]:
        g = f["geom"]
        if g not in tabs: continue
        tab = tabs[g]; hi0 = math.log(P["q_cap_mean_w_cm2"]); lo_edge = hi0 - float(tab["W"]) + RES / 2.0
        d0 = FL.run_dir  # placeholder to keep import used
        rows.append(dict(cell=f["arm"].replace("-", "") + "@" + g, run_id=f["run_id"], geom=g, dv=f["dv"], rp_f=f["rp_f"], alt_f=f["alt_f"], ra_f=f["ra_f"], x_lo=lo_edge, wmin=float(np.min(tab["w"][1:])), H=float(tab["H"][-1])))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights", required=True); ap.add_argument("--flights-src", required=True); ap.add_argument("--out", required=True); ap.add_argument("--h-des-km", type=float, default=250.0); a = ap.parse_args(); h = a.h_des_km * 1e3
    out = dict(slope=slope_check(h)); kap = {}
    cal = pathlib.Path(a.flights).parent / "rev18_calibrate.json"
    if cal.exists():
        for g in json.load(open(cal))["geometry"]: O.KAPPA[g["geom"]] = float(g["kappa_med"])
    out["model"] = {}
    for tok in ("ODY", "G1"):
        tab = O.load_tab(tok); out["model"][tok] = dict(no_term=model_plan(tab, False, h), with_term=model_plan(tab, True, h))
    # flights: need the flown ln q at the last counted pass and the table lower edge
    src = {r["run_id"]: r for r in json.load(open(a.flights_src)) if r.get("status") != "not_corrected"}; Fl = json.load(open(a.flights)); tabs = {g: O.load_tab(g) for g in ("ODY", "G1")}; rows = []
    for f in Fl["flights"]:
        g = f["geom"]
        if g not in tabs: continue
        r = src[f["run_id"]]; d = FL.run_dir(r["root"], r["phase"], r["run_id"]); Pp = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); i = f["n"] - 1
        tab = tabs[g]; lo_edge = math.log(P["q_cap_mean_w_cm2"]) - float(tab["W"]) + RES / 2.0; x_f = math.log(float(Pp.q_w_cm2.iloc[i])); Hf = float(Pp.h_local_m.iloc[i]); delta = max(0.0, x_f - lo_edge)
        # raising the final level to the lower tube edge by an early burn: periapsis +H*delta, cost at least delta * w_min (cheapest table weight)
        fb0 = fb(f["ra_f"], f["rp_f"], h, f["alt_f"]); fb1 = fb(f["ra_f"], f["rp_f"] + Hf * delta, h, f["alt_f"] + Hf * delta); wmin = float(np.min(tab["w"][1:]))
        rows.append(dict(cell=f["arm"].replace("-", "") + "@" + g, geom=g, dv=f["dv"], fb0=fb0, fb1=fb1, delta_ln=delta, gross=fb0 - fb1, cost_min=wmin * delta, net_bound=(fb0 - fb1) - wmin * delta, x_f=x_f, x_lo=lo_edge, wmin=wmin))
    D = pd.DataFrame(rows); D["total0"] = D.dv + D.fb0; D["total_bound"] = D.dv + D.fb0 - D.net_bound.clip(lower=0)
    C = D.groupby("cell").agg(n=("dv", "size"), dv=("dv", "median"), fb0=("fb0", "median"), delta_ln=("delta_ln", "median"), gross=("gross", "median"), cost_min=("cost_min", "median"), net_bound=("net_bound", "median"), total0=("total0", "median"), total_bound=("total_bound", "median")).reset_index()
    out["cells"] = C.to_dict("records"); out["overall"] = dict(net_bound_median=float(D.net_bound.median()), net_bound_p10=float(D.net_bound.quantile(.1)), net_bound_p90=float(D.net_bound.quantile(.9)), gross_median=float(D.gross.median()), delta_median=float(D.delta_ln.median()))
    json.dump(out, open(a.out, "w"), indent=1, default=float); print(json.dumps(out["slope"], indent=1)); print(json.dumps(out["model"], indent=1)); print(C.round(3).to_string()); print(out["overall"])


if __name__ == "__main__":
    main()
