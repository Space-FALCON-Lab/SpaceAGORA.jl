"""Revision 18 report data (analysis only): dV-versus-arrival-time curves per geometry (lambda sweep and a band family), D_min / D_max, model-world arm points, and the flight statistics tables. Writes JSON; prints markdown tables.
Usage: python3 -I planner_rev18_report.py --flights-json rev18_flights.json --out rev18_report.json"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev18_offline as O
import planner_lib as PLB
WMAX_ORBITS = 3          # ASSUMED named parameter W_max (final-orbit periods)
BAND = 0.10              # ASSUMED half width (ln q) of the band family used for the dense curve
LAMS = sorted(set([0.0] + [s * float(v) for s in (-1, 1) for v in list(np.logspace(-3.5, -1.3, 36)) + [0.1, 0.5]]))     # lambda_orb grid (m/s per orbit): ASSUMED sweep range


def band_family(tab, centers):
    inp = PLB.ol_inputs(dict(), tab, O.P, math.log(O.Q_FIRST)); lvl = np.array(inp["level"]); lo = np.array(inp["lo"]); hi = np.array(inp["hi"]); w = np.array(inp["w"]); rows = []
    for c in centers:
        b_lo = np.maximum(lo, c - BAND); b_hi = np.minimum(hi, c + BAND)
        if np.any(b_lo > b_hi + 1e-12): continue
        nu = O.dp_lin(b_lo - lvl, b_hi - lvl, PLB.tie_weights(w, 1e-6), np.zeros(len(lvl)))
        if nu is None: continue
        du = np.diff(np.concatenate(([0.0], nu))); x = np.concatenate(([math.log(O.Q_FIRST)], lvl + nu)); r = O.arrival(tab, x)
        rows.append(dict(center_ln_q=float(c), dv=float(np.sum(w * np.abs(du))), days=r["days"], N=r["N"]))
    return rows


def model_arms(tab, token):
    out = {}
    B = PLB.flat_P(O.CFG)
    for arm in ("OL", "OLR", "L0", "L1", "L2", "L3"):
        sw = PLB.arm_switches(O.CFG, arm); Pa = PLB.flat_P(O.CFG, arm)
        try:
            s = PLB.simulate(sw, tab, Pa, {}, math.log(O.Q_FIRST), fast=True)
        except Exception as e:
            out[arm] = dict(error=str(e)[:80]); continue
        x = np.asarray(s["lnq"], float); r = O.arrival(tab, x); out[arm] = dict(dv=s["dv"], burns=s["burns"], days=r["days"], N=r["N"])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights-json", required=True); ap.add_argument("--out", required=True); ap.add_argument("--kappa-json"); a = ap.parse_args(); res = {}
    if a.kappa_json:
        for g in json.load(open(a.kappa_json))["geometry"]: O.KAPPA[g["geom"]] = float(g["kappa_med"])
    res["_kappa"] = dict(O.KAPPA)
    for tok in ("ODY", "G1"):
        tab = O.load_tab(tok); tab["token"] = tok; r0 = O.plan_for_price(tab, 0.0); hi = O.plan_for_price(tab, 0.0, pinned="hi"); lo = O.plan_for_price(tab, 0.0, pinned="lo")
        lam = []
        for l in LAMS:
            r = O.plan_for_price(tab, l)
            if r: lam.append(dict(lam_orb=l, dv=r["dv"], days=r["days"], N=r["N"], burns=r["burns"]))
        inp = PLB.ol_inputs(dict(), tab, O.P, math.log(O.Q_FIRST)); centers = np.linspace(float(inp["lo"][0]) + BAND, float(inp["hi"][0]) - BAND, 25)
        res[tok] = dict(nominal=dict(dv=r0["dv"], days=r0["days"], N=r0["N"], burns=r0["burns"]), D_min=dict(dv=hi["dv"], days=hi["days"], N=hi["N"]), D_max=dict(dv=lo["dv"], days=lo["days"], N=lo["N"]), lam=lam, band=band_family(tab, centers),
                        model_arms=model_arms(tab, tok), q_req=float(tab["q_req"]), W=float(tab["W"]), T_nom_median_s=float(np.median(np.diff(tab["t"]))), T_first_s=float(np.diff(tab["t"])[0]), T_last_s=float(np.diff(tab["t"])[-1]))
    json.dump(res, open(a.out, "w"), indent=1, default=float); print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ("nominal", "D_min", "D_max", "model_arms")} for k, v in res.items()}, indent=1, default=float))


if __name__ == "__main__":
    main()
