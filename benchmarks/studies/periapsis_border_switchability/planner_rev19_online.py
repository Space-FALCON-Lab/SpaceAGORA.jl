"""Revision 19 (analysis only): noise-free model-world replay of the online replanner arms (OLR, L3, L2, L1) with and without the terminal-periapsis term. Same loop as planner_lib.simulate; the DP of each decision is dp_lin with a linear cost
c_T nu on the level at the arrival pass N (N from the arm's own no-term run). Total = campaign dV of the decisions up to pass N + first burn FB at h_des. Usage: python3 -I planner_rev19_online.py --out F [--h-des-km 250]"""
import argparse, json, math, pathlib, sys
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_analysis as _AN
import planner_rev18_offline as O
import planner_rev19_offline as R19
import planner_rev18_flights as FL
P = O.P; MU = O.MU


def sim(arm, tab, term_N=None, cT=0.0):
    sw = PLB.arm_switches(O.CFG, arm); Pa = PLB.flat_P(O.CFG, arm); n = int(tab["n"]); path = np.asarray(tab["path"], float); q0 = math.log(O.Q_FIRST); hi0 = math.log(Pa["q_cap_mean_w_cm2"])
    nu = 0.0; dvs = {}; nus = np.zeros(n); base = lambda p: q0 + float(path[min(p, len(path)) - 1] - path[0])
    for p in range(1, n + 1):
        if p >= 2:
            z0 = base(p) + nu; w0 = PLB.tget(tab["w"], p); inp = PLB.plan_inputs(sw, tab, Pa, p, z0, w0, 0.0)
            lvl = np.asarray(inp["level"], float); lin = np.zeros(len(lvl))
            if term_N is not None and 0 <= term_N - p < len(lvl): lin[term_N - p] = cT
            nuv = O.dp_lin(np.asarray(inp["lo"], float) - lvl, np.asarray(inp["hi"], float) - lvl, PLB.tie_weights(inp["w"], Pa["dp_tie_eps"]), lin)
            dv = -float(nuv[0]) * w0 if nuv is not None else 0.0
            if abs(dv) < Pa["dp_min_burn_mps"]: dv = 0.0
            if dv != 0.0: dvs[p] = dv; nu += -dv / w0
        nus[p - 1] = nu
    x = np.array([base(p) + nus[p - 1] for p in range(1, n + 1)]); return dict(dvs=dvs, nus=nus, x=x)


def evaluate(arm, tab, h_des_m, term):
    r0 = sim(arm, tab); a0 = O.arrival(tab, r0["x"]); N0 = a0["N"]
    if term:
        kk = min(N0, tab["n"]); rp_t = float(tab["rp"][kk - 1]); H = float(tab["H"][kk - 1]); ra_f = _AN.RA_TARGET
        dFB = (FL.first_burn_dv(ra_f, rp_t - 250.0, rp_t + 150e3) - FL.first_burn_dv(ra_f, rp_t + 250.0, rp_t + 150e3)) / 500.0; r = sim(arm, tab, N0, H * dFB)
    else:
        r = r0
    a = O.arrival(tab, r["x"]); N = a["N"]; kk = min(N, tab["n"]); nuN = float(r["nus"][min(N, tab["n"]) - 1]); rp_f = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nuN; alt_f = R19.H_REF - float(tab["H"][kk - 1]) * nuN
    dvc = float(sum(abs(v) for p, v in r["dvs"].items() if p <= N)); f = R19.fb(_AN.RA_TARGET, rp_f, h_des_m, alt_f)
    return dict(N=N, days=a["days"], dv_campaign=dvc, burns=int(sum(1 for p in r["dvs"] if p <= N)), nuN=nuN, alt_final_km=alt_f / 1e3, first_burn=f, total=dvc + f)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--h-des-km", type=float, default=250.0); ap.add_argument("--kappa-json"); a = ap.parse_args(); h = a.h_des_km * 1e3; out = {}
    if a.kappa_json:
        for g in json.load(open(a.kappa_json))["geometry"]: O.KAPPA[g["geom"]] = float(g["kappa_med"])
    for tok in ("ODY", "G1"):
        tab = O.load_tab(tok); out[tok] = {}
        for arm in ("OLR", "L3", "L2", "L1"):
            out[tok][arm] = dict(no_term=evaluate(arm, tab, h, False), term=evaluate(arm, tab, h, True))
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    for tok, d in out.items():
        for arm, v in d.items():
            n, t = v["no_term"], v["term"]; print(f"{tok} {arm}: total {n['total']:.3f} -> {t['total']:.3f} ({t['total'] - n['total']:+.3f}); campaign {n['dv_campaign']:.2f} -> {t['dv_campaign']:.2f}; alt {n['alt_final_km']:.1f} -> {t['alt_final_km']:.1f} km; days {n['days']:.1f} -> {t['days']:.1f}; burns {n['burns']} -> {t['burns']}; N {n['N']} -> {t['N']}")


if __name__ == "__main__":
    main()
