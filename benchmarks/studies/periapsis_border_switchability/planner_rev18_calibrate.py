"""Revision 18b (analysis only): diagnose and calibrate the in-model arrival time. Per flight: (1) period check T_obs / T(a), (2) drop-per-q ratio flown / model with the drop read at the pass index (old) and at the current apoapsis radius (new), (3) per-geometry
calibration factor kappa = sum(flown drops) / sum(model drops) over the flight, (4) validation of the arrival with kappa. Usage: python3 -I planner_rev18_calibrate.py --flights-src F --flights G --out H"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev18_offline as O
import planner_rev18_flights as FL


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights-src", required=True); ap.add_argument("--flights", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    src = {r["run_id"]: r for r in json.load(open(a.flights_src)) if r.get("status") != "not_corrected"}; J = json.load(open(a.flights))["flights"]; tabs = {g: O.load_tab(g) for g in ("ODY", "G1")}; rows = []
    for f in J:
        g = f["geom"]
        if g not in tabs: continue
        tab = tabs[g]; r = src[f["run_id"]]; d = FL.run_dir(r["root"], r["phase"], r["run_id"]); P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); n = f["n"]
        ra = P.ra_m.values.astype(float)[:n]; q = P.q_w_cm2.values.astype(float)[:n]; rp = P.rp_m.values.astype(float); t = P.time_s.values.astype(float)
        ra_prev = np.concatenate(([float(tab["ra"][0])], ra[:-1])); drop = ra_prev - ra; d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])]))
        m_idx = np.array([d0[min(k, tab["n"]) - 1] * q[k - 1] / tab["q_req"] for k in range(1, n + 1)]); m_st = np.interp(ra_prev, tab["ra"][::-1], d0[::-1]) * q / tab["q_req"]; ok = np.isfinite(drop) & (m_st > 0)
        Tobs = np.diff(t[:n]); Ta = 2 * np.pi * np.sqrt(((ra[:-1] + rp[:n - 1]) / 2) ** 3 / O.MU); okT = np.isfinite(Ta)
        kappa = float(np.sum(drop[ok]) / np.sum(m_st[ok])); x = np.log(q)
        m0 = O.arrival(tab, x, math.log(q[0]), kappa=1.0); m1 = O.arrival(tab, x, math.log(q[0]), kappa=kappa)
        rows.append(dict(cell=f["arm"].replace("-", "") + "@" + g, geom=g, flown_d=f["t_arr_d"], T_ratio=float(np.nanmedian(Tobs[okT[:len(Tobs)]] / Ta[okT][:len(Tobs)])) if okT.any() else float("nan"), ratio_idx=float(np.sum(drop[ok]) / np.sum(m_idx[ok])), ratio_state=kappa,
                         model_state_d=m0["days"], model_cal_d=m1["days"], H_ratio=float(np.nanmedian(P.h_local_m.values[:n] / np.interp(ra_prev, tab["ra"][::-1], np.array(tab["H"])[::-1])))))
    D = pd.DataFrame(rows); D["val_state"] = D.flown_d / D.model_state_d; D["val_cal"] = D.flown_d / D.model_cal_d
    g = D.groupby("cell").agg(n=("val_state", "size"), T_ratio=("T_ratio", "median"), ratio_idx=("ratio_idx", "median"), kappa=("ratio_state", "median"), kappa_min=("ratio_state", "min"), kappa_max=("ratio_state", "max"), H_ratio=("H_ratio", "median"),
                              val_state=("val_state", "median"), val_cal=("val_cal", "median"), val_cal_min=("val_cal", "min"), val_cal_max=("val_cal", "max"))
    geo = D.groupby("geom").agg(n=("ratio_state", "size"), kappa_med=("ratio_state", "median"), kappa_iqr_lo=("ratio_state", lambda s: s.quantile(.25)), kappa_iqr_hi=("ratio_state", lambda s: s.quantile(.75)), val_state=("val_state", "median"), val_cal=("val_cal", "median"), H_ratio=("H_ratio", "median"))
    out = dict(cells=g.reset_index().to_dict("records"), geometry=geo.reset_index().to_dict("records"), n=int(len(D)))
    json.dump(out, open(a.out, "w"), indent=1, default=float); print(geo.round(3).to_string()); print(g.round(3).to_string())


if __name__ == "__main__":
    main()
