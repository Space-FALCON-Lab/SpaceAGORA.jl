"""Revision 18 flight tables from rev18_flights.json (analysis only). Usage: python3 -I planner_rev18_tables.py --flights F --out-json G --out-md H"""
import argparse, json, math, pathlib
import numpy as np, pandas as pd
WM = 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights", required=True); ap.add_argument("--out-json", required=True); ap.add_argument("--out-md", required=True); a = ap.parse_args()
    J = json.load(open(a.flights)); df = pd.DataFrame(J["flights"]); H = J["H_DES_KM"]; df["cell"] = df.arm.str.replace("-", "") + "@" + df.geom
    df["Wmax_d"] = WM * df.T_final_s / 86400.0; df["resid_deg"] = df.omega_obs_deg - df.omega_j2_deg; df["resid_ratio"] = (df.resid_deg.abs() / df.omega_j2_deg.abs().clip(lower=1e-9))
    rows = []
    for c, g in df.groupby("cell"):
        if len(g) < 4: continue
        m = g.t_arr_d.median(); W = g.Wmax_d.median()
        rows.append(dict(cell=c, n=len(g), t_med=m, t_std=g.t_arr_d.std(), t_iqr=g.t_arr_d.quantile(.75) - g.t_arr_d.quantile(.25), t_range=g.t_arr_d.max() - g.t_arr_d.min(), T_final_h=g.T_final_s.median() / 3600, Wmax_d=W,
                         std_over_W=g.t_arr_d.std() / W, frac_within=float(((g.t_arr_d - m).abs() <= g.Wmax_d).mean()), j2_deg=g.omega_j2_deg.median(), obs_deg=g.omega_obs_deg.median(), resid_med=g.resid_deg.median(), resid_ratio_med=g.resid_ratio.median(),
                         omega_f_std=g.omega_final_deg.std(), arc_km=math.radians(g.omega_final_deg.std()) * g.ra_f.median() / 1e3, dv=g.dv.median(), rp_f=g.rp_f.median(), alt_f=g.alt_f.median(), **{f"fb_{h:g}": g[f"fb_{h:g}"].median() for h in H}))
    T = pd.DataFrame(rows).sort_values("cell")
    # rankings: by dv, and by dv + first burn, per h_des, within the Odyssey cells and within G1 cells
    rk = {}
    for geom in ("ODY", "G1", "GC"):
        S = T[T.cell.str.endswith("@" + geom)]
        if len(S) < 2: continue
        rk[geom] = {"dv": list(S.sort_values("dv").cell)}
        for h in H: rk[geom][f"dv+fb{h:g}"] = list(S.assign(tot=S.dv + S[f"fb_{h:g}"]).sort_values("tot").cell)
    out = dict(table=T.to_dict("records"), rankings=rk, overall=dict(first_burn_median_250=float(df["fb_250"].median()), first_burn_min_250=float(df["fb_250"].min()), first_burn_max_250=float(df["fb_250"].max()), cell_median_spread_250=float(T["fb_250"].max() - T["fb_250"].min()),
                    resid_ratio_median_all=float(df.resid_ratio.median()), std_over_W_median=float(T.std_over_W.median()), frac_within_max=float(T[T.n >= 8].frac_within.max()), omega_std_median_deg=float(T.omega_f_std.median()), T_final_median_s=float(df.T_final_s.median())))
    json.dump(out, open(a.out_json, "w"), indent=1, default=float)
    f = lambda x, d=2: f"{x:.{d}f}"
    L = ["| cell | n | median arrival (d) | std (d) | IQR (d) | range (d) | T_final (h) | W_max (d) | std / W_max | within +/-W_max of cell median | J2 apse change (deg) | observed (deg) | residual (deg) | |resid|/|J2| | std of final omega (deg) | arc at r_a (km) |", "|" + "---|" * 16]
    for r in T.itertuples():
        L.append(f"| {r.cell} | {r.n} | {f(r.t_med,1)} | {f(r.t_std)} | {f(r.t_iqr)} | {f(r.t_range)} | {f(r.T_final_h)} | {f(r.Wmax_d,3)} | {f(r.std_over_W,1)} | {f(r.frac_within,2)} | {f(r.j2_deg,1)} | {f(r.obs_deg,1)} | {f(r.resid_med,1)} | {f(r.resid_ratio_med,2)} | {f(r.omega_f_std,2)} | {f(r.arc_km,0)} |")
    L += ["", "First burn (periapsis raise at the closing apoapsis to altitude h_des; median over the flights of the cell; m/s) and the dV to target (corrected):", "", "| cell | dV to target | final rp altitude (km) | " + " | ".join(f"first burn {h:g} km" for h in H) + " | total at 250 km |", "|" + "---|" * (len(H) + 4)]
    for r in T.itertuples():
        L.append(f"| {r.cell} | {f(r.dv,2)} | {f(r.alt_f/1e3,1)} | " + " | ".join(f(getattr(r, f"fb_{h:g}".replace('.', '_')) if False else T.loc[T.cell == r.cell, f"fb_{h:g}"].iloc[0]) for h in H) + f" | {f(r.dv + T.loc[T.cell == r.cell, 'fb_250'].iloc[0],2)} |")
    L += ["", "Rankings (cell medians, lowest first):"]
    for g, d in rk.items():
        for k, v in d.items(): L.append(f"- {g} {k}: " + " < ".join(v))
    pathlib.Path(a.out_md).write_text("\n".join(L)); print("\n".join(L)); print(json.dumps(out["overall"], indent=1))


if __name__ == "__main__":
    main()
