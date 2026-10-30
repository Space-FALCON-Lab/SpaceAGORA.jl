"""Revision 18 flight statistics (analysis only; reads the 366 flights of dv_correction/flights.json): closing-apoapsis time spread per arm@geometry against W_max = 3 final-orbit periods, apse-angle drift against the secular J2 prediction,
the first-burn dV (periapsis raise at the closing apoapsis to a named ASSUMED periapsis altitude) and its effect on the arm ranking, and the flight check of the apoapsis-burn lever and of dT = 3 T sqrt((1-e)/(1+e)) dv / V_c.
Usage: python3 -I planner_rev18_flights.py --flights F.json --out F2.json [--h-des-km 250]"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_lib as PLB
import geometry_rank as G
MU = A.MU; RA_T = A.RA_TARGET; PL = G.K.planet("Mars", 20); J2 = PL["J"][2]; RJ = PL["R"]
TOK = {"ODY": "ODY", "G5863": "G1", "G0975": "GC", "nan": "ODY"}
H_DES_KM = [120.0, 200.0, 250.0, 300.0, 400.0]      # ASSUMED sensitivity set for the desired periapsis altitude; the base value is a named parameter (default 250 km)


def run_dir(root, phase, rid):
    root = pathlib.Path(root)
    for c in (root / phase / rid, root / "runs" / rid, root / "mars" / "runs" / rid, root / "dev" / rid, root / "sens" / rid, root / "gen1" / rid):
        if (c / "passes.csv").exists(): return c
    return None


def first_burn_dv(ra, rp, rp_des):
    """Tangential periapsis-raise burn at apoapsis: dv = v_a(a') - v_a(a), v_a = sqrt(mu (2/ra - 1/a)), a = (ra + rp)/2."""
    va = lambda r_p: math.sqrt(MU * (2 / ra - 2 / (ra + r_p)))
    return va(rp_des) - va(rp)


def flight(d, row, h_des_m):
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= RA_T))[0]
    n = int(hit[0] + 1) if len(hit) else len(P); i = n - 1; a = P.a_m.values.astype(float); e = P.e.values.astype(float); inc = P.i_rad.values.astype(float)
    Tn = 2 * math.pi * math.sqrt(a[i] ** 3 / MU); t_arr = float(P.time_s.iloc[i]) + 0.5 * Tn
    om = np.unwrap(P.omega_rad.values.astype(float)); pred = sum(PLB_apsidal(a[k], e[k], inc[k]) for k in range(i)); obs = float(om[i] - om[0])
    ra_f = float(P.ra_m.iloc[i]); rp_f = float(P.rp_m.iloc[i]); alt_f = float(P.alt_m.iloc[i]); out = dict(run_id=row["run_id"], arm=row["arm"], geom=TOK[str(row["geometry"])], phase=row["summary"].split("_")[1], n=n, t_arr_d=t_arr / 86400, T_final_s=Tn,
        omega_obs_deg=math.degrees(obs), omega_j2_deg=math.degrees(pred), ra_f=ra_f, rp_f=rp_f, alt_f=alt_f, e_f=float(e[i]), a_f=float(a[i]), dv=row["dv_corrected"], omega_final_deg=math.degrees(float(om[i])))
    for h in H_DES_KM + [h_des_m / 1e3]:
        rp_des = rp_f - alt_f + h * 1e3; out[f"fb_{h:g}"] = first_burn_dv(ra_f, rp_f, rp_des)
    return out, P


def PLB_apsidal(a, e, i):
    p = a * (1 - e * e); return 3 * math.pi * J2 * (RJ / p) ** 2 * (2 - 2.5 * math.sin(i) ** 2)


def lever_checks(rows_dirs):
    """Regress the pass-to-pass change of the periapsis radius on the derived apoapsis lever (burn before pass k+1 = decision target_pass k+1, elements of the orbit at its closing apoapsis (ra, rp of pass k-1)),
    and the period change on the derived timing lever dT with a drag regressor, overall and by orbit size a (bins) and for isolated burns (no burn on the neighboring orbits)."""
    rec = []
    for d, run in rows_dirs:
        P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); D = pd.read_csv(d / "decisions.csv")
        dv = {int(r.target_pass): float(r.delta_v_mps) for r in D.itertuples()}; rp = P.rp_m.values.astype(float); ra = P.ra_m.values.astype(float); t = P.time_s.values.astype(float); a_m = P.a_m.values.astype(float)
        for k in range(3, len(P) - 3):
            b = dv.get(k + 1, 0.0)
            if not (np.isfinite(ra[k - 1]) and np.isfinite(ra[k - 2]) and np.isfinite(ra[k + 1])): continue
            a_apo = 0.5 * (ra[k - 1] + rp[k - 1]); e_apo = (ra[k - 1] - rp[k - 1]) / (ra[k - 1] + rp[k - 1]); T = 2 * math.pi * math.sqrt(a_apo ** 3 / MU)
            rec.append(dict(b=b, a=a_apo, x_rp=b * PLB.lever_m_per_mps(MU, a_apo, e_apo), y_rp=rp[k] - rp[k - 1], x_T=b * lever_T(a_apo, e_apo), y_T=(t[k + 2] - t[k + 1]) - (t[k - 1] - t[k - 2]),
                            z_drag=1.5 * T / a_apo * (-(ra[k - 2] - ra[k + 1]) / 2.0), iso=(dv.get(k, 0.0) == 0.0 and dv.get(k + 2, 0.0) == 0.0)))
    R = pd.DataFrame(rec); out = dict(n_passes=int(len(R)), n_burns=int((R.b != 0).sum()))
    def fit(S, tag):
        if (S.b != 0).sum() < 8: return
        c = np.polyfit(S.x_rp, S.y_rp, 1); r = S.y_rp - np.polyval(c, S.x_rp); se = math.sqrt(np.sum(r ** 2) / (len(S) - 2) / np.sum((S.x_rp - S.x_rp.mean()) ** 2))
        M = np.column_stack([S.x_T, S.z_drag, np.ones(len(S))]); co = np.linalg.lstsq(M, S.y_T.values, rcond=None)[0]; rs = S.y_T.values - M @ co; cov = np.linalg.inv(M.T @ M) * np.sum(rs ** 2) / (len(S) - 3)
        out[tag] = dict(n=int(len(S)), n_burns=int((S.b != 0).sum()), rp_slope=float(c[0]), rp_slope_se=float(se), T_burn_coef=float(co[0]), T_burn_se=float(math.sqrt(cov[0, 0])), T_drag_coef=float(co[1]))
    fit(R, "all")
    for lo, hi in ((0, 8e6), (8e6, 1.2e7), (1.2e7, 2e7)):
        S = R[(R.a >= lo) & (R.a < hi)]; fit(S, f"a_{lo:.0e}_{hi:.0e}"); fit(S[S.iso | (S.b == 0)], f"a_{lo:.0e}_{hi:.0e}_isolated")
    return out


def lever_T(a, e):
    return 3 * (2 * math.pi * math.sqrt(a ** 3 / MU)) * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights", required=True); ap.add_argument("--out", required=True); ap.add_argument("--h-des-km", type=float, default=250.0); a = ap.parse_args()
    R = json.load(open(a.flights)); rows = []; dirs = []
    for r in R:
        if r.get("status") == "not_corrected": continue
        d = run_dir(r["root"], r["phase"], r["run_id"])
        if d is None: continue
        o, P = flight(d, r, a.h_des_km * 1e3); rows.append(o)
        if r["root"].endswith("planner") and r["phase"] in ("held", "dev", "devT", "devR", "gen1", "sens"): dirs.append((d, r["run_id"]))
    df = pd.DataFrame(rows); chk = lever_checks(dirs)
    json.dump(dict(flights=rows, lever_checks=chk, h_des_km=a.h_des_km, H_DES_KM=H_DES_KM), open(a.out, "w"), indent=1, default=float)
    print(len(df), "flights;", chk)


if __name__ == "__main__":
    main()
