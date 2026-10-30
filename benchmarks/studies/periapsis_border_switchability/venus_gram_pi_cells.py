"""Venus with GRAM truth, part (c): Pi values per cell at the first truth pass (venus_gram_plan.md section 2c: H = the GRAM local scale height at the pass, companion with the cell-start single-fit H) (definitions of pi_policy_interface_spec.md section 1.2), group ranges and the hand-prediction check (venus_gram_plan.md section 6.4). Python only (numpy, pandas, scipy).
Usage: python3 -I venus_pi_cells.py [--results-root DIR] [--gravity-dir DIR] [--nsum 200] -> <root>/venus/pi/{pi_cells.csv, pi_groups.csv, hand_check.csv, kaula_constants.csv}
"""
import argparse, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import apoapsis_pi_lib as L
import kaula_lib as kl

VC = tomllib.load(open(HERE / "venus_generality_config.toml", "rb"))   # inherited band, grid (the GRAM config only overrides)
Q_LO, Q_HI, Q_HARD = (VC["band"][k] for k in ("q_low_w_cm2", "q_high_w_cm2", "q_hard_w_cm2"))
HORIZON_D = VC["grid"]["horizon_days"]; SPIN = 2.99e-7     # SOURCED: Venus() struct spin magnitude
DEFAULT_GDIR = HERE.parents[2] / "data" / "planets" / "gravity"
NKEY = (2, 3, 5)


def sigma_n(df, nmax):
    s = np.full(nmax + 1, np.nan)
    for n in range(2, nmax + 1):
        g = df[(df.degree == n) & (df.order >= 1)]; s[n] = np.sqrt((g.C ** 2 + g.S ** 2).sum() / (2 * n))
    return s


def lever_fd(a, e, dv=1e-3, mu=1.0):
    ra = a * (1 + e); va = np.sqrt(mu * (2 / ra - 1 / a)); a2 = 1 / (2 / ra - (va + dv) ** 2 / mu); return (2 * a2 - ra) - (2 * a - ra), dv


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--gravity-dir", default=str(DEFAULT_GDIR)); ap.add_argument("--nsum", type=int, default=200)
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); v = root / "venus"; out = v / "pi"; out.mkdir(exist_ok=True)
    pc = pd.read_csv(v / "planet_constants.csv").iloc[0]; bd = pd.read_csv(v / "body_distances.csv").set_index("body")
    mu, R, beta = float(pc.mu_m3s2), float(pc.Re_m), float(pc.beta_kg_m2); J2, J3 = float(pc.J2), float(pc.J3)
    muS, rS = float(bd.loc["Sun", "mu_b_m3s2"]), float(bd.loc["Sun", "r_b_m"])
    gdf, Rf, muf, _ = kl.read_gravity("Venus", a.gravity_dir); assert abs(Rf - R) < 1 and abs(muf - mu) < 1e3
    n_, A_ = kl.degree_amplitude(gdf); nz = n_[np.isfinite(A_) & (A_ > 0)]; K, scat, nused, slope = kl.fit_K(n_, A_, 10, int(nz.max()))
    Jf = kl.zonal_J(gdf); nfile = max(Jf); sig = sigma_n(gdf, 20)
    pd.DataFrame([dict(planet="Venus", K=K, scatter_dex=scat, n_used=nused, free_slope=slope, nmax_file=int(nz.max()), J2_file=Jf[2], J3_file=Jf[3])]).to_csv(out / "kaula_constants.csv", index=False)
    ns = np.arange(3, a.nsum + 1); Jk = np.sqrt(2 * ns + 1) * K / ns ** 2; Jh = np.array([abs(Jf[n]) if n <= nfile else Jk[n - 3] for n in ns])
    cells = pd.read_csv(v / "cells.csv"); Hfit_of = cells.set_index("tag").fit_H_m.to_dict(); hw = cells.set_index("tag")[["hloc_min_m", "hloc_max_m"]].to_dict("index"); rs = pd.read_csv(v / "truth_run_summary.csv").set_index("tag"); gcache = {}; rows = []
    def gmax(e, inc):
        key = (round(inc, 4), round(e, 3))
        if key not in gcache: gcache[key] = np.array([0.0, 0.0, 0.0] + [kl.gamma_max(n, key[1], key[0], nw=64, nq=4096 if n < 120 else 8192) for n in range(3, a.nsum + 1)])
        return gcache[key]
    for c in cells.itertuples():
        T = pd.read_csv(v / "truth" / c.tag / "passes_3b.csv"); t1 = T.iloc[0]
        status = "completed" if rs.loc[c.tag, "retcode"] == "Success" else "decayed"
        H = float(t1.h_local_m); Hf = float(Hfit_of[c.tag]); aa, e, rp, Q, rho, vr, lat = float(t1.a_m), float(t1.e), float(t1.radius_m), float(t1.q_w_cm2) * 1e4, float(t1.rho), float(t1.vrel_mps), float(t1.lat_rad)
        inc, om = np.radians(c.inc), np.radians(c.omega); p = aa * (1 - e * e); Vc = np.sqrt(mu / rp); N = np.sqrt(mu / aa ** 3); Tp = 2 * np.pi * np.sqrt(aa ** 3 / mu)
        fr = L.f_rot(rp, aa, vr, lat, SPIN, mu); PiE = L.pi_E(aa, e, rp, Q, H, beta, mu, fr); PiEf = L.pi_E(aa, e, rp, Q, Hf, beta, mu, fr); lev = 4 * aa * np.sqrt((1 - e) / (1 + e)); LH = lev / H
        fdv, fdd = lever_fd(aa, e, mu=mu); si = np.sin(inc)
        Pi2, Pi3 = J2 * (R / p) ** 2, abs(J3) * (R / p) ** 3; Pi3B = muS / mu * (aa / rS) ** 3
        G3 = 0.375 * (1 - e * e) * si * (5 * si ** 2 - 4) * np.cos(om); G4 = (15 / 32) * e * (1 - e * e) * si ** 2 * (7 * si ** 2 - 6) * np.sin(2 * om)
        Th = min(HORIZON_D, float(T.time_s.iloc[-1]) / 86400.0) * 86400.0
        _, wdot2, _ = kl.zonal_rates(2, e, inc, om, J2 * (R / p) ** 2, N); Teff = Th if abs(wdot2) * Th < 2 else 2.0 / abs(wdot2)
        Gm = gmax(e, inc)[3:]; pfac = (R / p) ** ns; Dh = aa * N * Jh * pfac * Gm * Teff / 1e3
        row = dict(tag=c.tag, inc_deg=c.inc, raan_deg=c.raan, omega_deg=c.omega, status=status, t_end_days=float(T.time_s.iloc[-1]) / 86400.0,
                   alt_p1_km=float(t1.alt_m) / 1e3, a_km=aa / 1e3, e=e, rp_km=rp / 1e3, T_h=Tp / 3600, Q_p1_w_cm2=Q / 1e4, vrel_p1_mps=vr, f_rot=fr,
                   Pi_q=L.pi_q(rp, Q, beta, mu), Pi_H=H / rp, Pi_H_fit=Hf / rp, Pi_E=PiE, Pi_E_fit=PiEf, dra_per_pass_fit_km=2 * L.delta_a_pass(aa, e, rp, Q, Hf, beta, mu, fr) / 1e3, L_H_fit=lev / Hf, L_H_c_fit=lev / Hf / np.sqrt(1 - e),
                   H_local_p1_km=H / 1e3, H_fit_km=Hf / 1e3, H_local_passes_min_km=float(T.h_local_m.min()) / 1e3, H_local_passes_max_km=float(T.h_local_m.max()) / 1e3, H_local_window_min_km=hw[c.tag]["hloc_min_m"] / 1e3, H_local_window_max_km=hw[c.tag]["hloc_max_m"] / 1e3, Pi_drag=0.5 * rho * rp * (1 + e) / beta, k_Q=L.k_Q(rp, Q, H, beta, mu, fr),
                   dra_per_pass_km=2 * L.delta_a_pass(aa, e, rp, Q, H, beta, mu, fr) / 1e3, lever_km_per_mps=lev / np.sqrt(mu / aa) / 1e3,
                   lever_fd_rel_err=abs(fdv / fdd / (lev / np.sqrt(mu / aa)) - 1) if False else abs((lever_fd(aa, e, mu=mu)[0] / lever_fd(aa, e, mu=mu)[1]) / (lev / np.sqrt(mu / aa)) - 1),
                   L_H=LH, L_H_c=LH / np.sqrt(1 - e), Pi_J2=Pi2, Pi_J3=Pi3, Pi_3B=Pi3B, Gamma_3B_J2=Pi3B / Pi2, Gamma_3=G3, Gamma_4=G4, Gamma_3_max_over_omega=0.375 * (1 - e * e) * si * abs(5 * si ** 2 - 4), Gamma_4_max_over_omega=(15 / 32) * e * (1 - e * e) * si ** 2 * abs(7 * si ** 2 - 6),
                   drift_J3_per_period_over_H=2 * np.pi * aa * Pi3 * abs(G3) / H, drift_3B_max_per_period_over_H=2 * np.pi * aa * 1.875 * e * np.sqrt(1 - e * e) * Pi3B / H,
                   band_width_km=float(t1.h_qlow_m - t1.h_qhigh_m) / 1e3, band_width_over_Hloc=float(t1.h_qlow_m - t1.h_qhigh_m) / H, w_Q=np.log(Q_HI / Q_LO), Qhard_over_Qhi=Q_HARD / Q_HI, Teff_days=Teff / 86400)
        for Nin in NKEY:
            s = ns > Nin; row[f"PiK_wc_N{Nin}_km"] = Dh[s].sum(); row[f"PiK_rss_N{Nin}_km"] = np.sqrt((Dh[s] ** 2).sum())
        for n in (2, 3, 4, 6, 10, 20): row[f"Pi_T_{n}"] = sig[n] * (R / rp) ** n
        # commensurability: nu at pass 1 (period of the first full pass interval from the truth) and at the last interval; m nu relative to integers
        tt = T.time_s.values; nu0 = SPIN * Tp / (2 * np.pi); nuE = SPIN * np.diff(tt)[-1] / (2 * np.pi) if len(tt) > 1 else np.nan
        row.update(nu_0=nu0, nu_end=nuE, **{f"delta_{m}": abs(m * nu0 - round(m * nu0)) for m in (2, 3, 4)}, nu_crosses_integer=any(np.floor(m * nu0) != np.floor(m * nuE) for m in (2, 3, 4)))
        rows.append(row)
    P = pd.DataFrame(rows); P.to_csv(out / "pi_cells.csv", index=False)
    cols = ["alt_p1_km", "a_km", "e", "rp_km", "T_h", "Pi_q", "Pi_H", "Pi_H_fit", "Pi_E", "Pi_E_fit", "dra_per_pass_fit_km", "L_H_fit", "L_H_c_fit", "H_local_p1_km", "H_fit_km", "H_local_passes_min_km", "H_local_passes_max_km", "H_local_window_min_km", "H_local_window_max_km", "Pi_drag", "k_Q", "dra_per_pass_km", "lever_km_per_mps", "L_H", "L_H_c", "Pi_J2", "Pi_J3", "Pi_3B", "Gamma_3B_J2", "Gamma_3", "Gamma_4", "Gamma_3_max_over_omega", "Gamma_4_max_over_omega",
            "drift_J3_per_period_over_H", "drift_3B_max_per_period_over_H", "PiK_rss_N3_km", "PiK_rss_N5_km", "PiK_wc_N3_km", "PiK_wc_N5_km", "band_width_km", "band_width_over_Hloc", "Pi_T_2", "Pi_T_3", "Pi_T_4", "Pi_T_6", "Pi_T_10", "Pi_T_20", "nu_0", "nu_end", "delta_2", "delta_3", "delta_4"]
    groups = [("all", P)] + [(f"inc {i}", P[P.inc_deg == i]) for i in sorted(P.inc_deg.unique())] + [(f"omega {w:g}", P[P.omega_deg == w]) for w in sorted(P.omega_deg.unique())] + [("completed", P[P.status == "completed"]), ("decayed", P[P.status == "decayed"])]
    G = pd.DataFrame([dict(group=g, n=len(d), quantity=q, min=d[q].min(), median=d[q].median(), max=d[q].max()) for g, d in groups if len(d) for q in cols]); G.to_csv(out / "pi_groups.csv", index=False)
    # hand check (venus_gram_plan.md 6.4): hand values from the plan, medians over cells within 25%
    hand = [("alt_p1_km", 135.15, 3.0, "abs"), ("a_km", 39419.0, 0.25, "rel"), ("e", 0.84305, 0.25, "rel"), ("T_h", 23.97, 0.25, "rel"), ("Pi_H", 6.2e-4, .25, "rel"), ("Pi_H_fit", 7.3e-4, .25, "rel"), ("Pi_q", 1.33e-3, .25, "rel"),
            ("Pi_E", 1.16e-3, .25, "rel"), ("dra_per_pass_km", 91.3, .25, "rel"), ("lever_km_per_mps", 16.03, .25, "rel"), ("L_H", 1.19e4, .25, "rel"), ("L_H_fit", 1.02e4, .25, "rel"),
            ("Pi_J2", 1.24e-6, .25, "rel"), ("Pi_J3", 3.15e-7, .25, "rel"), ("Pi_3B", 1.94e-5, .25, "rel"), ("Gamma_3B_J2", 15.6, .25, "rel"), ("nu_0", 4.1e-3, .25, "rel")]
    H_ = []
    for q, hv, tol, kind in hand:
        col = P[q].abs() if q == "dra_per_pass_km" else P[q]    # the plan states the apoapsis change per pass as a magnitude (86.9 km): sign is a decrease
        m = float(col.median()); err = abs(m - hv) if kind == "abs" else abs(m / hv - 1)
        H_.append(dict(quantity=q, hand=hv, median_over_cells=m, min=float(col.min()), max=float(col.max()), tolerance=tol, kind=kind, deviation=err, passed=bool(err <= tol)))
    pd.DataFrame(H_).to_csv(out / "hand_check.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print("Kaula K = %.4g, scatter %.3f dex (n=%d, free slope %.3f); lever finite-difference max rel err %.2e" % (K, scat, nused, slope, P.lever_fd_rel_err.max()))
    print(pd.DataFrame(H_).to_string(index=False)); print(G[G.group == "all"].to_string(index=False))
    print(G[G.quantity.isin(["Gamma_3", "Gamma_4", "Pi_3B", "PiK_rss_N5_km", "L_H", "Pi_q", "Pi_E", "Pi_H"])].pivot(index="quantity", columns="group", values="median").to_string())


if __name__ == "__main__":
    main()
