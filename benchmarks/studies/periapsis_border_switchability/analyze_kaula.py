"""Predictive test of the Kaula harmonic parameter against the pre-registered criteria (kaula_preregistration.md).
Usage: python3 -I analyze_kaula.py [--results-root DIR]. Needs kaula_calibrate.py and kaula_cells.py (both planets, both stages) first.
Writes <root>/kaula/{criteria.csv, cell_table.csv, rule_confusion.csv, tau_sweep_posthoc.csv, j4_sizing.csv, summary.txt, fig_*.png}.
"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import kaula_lib as kl

EPS = 1.0            # km, T5b of mars_generality_plan.md
FLOOR = 0.01         # km, short-period offset floor for the bound criteria


def auc(score, label):
    s1, s0 = score[label], score[~label]
    return float(mannwhitneyu(s1, s0).statistic / (len(s1) * len(s0))) if len(s1) and len(s0) else np.nan


def conf(flag, lab):
    tp, fp, fn, tn = int((flag & lab).sum()), int((flag & ~lab).sum()), int((~flag & lab).sum()), int((~flag & ~lab).sum())
    return dict(TP=tp, FP=fp, FN=fn, TN=tn, recall=tp / (tp + fn) if tp + fn else np.nan, fpr=fp / (fp + tn) if fp + tn else np.nan)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); k = root / "kaula"
    D = {}
    for pl in ("Earth", "Mars"):
        A = pd.read_csv(k / f"cells_{pl}_analytic.csv"); C = pd.read_csv(k / f"cells_{pl}_cowell.csv"); D[pl] = A.merge(C, on="tag")
    J = pd.read_csv(root / "mars" / "j3_exposure.csv").set_index("tag")
    D["Mars"]["rms_frozen_km"] = D["Mars"].tag.map(J.rms_frozen_km); D["Mars"]["rms_j3diag_km"] = D["Mars"].tag.map(J.rms_j3diag_km)
    tau = EPS / D["Earth"].D_3B_km.median()
    both = pd.concat(D.values(), ignore_index=True)
    for n in (3, 4):
        both[f"R1_{n}"] = both[f"rho_{n}"] > tau; both[f"R2_{n}"] = both[f"D_{n}_km"] > EPS; both[f"Lnum_{n}"] = both[f"E{n}_rms_km"] > EPS
    both["Ltruth_3"] = both.rms_frozen_km > EPS
    # POST HOC diagnostic (not pre-registered): D_3 evaluated with the elements at the end of the horizon instead of pass 1
    import kaula_cells as kc
    ends = []
    for pl in ("Earth", "Mars"):
        L_, P_, pc_, hist_, _ = kc.paths(root, pl); gdf_, R_, mu_, _ = kl.read_gravity(pl); J3_ = kl.zonal_J(gdf_)[3]
        for c in both[both.planet == pl].itertuples():
            T_ = pd.read_csv(hist_(c.tag, c.inc_deg), usecols=["time_s", "a_m", "e"]); T_ = T_[T_.time_s <= c.t_end_days * 86400 + 1].iloc[-1]
            ae, ee = float(T_.a_m), float(T_.e); pe = ae * (1 - ee * ee); Ne = np.sqrt(mu_ / ae ** 3)
            ends.append(ae * Ne * abs(J3_) * (R_ / pe) ** 3 * kl.gamma_max(3, ee, np.radians(c.inc_deg), nw=64) * c.Teff_days * 86400 / 1e3)
    both["D_3_end_km"] = ends
    both.to_csv(k / "cell_table.csv", index=False)
    E, M = both[both.planet == "Earth"], both[both.planet == "Mars"]
    rows, cr = [], []
    def add(i, stmt, val, thr, ok): cr.append(dict(id=i, statement=stmt, value=val, threshold=thr, verdict=ok))
    c1 = conf(M.R1_3, M.Ltruth_3); add("C1", "R1 vs L_truth(Mars J3), recall", f"{c1['recall']:.3f} (TP {c1['TP']}, FN {c1['FN']})", ">= 0.95", "PASS" if c1["recall"] >= 0.95 else "FAIL")
    crit = M[M.inc_deg.isin([63.4, 116.6])]; nfl = int(crit.R1_3.sum()); add("C2", "R1 flags at critical inclinations", f"{nfl} of {len(crit)}", "0", "PASS" if nfl == 0 else "FAIL")
    c3 = conf(E.R1_3, E.Lnum_3); c3b = conf(E.R1_4, E.Lnum_4)
    add("C3", "Earth R1 n=3 / n=4 vs L_num: needed cells; flagged (FP)", f"positives {int(E.Lnum_3.sum())}/{int(E.Lnum_4.sum())}; flagged {int(E.R1_3.sum())}/{int(E.R1_4.sum())}", "no missed need", "PASS" if c3["FN"] == 0 and c3b["FN"] == 0 else "FAIL")
    c4 = conf(both.R1_3, both.Lnum_3); add("C4", "R1 vs L_num(3) pooled, recall", f"{c4['recall']:.3f} (TP {c4['TP']}, FN {c4['FN']}, FP {c4['FP']})", ">= 0.95", "PASS" if c4["recall"] >= 0.95 else "FAIL")
    c5 = conf(both.R1_4, both.Lnum_4); add("C5", "R1 vs L_num(4) pooled, recall", "undefined (no positives)" if c5["TP"] + c5["FN"] == 0 else f"{c5['recall']:.3f} (TP {c5['TP']}, FN {c5['FN']})", ">= 0.95", "n/a" if c5["TP"] + c5["FN"] == 0 else ("PASS" if c5["recall"] >= 0.95 else "FAIL"))
    rho6 = spearmanr(M.Hrms_3_km, M.rms_frozen_km).correlation; add("C6", "Spearman Hrms_3 vs frozen RMS error, Mars", f"{rho6:.3f}", ">= 0.7 (assumed)", "PASS" if rho6 >= 0.7 else "FAIL")
    pos = M[M.Ltruth_3]; r7 = float(np.median(pos.Hrms_3_km / pos.rms_frozen_km)); add("C7", "median Hrms_3 / frozen RMS, needed cells", f"{r7:.3f} (q1-q3 {np.quantile(pos.Hrms_3_km / pos.rms_frozen_km, .25):.2f}-{np.quantile(pos.Hrms_3_km / pos.rms_frozen_km, .75):.2f})", "[0.5, 2] (assumed)", "PASS" if 0.5 <= r7 <= 2 else "FAIL")
    for n, cid in ((3, "C8a"), (4, "C8b")):
        s = both[both[f"E{n}_max_km"] > FLOOR]; f = float((s[f"D_{n}_km"] >= s[f"E{n}_max_km"]).mean()); add(cid, f"D_{n} >= Cowell E_{n},max, pooled (n cells {len(s)})", f"{f:.3f}", ">= 0.95", "PASS" if f >= 0.95 else "FAIL")
    s = both[both.Eall_max_km > FLOOR]; f9 = float((s["PiK_wc_N3_to12_km"] >= s.Eall_max_km).mean()); add("C9", f"Pi_K(3; n=4..12) >= Cowell E_all,max (n cells {len(s)})", f"{f9:.3f}", ">= 0.95", "PASS" if f9 >= 0.95 else "FAIL")
    r10 = spearmanr(s["PiK_wc_N3_to12_km"], s.Eall_max_km).correlation; add("C10", "Spearman Pi_K(3; 4..12) vs Cowell E_all,max, pooled", f"{r10:.3f}", ">= 0.7 (assumed)", "PASS" if r10 >= 0.7 else "FAIL")
    ns = both[~both.R2_4]; f11 = float((ns.E4_rms_km <= EPS).mean()); add("C11", "cells with D_4 <= eps have Cowell E_4,rms <= eps", f"{f11:.3f} of {len(ns)}", "1.0", "PASS" if f11 == 1.0 else "FAIL")
    CR = pd.DataFrame(cr); CR.to_csv(k / "criteria.csv", index=False)
    # confusion tables (reported)
    ct = []
    for pl, df in (("Earth", E), ("Mars", M)):
        for rule in ("R1", "R2"):
            for lab in ("Lnum_3", "Lnum_4") + (("Ltruth_3",) if pl == "Mars" else ()):
                n = 3 if lab.endswith("3") else 4; ct.append(dict(planet=pl, rule=f"{rule}_{n}", label=lab, **conf(df[f"{rule}_{n}"], df[lab])))
    CT = pd.DataFrame(ct); CT.to_csv(k / "rule_confusion.csv", index=False)
    byinc = M.groupby("inc_deg").apply(lambda g: pd.Series(dict(n=len(g), R1_flag=int(g.R1_3.sum()), R2_flag=int(g.R2_3.sum()), Ltruth=int(g.Ltruth_3.sum()), Lnum3=int(g.Lnum_3.sum()),
              rho3_med=g.rho_3.median(), D3_med=g.D_3_km.median(), Hrms3_med=g.Hrms_3_km.median(), E3rms_med=g.E3_rms_km.median(), frozen_med=g.rms_frozen_km.median(), D4_med=g.D_4_km.median(), E4rms_med=g.E4_rms_km.median(), PiK3_wc_med=g.PiK_wc_N3_km.median())))
    byincE = E.groupby("inc_deg").apply(lambda g: pd.Series(dict(n=len(g), R1_flag=int(g.R1_3.sum()), R2_flag=int(g.R2_3.sum()), rho3_med=g.rho_3.median(), D3_med=g.D_3_km.median(), Hrms3_med=g.Hrms_3_km.median(), E3rms_med=g.E3_rms_km.median(), est_med=g.est_rms_err_km.median(), D4_med=g.D_4_km.median(), E4rms_med=g.E4_rms_km.median(), PiK3_wc_med=g.PiK_wc_N3_km.median())))
    # extra rank statistics
    ex = {}
    for nm, x, y in (("Mars: Spearman D_3 vs frozen RMS", M.D_3_km, M.rms_frozen_km), ("Mars: Spearman rho_3 vs frozen RMS", M.rho_3, M.rms_frozen_km), ("Mars: Spearman Cowell E3_rms vs frozen RMS", M.E3_rms_km, M.rms_frozen_km),
                     ("Mars: Spearman Hrms_3 vs Cowell E3_rms", M.Hrms_3_km, M.E3_rms_km), ("Earth: Spearman Hrms_3 vs v2 RMS error", E.Hrms_3_km, E.est_rms_err_km), ("Earth: Spearman D_3 vs v2 RMS error", E.D_3_km, E.est_rms_err_km),
                     ("Earth: Spearman Hrms_3 vs Cowell E3_rms", E.Hrms_3_km, E.E3_rms_km), ("pooled: Spearman Pi_K(2) vs error (Earth v2, Mars frozen)", pd.concat([E.PiK_wc_N2_km, M.PiK_wc_N2_km]), pd.concat([E.est_rms_err_km, M.rms_frozen_km])),
                     ("pooled: Spearman Hrms_3 vs error (Earth v2, Mars frozen)", pd.concat([E.Hrms_3_km, M.Hrms_3_km]), pd.concat([E.est_rms_err_km, M.rms_frozen_km])),
                     ("Mars: Spearman Hrms_3 vs frozen RMS minus J3-diag RMS", M.Hrms_3_km, M.rms_frozen_km - M.rms_j3diag_km)):
        m = np.isfinite(x) & np.isfinite(y); ex[nm] = float(spearmanr(x[m], y[m]).correlation)
    ex["Mars AUC rho_3 vs L_truth"] = auc(M.rho_3.values, M.Ltruth_3.values); ex["Mars AUC D_3 vs L_truth"] = auc(M.D_3_km.values, M.Ltruth_3.values); ex["Mars AUC Hrms_3 vs L_truth"] = auc(M.Hrms_3_km.values, M.Ltruth_3.values)
    # median ratios of predictions to Cowell numerics (validation of the geometry)
    for n in (3, 4):
        for pl, df in (("Earth", E), ("Mars", M)):
            s = df[df[f"E{n}_max_km"] > FLOOR]; ex[f"{pl}: median D_{n} / Cowell E_{n},max (n={len(s)})"] = float((s[f"D_{n}_km"] / s[f"E{n}_max_km"]).median())
    for pl, df in (("Earth", E), ("Mars", M)):
        df = both[both.planet == pl]; q = df[df.D_3_km > 0.5]
        if len(q): ex[f"POST HOC {pl}: median Hmax_3 / D_3 (pass-1 elements), n={len(q)}"] = float((q.Hmax_3_km / q.D_3_km).median()); ex[f"POST HOC {pl}: median Hmax_3 / D_3_end (end elements)"] = float((q.Hmax_3_km / q.D_3_end_km).median()); ex[f"POST HOC {pl}: fraction Hmax_3 <= D_3_end"] = float((q.Hmax_3_km <= q.D_3_end_km).mean()); ex[f"POST HOC {pl}: fraction Hmax_3 <= D_3"] = float((q.Hmax_3_km <= q.D_3_km).mean())
    s_ = both[both.E3_max_km > 0.05]; ex["POST HOC C8a rescored with floor 0.05 km: fraction D_3 >= E3_max (n=%d)" % len(s_)] = float((s_.D_3_km >= s_.E3_max_km).mean())
    # post hoc tau sweep
    sw = []
    for t in np.logspace(-4, 0, 25):
        c = conf(M.rho_3 > t, M.Ltruth_3); sw.append(dict(tau=t, recall=c["recall"], fpr=c["fpr"], FN=c["FN"], FP=c["FP"], earth_flagged=int((E.rho_3 > t).sum())))
    pd.DataFrame(sw).to_csv(k / "tau_sweep_posthoc.csv", index=False)
    # mis-classified Mars cells under R1
    mis = M[(M.R1_3 != M.Ltruth_3)][["tag", "inc_deg", "omega_deg", "status", "rho_3", "D_3_km", "Hrms_3_km", "E3_rms_km", "rms_frozen_km", "rms_j3diag_km"]]
    # J4 / J2^2 sizing on representative cells (first cell of each inclination, omega = 0 epoch 1)
    sz = []
    for pl, df in (("Earth", E), ("Mars", M)):
        gdf, R, mu, _ = kl.read_gravity(pl); Jf = kl.zonal_J(gdf)
        for inc, g in df.groupby("inc_deg"):
            c = g.iloc[0]; e, a = c.e, c.a_m; p = a * (1 - e * e); N = np.sqrt(mu / a ** 3); i = np.radians(inc); row = dict(planet=pl, inc_deg=inc, e=e)
            w2 = kl.zonal_rates(2, e, i, 0.3, Jf[2] * (R / p) ** 2, N)[1]; row["wdot_J2_rad_per_day"] = w2 * 86400
            for n in (3, 4, 5):
                ed, wd, Od = kl.zonal_rates(n, e, i, 0.3, Jf[n] * (R / p) ** n, N); row[f"wdot{n}_over_wdotJ2"] = wd / w2 if abs(w2) > 0 else np.nan
            row["Pi_J2"] = abs(Jf[2]) * (R / p) ** 2; row["J2sq_phase_err_rad_60d_assumed_coeff1"] = abs(w2) * row["Pi_J2"] * 60 * 86400
            row["D_3_med_km"] = g.D_3_km.median(); row["D_4_med_km"] = g.D_4_km.median(); row["D_4_max_km"] = g.D_4_km.max(); row["E4_rms_med_km"] = g.E4_rms_km.median(); row["E4_rms_max_km"] = g.E4_rms_km.max(); row["Pi_K3_wc_med_km"] = g.PiK_wc_N3_km.median(); row["Pi_K3_wc_max_km"] = g.PiK_wc_N3_km.max()
            row["J2sq_effect_on_J3_km_bound"] = row["J2sq_phase_err_rad_60d_assumed_coeff1"] * g.D_3_km.median(); sz.append(row)
    pd.DataFrame(sz).to_csv(k / "j4_sizing.csv", index=False)
    with open(k / "summary.txt", "w") as f:
        pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
        print(f"tau (Earth, eps/median D_3B) = {tau:.4e}; Earth median D_3B = {D['Earth'].D_3B_km.median():.2f} km; Mars median D_3B = {D['Mars'].D_3B_km.median():.3f} km", file=f)
        print(CR.to_string(), file=f); print(CT.to_string(), file=f); print("\nMars by inclination\n", byinc.round(4).T.to_string(), file=f); print("\nEarth by inclination\n", byincE.round(4).T.to_string(), file=f)
        print("\nextra statistics", file=f); [print(f"  {a_}: {b:.4g}", file=f) for a_, b in ex.items()]
        print("\nR1 misclassified Mars cells vs L_truth\n", mis.round(4).to_string(), file=f); print("\nJ4 / J2^2 sizing\n", pd.DataFrame(sz).round(5).T.to_string(), file=f)
        print("\ntau sweep (post hoc)\n", pd.DataFrame(sw).round(4).to_string(), file=f)
    print(open(k / "summary.txt").read())
    # figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4)); col = {23.4: "#0b6e99", 25.2: "#c4612b", 50.0: "#6a8f2a", 63.4: "#777", 93.522: "#a23b8f", 116.6: "#444"}
    for inc, g in M.groupby("inc_deg"): ax[0].loglog(g.Hrms_3_km.clip(lower=1e-4), g.rms_frozen_km, ".", c=col[inc], label=f"Mars {inc}")
    ax[0].plot([1e-4, 20], [1e-4, 20], "k:"); ax[0].axhline(EPS, c="r", lw=.5); ax[0].set_xlabel("predicted J3 effect along truth history, RMS (km; clipped at 1e-4)"); ax[0].set_ylabel("Stage C frozen-estimator RMS error (km)"); ax[0].legend(fontsize=7); ax[0].set_title("Mars: omitted J3")
    for pl, df, m_ in (("Earth", E, "o"), ("Mars", M, "s")): ax[1].loglog(df.rho_3.clip(lower=1e-7), df.D_3_km.clip(lower=1e-4), m_, ms=3, alpha=.5, label=pl)
    ax[1].axvline(tau, c="r", lw=.8); ax[1].set_xlabel("rho_3 = rate ratio to third body"); ax[1].set_ylabel("D_3 worst case (km)"); ax[1].legend(); ax[1].set_title("R1 threshold (red)")
    for pl, df, m_ in (("Earth", both[both.planet == "Earth"], "o"), ("Mars", both[both.planet == "Mars"], "s")): ax[2].loglog(df.PiK_wc_N3_to12_km.clip(lower=1e-4), df.Eall_max_km.clip(lower=1e-4), m_, ms=3, alpha=.5, label=pl)
    ax[2].plot([1e-3, 5], [1e-3, 5], "k:"); ax[2].set_xlabel("Pi_K(3), n = 4..12, worst case (km)"); ax[2].set_ylabel("Cowell max |drp|, J2..J12 vs J2..J3 (km)"); ax[2].legend(); ax[2].set_title("Truncation bound")
    for x in ax: x.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(k / "fig_kaula_predictive.png", dpi=130)


if __name__ == "__main__":
    main()
