"""Venus with GRAM truth: the GRAM-specific criteria G1-G5 (venus_gram_plan.md section 6.1), the UNFROZEN diagnostics U1/U2 (fit-induced density error per pass,
Pi-rule site error), and the side-by-side comparison with the exponential-atmosphere Venus run (same metrics). Python only.
Usage: python3 -I venus_gram_diagnostics.py --results-root GRAM_ROOT [--single-root SINGLE_ROOT] [--exp-root EXP_ROOT] -> <GRAM_ROOT>/venus/gram/{criteria.csv, fit_error_by_bin.csv, site_error.csv, side_by_side.csv}
"""
import argparse, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
VC = tomllib.load(open(HERE / "venus_generality_config.toml", "rb"))
BINS = VC["analysis"]["lookahead_bin_edges_days"]; LAB = [f"{BINS[i]:g}-{BINS[i + 1]:g}" for i in range(len(BINS) - 1)]; DAY = 86400.0
Q_LO, Q_HI = VC["band"]["q_low_w_cm2"], VC["band"]["q_high_w_cm2"]; EXPECT_EPOCH = "2014-5-19 4:0:0.0"


def rms(x): x = np.asarray(x, float); return float(np.sqrt(np.mean(x ** 2))) if len(x) else float("nan")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--single-root"); ap.add_argument("--exp-root")
    a = ap.parse_args(); v = pathlib.Path(a.results_root) / "venus"; out = v / "gram"; out.mkdir(exist_ok=True)
    cells = pd.read_csv(v / "cells.csv"); rs = pd.read_csv(v / "truth_run_summary.csv"); mc = tomllib.load(open(v / "truth_mc.toml", "rb"))
    crit = []
    def C(i, stmt, value, thr, ok, note=""): crit.append(dict(id=i, statement=stmt, value=value, threshold=thr, passed=bool(ok), note=note))
    C("G1", "GRAM core epoch of every truth cell equals the run epoch", int((rs.gram_core_epoch == EXPECT_EPOCH).sum()), "32 of 32", (rs.gram_core_epoch == EXPECT_EPOCH).all(), f"epochs seen: {sorted(rs.gram_core_epoch.unique())}")
    C("G2", "route process and no failed cell", f"{mc['mc_route']}, failed {mc['mc_failed']}", "process, 0", mc["mc_route"] == "process" and mc["mc_failed"] == 0, f"workers {mc['mc_workers']}, density_family {mc.get('mc_density_family')}")
    cal = {t: pd.read_csv(v / "truth" / t / "calibration.csv") for t in cells.tag}
    conv = [abs(c.q_w_cm2.iloc[-1] / c.q_target.iloc[-1] - 1) <= 1e-3 and len(c) <= 8 for c in cal.values()]
    C("G3", "calibration converged (|Q/Q_target - 1| <= 1e-3) within 8 iterations", int(sum(conv)), "32 of 32", all(conv), "iterations per cell: " + str(pd.Series([len(c) for c in cal.values()]).value_counts().sort_index().to_dict()))
    q1 = pd.Series({t: c.q_w_cm2.iloc[0] for t, c in cal.items()}); inb = ((q1 >= Q_LO) & (q1 <= Q_HI))
    C("G4", "uncalibrated first pass at the Pi-rule altitude has Q in [0.14, 0.32] W/cm2", int(inb.sum()), ">= 29 of 32", inb.sum() >= 29, f"Q at the Pi-rule altitude: min {q1.min():.4f}, median {q1.median():.4f}, max {q1.max():.4f}; ratio to target median {(q1 / np.sqrt(Q_LO * Q_HI)).median():.3f}")
    C("G5", "primary (piecewise) fit: max |ln(fit/GRAM)| over the window", float(cells.pw_max_abs_lnrho.max()), "<= 0.05", cells.pw_max_abs_lnrho.max() <= 0.05, f"RMS median {cells.pw_rms_lnrho.median():.4f}; single fit RMS median {cells.fit_rms_lnrho.median():.3f}, max abs median {cells.fit_max_abs_lnrho.median():.3f}")
    pd.DataFrame(crit).to_csv(out / "criteria.csv", index=False)
    # Pi-rule site and altitude error (why G4 may miss): guess vs calibrated seed, guess site vs first-pass site
    rows = []
    for c in cells.itertuples():
        T = pd.read_csv(v / "truth" / c.tag / "passes_3b.csv"); s = cal[c.tag].seed_alt_m.iloc[-1]; t1 = T.iloc[0]
        dlon = (t1.lon_rad - c.site_lon_rad + np.pi) % (2 * np.pi) - np.pi
        rows.append(dict(tag=c.tag, h_guess_km=c.h_guess_m / 1e3, seed_calibrated_km=s / 1e3, seed_minus_guess_km=(s - c.h_guess_m) / 1e3, q_at_guess=q1[c.tag], dlat_deg=np.degrees(t1.lat_rad - c.site_lat_rad), dlon_deg=np.degrees(dlon),
                         dt_s=float(t1.time_s - c.site_t_s), site_lon_deg=np.degrees(c.site_lon_rad), site_lat_deg=np.degrees(c.site_lat_rad), inc=c.inc, omega=c.omega))
    SE = pd.DataFrame(rows); SE.to_csv(out / "site_error.csv", index=False)
    # U1/U2: fit-induced density error per truth pass (completed cells), by lookahead bin
    recs = []
    done = set(rs[rs.retcode == "Success"].tag)
    for t in cells.tag:
        T = pd.read_csv(v / "truth" / t / "passes_3b.csv"); T["tag"] = t; T["t_days"] = T.time_s / DAY; T["status"] = "completed" if t in done else "decayed"; recs.append(T)
    P = pd.concat(recs); P["bin"] = pd.cut(P.t_days, BINS, right=True, include_lowest=True, labels=LAB)
    for nm, col in (("piecewise", "rho_pw"), ("single", "rho_fit")):
        P[f"err_{nm}"] = P[col] / P.rho - 1; P[f"dh_{nm}_m"] = P.h_local_m * np.log(P[col] / P.rho)   # altitude shift that would equate fit and GRAM density (local scale height)
    fe = []
    for st, D in (("completed", P[P.status == "completed"]), ("all", P)):
        for b, g in D.groupby("bin", observed=False):
            if not len(g): continue
            fe.append(dict(set=st, bin=b, n=len(g), rms_rel_piecewise=rms(g.err_piecewise), median_rel_piecewise=float(g.err_piecewise.median()), rms_rel_single=rms(g.err_single), median_rel_single=float(g.err_single.median()),
                           rms_dh_piecewise_km=rms(g.dh_piecewise_m) / 1e3, rms_dh_single_km=rms(g.dh_single_m) / 1e3, h_local_median_km=float(g.h_local_m.median()) / 1e3, h_local_min_km=float(g.h_local_m.min()) / 1e3, h_local_max_km=float(g.h_local_m.max()) / 1e3))
    FE = pd.DataFrame(fe); FE.to_csv(out / "fit_error_by_bin.csv", index=False)
    p1 = P[P.pass_ == 1] if "pass_" in P else P[P["pass"] == 1]
    u2 = dict(pass1_rms_rel_piecewise=rms(p1.err_piecewise), pass1_rms_rel_single=rms(p1.err_single), later_rms_rel_piecewise_completed=rms(P[(P["pass"] > 1) & (P.status == "completed")].err_piecewise),
              later_rms_rel_single_completed=rms(P[(P["pass"] > 1) & (P.status == "completed")].err_single))
    pd.DataFrame([u2]).to_csv(out / "u2_pass1_vs_later.csv", index=False)
    # U3: relative peak-Q error of v5 (truth seed) by truth Q level, and fit-density error inside / outside the fit window (UNFROZEN; the exponential run's diagnostic 1, extended)
    hg = cells.set_index("tag").h_guess_m.to_dict(); fw = VC_FIT = tomllib.load(open(HERE / "venus_gram_config.toml", "rb"))["fit"]
    P["in_window"] = [(h >= hg[t] - fw["window_below_m"]) and (h <= hg[t] + fw["window_above_m"]) for t, h in zip(P.tag, P.alt_m)]
    qe = []
    for t in cells.tag:
        pr = pd.read_csv(v / "pred_v5" / f"{t}.csv"); tr = P[P.tag == t]; n = min(len(pr), len(tr))
        d = tr.iloc[:n].copy(); d["q_err_v5"] = pr.q_w_cm2.values[:n] / d.q_w_cm2.values - 1; qe.append(d)
    Q = pd.concat(qe); Qc = Q[Q.status == "completed"]; u3 = []
    for nm, m in (("Q in band [0.14, 0.32]", (Qc.q_w_cm2 >= Q_LO) & (Qc.q_w_cm2 <= Q_HI)), ("Q below 0.01", Qc.q_w_cm2 < 0.01), ("Q above 0.54 (hard)", Qc.q_w_cm2 > VC["band"]["q_hard_w_cm2"]), ("pass altitude inside fit window", Qc.in_window), ("pass altitude outside fit window", ~Qc.in_window)):
        g = Qc[m]; u3.append(dict(subset=nm, n=len(g), share=len(g) / len(Qc), rms_q_err_v5=rms(g.q_err_v5), median_q_err_v5=float(g.q_err_v5.median()) if len(g) else np.nan, rms_fit_density_err_piecewise=rms(g.err_piecewise), rms_fit_density_err_single=rms(g.err_single)))
    pd.DataFrame(u3).to_csv(out / "u3_by_q_level.csv", index=False)
    # side by side with the exponential Venus run
    sb = []
    def load(root):
        r = pathlib.Path(root) / "venus"; b = pd.read_csv(r / "analysis" / "bin_errors.csv"); b = b[(b.seed == "truth_seed") & (b.set == "completed")]
        ver = pd.read_csv(r / "analysis" / "verdict.csv"); apf = r / "apoapsis_pi" / "criteria.csv"; ap_ = pd.read_csv(apf).set_index("id") if apf.exists() else pd.DataFrame(index=[]); rsx = pd.read_csv(r / "truth_run_summary.csv")
        return b, ver, ap_, rsx
    roots = {"gram_piecewise": a.results_root, "gram_single": a.single_root, "exponential_H3.5km": a.exp_root}
    for nm, root in roots.items():
        if not root: continue
        b, ver, ap_, rsx = load(root)
        for est in ("v3", "v5"):
            for r in b[b.estimator == est].itertuples():
                sb.append(dict(world=nm, estimator=est, bin=r.bin, rms_alt_km=r.rms_alt_km, rms_q_rel=r.rms_q_rel, rms_dt_s=r.rms_dt_s))
        sb.append(dict(world=nm, estimator="-", bin="decayed cells / completed", rms_alt_km=int((rsx.retcode != "Success").sum()), rms_q_rel=int((rsx.retcode == "Success").sum())))
        if "D1V" in ap_.index:
            sb.append(dict(world=nm, estimator="-", bin="D1V median, D2V within 5% (H as primary of that run)", rms_alt_km=float(ap_.loc["D1V", "value"]), rms_q_rel=float(ap_.loc["D2V", "value"])))
    SB = pd.DataFrame(sb); SB.to_csv(out / "side_by_side.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 120)
    print(pd.DataFrame(crit).to_string(index=False)); print(SE.describe().round(4).to_string()); print(FE.round(4).to_string(index=False)); print(pd.DataFrame([u2]).to_string(index=False)); print(pd.DataFrame(u3).round(4).to_string(index=False))
    print(SB.pivot_table(index=["estimator", "bin"], columns="world", values="rms_q_rel", aggfunc="first").to_string())


if __name__ == "__main__":
    main()
