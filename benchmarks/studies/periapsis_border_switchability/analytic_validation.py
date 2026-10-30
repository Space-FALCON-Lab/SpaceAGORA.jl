"""WP-A3: offline validation of the analytic periapsis-altitude terms (plan analytic_framework_and_transfer_plan.md 3.1, 4; Revisions 1-5). Python only, no flights.

Usage: python3 -I analytic_validation.py [--terms T2,T3,T4a,T4b,T4c,T5,T6,T8,T9,VSUM,PI] [--out analytic_results.json]
Data (all on disk): Stage O dev flights (stageo_data/dev/runs), Stage K pair tables and Mars truth ablation tables (snapshot in analytic_data/snapshot), apoapsis-Pi dev pairs.
Pass/fail thresholds are the plan's ASSUMED accuracy targets (section 3.1), frozen before the first run of this script (they are the literals in THRESH).
Split: Stage O seeds 1-6 development, 7-12 held out; Stage K "dev"/"heldout" cell split. Fitted constants: none (everything is predicted), except where a line says "fitted".
Burn pairs (a nonzero commanded dv between the two passes) are excluded from every natural-drift test.
"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analytic_lib as L, analytic_pairs as PP, apoapsis_pi_lib as AP, stagek_lib as SK, tesseral_lib as TL

THRESH = dict(T2_pass_frac=0.95, T2_pass_abs_m=5.0, T2_pass_rel=0.05, T2_cum_m=50.0, T4a_EV=0.90, T4b_shrink=0.50, T4c_dk=2, T4c_jump_factor=2.0, T5_lo=0.97, T5_hi=1.03,
              T6_quad=0.01, T6_laplace=0.005, T8_EV=0.80, T8_mars_reduction=0.30, VSUM_EV=0.85)
EV = lambda y, p: float(1 - ((y - p) ** 2).sum() / (y ** 2).sum())
RES = {}


def verdict(ok): return "PASS" if ok else "FAIL"


def stageo_natural():
    df = PP.stageo_pairs(); n = df[~df.burn].copy(); n["held"] = n.seed > 6; return n


# ---------------------------------------------------------------- T2 ----------------------------------------------------------------------
def t2():
    n = stageo_natural(); T = THRESH; out = {}
    n["g_obsdw"] = [L.dh_geom(r.inc, r.om, r.o_dom, "mars", r.o_di) for r in n.itertuples()]
    n["g_pred"] = [L.dh_geom(r.inc, r.om, r.z_dom + r.t_dom, "mars", r.t_di) for r in n.itertuples()]
    n["g_pred_J2"] = [L.dh_geom(r.inc, r.om, r.z_dom, "mars") for r in n.itertuples()]
    n["g_first"] = [L.dh_geom(r.inc, r.om, r.o_dom, "mars", r.o_di, first_order=True) for r in n.itertuples()]
    ok = lambda err, truth: np.mean(err <= np.maximum(T["T2_pass_abs_m"], T["T2_pass_rel"] * truth.abs()))
    for nm in ("g_obsdw", "g_first", "g_pred", "g_pred_J2"):
        err = (n[nm] - n.o_dg).abs(); out[nm] = dict(frac_within=float(ok(err, n.o_dg)), median_err_m=float(err.median()), p95_err_m=float(err.quantile(.95)), max_err_m=float(err.max()))
    out["truth_rms_m"] = float((n.o_dg ** 2).mean() ** .5)
    g = n.groupby("run"); cum = pd.DataFrame(dict(truth=g.o_dg.sum(), obsdw=g.g_obsdw.sum(), pred=g.g_pred.sum(), predJ2=g.g_pred_J2.sum(), n=g.size()))
    for c in ("obsdw", "pred", "predJ2"): cum[f"err_{c}"] = (cum[c] - cum.truth).abs()
    out["cumulative_natural_m"] = {c: dict(median=float(cum[f"err_{c}"].median()), max=float(cum[f"err_{c}"].max())) for c in ("obsdw", "pred", "predJ2")}
    out["cumulative_truth_median_m"] = float(cum.truth.median())
    out["dom_ratio_by_e"] = {str(k): float(v) for k, v in (n.o_dom / (n.z_dom + n.t_dom)).groupby(pd.cut(n.e, [0.3, 0.5, 0.6, 0.7, 0.8]), observed=True).median().items()}
    res = n.o_dom - (n.z_dom + n.t_dom); lam = np.array([AP.lam_rp(r.e, -r.o_da / r.a) for r in n.itertuples()])
    out["dom_residual_corr_with_Lambda"] = float(np.corrcoef(res, lam)[0, 1]); out["dom_residual_corr_with_1_over_e"] = float(np.corrcoef(res, 1 / n.e)[0, 1])
    out["geodetic_offset_first_order_m"] = [float((g_ := n.g_first - n.o_dg).min()), float(g_.max())]
    out["verdict_form"] = verdict(out["g_obsdw"]["frac_within"] >= T["T2_pass_frac"])
    out["verdict_perpass_predicted_dom"] = verdict(out["g_pred"]["frac_within"] >= T["T2_pass_frac"])
    out["verdict_cumulative_predicted_dom"] = verdict(out["cumulative_natural_m"]["pred"]["max"] <= T["T2_cum_m"])
    RES["T2"] = out


# ---------------------------------------------------------------- T3 (existing; numerical-map cross-check) ---------------------------------
def t3():
    n = stageo_natural(); P = SK.planet("Mars"); s = n[n.k % 7 == 0].head(250).copy(); out = {}
    num = []
    for r in s.itertuples():
        a1, e1, i1, O1, w1, T = SK.zonal_map(P, r.a, r.e, r.inc, r.raan, r.om, nmax=20, rtol=1e-10, atol=1e-6); num.append(a1 * (1 - e1) - r.a * (1 - r.e))
    s["z_num"] = num
    out["n"] = int(len(s)); out["rms_numeric_m"] = float((s.z_num ** 2).mean() ** .5); out["rms_analytic_m"] = float((s.z_drp ** 2).mean() ** .5)
    out["rms_diff_m"] = float(((s.z_num - s.z_drp) ** 2).mean() ** .5); out["median_ratio_analytic_to_numeric"] = float((s.z_drp / s.z_num).median())
    out["corr"] = float(np.corrcoef(s.z_num, s.z_drp)[0, 1])
    RES["T3"] = out


# ---------------------------------------------------------------- T4a ---------------------------------------------------------------------
def t4a():
    n = stageo_natural(); out = {}; H = L.PLANETS["mars"]["H"]
    n["drag_rp"] = [L.drag_drop_rp_m(r.a, r.e, r.rp, r.q_eff * 1e4, H, r.beta, PP.MU_M, r.frot) for r in n.itertuples()]
    n["drag_rp_Qs"] = [L.drag_drop_rp_m(r.a, r.e, r.rp, 0.5 * r.rho_peak * r.vp ** 3, H, r.beta, PP.MU_M, r.frot) for r in n.itertuples()]
    cand = {"zonal": n.z_drp, "tesseral": n.t_drp, "zonal+tesseral": n.z_drp + n.t_drp, "zonal+tesseral+drag(q_eff)": n.z_drp + n.t_drp + n.drag_rp,
            "zonal+tesseral+drag(Q_s)": n.z_drp + n.t_drp + n.drag_rp_Qs}
    for nm, p in cand.items():
        out[nm] = {sp: EV(n.o_dr[m], p[m]) for sp, m in (("dev", ~n.held), ("held", n.held), ("all", n.held | ~n.held))}
    out["obs_dr_rms_m"] = float((n.o_dr ** 2).mean() ** .5); out["resid_rms_m"] = float(((n.o_dr - cand["zonal+tesseral+drag(q_eff)"]) ** 2).mean() ** .5)
    for key in ("a",):
        by = {}
        for lo, hi in ((3e6, 8e6), (8e6, 1.2e7), (1.2e7, 1.6e7), (1.6e7, 3e7)):
            m = (n.a > lo) & (n.a <= hi) & n.held
            if m.sum() > 30: by[f"{lo/1e6:.0f}-{hi/1e6:.0f}e6 m"] = dict(n=int(m.sum()), EV=EV(n.o_dr[m], cand["zonal+tesseral+drag(q_eff)"][m]))
        out["EV_by_a_bin_held"] = by
    kj = json.load(open(PP.SNAP / "stagek" / "k_tess_results.json")); out["stagek_mars_heldout_EV_Z_T_D"] = kj["EV_heldout"]["dr_ZTD"]
    out["verdict"] = verdict(out["zonal+tesseral+drag(q_eff)"]["held"] >= THRESH["T4a_EV"]); RES["T4a"] = out


# ---------------------------------------------------------------- T4b ---------------------------------------------------------------------
def hw(v): return float((np.percentile(v, 95) - np.percentile(v, 5)) / 2)


def t4b():
    df = pd.read_csv(PP.CACHE / "mars_truth_pairs.csv"); df["T"] = df.t_dra + df.t_drp; out = {"bins": {}}; ok_all = True
    df["keep"] = df.pred.abs() >= 5000; df["res"] = df.obs - df.pred       # pred: drag of pass k (apoapsis_pi_check convention, verified: regression coefficient 1.05 on pass k)
    for lo, hi in ((3.5e6, 6e6), (6e6, 9e6), (9e6, 1.2e7), (1.2e7, 2e7)):
        row = {}
        for arm in ("3b", "zonal"):
            s = df[(df.arm == arm) & df.keep & (df.a > lo) & (df.a <= hi)]
            row[arm] = dict(n=int(len(s)), raw_halfwidth_km=hw(s.res) / 1e3, corrected_halfwidth_km=hw(s.res - s["T"]) / 1e3, predicted_T_halfwidth_km=hw(s["T"]) / 1e3, corr_res_T=float(np.corrcoef(s.res, s["T"])[0, 1]))
        shrink = 1 - row["3b"]["corrected_halfwidth_km"] / row["3b"]["raw_halfwidth_km"]; row["3b"]["shrink"] = shrink; ok_all &= shrink >= THRESH["T4b_shrink"]
        fl, rw, cr = row["zonal"]["raw_halfwidth_km"], row["3b"]["raw_halfwidth_km"], row["3b"]["corrected_halfwidth_km"]       # floor = drag-model spread of the zonal arm (no tesseral)
        row["3b"]["floor_km"] = fl; row["3b"]["excess_variance_removed"] = float(1 - max(cr ** 2 - fl ** 2, 0.0) / max(rw ** 2 - fl ** 2, 1e-12)); row["3b"]["predicted_over_excess_halfwidth"] = row["3b"]["predicted_T_halfwidth_km"] / math.sqrt(max(rw ** 2 - fl ** 2, 1e-12))
        # tesseral-only observed part: excess of the 3b spread over the zonal-arm spread (same atmosphere model, same cells)
        s3 = df[(df.arm == "3b") & df.keep & (df.a > lo) & (df.a <= hi)]; row["Pi_T_a"] = float(hw(s3["T"]) / s3.pred.abs().median())   # plan 3.2: tesseral energy noise / per-pass drop (a feature, not an independent group)
        row["excess_halfwidth_km"] = math.sqrt(max(row["3b"]["raw_halfwidth_km"] ** 2 - row["zonal"]["raw_halfwidth_km"] ** 2, 0.0)); out["bins"][f"{lo/1e6:g}-{hi/1e6:g}e6 m"] = row
    s = df[(df.arm == "3b") & df.keep]; x = s["T"].values; y = s.res.values
    out["pooled_3b"] = dict(n=int(len(s)), corr=float(np.corrcoef(x, y)[0, 1]), slope=float((x * y).sum() / (x * x).sum()), EV=EV(y, x))
    jac = df[(df.arm == "3b")]; out["jacobi_apoapsis_arc"] = dict(median_ratio_dE_over_w_dhz=float((jac.t_dE_apo / (L.PLANETS["mars"]["w_rot"] * jac.t_dhz_apo)).median()),
                                                                  p05_p95=[float(v) for v in np.percentile(jac.t_dE_apo / (L.PLANETS["mars"]["w_rot"] * jac.t_dhz_apo), [5, 95])])
    x = (L.PLANETS["mars"]["w_rot"] * jac.t_dhz_apo).values; y = jac.t_dE_apo.values; out["jacobi_apoapsis_arc"].update(slope_dE_vs_w_dhz=float((x * y).sum() / (x * x).sum()), corr=float(np.corrcoef(x, y)[0, 1]))
    out["verdict"] = verdict(ok_all); out["note"] = "As pre-registered (corrected half-width <= 50% of raw) the target is not reachable: the drag-model spread of the zonal arm (no tesseral) is 70-95% of the raw 3b spread; see excess_variance_removed"; RES["T4b"] = out


# ---------------------------------------------------------------- T4c ---------------------------------------------------------------------
def t4c():
    """(A) nu formula vs the measured periapsis-to-periapsis rotation w_rot dt / 2 pi on the Stage O tables; (B) calendar: the pass at which m nu crosses an integer, predicted
    from the constant-Q decay law started at pass 1 of each Stage O run (T5 map, K fitted? no: theory K, oracle mean Q of the run) against the observed crossing pass;
    (C) the Fresnel jump |Z|/sqrt(m |nu'|) against the numerically swept sum of per-pass kicks of the order-m field alone (synthetic sweep on the Kepler/T5 trajectory)."""
    n = stageo_natural(); W = L.PLANETS["mars"]["w_rot"]; mu = PP.MU_M; out = {}
    nu_obs = W * n.dt / (2 * math.pi); J2 = SK.planet("Mars")["J"][2]; R0 = SK.planet("Mars")["R"]; nu_f0 = np.array([L.nu_commens(a, mu, W) for a in n.a])
    nu_f = np.array([L.nu_from_osc(r.a, r.e, r.om, r.a + r.o_da, r.e + r.o_de, r.om + r.o_dom, r.inc, mu, W, J2, R0) for r in n.itertuples()])
    out["A_nu_event_a_osculating_median_rel_err"] = float(((nu_f0 - nu_obs) / nu_obs).abs().median())
    out["A_nu_formula_median_rel_err_mean_a_midpoint"] = float(((nu_f - nu_obs) / nu_obs).abs().median()); out["A_nu_formula_p95_rel_err_mean_a_midpoint"] = float(((nu_f - nu_obs) / nu_obs).abs().quantile(.95))
    rows = []
    for run, g in PP.stageo_pairs().groupby("run"):
        g = g.sort_values("k"); nat = g[~g.burn]
        if len(g) < 60 or g.arm.iloc[0] != "FO": continue          # F-O: few large burns; use runs whose a(k) decays smoothly enough
        a0, e0, rp = g.a.iloc[0], g.e.iloc[0], g.rp.iloc[0]; beta = g.beta.iloc[0]; Qm = float(g.q_eff.mean()) * 1e4
        # predicted a(k) from the constant-Q law (T5 map) started at pass 1 with the run's mean heat rate
        nu_meas = W * g.dt.values / (2 * math.pi)
        # in-flight calibration of the decay law: ratio of observed to predicted 2 delta a over the first 15 natural pairs (constant afterwards)
        pr = [2 * AP.delta_a_pass(r.a, r.e, r.rp, r.q_eff * 1e4, L.PLANETS["mars"]["H"], r.beta, mu, r.frot) for r in nat.head(15).itertuples()]; cal_ratio = float(np.median(2 * nat.head(15).o_da.values / np.array(pr)))
        for variant, scale in (("theory", 1.0), ("calibrated15", cal_ratio)):
            ek, ak, wk = [e0], [a0], [g.om.iloc[0]]; inc0 = g.inc.iloc[0]
            for _ in range(len(g) - 1):
                da = scale * AP.delta_a_pass(ak[-1], ek[-1], rp, Qm, L.PLANETS["mars"]["H"], beta, mu, 1.0); ak.append(ak[-1] + da); ek.append(1 - rp / ak[-1])
                wk.append(wk[-1] + L.zonal_rates_per_pass(SK.planet("Mars"), ak[-2], ek[-2], inc0, wk[-1])[1])
            nu_pred = np.array([L.nu_from_osc(ak[i], ek[i], wk[i], ak[min(i + 1, len(ak) - 1)], ek[min(i + 1, len(ak) - 1)], wk[min(i + 1, len(ak) - 1)], inc0, mu, W, J2, R0) for i in range(len(ak))])
            for m in (1, 2, 3, 4):
                for j in range(1, 4 * m):
                    xp = m * nu_pred; xm = m * nu_meas
                    cp = np.where(np.diff(np.sign(xp - j)) != 0)[0]; cm = np.where(np.diff(np.sign(xm - j)) != 0)[0]
                    if len(cp) and len(cm): rows.append(dict(variant=variant, run=run, m=m, j=j, k_pred=float(cp[0] + (j - xp[cp[0]]) / (xp[cp[0] + 1] - xp[cp[0]])), k_obs=float(cm[0] + (j - xm[cm[0]]) / (xm[cm[0] + 1] - xm[cm[0]]))))
    cal_all = pd.DataFrame(rows); cal_all["dk"] = cal_all.k_pred - cal_all.k_obs
    # rolling look-ahead: at pass k_obs - Lh the keeper knows the past; decay law with the ratio calibrated on the last 10 passes, heat rate = mean of the last 10 q_eff
    roll = []
    for run, g in PP.stageo_pairs().groupby("run"):
        g = g.sort_values("k").reset_index(drop=True); nu_meas = W * g.dt.values / (2 * math.pi)
        for m in (1, 2, 3, 4):
            for j in range(1, 4 * m):
                xm = m * nu_meas; cm = np.where(np.diff(np.sign(xm - j)) != 0)[0]
                if not len(cm): continue
                k_obs = float(cm[0] + (j - xm[cm[0]]) / (xm[cm[0] + 1] - xm[cm[0]]))
                for Lh in (3, 5, 10, 20, 40):
                    k0 = int(k_obs) - Lh
                    if k0 < 12: continue
                    past = g.iloc[k0 - 10:k0]; pr = np.array([2 * AP.delta_a_pass(r.a, r.e, r.rp, r.q_eff * 1e4, L.PLANETS["mars"]["H"], r.beta, mu, r.frot) for r in past.itertuples()])
                    ratio = float(np.median(2 * past.o_da.values / pr)); ak = g.a.iloc[k0]; rp = g.rp.iloc[k0]; Qm = float(past.q_eff.mean()) * 1e4; beta = g.beta.iloc[k0]; seq = [ak]; wq = [g.om.iloc[k0]]; inc0 = g.inc.iloc[k0]
                    for _ in range(Lh + 60):
                        wq.append(wq[-1] + L.zonal_rates_per_pass(SK.planet("Mars"), ak, 1 - rp / ak, inc0, wq[-1])[1]); ak = ak + ratio * AP.delta_a_pass(ak, 1 - rp / ak, rp, Qm, L.PLANETS["mars"]["H"], beta, mu, 1.0); seq.append(ak)
                    es = [1 - rp / a for a in seq]
                    xp = m * np.array([L.nu_from_osc(seq[i], es[i], wq[i], seq[min(i + 1, len(seq) - 1)], es[min(i + 1, len(seq) - 1)], wq[min(i + 1, len(seq) - 1)], inc0, mu, W, J2, R0) for i in range(len(seq))]); cp = np.where(np.diff(np.sign(xp - j)) != 0)[0]
                    if len(cp): roll.append(dict(Lh=Lh, m=m, dk=k0 + cp[0] + (j - xp[cp[0]]) / (xp[cp[0] + 1] - xp[cp[0]]) - k_obs))
    roll = pd.DataFrame(roll)
    out["B_rolling_lookahead"] = {int(Lh): dict(n=int(len(g)), frac_within_2=float((g.dk.abs() <= 2).mean()), median_abs_dk=float(g.dk.abs().median()), p90_abs_dk=float(g.dk.abs().quantile(.9))) for Lh, g in roll.groupby("Lh")}
    out["B_calendar"] = {}
    for variant, cal in cal_all.groupby("variant"):
        out["B_calendar"][variant] = dict(n_crossings=int(len(cal)), median_abs_dk=float(cal.dk.abs().median()), frac_within_2=float((cal.dk.abs() <= THRESH["T4c_dk"]).mean()), p95_abs_dk=float(cal.dk.abs().quantile(.95)),
                                          bias_dk=float(cal.dk.median()), by_m={int(m): float((g.dk.abs() <= 2).mean()) for m, g in cal.groupby("m")})
    # (C) synthetic swept sum of order-m kicks against |Z| / sqrt(m |nu'|)
    cs = TL.coefficients("Mars", 20); jr = []
    rp = 3.5045e6; beta = 20.24; H = L.PLANETS["mars"]["H"]; inc, raan, om = 1.6322, 0.49, 1.9
    for m in (2, 3):
        Ct, St = TL.select(cs["C"], cs["S"], lambda nn, mm_: mm_ == m)
        for j in (1, 2) if m == 3 else (1,):
            a_c = (j / (m * W)) ** (2 / 3) * mu ** (1 / 3); e_c = 1 - rp / a_c
            if not (0.1 < e_c < 0.9): continue
            # sweep: a(k) from the T5 law at heat rate chosen so that nu' is moderate; kick sum over passes with theta advancing by 2 pi nu_k
            Qm = 0.2e4; ak = a_c; ks = np.arange(-120, 121); aseq = {}
            ak = a_c
            for k in range(0, 121): aseq[k] = ak; ak = ak + AP.delta_a_pass(ak, 1 - rp / ak, rp, Qm, H, beta, mu, 1.0)
            # backward sweep (upstream of the crossing): invert the map approximately by stepping with the (negative of the) local decrement
            ak = a_c
            for k in range(0, -121, -1): aseq[k] = ak; ak = ak - AP.delta_a_pass(ak, 1 - rp / ak, rp, Qm, H, beta, mu, 1.0)
            th = 0.0; kicks = []; nup = (L.nu_commens(aseq[1], mu, W) - L.nu_commens(aseq[-1], mu, W)) / 2
            thetas = {0: 0.3}
            for k in range(1, 121): thetas[k] = thetas[k - 1] + 2 * math.pi * L.nu_commens(0.5 * (aseq[k] + aseq[k - 1]), mu, W)
            for k in range(-1, -121, -1): thetas[k] = thetas[k + 1] - 2 * math.pi * L.nu_commens(0.5 * (aseq[k] + aseq[k + 1]), mu, W)
            for k in ks:
                a_ = aseq[int(k)]; kicks.append(L.tess_arc(cs, Ct, St, a_, 1 - rp / a_, inc, raan, om, thetas[int(k)] % (2 * math.pi), W, 64, f_lo=-math.pi, f_hi=math.pi)["rp"])
            cum = np.cumsum(kicks)
            # amplitude Z at exact resonance: kick at phases th0 and th0 + pi/(2 m)
            z1 = L.tess_arc(cs, Ct, St, a_c, e_c, inc, raan, om, 0.3, W, 64, f_lo=-math.pi, f_hi=math.pi)["rp"]
            z2 = L.tess_arc(cs, Ct, St, a_c, e_c, inc, raan, om, 0.3 + math.pi / (2 * m), W, 64, f_lo=-math.pi, f_hi=math.pi)["rp"]
            Zres = math.hypot(z1, z2); pred = Zres / math.sqrt(m * abs(nup))
            psi0 = math.atan2(-z2, z1); sgn = -1.0 if nup < 0 else 1.0                         # kick(theta) = R cos(m theta - delta); at the crossing theta = 0.3: psi0 = m 0.3 - delta
            pred_signed = Zres * math.cos(psi0 + sgn * math.pi / 4) / math.sqrt(m * abs(nup))
            # observed swept jump: peak-to-peak envelope of the cumulative sum about its mean ramp (Fresnel integral amplitude ~ jump)
            half = len(ks) // 2; before = cum[half - 60:half - 40].mean(); after = cum[half + 40:half + 60].mean()
            jr.append(dict(m=m, j=j, nu_dot_per_pass=float(nup), Z_res_m=Zres, jump_pred_m=pred, jump_swept_m=float(after - before), jump_pred_signed_m=float(pred_signed), ratio=float(abs(after - before) / pred), ratio_signed=float((after - before) / pred_signed) if abs(pred_signed) > 1e-9 else float("nan")))
    out["C_fresnel_vs_swept"] = jr
    out["verdict_calendar_theory"] = verdict(out["B_calendar"]["theory"]["frac_within_2"] >= 0.9); out["verdict_calendar_calibrated15"] = verdict(out["B_calendar"]["calibrated15"]["frac_within_2"] >= 0.9)
    out["verdict_fresnel_phase_aware"] = verdict(all(abs(r["jump_swept_m"] - r["jump_pred_signed_m"]) <= 0.2 * r["jump_pred_m"] for r in jr)) if jr else "N/A"
    out["verdict_calendar_rolling_L10"] = verdict(out["B_rolling_lookahead"].get(10, {}).get("frac_within_2", 0) >= 0.9)
    RES["T4c"] = out


# ---------------------------------------------------------------- T5 ---------------------------------------------------------------------
def t5():
    n = stageo_natural(); out = {}; H = L.PLANETS["mars"]["H"]
    for nm, q in (("q_mean", n.q_mean), ("q_eff", n.q_eff), ("Q_s_inertial", 0.5 * n.rho_peak * n.vp ** 3 / 1e4), ("Q_s_vrel", 0.5 * n.rho_peak * n.vrel ** 3 / 1e4)):
        pred = np.array([L.drag_drop_ra_m(r.a, r.e, r.rp, qq * 1e4, H, r.beta, PP.MU_M, r.frot) for r, qq in zip(n.itertuples(), q)]); ratio = 2 * n.o_da.values / pred
        row = {}
        for lo, hi in ((3e6, 1e7), (1e7, 1.3e7), (1.3e7, 1.6e7), (1.6e7, 3e7)):
            m = ((n.a > lo) & (n.a <= hi)).values
            if m.sum() > 30: row[f"{lo/1e6:g}-{hi/1e6:g}e6 m"] = float(np.median(ratio[m]))
        out[nm] = dict(median=float(np.median(ratio)), mean_obs_over_mean_pred=float(n.o_da.mul(2).sum() / pred.sum()), by_a_bin=row,
                       pass_in_all_bins=bool(all(THRESH["T5_lo"] <= v <= THRESH["T5_hi"] for v in row.values())))
    sig = float(np.log(n.q_eff / n.q_mean).std()); out["sigma_ln_qeff_over_qmean"] = sig; out["exp_sigma2_over_2"] = L.drag_density_factor(sig)
    dp = pd.read_csv(PP.SNAP / "apoapsis_pi" / "dev_pairs.csv"); m = dp[(dp.planet == "Mars") & dp.keep]
    out["meanfield_truth_Mars"] = {f: float(g.ratio.median()) for f, g in m.groupby("file")}
    out["meanfield_truth_Mars_by_c_bin"] = {str(k): float(v) for k, v in m[m.file == "passes_3b.csv"].groupby(pd.cut(m[m.file == "passes_3b.csv"].c, [0, 100, 300, 1000, 3000]), observed=True).ratio.median().items()}
    out["verdict_meanfield"] = verdict(THRESH["T5_lo"] <= out["meanfield_truth_Mars"]["passes_3b.csv"] <= THRESH["T5_hi"])
    out["verdict_perturbed_qeff"] = verdict(out["q_eff"]["pass_in_all_bins"]); out["verdict_perturbed_qmean"] = verdict(out["q_mean"]["pass_in_all_bins"]); RES["T5"] = out


# ---------------------------------------------------------------- T6 ---------------------------------------------------------------------
def t6():
    H = 6308.278108; mu = 4.2828e13; beta = 20.24; rho = 1e-9; out = {"grid": []}; bad_quad = bad_lap = 0; worst_quad = worst_lap = 0.0
    for e in (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.78):
        for c in (2, 5, 10, 25, 50, 100, 300, 1000):
            a = c * H / e
            if not (3.0e6 < a < 5e7): continue
            rp = a * (1 - e); N = 12 if e <= 0.3 else 40
            da, de, drp = L.bessel_drag_pass(a, e, H, rho, beta, mu, N=N); q = AP.gauss_pass_quad(a, e, H, rho, beta, mu)
            Q = 0.5 * rho * (mu * (1 + e) / rp) ** 1.5; lap_da = AP.delta_a_pass(a, e, rp, Q, H, beta, mu, 1.0); lap_rp = -AP.lam_rp(e, AP.pi_E(a, e, rp, Q, H, beta, mu, 1.0)) * H
            row = dict(e=e, c=c, a_km=a / 1e3, quad_err_da=abs(da / q[0] - 1), quad_err_de=abs(de / q[1] - 1), quad_err_drp=abs(drp / q[2] - 1), laplace_err_da=abs(da / lap_da - 1), laplace_err_drp=abs(drp / lap_rp - 1))
            out["grid"].append(row)
            if c >= 2: worst_quad = max(worst_quad, row["quad_err_da"], row["quad_err_de"]); bad_quad += max(row["quad_err_da"], row["quad_err_de"]) > THRESH["T6_quad"]
            if c >= 25: worst_lap = max(worst_lap, row["laplace_err_da"]); bad_lap += row["laplace_err_da"] > THRESH["T6_laplace"]
    g = pd.DataFrame(out["grid"]); out["worst_quad_err_da_de_c>=2"] = worst_quad; out["worst_laplace_err_da_c>=25"] = worst_lap
    out["laplace_err_da_by_c"] = {str(int(c)): float(v) for c, v in g.groupby("c").laplace_err_da.max().items()}
    out["laplace_err_drp_by_c"] = {str(int(c)): float(v) for c, v in g.groupby("c").laplace_err_drp.median().items()}
    out["quad_err_drp_note"] = "drp = (1-e) da - a de cancels at large e; relative error of drp is not gated (da and de are)"
    out["verdict_quadrature"] = verdict(bad_quad == 0); out["verdict_laplace"] = verdict(bad_lap == 0); RES["T6"] = out


# ---------------------------------------------------------------- T8 ---------------------------------------------------------------------
def t8():
    out = {}
    for pl in ("venus", "earth"):
        d = PP.tide_pairs(pl); t = PP.tess_pairs_other(pl); row = {}
        for sp in ("dev", "heldout"):
            m = d.split == sp
            row[sp] = dict(EV_zonal_drag=EV(d.dr[m], d.dr_M[m]), EV_plus_T8_direct=EV(d.dr[m], d.dr_M[m] + d.dr_T8[m]), EV_plus_T8_averaged=EV(d.dr[m], d.dr_M[m] + d.dr_T8avg[m]),
                           EV_plus_T8_direct_plus_tesseral=EV(d.dr[m], d.dr_M[m] + d.dr_T8[m] + t.dr_T[m]))
        row["tesseral_rms_m"] = float((t.dr_T ** 2).mean() ** .5); row["obs_dr_rms_m"] = float((d.dr ** 2).mean() ** .5); row["T8_rms_m"] = float((d.dr_T8 ** 2).mean() ** .5)
        row["averaged_over_direct_median"] = float((d.dr_T8avg / d.dr_T8).median()); out[pl] = row
    d = PP.tide_pairs("mars"); m = d.split == "heldout"; res0 = d.dr - d.dr_ZTD if "dr_ZTD" in d else None
    res0 = (d.dr - d.dr_ZTD)[m]; res1 = (d.dr - d.dr_ZTD - d.dr_T8)[m]
    out["mars"] = dict(heldout_resid_rms_ZTD_m=float((res0 ** 2).mean() ** .5), heldout_resid_rms_ZTD_T8_m=float((res1 ** 2).mean() ** .5), reduction_rms=float(1 - (res1 ** 2).mean() ** .5 / (res0 ** 2).mean() ** .5),
                       reduction_sumsq=float(1 - (res1 ** 2).sum() / (res0 ** 2).sum()), T8_rms_m=float((d.dr_T8 ** 2).mean() ** .5), EV_ZTD=EV(d.dr[m], d.dr_ZTD[m]), EV_ZTD_T8=EV(d.dr[m], (d.dr_ZTD + d.dr_T8)[m]))
    for pl in ("venus", "earth"): out[f"verdict_{pl}"] = verdict(out[pl]["heldout"]["EV_plus_T8_direct"] >= THRESH["T8_EV"]); out[f"verdict_{pl}_averaged"] = verdict(out[pl]["heldout"]["EV_plus_T8_averaged"] >= THRESH["T8_EV"])
    out["verdict_mars"] = verdict(out["mars"]["reduction_rms"] >= THRESH["T8_mars_reduction"]); RES["T8"] = out


# ---------------------------------------------------------------- T9 ---------------------------------------------------------------------
def t9():
    import tesseral_lib as TLm
    mu = PP.MU_M; area, mass = 2.2 * 1.7 + 2 * 1.945 * 1.7, 461.0; a_srp = L.srp_accel(area, mass, 1.3, 1.524 * L.AU); out = dict(a_srp=a_srp)
    for lab, ra in (("start", 28533e3), ("7408 km", 7408e3)):
        rp = 3501e3; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp); bound = L.srp_drp_bound(a, e, mu, a_srp)
        out[lab] = dict(a_km=a / 1e3, e=e, bound_m=bound)
    # exact numeric of the constant-force quadrature (own small routine so the force can be fixed in the orbit frame)
    def const_force_kick(a, e, F_pqw, nf=128, f_lo=0.0, f_hi=2 * math.pi):
        x, w = np.polynomial.legendre.leggauss(nf); half = 0.5 * (f_hi - f_lo); f = f_lo + half * (x + 1); w = half * w
        p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f))
        # unit vectors in the orbit plane: radial rhat = (cos f, sin f), transverse that = (-sin f, cos f); F in (P, Q) components
        aR = F_pqw[0] * np.cos(f) + F_pqw[1] * np.sin(f); aS = -F_pqw[0] * np.sin(f) + F_pqw[1] * np.cos(f)
        rates = TL.gauss_rates(a, e, 0.9, 0.0, f, rr, h, aR, aS, 0 * f, mu)    # in-plane force: aW = 0 (the in-plane bound is the maximum)
        return dict(zip(TL.KICK_NAMES, rates @ (w * rr ** 2 / h)))
    for lab, ra in (("start", 28533e3), ("7408 km", 7408e3)):
        rp = 3501e3; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp); mx = 0.0; mda = 0.0; mda_half = 0.0
        for psi in np.linspace(0, 2 * math.pi, 361):
            k = const_force_kick(a, e, a_srp * np.array([math.cos(psi), math.sin(psi)])); mx = max(mx, abs(k["rp"])); mda = max(mda, abs(k["a"]))
            kh = const_force_kick(a, e, a_srp * np.array([math.cos(psi), math.sin(psi)]), f_lo=-math.pi, f_hi=0.0); mda_half = max(mda_half, abs(kh["a"]))
        out[lab].update(max_drp_numeric_m=mx, max_abs_da_full_orbit_m=mda, max_abs_da_half_orbit_m=mda_half, bound_holds=bool(mx <= out[lab]["bound_m"] * 1.0001), bound_tightness=mx / out[lab]["bound_m"],
                        plan_energy_bound_m=4 * math.pi * a ** 3 * a_srp / mu)
    out["drag_da_start_m"] = float(-2 * AP.delta_a_pass(0.5 * (28533e3 + 3501e3), 0.78, 3501e3, 0.224e4, 6308.278, 20.24, mu, 1.0) / 2)
    out["verdict"] = verdict(all(out[l]["bound_holds"] for l in ("start", "7408 km"))); RES["T9"] = out


# ---------------------------------------------------------------- V-SUM -------------------------------------------------------------------
def stageo_tide():
    f = PP.CACHE / "stageo_tide.csv"
    if f.exists(): return pd.read_csv(f)
    import analytic_ephem as E
    from scipy.interpolate import CubicSpline
    df = PP.stageo_pairs(); et0 = E.utc_to_et(2001, 11, 6, 19, 0, 32); rows = []
    for r in df.itertuples():
        T = r.dt; ts = np.linspace(0, T, 5); pos = np.array([E.third_body_in_frame("mars", 10, et0 + r.t + x) for x in ts]).T; cs_ = CubicSpline(ts, pos, axis=1)
        rows.append(L.tide_arc_direct(PP.MU_M, E.GM[10], r.a, r.e, r.inc, r.raan, r.om, lambda tau: cs_(np.minimum(tau, T)), 32)["rp"])
    df["t8_rp"] = rows; df[["run", "k", "t8_rp"]].to_csv(f, index=False); return df[["run", "k", "t8_rp"]]


def vsum():
    n = stageo_natural().merge(stageo_tide(), on=["run", "k"]); H = L.PLANETS["mars"]["H"]; out = {}
    n["drag_rp"] = [L.drag_drop_rp_m(r.a, r.e, r.rp, r.q_eff * 1e4, H, r.beta, PP.MU_M, r.frot) for r in n.itertuples()]
    n["g_pred"] = [L.dh_geom(r.inc, r.om, r.z_dom + r.t_dom, "mars", r.t_di) for r in n.itertuples()]
    comp = dict(T2=n.g_pred, T3=n.z_drp, T4a=n.t_drp, T5=n.drag_rp, T8=n.t8_rp); pred = sum(comp.values())
    out["components_rms_m"] = {k: float((v ** 2).mean() ** .5) for k, v in comp.items()}; out["obs_dalt_rms_m"] = float((n.o_dalt ** 2).mean() ** .5)
    for sp, m in (("dev", ~n.held), ("held", n.held)): out[f"EV_{sp}"] = EV(n.o_dalt[m], pred[m])
    out["EV_held_without"] = {k: EV(n.o_dalt[n.held], (pred - v)[n.held]) for k, v in comp.items()}
    r = (n.o_dalt - n.t_drp)[n.held]; sm = (pred - comp["T4a"])[n.held]; out["after_removing_T4a"] = dict(resid_rms_m=float((r ** 2).mean() ** .5), EV_by_T2_T3_T5_T8=EV(r, sm), EV_by_T2=EV(r, comp["T2"][n.held]), EV_by_T8=EV(r, comp["T8"][n.held]), EV_by_T5=EV(r, comp["T5"][n.held]), EV_by_T3=EV(r, comp["T3"][n.held]))
    out["EV_held_by_arm"] = {a: EV(n.o_dalt[n.held & (n.arm == a)], pred[n.held & (n.arm == a)]) for a in n.arm.unique()}
    cumerr = (n.assign(p=pred).groupby("run").p.sum() - n.groupby("run").o_dalt.sum()).abs(); out["cumulative_error_m"] = dict(median=float(cumerr.median()), max=float(cumerr.max()))
    out["cumulative_truth_median_m"] = float(n.groupby("run").o_dalt.sum().median())
    out["verdict"] = verdict(out["EV_held"] >= THRESH["VSUM_EV"]); RES["VSUM"] = out


# ---------------------------------------------------------------- Pi groups ---------------------------------------------------------------
def pi():
    out = {}; mu = PP.MU_M; a = 0.5 * (28533e3 + 3501e3); e = (28533e3 - 3501e3) / (28533e3 + 3501e3); rp = 3501e3
    out["Pi_f"] = {pl: L.pi_f(pl) for pl in L.PLANETS}
    out["Mars_Odyssey_start"] = dict(Pi_omega=L.pi_omega(rp, mu, L.PLANETS["mars"]["w_rot"]), nu=L.nu_of_e(rp, e, mu, L.PLANETS["mars"]["w_rot"]), Pi_3B_Sun=L.pi_3b(1.32712440041e20, mu, a, 1.524 * L.AU),
                                      Pi_SRP=L.pi_srp(L.srp_accel(2.2 * 1.7 + 2 * 1.945 * 1.7, 461.0, 1.3, 1.524 * L.AU), a, mu), Pi_J2=1.958616e-3 * (3394200.0 / (a * (1 - e * e))) ** 2)
    out["Pi_f_note"] = "Pi_f sin^2(i) sin^2(omega) H is the geometry term of the periapsis height; its per-pass rate is Pi_f H sin^2(i) sin(2 omega) d omega"
    RES["PI"] = out


ALL = dict(T2=t2, T3=t3, T4a=t4a, T4b=t4b, T4c=t4c, T5=t5, T6=t6, T8=t8, T9=t9, VSUM=vsum, PI=pi)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--terms", default=",".join(ALL)); ap.add_argument("--out", default=str(PP.TASKS / "analytic_data" / "analytic_results.json")); a = ap.parse_args()
    for t in a.terms.split(","):
        print("running", t, flush=True); ALL[t]()
    p = pathlib.Path(a.out); old = json.loads(p.read_text()) if p.exists() else {}; old.update(RES); p.write_text(json.dumps(old, indent=1, default=float)); print(json.dumps(RES, indent=1, default=float))


if __name__ == "__main__":
    main()
