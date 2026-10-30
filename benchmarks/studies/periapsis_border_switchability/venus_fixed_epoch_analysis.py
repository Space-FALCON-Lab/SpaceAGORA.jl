"""Analysis half of venus_fixed_epoch.py: parameter table, closed-form and fitted predictors, out-of-sample tests, ablation, failure cases, planner. Offline; Python only.
Reads OUT/truth.csv and OUT/truth_paths.npz (written by `venus_fixed_epoch.py --truth`), writes OUT/analysis/*.csv, *.json, report.txt and figures."""
import json, math, pathlib, sys, time
import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V
import venus_fixed_epoch as E

GROUPS = {
    "sun_scale": ["Pi_3B0", "A3"],
    "sun_geom0": ["bP0", "bQ0", "bW0", "X3_0", "cos2beta0", "sin2psi0", "cos2psi0"],
    "sun_travel": ["dpsi", "bW_mid", "bW_end", "sin2psi_mid", "cos2psi_mid"],
    "sun_integrals": ["N3_cf", "R3_cf", "N_cf", "R_cf", "Rraw_cf", "maxabs_cf", "Sun_signchg", "floor_cf"],
    "j2": ["Pi_J2", "dom_J2", "Gamma_3B_J2"],
    "zonal_J3plus": ["X_zon0"],
    "tesseral": ["X_T0"],
    "drag_surface": ["Lam0", "Pi_H", "Pi_f"],
    "ratios": ["rho_zon", "rho_T", "rho_drag"],
    "elements": ["cosi", "sini", "cos_w", "sin_w", "cos2w", "sin2w", "cos_O", "sin_O"],
}
ALL = [c for g in GROUPS.values() for c in g]
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.append(s)


def metrics(df, pred, tag=""):
    """Floor-prediction metrics (pred in m/s)."""
    y = df.floor.values; p = np.asarray(pred, float); ly, lp = E.fy(y), E.fy(np.maximum(p, 0))
    out = dict(n=len(y), R2_log=E.r2(ly, lp), MAE_mps=E.mae(y, p), medAE_mps=float(np.median(np.abs(y - p))), spearman=float(stats.spearmanr(y, p)[0]))
    w = []
    for ep, g in df.assign(_p=p).groupby("epoch_off"):
        if len(g) > 5 and g.floor.std() > 0 and g._p.std() > 0: w.append(stats.spearmanr(g.floor, g._p)[0])
    out["spearman_within_epoch"] = float(np.mean(w)) if w else float("nan")
    for th in (0.1, 1.0):
        pos = y < th; pr = p < th
        out["recall_lt%g" % th] = float((pr & pos).sum() / max(pos.sum(), 1)); out["precision_lt%g" % th] = float((pr & pos).sum() / max(pr.sum(), 1)); out["n_true_lt%g" % th] = int(pos.sum())
    return out


def folds(df, kind, seed=0):
    """Yield (train_idx, test_idx) positions. kind: geom (5-fold by geometry key), epoch (leave one epoch out), inc (hold out inclinations 50, 85, 110)."""
    if kind == "geom":
        keys = df.inc.round(3).astype(str) + "_" + df.omega0.astype(str) + "_" + df.raan0.astype(str); u = keys.unique(); rng = np.random.default_rng(seed); f = dict(zip(u, rng.integers(0, 5, len(u)))); fid = keys.map(f).values
        for k in range(5): yield np.where(fid != k)[0], np.where(fid == k)[0]
    elif kind == "epoch":
        for ep in sorted(df.epoch_off.unique()): yield np.where(df.epoch_off.values != ep)[0], np.where(df.epoch_off.values == ep)[0]
    elif kind == "inc":
        te = df.inc.isin([50.0, 85.0, 110.0]).values; yield np.where(~te)[0], np.where(te)[0]


def cv_gbm(df, cols, target, kind, n_trees=100, seed=0):
    X = df[cols].values; y = target(df); pred = np.full(len(df), np.nan)
    for tr, te in folds(df, kind, seed):
        pred[te] = E.GBM(n_trees).fit(X[tr], y[tr]).predict(X[te])
    return pred


def cv_baseline_inc(df, kind):
    """Inclination-class lookup: mean log floor by inclination from the training rows (nearest training inclination for an unseen one)."""
    ly = E.fy(df.floor.values); pred = np.full(len(df), np.nan)
    for tr, te in folds(df, kind):
        m = pd.Series(ly[tr]).groupby(df.inc.values[tr]).mean(); inc = np.array(m.index)
        pred[te] = [m.values[np.argmin(np.abs(inc - i))] for i in df.inc.values[te]]
    return E.fy_inv(pred)


def main(out_dir):
    out = pathlib.Path(out_dir); an = out / "analysis"; an.mkdir(exist_ok=True); t0 = time.time()
    ctx, info = V.make_ctx(); nom = E.Nominal(ctx); H = nom.H[0]
    df = pd.read_csv(out / "truth.csv"); P = np.load(out / "truth_paths.npz", allow_pickle=False); assert (P["ids"] == df["id"].values.astype(str)).all(); arr = P["arr"]
    say("candidates", len(df), "sets", df.set.value_counts().to_dict(), "nominal: q_req %.4f W_eff %.4f n %d H %.0f m" % (nom.q_req, nom.W_eff, nom.n, H))
    ft = E.build_feature_table(ctx, nom, df); ft.to_csv(an / "parameters.csv"); say("parameter table built, %.0f s" % (time.time() - t0))
    D = pd.concat([df, ft], axis=1); D["net3_truth"] = D["net_third_body"]; D.to_csv(an / "truth_and_parameters.csv", index=False)
    gr = D[D.set == "grid"].reset_index(drop=True); rg = D[D.set == "ridge"].reset_index(drop=True); vx = D[D.set == "vex"].reset_index(drop=True)
    res = {}

    # ---- 0. how good is the closed-form path (nominal schedule, quadrupole) against the model's own path
    idx = D.index[D.set == "grid"]; sun_truth = arr[idx.values, 2, :]; sched = arr[idx.values, 3:6, :]
    # closed-form Sun path per row at its own epoch
    cfp = np.zeros((len(gr), nom.n))
    for ep, g in gr.groupby("epoch_off"):
        cf = E.closed_form(ctx, nom, ep, np.radians(g.inc.values), np.radians(g.omega0.values), np.radians(g.raan0.values)); cfp[g.index.values] = cf["path3"]
    tr = sun_truth[:, :nom.n]; err = cfp - tr; nlen = (gr.n_passes.values)
    say("passes per candidate: min %d max %d (nominal %d)" % (nlen.min(), nlen.max(), nom.n))
    ok = nlen == nom.n
    e_rms = np.sqrt(np.nanmean(err[ok] ** 2, 1)); amp = np.nanmax(np.abs(tr[ok]), 1)
    res["path_sun_vs_model"] = dict(frac_same_npasses=float(ok.mean()), rms_err_ln_q_median=float(np.median(e_rms)), rms_err_ln_q_p95=float(np.percentile(e_rms, 95)), amp_max_abs_ln_q_median=float(np.median(amp)),
                                    max_abs_err_ln_q=float(np.nanmax(np.abs(err[ok]))), end_err_ln_q_median=float(np.median(np.abs(err[ok][:, -1]))))
    say("closed-form Sun path vs model Sun path:", json.dumps(res["path_sun_vs_model"]))
    sched_dev = dict(q_req_rel_dev_max=float(np.abs(gr.q_req / nom.q_req - 1).max()), W_eff_abs_dev_max=float(np.abs(gr.W_eff - nom.W_eff).max()))
    res["nominal_schedule_vs_candidates"] = sched_dev; say("nominal schedule vs per-candidate:", json.dumps(sched_dev))

    # ---- 1. closed-form (zero-fit) predictors
    def reg(y, p, scale=1.0):
        return dict(R2=E.r2(y, p), pearson=float(np.corrcoef(y, p)[0, 1]), slope_through_origin=float((y * p).sum() / (p * p).sum()), MAE=E.mae(y, p) * scale, medAE=float(np.median(np.abs(y - p))) * scale)
    for nm, sub in (("grid", gr), ("ridge_scan", rg)):
        r = {}
        r["net_total_N_ln_q_vs_N_cf"] = reg(sub.N.values, sub.N_cf.values); r["net_total_N_km"] = dict(MAE_km=E.mae(sub.N.values, sub.N_cf.values) * H / 1e3, medAE_km=float(np.median(np.abs(sub.N - sub.N_cf))) * H / 1e3)
        r["net_sun_channel_vs_N3_cf"] = reg(sub.net3_truth.values, sub.N3_cf.values)
        r["range_R_vs_R_cf"] = reg(sub.R.values, sub.R_cf.values); r["range_km"] = dict(MAE_km=E.mae(sub.R.values, sub.R_cf.values) * H / 1e3, medAE_km=float(np.median(np.abs(sub.R - sub.R_cf))) * H / 1e3)
        r["floor_closed_form_DP"] = metrics(sub, sub.floor_cf.values)
        res["closed_form_" + nm] = r
    say("closed form on the grid:", json.dumps(res["closed_form_grid"], indent=1)); say("closed form on the ridge scan (omega0 every 10 deg, never in the grid):", json.dumps(res["closed_form_ridge_scan"], indent=1))
    by_ep = []
    for ep, g in gr.groupby("epoch_off"):
        by_ep.append(dict(epoch_off=ep, R2_N=E.r2(g.N.values, g.N_cf.values), R2_R=E.r2(g.R.values, g.R_cf.values), floor_MAE=E.mae(g.floor.values, g.floor_cf.values), floor_spearman=float(stats.spearmanr(g.floor, g.floor_cf)[0]), true_best=float(g.floor.min()), cf_best=float(g.floor_cf.min())))
    pd.DataFrame(by_ep).to_csv(an / "closed_form_by_epoch.csv", index=False); say(pd.DataFrame(by_ep).round(3).to_string())

    # ---- 2. fitted predictors (parameters -> floor, N, R), out-of-sample splits
    tg = {"floor": lambda d: E.fy(d.floor.values), "N": lambda d: d.N.values, "R": lambda d: d.R.values}
    fit_tab = []
    for kind in ("geom", "epoch", "inc"):
        pf = E.fy_inv(cv_gbm(gr, ALL, tg["floor"], kind)); te = ~np.isnan(pf)
        sub = gr[te]; m = metrics(sub, pf[te]); m.update(split=kind, model="GBM(all parameters)"); fit_tab.append(m)
        m = metrics(sub, sub.floor_cf.values); m.update(split=kind, model="closed form + DP (no fit)"); fit_tab.append(m)
        m = metrics(sub, cv_baseline_inc(gr, kind)[te]); m.update(split=kind, model="inclination-class mean"); fit_tab.append(m)
        m = metrics(sub, np.full(te.sum(), 0)); m.update(split=kind, model="constant 0"); fit_tab.append(m)
        for t in ("N", "R"):
            pt = cv_gbm(gr, ALL, tg[t], kind); y = gr[t].values[te]; cfc = gr[t + "_cf"].values[te]
            fit_tab.append(dict(split=kind, model="GBM(all parameters)", target=t, n=int(te.sum()), R2=E.r2(y, pt[te]), MAE_ln=E.mae(y, pt[te]), MAE_km=E.mae(y, pt[te]) * H / 1e3))
            fit_tab.append(dict(split=kind, model="closed form (no fit)", target=t, n=int(te.sum()), R2=E.r2(y, cfc), MAE_ln=E.mae(y, cfc), MAE_km=E.mae(y, cfc) * H / 1e3))
        say("split", kind, "done, %.0f s" % (time.time() - t0))
    ft_df = pd.DataFrame(fit_tab); ft_df.to_csv(an / "fitted_predictors_oos.csv", index=False); say(ft_df.round(3).to_string())

    # ---- 3. ablation (floor target, 5-fold by geometry; and leave-one-epoch-out for the key rows)
    ab = []
    base_p = E.fy_inv(cv_gbm(gr, ALL, tg["floor"], "geom")); b = metrics(gr, base_p); ab.append(dict(config="ALL", **b))
    for g, cols in GROUPS.items():
        rest = [c for c in ALL if c not in cols]
        p1 = E.fy_inv(cv_gbm(gr, rest, tg["floor"], "geom")); m1 = metrics(gr, p1); ab.append(dict(config="drop " + g, **m1))
        p2 = E.fy_inv(cv_gbm(gr, cols, tg["floor"], "geom")); m2 = metrics(gr, p2); ab.append(dict(config="only " + g, **m2))
        say("ablation", g, "drop R2_log %.3f MAE %.3f | only R2_log %.3f MAE %.3f" % (m1["R2_log"], m1["MAE_mps"], m2["R2_log"], m2["MAE_mps"]))
    for nm, cols in (("only start-state Sun geometry + scale (no campaign integral, no travel)", GROUPS["sun_scale"] + GROUPS["sun_geom0"]),
                     ("start-state geometry + travel (no closed-form integrals)", GROUPS["sun_scale"] + GROUPS["sun_geom0"] + GROUPS["sun_travel"]),
                     ("Sun only (scale + geometry + travel + integrals)", GROUPS["sun_scale"] + GROUPS["sun_geom0"] + GROUPS["sun_travel"] + GROUPS["sun_integrals"]),
                     ("non-Sun parameters only (j2, zonal, tesseral, drag, ratios, elements)", GROUPS["j2"] + GROUPS["zonal_J3plus"] + GROUPS["tesseral"] + GROUPS["drag_surface"] + GROUPS["ratios"] + GROUPS["elements"]),
                     ("floor_cf alone", ["floor_cf"]), ("R_cf alone", ["R_cf"]), ("R_cf + N_cf", ["R_cf", "N_cf"])):
        p = E.fy_inv(cv_gbm(gr, cols, tg["floor"], "geom")); m = metrics(gr, p); ab.append(dict(config=nm, **m)); say("ablation", nm, "R2_log %.3f MAE %.3f" % (m["R2_log"], m["MAE_mps"]))
        pe = E.fy_inv(cv_gbm(gr, cols, tg["floor"], "epoch")); me = metrics(gr, pe); ab[-1].update({"loeo_" + k: v for k, v in me.items() if k in ("R2_log", "MAE_mps", "spearman")})
    # the channel-share fact behind it: share of range^2 by channel in the model truth
    sh = {k: float(np.nanmean(gr["rng_" + k] ** 2) / np.nansum([np.nanmean(gr["rng_" + c] ** 2) for c in V.CH if "rng_" + c in gr])) for k in V.CH if "rng_" + k in gr}
    res["range_sq_share_by_channel_mean_over_grid"] = sh; say("share of mean range^2 by channel:", json.dumps({k: round(v, 6) for k, v in sh.items()}))
    pd.DataFrame(ab).to_csv(an / "ablation.csv", index=False)

    # ---- 4. failures: polar ridge, near-zero classification, worst misses
    fails = {}
    pol = gr[(gr.inc - 89.876).abs() < 0.01]; fails["polar_by_omega0"] = []
    for ow in sorted(pol.omega0.unique()):
        g = pol[pol.omega0 == ow]; fails["polar_by_omega0"].append(dict(omega0=ow, true_median=float(g.floor.median()), cf_median=float(g.floor_cf.median()), true_min=float(g.floor.min()), cf_min=float(g.floor_cf.min())))
    say(pd.DataFrame(fails["polar_by_omega0"]).round(2).to_string())
    onr = pol[(pol.omega0 % 90) == 0]; offr = pol[(pol.omega0 % 90) != 0]
    fails["ridge_pattern"] = dict(true_median_on_90deg_multiples=float(onr.floor.median()), true_median_other=float(offr.floor.median()), cf_median_on=float(onr.floor_cf.median()), cf_median_other=float(offr.floor_cf.median()))
    say("polar apse-line pattern:", json.dumps(fails["ridge_pattern"]))
    # ridge scan: omega0 every 10 deg
    rs = []
    for ow in sorted(rg.omega0.unique()):
        g = rg[rg.omega0 == ow]; rs.append(dict(omega0=ow, true_median=float(g.floor.median()), cf_median=float(g.floor_cf.median()), true_min=float(g.floor.min()), cf_min=float(g.floor_cf.min()), frac_true_lt1=float((g.floor < 1).mean()), frac_cf_lt1=float((g.floor_cf < 1).mean())))
    pd.DataFrame(rs).to_csv(an / "ridge_scan_by_omega0.csv", index=False); say(pd.DataFrame(rs).round(2).to_string())
    # GBM fit on the whole grid, tested on the ridge scan
    gb = E.GBM(100).fit(gr[ALL].values, E.fy(gr.floor.values)); pr = E.fy_inv(gb.predict(rg[ALL].values)); fails["gbm_on_ridge_scan"] = metrics(rg, pr); fails["cf_on_ridge_scan"] = metrics(rg, rg.floor_cf.values)
    say("ridge scan (unseen omega0, Omega0 phases): GBM", json.dumps(fails["gbm_on_ridge_scan"]), "\n closed form", json.dumps(fails["cf_on_ridge_scan"]))
    # confusion at the zero-floor boundary (closed form)
    conf = {}
    for th in (0.1, 0.5, 1.0):
        t = gr.floor < th; p = gr.floor_cf < th; conf["lt%g" % th] = dict(tp=int((t & p).sum()), fp=int((~t & p).sum()), fn=int((t & ~p).sum()), tn=int((~t & ~p).sum()))
    fails["confusion_closed_form_grid"] = conf; say("confusion (closed form vs truth):", json.dumps(conf))
    gr["abs_err_cf"] = (gr.floor - gr.floor_cf).abs(); worst = gr.sort_values("abs_err_cf", ascending=False).head(15)[["id", "epoch_off", "inc", "omega0", "raan0", "floor", "floor_cf", "R", "R_cf", "N", "N_cf"]]
    worst.to_csv(an / "closed_form_worst_misses.csv", index=False); say(worst.round(3).to_string())
    # near-boundary: error of R_cf relative to W_eff
    near = gr[(gr.R - nom.W_eff).abs() < 0.15]; fails["near_boundary_|R-W_eff|<0.15"] = dict(n=int(len(near)), floor_true_median=float(near.floor.median()), floor_cf_median=float(near.floor_cf.median()), R_rel_err_median=float(np.median(np.abs(near.R_cf - near.R) / near.R)))
    say("near the boundary R ~ W_eff:", json.dumps(fails["near_boundary_|R-W_eff|<0.15"]))
    # epoch/node sensitivity (truth and closed form) at polar, omega0 = 0 : floor vs Omega0, per epoch
    en = rg[(rg.omega0 == 0.0)].pivot_table(index="epoch_off", columns="raan0", values="floor"); enc = rg[(rg.omega0 == 0.0)].pivot_table(index="epoch_off", columns="raan0", values="floor_cf")
    say("polar omega0=0, truth floor (rows headline epochs, columns Omega0)\n" + en.round(2).to_string()); say("closed-form floor\n" + enc.round(2).to_string())
    res["failures"] = fails

    # ---- 5. planner
    plan = run_planner(ctx, nom, gr, vx, an, say)
    res["planner"] = plan
    (an / "results.json").write_text(json.dumps(res, indent=1, default=float)); (an / "report.txt").write_text("\n".join(LOG)); say("done %.0f s" % (time.time() - t0))
    figures(gr, rg, an)


def run_planner(ctx, nom, gr, vx, an, say):
    """For each fixed epoch pick (i, omega0, Omega0) with parameter-only predictors; look the truth up in the full-model grid. No access to the truth of the epoch being planned (the GBM is trained leave-one-epoch-out)."""
    eps = sorted(gr.epoch_off.unique()); rows = []; rng = np.random.default_rng(0)
    gbm_pred = np.full(len(gr), np.nan)
    for ep in eps:
        te = (gr.epoch_off == ep).values; gbm_pred[te] = E.fy_inv(E.GBM(100).fit(gr.loc[~te, ALL].values, E.fy(gr.loc[~te, "floor"].values)).predict(gr.loc[te, ALL].values))
    gr = gr.assign(floor_gbm=gbm_pred)
    for ep in eps:
        g = gr[gr.epoch_off == ep]; best = g.loc[g.floor.idxmin()]; v = vx[vx.epoch_off == ep].iloc[0]
        pb = g.sort_values(["floor_cf", "R_cf"]).iloc[0]; pg = g.sort_values(["floor_gbm", "R_cf"]).iloc[0]
        top10 = g.sort_values(["floor_cf", "R_cf"]).head(10); pv = top10.loc[top10.floor.idxmin()]
        rnd = rng.choice(g.floor.values, size=(2000, 10)); rnd1 = rng.choice(g.floor.values, 2000)
        rows.append(dict(epoch_off=ep, dp_best=best.floor, dp_best_geom="i%.1f/w%.0f/O%.0f" % (best.inc, best.omega0, best.raan0), plan_cf=pb.floor, plan_cf_pred=pb.floor_cf, plan_cf_geom="i%.1f/w%.0f/O%.0f" % (pb.inc, pb.omega0, pb.raan0),
                         plan_gbm=pg.floor, plan_gbm_pred=pg.floor_gbm, plan_gbm_geom="i%.1f/w%.0f/O%.0f" % (pg.inc, pg.omega0, pg.raan0), plan_cf_verify10=pv.floor, vex=v.floor, random_mean=g.floor.mean(), random_median=g.floor.median(),
                         random_best_of_10_mean=float(rnd.min(1).mean()), random_p_lt01=float((g.floor < 0.1).mean()), random1_p_lt01=float((rnd1 < 0.1).mean()), pool_n=len(g),
                         rank_pct_cf=float((g.floor < pb.floor).mean()), rank_pct_gbm=float((g.floor < pg.floor).mean()), n_pool_lt01=int((g.floor < 0.1).sum())))
    T = pd.DataFrame(rows); T.to_csv(an / "planner_grid_pool.csv", index=False); say(T.round(3).to_string())
    s = {}
    for nm in ("plan_cf", "plan_gbm", "plan_cf_verify10", "vex", "random_mean", "random_best_of_10_mean"):
        reg = T[nm] - T.dp_best; s[nm] = dict(mean_floor=float(T[nm].mean()), median_floor=float(T[nm].median()), mean_regret=float(reg.mean()), max_regret=float(reg.max()), frac_within_0p1=float((reg <= 0.1).mean()), frac_lt0p1=float((T[nm] < 0.1).mean()), frac_lt1=float((T[nm] < 1).mean()))
    s["dp_best"] = dict(mean_floor=float(T.dp_best.mean()), median_floor=float(T.dp_best.median()), frac_lt0p1=float((T.dp_best < 0.1).mean()), frac_lt1=float((T.dp_best < 1).mean()), n_epochs=len(T))
    say("planner summary (grid pool, 12 epochs):", json.dumps(s, indent=1))
    # fine pool, closed form only: i 20..160 step 5, omega0 and Omega0 every 5 deg; verify the best 10 picks per epoch with the full model
    fine = []
    incs = np.arange(20.0, 160.1, 5.0); oms = np.arange(0.0, 360.0, 5.0); Oms = np.arange(0.0, 360.0, 5.0)
    II, WW, OO = np.meshgrid(incs, oms, Oms, indexing="ij"); II, WW, OO = II.ravel(), WW.ravel(), OO.ravel()
    picks = []
    for ep in eps:
        Rcf = np.empty(len(II)); N_ = np.empty(len(II))
        for s0 in range(0, len(II), 20000):
            sl = slice(s0, s0 + 20000); F, _ = E.features(ctx, nom, ep, np.radians(II[sl]), np.radians(WW[sl]), np.radians(OO[sl])); Rcf[sl] = F["R_cf"]; N_[sl] = F["N_cf"]
        cand = np.argsort(Rcf)[:3000]
        F, _ = E.features(ctx, nom, ep, np.radians(II[cand]), np.radians(WW[cand]), np.radians(OO[cand])); fl = np.array([E.dp_floor(s, nom.w, nom.W_eff) for s in F["_slow"]])
        order = np.lexsort((F["R_cf"], fl))[:10]
        for rk, j in enumerate(order): picks.append(dict(epoch_off=ep, rank=rk, inc=II[cand[j]], omega0=WW[cand[j]], raan0=OO[cand[j]], floor_cf=fl[j], R_cf=F["R_cf"][j]))
        say("fine pool epoch %g: %d candidates, closed-form best floors %s" % (ep, len(II), np.round(fl[order[:3]], 3)))
    pk = pd.DataFrame(picks); pk["id"] = ["P%05d" % i for i in range(len(pk))]; pk["set"] = "fine_pick"
    import multiprocessing as mp
    E._W["ctx"] = ctx; cl = [dict(id=r.id, set="fine_pick", inc=float(r.inc), epoch_off=float(r.epoch_off), omega=float(r.omega0), raan=float(r.raan0)) for r in pk.itertuples()]
    with mp.get_context("fork").Pool(6) as pool: rr = pool.map(E.eval_truth, cl, chunksize=4)
    pk["floor_true"] = [r["floor"] for r, _ in rr]; pk["R_true"] = [r["R"] for r, _ in rr]; pk.to_csv(an / "planner_fine_pool_picks.csv", index=False)
    t0p = pk[pk["rank"] == 0].set_index("epoch_off"); b10 = pk.groupby("epoch_off").floor_true.min()
    FT = pd.DataFrame(dict(epoch_off=eps, dp_best_grid=T.dp_best.values, fine_pick_top1_true=[t0p.floor_true[e] for e in eps], fine_pick_top1_pred=[t0p.floor_cf[e] for e in eps], fine_pick_best_of_10_true=[b10[e] for e in eps],
                           fine_top1_geom=["i%.0f/w%.0f/O%.0f" % (t0p.inc[e], t0p.omega0[e], t0p.raan0[e]) for e in eps]))
    FT.to_csv(an / "planner_fine_pool_summary.csv", index=False); say(FT.round(3).to_string())
    s["fine_pool"] = dict(top1_mean_floor=float(FT.fine_pick_top1_true.mean()), top1_frac_lt0p1=float((FT.fine_pick_top1_true < 0.1).mean()), top1_mean_regret_vs_grid_best=float((FT.fine_pick_top1_true - FT.dp_best_grid).mean()),
                          best10_mean_floor=float(FT.fine_pick_best_of_10_true.mean()), best10_frac_lt0p1=float((FT.fine_pick_best_of_10_true < 0.1).mean()), top1_pred_vs_true_MAE=float((FT.fine_pick_top1_true - FT.fine_pick_top1_pred).abs().mean()))
    say(json.dumps(s["fine_pool"], indent=1)); return s


def figures(gr, rg, an):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].scatter(gr.N_cf, gr.N, s=3, alpha=.3); lim = [gr.N.min(), gr.N.max()]; ax[0].plot(lim, lim, "k-", lw=.8); ax[0].set_xlabel("closed-form net drift N (ln q)"); ax[0].set_ylabel("ranking model N (ln q)")
    ax[1].scatter(gr.R_cf, gr.R, s=3, alpha=.3); lim = [0, gr.R.max()]; ax[1].plot(lim, lim, "k-", lw=.8); ax[1].set_xlabel("closed-form range R (ln q)"); ax[1].set_ylabel("ranking model R (ln q)")
    ax[2].scatter(gr.floor_cf + .01, gr.floor + .01, s=3, alpha=.3); ax[2].set_xscale("log"); ax[2].set_yscale("log"); ax[2].plot([.01, 30], [.01, 30], "k-", lw=.8); ax[2].set_xlabel("floor from closed-form path (m/s, +0.01)"); ax[2].set_ylabel("DP floor, ranking model (m/s, +0.01)")
    fig.tight_layout(); fig.savefig(an / "fig_closed_form_vs_model.png", dpi=110); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); s = rg.groupby("omega0")[["floor", "floor_cf"]].median(); ax.plot(s.index, s.floor, "o-", label="model (median over Omega0, 4 epochs)"); ax.plot(s.index, s.floor_cf, "s--", label="closed form"); ax.set_xlabel("omega0 (deg), i = 89.876"); ax.set_ylabel("floor (m/s)"); ax.legend(); fig.tight_layout(); fig.savefig(an / "fig_polar_apse_line.png", dpi=110); plt.close(fig)
