"""Stage J analysis (plan: tasks/2026-10-15_periapsis-region-border-targeting/stagej_plan.md, frozen before any statistic).
python3 -B stagej_analysis.py --phase selftest | run --f20 <root with stagef20/> --stagei <root with stagei/> [--out file.json]
Offline, Python only. Predictable (feature) versus white density error, controller reaction test, buffer implication."""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import density_buffer_analysis as dba
PLANETS = ["mars", "earth", "venus"]
WB, WA = 15000.0, 25000.0          # in-window rule (stagef20_config.toml)
LAMS = [1.0, 10.0, 100.0, 1000.0]; NB = 10000; SEED = 20261022
W_OLD, W_PAPER = float(np.log(0.32 / 0.14)), float(np.log(5.0))
FEAT = ["dh", "sinlat", "c1lon", "s1lon", "c2lon", "s2lon", "c1lst", "s1lst", "c2lst", "s2lst", "tday"]


def feats(P, hg):
    lon, lst = P.lon_rad.values, 2 * np.pi * P.lst_h.values / 24.0
    return np.column_stack([(P.alt_m.values - hg) / P.h_local_m.values, np.sin(P.lat_rad.values), np.cos(lon), np.sin(lon), np.cos(2 * lon), np.sin(2 * lon),
                            np.cos(lst), np.sin(lst), np.cos(2 * lst), np.sin(2 * lst), P.time_s.values / 86400.0 / 30.0])


def load(root, sub, pl, arm):
    d = pathlib.Path(root) / sub / pl; cells = pd.read_csv(d / "cells_base.csv"); out = []
    for t, sp in zip(cells.tag, cells.split):
        f = d / "runs" / f"{arm}_{t}" / "passes.csv"; g = d / "cal" / t / "guess.csv"
        if not f.exists() or not g.exists(): continue
        P = pd.read_csv(f); hg = float(pd.read_csv(g).h_guess_m[0])
        if len(P) == 0: continue
        P["tag"] = t; P["split"] = "dev" if sp == "dev" else "heldout"; P["hg"] = hg
        P["inw"] = (P.alt_m >= hg - WB) & (P.alt_m <= hg + WA)
        P["y"] = np.log(P.rho_eff / P.rho_pw); P["z"] = np.log(P.r_pert); P["ysys"] = P.y - P.z
        out.append(P)
    return pd.concat(out, ignore_index=True)


def design(P):
    X = np.vstack([feats(g, g.hg.iloc[0]) for _, g in P.groupby("tag", sort=False)]); return X


def ridge_fit(X, y, lam):
    mu, sd = X.mean(0), X.std(0) + 1e-12; Z = (X - mu) / sd; ym = y.mean()
    b = np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ (y - ym)); return dict(mu=mu, sd=sd, b=b, ym=ym)


def ridge_pred(m, X): return m["ym"] + ((X - m["mu"]) / m["sd"]) @ m["b"]


def choose_lam(X, y, tags):
    ut = sorted(set(tags)); fold = {t: i % 6 for i, t in enumerate(ut)}; f = np.array([fold[t] for t in tags]); best = None
    for lam in LAMS:
        se = 0.0
        for k in range(6):
            tr = f != k
            if tr.sum() < 20 or (~tr).sum() == 0: continue
            se += np.sum((y[~tr] - ridge_pred(ridge_fit(X[tr], y[tr], lam), X[~tr])) ** 2)
        if best is None or se < best[1]: best = (lam, se)
    return best[0]


def lag1_pairs(D, col):
    """D with tag, pass, inw rows sorted; consecutive in-window passes; per-cell mean removed (dba.lag1 convention). returns (a, b, cell index)"""
    a, b, c = [], [], []
    for t, G in D.groupby("tag", sort=False):
        v = G[col].values - G[col].mean(); ps = G["pass"].values
        for j in range(len(G) - 1):
            if ps[j + 1] - ps[j] == 1: a.append(v[j]); b.append(v[j + 1]); c.append(t)
    return np.array(a), np.array(b), np.array(c)


def corr(a, b): return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 else np.nan


def cboot(tags, fn, rng):
    """cluster bootstrap over unique tags; fn(list of tags with repeats) -> value."""
    ut = np.array(sorted(set(tags))); v = []
    for _ in range(NB):
        s = rng.choice(ut, len(ut), replace=True); r = fn(s)
        if np.isfinite(r): v.append(r)
    return [float(np.percentile(v, 5)), float(np.percentile(v, 95))] if v else [np.nan, np.nan]


def spread(e): return float(np.quantile(e, 0.95) - np.quantile(e, 0.05))


def ctrl_errors(P, pred):
    """one-pass-ahead errors of the controller predictor on consecutive in-window runs (Stage F build)."""
    out = []
    for t, T in P.groupby("tag", sort=False):
        T = T[T.inw]; rid = np.cumsum(np.r_[True, np.diff(T["pass"].values) != 1])
        for r in np.unique(rid):
            R = T[rid == r]
            if len(R) < 2: continue
            pr = dba.predict_run(R.y.values, np.unwrap(R.lon_rad.values))[pred]
            out += [dict(tag=t, e=R.y.values[j] - pr[j]) for j in range(1, len(R))]
    return pd.DataFrame(out)


def part_a(root, sub, rng, boot=True):
    data = {pl: load(root, sub, pl, "fixed_pert") for pl in PLANETS}
    for pl in PLANETS:
        P = data[pl]; P = P[P.inw].reset_index(drop=True); data[pl] = P
    Xd, yd, td, Xh, yh, th, mods = {}, {}, {}, {}, {}, {}, {}
    for pl, P in data.items():
        D, H = P[P.split == "dev"], P[P.split == "heldout"]
        Xd[pl], yd[pl], td[pl] = design(D), D.y.values, D.tag.values; Xh[pl], yh[pl], th[pl] = design(H), H.y.values, H.tag.values
    res = {}
    # pooled model B (planet intercepts via per-planet centering of y and shared Z features standardized per planet)
    for pl in PLANETS:
        P = data[pl]; H = P[P.split == "heldout"].reset_index(drop=True); D = P[P.split == "dev"]
        lam = choose_lam(Xd[pl], yd[pl], td[pl]); m = ridge_fit(Xd[pl], yd[pl], lam); yhat = ridge_pred(m, Xh[pl]); e = yh[pl] - yhat
        sst = np.sum((yh[pl] - yd[pl].mean()) ** 2); r2 = lambda yy, hh, yb: 1 - np.sum((yy - hh) ** 2) / np.sum((yy - yb) ** 2)
        H["e"] = e; H["yfeat"] = yhat
        a, b, c = lag1_pairs(H, "e"); a0, b0, _ = lag1_pairs(H, "y")
        ysys = H.ysys.values; r2s = 1 - np.sum((ysys - (yhat - yd[pl].mean() + D.ysys.mean())) ** 2) / np.sum((ysys - D.ysys.mean()) ** 2)
        # controller predictor from adapt_params
        ap = pd.read_csv(pathlib.Path(root) / sub / pl / "adapt_params.csv").iloc[0]; ce = ctrl_errors(P[P.split == "heldout"].assign(inw=True), ap.predictor)
        static = yh[pl] - yd[pl].mean()
        hl = float(np.median(H.h_local_m)); sp_e, sp_c, sp_s = spread(e), spread(ce.e.values), spread(static)
        out = dict(planet=pl, n_dev=len(D), n_held=len(H), n_held_cells=int(H.tag.nunique()), lam=lam, r2=float(r2(yh[pl], yhat, yd[pl].mean())),
                   r2_vs_ysys=float(r2s), ub_r2=float(np.var(ysys) / np.var(H.y.values)), var_y=float(np.var(H.y.values)), var_ysys=float(np.var(ysys)), var_z=float(np.var(H.z.values)),
                   lag1_e=corr(a, b), lag1_y=corr(a0, b0), lag1_ysys=corr(*lag1_pairs(H, "ysys")[:2]), q95_e=float(np.quantile(e, .95)), q05_e=float(np.quantile(e, .05)),
                   spread_e=sp_e, spread_ctrl=sp_c, spread_static=sp_s, ctrl_pred=ap.predictor, ctrl_ehi_minus_elo_dev=float(ap.e_hi - ap.e_lo), n_ctrl_scored=len(ce),
                   h_local_median_m=hl, ratio_ctrl=sp_e / sp_c, ratio_static=sp_e / sp_s, coef=dict(zip(FEAT, map(float, m["b"]))))
        if boot:
            Hs = H.reset_index(drop=True); idx = {t: np.where(Hs.tag.values == t)[0] for t in Hs.tag.unique()}; cidx = {t: np.where(ce.tag.values == t)[0] for t in ce.tag.unique()}
            yb = float(yd[pl].mean()); yv, hv, ev, tv = Hs.y.values, yhat, e, Hs.tag.values
            def pick(s, d): return np.concatenate([d[t] for t in s if t in d])
            def f_r2(s): i = pick(s, idx); return 1 - np.sum((yv[i] - hv[i]) ** 2) / np.sum((yv[i] - yb) ** 2)
            def f_ratio_c(s): i = pick(s, idx); j = pick(s, cidx); return spread(ev[i]) / spread(ce.e.values[j])
            def f_ratio_s(s): i = pick(s, idx); return spread(ev[i]) / spread(yv[i] - yb)
            def f_lag(s):
                g = Hs.set_index("tag", drop=False); aa, bb = [], []
                for t in s:
                    A, B, _ = lag1_pairs(Hs[Hs.tag == t], "e"); aa += list(A); bb += list(B)
                return corr(np.array(aa), np.array(bb))
            out.update(r2_ci=cboot(Hs.tag.values, f_r2, rng), ratio_ctrl_ci=cboot(Hs.tag.values, f_ratio_c, rng), ratio_static_ci=cboot(Hs.tag.values, f_ratio_s, rng))
            out["lag1_e_ci"] = cboot(Hs.tag.values, f_lag, rng) if NB <= 2000 else _lag_fast(Hs, rng)
        # buffer implication
        for nm, W in (("old", W_OLD), ("paper", W_PAPER)):
            out[f"band_{nm}_ln_left_features"] = W - sp_e; out[f"band_{nm}_ln_left_ctrl"] = W - sp_c; out[f"band_{nm}_corridor_m_features"] = (W - sp_e) * hl; out[f"band_{nm}_corridor_m_ctrl"] = (W - sp_c) * hl
        out["margin_m_features"] = sp_e * hl; out["margin_m_ctrl"] = sp_c * hl; out["fixed_1000m_in_H"] = 1000.0 / hl
        out["J1"] = bool(out["r2"] > 0.2 and out.get("r2_ci", [0, 0])[0] > 0); out["J2"] = bool(abs(out["lag1_e"]) < 0.1)
        out["J3_ctrl"] = bool(out["ratio_ctrl"] <= 0.8); out["J3_static"] = bool(out["ratio_static"] <= 0.8)
        res[pl] = out; mods[pl] = m
    # B pooled (planet intercepts) and C leave-one-planet-out
    def pooled(train_pl, target, recenter):
        X = np.vstack([Xd[p] for p in train_pl]); y = np.concatenate([yd[p] - yd[p].mean() for p in train_pl]); lam = choose_lam(X, y, np.concatenate([td[p] for p in train_pl]))
        m = ridge_fit(X, y, lam); return ridge_pred(m, Xh[target]) + yd[target].mean()
    for pl in PLANETS:
        for nm, tr in (("B_pooled", PLANETS), ("C_transfer", [p for p in PLANETS if p != pl])):
            yhat = pooled(tr, pl, True); e = yh[pl] - yhat
            res[pl][nm] = dict(r2=float(1 - np.sum(e ** 2) / np.sum((yh[pl] - yd[pl].mean()) ** 2)), spread_e=spread(e))
        res[pl]["J_gen"] = bool(res[pl]["C_transfer"]["r2"] > 0)
    return res, mods, data


def _lag_fast(Hs, rng):
    pairs = {t: lag1_pairs(Hs[Hs.tag == t], "e")[:2] for t in Hs.tag.unique()}; ut = np.array(sorted(pairs)); v = []
    for _ in range(NB):
        s = rng.choice(ut, len(ut)); a = np.concatenate([pairs[t][0] for t in s]); b = np.concatenate([pairs[t][1] for t in s]); r = corr(a, b)
        if np.isfinite(r): v.append(r)
    return [float(np.percentile(v, 5)), float(np.percentile(v, 95))]


def demean(df, cols, by):
    return df[cols] - df.groupby(by)[cols].transform("mean")


def part_b(root, sub, pl, arm, mod, rng):
    d = pathlib.Path(root) / sub / pl / "runs"; cells = pd.read_csv(d.parent / "cells_base.csv"); rows = []
    for t in cells.tag:
        fd, fp = d / f"{arm}_{t}" / "decisions.csv", d / f"{arm}_{t}" / "passes.csv"
        if not (fd.exists() and fp.exists()): continue
        D, P = pd.read_csv(fd), pd.read_csv(fp).set_index("pass"); hg = float(pd.read_csv(d.parent / "cal" / t / "guess.csv").h_guess_m[0])
        for _, r in D.iterrows():
            k = int(r.target_pass)
            if (k - 1) not in P.index or k not in P.index or not np.isfinite(r.yhat): continue
            x = P.loc[[k]].reset_index(); x["hg"] = hg
            rows.append(dict(tag=t, dv=r.delta_v_mps, yhat=r.yhat, target=(r.lower_m + r.upper_m) / 2, alt_pre=P.alt_m[k - 1], alt_post=P.alt_m[k], collapsed=bool(r.collapsed),
                             yfeat=float(ridge_pred(mod, feats(x, hg))[0]) if np.isfinite(x.alt_m[0]) else np.nan))
    R = pd.DataFrame(rows)
    if len(R) < 30: return dict(arm=arm, n=len(R), note="too few decisions")
    R["Tbar"] = R.groupby("tag").target.transform("mean"); R["R"] = R.target - R.Tbar; R["Pos"] = R.Tbar - R.alt_pre; R["dalt"] = R.alt_post - R.alt_pre
    Z = demean(R, ["dv", "yhat", "R", "Pos", "dalt", "yfeat"], R.tag)
    def slope(y, x): return float((x * y).sum() / (x * x).sum())
    def stats(idx_tags):
        Zs = pd.concat([Z[R.tag == t] for t in idx_tags]); s1 = slope(Zs.dv, Zs.yhat)
        A = np.column_stack([Zs.R, Zs.Pos]); b = np.linalg.lstsq(A, Zs.dv.values, rcond=None)[0]; v = A @ b
        sR = np.cov(v, b[0] * Zs.R)[0, 1] / np.var(v, ddof=1) if np.var(v) > 0 else np.nan
        return s1, b, sR, v, Zs
    s1, b, sR, v, Zs = stats(sorted(set(R.tag)))
    sl_ci = cboot(R.tag.values, lambda s: stats(s)[0], rng); sh_ci = cboot(R.tag.values, lambda s: stats(s)[2], rng)
    vnoR = Zs.Pos.values * b[1]; kappa = slope(Zs.dalt, Zs.dv)
    ok = R.yfeat.notna(); r2y = float(1 - np.sum((Z.yhat[ok] - slope(Z.yhat[ok], Z.yfeat[ok]) * Z.yfeat[ok]) ** 2) / np.sum(Z.yhat[ok] ** 2))
    return dict(arm=arm, n=len(R), n_cells=int(R.tag.nunique()), frac_collapsed=float(R.collapsed.mean()), slope_dv_on_yhat=s1, slope_ci=sl_ci, corr_dv_yhat=float(np.corrcoef(Zs.dv, Zs.yhat)[0, 1]),
                b_R=float(b[0]), b_Pos=float(b[1]), share_R=float(sR), share_R_ci=sh_ci, share_Pos=float(1 - sR), sum_abs_vhat=float(np.abs(v).sum()), sum_abs_vhat_noR=float(np.abs(vnoR).sum()),
                kappa_m_per_mps=kappa, yhat_R2_from_features=r2y, yhat_unexplained_frac=1 - r2y, mean_abs_dv=float(R.dv.abs().mean()), J4=bool(sl_ci[0] > 0 and sR > 0))


def drift(root, sub, pl, kappa):
    d = pathlib.Path(root) / sub / pl / "runs"; out = {}
    for arm in ("open_pert", "open_mean"):
        v = []
        for f in sorted(d.glob(f"{arm}_*/passes.csv")):
            P = pd.read_csv(f)
            if len(P) > 5: v.append((np.diff(P.alt_m.values).mean(), len(P)))
        if v: out[arm] = dict(n=len(v), mean_drift_m_per_pass=float(np.mean([a for a, _ in v])), median_abs_drift_m_per_pass=float(np.median([abs(a) for a, _ in v])), mean_passes=float(np.mean([n for _, n in v])),
                             dv_drift_per_run_mps=float(np.mean([abs(a) * n for a, n in v]) / kappa) if kappa and kappa > 0 else None)
    return out


def selftest():
    rng = np.random.default_rng(1); n = 4000; lon = rng.uniform(-np.pi, np.pi, n); z = rng.normal(0, 0.3, n)
    P = pd.DataFrame(dict(alt_m=1e5 + rng.normal(0, 100, n), h_local_m=6000.0, lat_rad=0.1, lon_rad=lon, lst_h=rng.uniform(0, 24, n), time_s=np.arange(n) * 3600.0))
    X = feats(P, 1e5); y = 0.5 * np.cos(lon) + z; tags = np.array([f"c{i // 100}" for i in range(n)]); dev = np.arange(n) < 2800
    m = ridge_fit(X[dev], y[dev], choose_lam(X[dev], y[dev], tags[dev])); e = y[~dev] - ridge_pred(m, X[~dev])
    r2 = 1 - np.sum(e ** 2) / np.sum((y[~dev] - y[dev].mean()) ** 2); assert 0.45 < r2 < 0.65, r2     # truth: 0.125/(0.125+0.09)=0.58
    assert abs(np.corrcoef(e[:-1], e[1:])[0, 1]) < 0.06 and spread(e) < spread(y[~dev] - y[dev].mean())
    print("selftest ok, r2 %.3f" % r2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", default="run"); ap.add_argument("--f20"); ap.add_argument("--stagei"); ap.add_argument("--out"); a = ap.parse_args()
    if a.phase == "selftest": selftest(); sys.exit()
    rng = np.random.default_rng(SEED); R = {}
    for name, root, sub in (("F20", a.f20, "stagef20"), ("I", a.stagei, "stagei")):
        res, mods, _ = part_a(root, sub, rng); R[name] = dict(a=res, b={}, drift={})
        arms = ["adaptive_pert", "adaptive_asym_pert", "adaptive_h025_pert"]
        for pl in PLANETS:
            R[name]["b"][pl] = {arm: part_b(root, sub, pl, arm, mods[pl], rng) for arm in arms}
        for pl in PLANETS:
            k = R[name]["b"][pl]["adaptive_pert"].get("kappa_m_per_mps"); R[name]["drift"][pl] = drift(a.f20, "stagef20", pl, k)
    json.dump(R, open(a.out, "w"), indent=1, default=float) if a.out else None
    print(json.dumps(R, indent=1, default=float))
