"""Stage J UNFROZEN supplement (written AFTER seeing stagej_analysis.py results; labeled UNFROZEN, changes no verdict).
Uses the logged perturbation z = ln r_pert to split the controller's predictor yhat exactly (linear predictors) into a mean-field part and a white-perturbation part,
and the white-only buffer spread(z). python3 -B stagej_unfrozen.py --f20 ROOT --stagei ROOT [--out f.json]"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import density_buffer_analysis as dba, stagej_analysis as sj
sj.NB = 4000; MU = dict(mars=4.282837e13, earth=3.986004e14, venus=3.24859e14)   # ASSUMED standard gravitational parameters (m3/s2), not read from the repo
PLANETS = sj.PLANETS


def spread_z(root, sub, pl, rng):
    P = sj.load(root, sub, pl, "fixed_pert"); P = P[P.inw]; H = P[P.split == "heldout"]; D = P[P.split == "dev"]
    ap = pd.read_csv(pathlib.Path(root) / sub / pl / "adapt_params.csv").iloc[0]; hl = float(H.h_local_m.median())
    cm = H.groupby("tag").ysys.mean(); vb = float(np.var(cm.reindex(H.tag).values)); vt = float(np.var(H.ysys)); sz = sj.spread(H.z.values); sc = float(ap.e_hi - ap.e_lo)
    out = dict(spread_z=sz, std_z=float(H.z.std()), var_ysys_total=vt, var_ysys_between_cells=vb, frac_ysys_between_cells=vb / vt, spread_ctrl_dev=sc, ratio_z_over_ctrl=sz / sc, h_local_m=hl, margin_m_z=sz * hl, margin_m_ctrl=sc * hl)
    for nm, W in (("old", sj.W_OLD), ("paper", sj.W_PAPER)): out[f"ln_left_z_{nm}"] = W - sz; out[f"corridor_m_z_{nm}"] = (W - sz) * hl; out[f"ln_left_ctrl_{nm}"] = W - sc
    # predictor noise injection: apply controller predictor to z alone (held-out, consecutive runs)
    e = []
    for t, T in H.groupby("tag", sort=False):
        pr = dba.predict_run(T.z.values, np.unwrap(T.lon_rad.values))[ap.predictor]; e += list(pr[1:])
    out["std_predictor_noise_ln"] = float(np.nanstd(e)); out["predictor"] = ap.predictor; out["predictor_noise_m"] = out["std_predictor_noise_ln"] * hl
    return out


def decompose(root, sub, pl, arm, rng):
    d = pathlib.Path(root) / sub / pl / "runs"; cells = pd.read_csv(d.parent / "cells_base.csv"); ap = pd.read_csv(d.parent / "adapt_params.csv").iloc[0]; rows = []; chk = []
    for t in cells.tag:
        fd, fm, fp = d / f"{arm}_{t}" / "decisions.csv", d / f"{arm}_{t}" / "measured.csv", d / f"{arm}_{t}" / "passes.csv"
        if not (fd.exists() and fm.exists() and fp.exists()): continue
        D, M, P = pd.read_csv(fd), pd.read_csv(fm), pd.read_csv(fp)
        y, lon = M.y.values, np.unwrap(M.lon_rad.values); z = np.log(M.r.values)
        py, pz = dba.predict_run(y, lon)[ap.predictor], dba.predict_run(z, lon)[ap.predictor]
        for _, r in D.iterrows():
            n = int(r.n_meas)
            if n < 2 or n >= len(M) + 1 or not np.isfinite(r.yhat): continue
            j = n  # target index in the predict_run sense = n (passes < n used)
            if j >= len(y): continue
            chk.append(abs(py[j] - r.yhat)); rows.append(dict(tag=t, dv=r.delta_v_mps, yhat=r.yhat, yz=pz[j], ysys_part=py[j] - pz[j], collapsed=bool(r.collapsed)))
    R = pd.DataFrame(rows)
    if len(R) < 30: return dict(arm=arm, n=len(R))
    Z = sj.demean(R, ["dv", "yhat", "yz", "ysys_part"], R.tag)
    def fit(s):
        Zs = pd.concat([Z[R.tag == t] for t in s]); A = np.column_stack([Zs.ysys_part, Zs.yz]); b = np.linalg.lstsq(A, Zs.dv.values, rcond=None)[0]; v = A @ b
        return b, np.cov(v, b[1] * Zs.yz)[0, 1] / np.var(v, ddof=1)
    tags = sorted(set(R.tag)); b, shz = fit(tags)
    return dict(arm=arm, n=len(R), frac_collapsed=float(R.collapsed.mean()), yhat_reproduction_max_abs_err=float(np.max(chk)), yhat_reproduction_median_err=float(np.median(chk)),
                var_within_cell_yhat=float(Z.yhat.var()), frac_var_yhat_from_z=float(Z.yz.var() / Z.yhat.var()), b_sys=float(b[0]), b_z=float(b[1]), b_z_ci=sj.cboot(R.tag.values, lambda s: fit(s)[0][1], rng),
                b_sys_ci=sj.cboot(R.tag.values, lambda s: fit(s)[0][0], rng), share_dv_command_from_z=float(shz), share_ci=sj.cboot(R.tag.values, lambda s: fit(s)[1], rng),
                corr_dv_yz=float(np.corrcoef(Z.dv, Z.yz)[0, 1]), corr_dv_ysyspart=float(np.corrcoef(Z.dv, Z.ysys_part)[0, 1]))


def drift_theory(root, sub, pl):
    d = pathlib.Path(root) / sub / pl / "runs"; out = {}; mu = MU[pl]
    for f in sorted(d.glob("adaptive_pert_*/passes.csv"))[:1]: P = pd.read_csv(f); a, e = P.a_m.median(), P.e.median()
    va = np.sqrt(mu * (1 - e) / (a * (1 + e))); kappa = 2 * 2 * a ** 2 * va / mu   # dr_p = 2 da, da = 2 a^2 v dv / mu at fixed r_a (ASSUMED impulsive at apoapsis)
    tot = []
    for f in sorted(d.glob("adaptive_pert_*/decisions.csv")): tot.append(pd.read_csv(f).delta_v_mps.abs().sum())
    return dict(kappa_theory_m_per_mps=float(kappa), a_m=float(a), e=float(e), mean_total_dv_per_run_adaptive=float(np.mean(tot)) if tot else None)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--f20"); ap.add_argument("--stagei"); ap.add_argument("--out"); a = ap.parse_args(); rng = np.random.default_rng(sj.SEED); R = {}
    for name, root, sub in (("F20", a.f20, "stagef20"), ("I", a.stagei, "stagei")):
        R[name] = {pl: dict(z=spread_z(root, sub, pl, rng), dec={arm: decompose(root, sub, pl, arm, rng) for arm in ("adaptive_pert", "adaptive_asym_pert", "adaptive_h025_pert")}, kappa=drift_theory(root, sub, pl)) for pl in PLANETS}
        for pl in PLANETS: R[name][pl]["drift_dv"] = sj.drift(a.f20, "stagef20", pl, R[name][pl]["kappa"]["kappa_theory_m_per_mps"])   # open-arm drift always from Stage F20 (band-independent orbit property)
    if a.out: json.dump(R, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))
