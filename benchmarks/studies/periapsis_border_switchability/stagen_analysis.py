"""Stage N analysis (plan: tasks/2026-10-15_periapsis-region-border-targeting/stagen_plan.md, frozen before any statistic).
python3 -B stagen_analysis.py --phase selftest | run --stagei ROOT --f20 ROOT [--out f.json]
Offline Python only. Decomposition, buffer ceiling, chance-constrained edge rule replay."""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
from scipy.stats import norm, beta
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import stagej_analysis as sj
PLANETS = sj.PLANETS; NB = 4000; SEED = 20261023; EPS = [0.001, 0.01, 0.05]; RHOS = [0, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0]
QMAX, QLOW, QHARD = 0.25, 0.05, 0.45; LN5 = float(np.log(5.0)); LNH = float(np.log(QHARD / QMAX))
ARMS = dict(fixed="fixed_pert", F="adaptive_pert", G="adaptive_asym_pert", H="adaptive_h025_pert")


def spread(e): return sj.spread(np.asarray(e))


def boot(tags, vals, fn, rng, nb=NB):
    """cluster bootstrap: vals dict of arrays aligned with tags; fn(dict of concatenated arrays)."""
    ut = np.array(sorted(set(tags))); idx = {t: np.where(tags == t)[0] for t in ut}; out = []
    for _ in range(nb):
        i = np.concatenate([idx[t] for t in rng.choice(ut, len(ut))]); r = fn({k: v[i] for k, v in vals.items()})
        if np.isfinite(r): out.append(r)
    return [float(np.percentile(out, 5)), float(np.percentile(out, 95))]


def decompose(P):
    """P in-window passes (all cells). returns arrays and shares."""
    b = P.groupby("tag").ysys.transform("mean").values; d = P.ysys.values - b; z = P.z.values; y = P.y.values
    return dict(y=y, b=b, d=d, z=z)


def shares(v):
    vy = np.var(v["y"]); s = dict(between=np.var(v["b"]) / vy, drift=np.var(v["d"]) / vy, white=np.var(v["z"]) / vy)
    s["cross"] = 1 - s["between"] - s["drift"] - s["white"]; return s


def kalman(P, r, q, m0, v0):
    """per-cell local-level filter over in-window passes; returns arrays m (predictive mean), s2 (predictive variance), e (error), tag, in pass order."""
    M, S, E, T = [], [], [], []
    for t, G in P.groupby("tag", sort=False):
        y, ps = G.y.values, G["pass"].values; m, p = m0, v0
        for k in range(len(y)):
            if k > 0: p += q * (ps[k] - ps[k - 1])
            s2 = p + r; M.append(m); S.append(s2); E.append(y[k] - m); T.append(t)
            K = p / s2; m = m + K * (y[k] - m); p = (1 - K) * p
    return np.array(M), np.array(S), np.array(E), np.array(T)


def fit_kalman(D):
    """dev-only: r, m0, v0 observable stand-ins, rho by dev NLL."""
    dy = np.concatenate([np.diff(G.y.values)[np.diff(G["pass"].values) == 1] for _, G in D.groupby("tag")]); r = 0.5 * np.var(dy)
    m0 = float(D.y.mean()); v0 = float(np.var(D.groupby("tag").y.mean())); best = None
    for rho in RHOS:
        _, S, E, _ = kalman(D, r, rho * r, m0, v0); nll = float(np.sum(0.5 * np.log(S) + 0.5 * E ** 2 / S))
        if best is None or nll < best[1]: best = (rho, nll)
    return dict(r=float(r), m0=m0, v0=v0, rho=best[0], q=best[0] * r)


def cp(k, n, a=0.10):
    lo = beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0; return [float(lo), float(hi)]


def replay(Hm, Hs2, He, tags, tq, rng, eps):
    """tq: multiplier t. c = sqrt(s2) t; counterfactual Q_aim = QMAX exp(-c), Q = Q_aim exp(e)."""
    c = np.sqrt(Hs2) * tq; over = He > c; hard = He > c + LNH; low = He < c - LN5; qa = QMAX * np.exp(-c); qr = qa * np.exp(He); n = len(He)
    v = dict(o=over.astype(float), h=hard.astype(float))
    return dict(eps=eps, n=n, over=float(over.mean()), over_k=int(over.sum()), over_ci_binom=cp(int(over.sum()), n), over_ci_cluster=boot(tags, v, lambda d: d["o"].mean(), rng),
               hard=float(hard.mean()), hard_k=int(hard.sum()), low=float(low.mean()), mean_Qaim=float(qa.mean()), mean_Qreal=float(qr.mean()), median_Qaim=float(np.median(qa)),
               calibrated=bool(cp(int(over.sum()), n)[0] <= eps <= cp(int(over.sum()), n)[1]))


def flown(root, sub, pl):
    out = {}
    for k, arm in ARMS.items():
        P = sj.load(root, sub, pl, arm); P = P[P.split == "heldout"]; q = P.q_eff_w_cm2.values
        out[k] = dict(n=len(q), over=float((q > QMAX).mean()), hard=float((q > QHARD).mean()), low=float((q < QLOW).mean()), mean_Q=float(q.mean()))
    for k in ("F", "G", "H"): out[k]["meanQ_ratio_fixed"] = out[k]["mean_Q"] / out["fixed"]["mean_Q"]
    out["fixed"]["meanQ_ratio_fixed"] = 1.0; return out


def planet(root, sub, pl, rng):
    P = sj.load(root, sub, pl, "fixed_pert"); P = P[P.inw].reset_index(drop=True); D = P[P.split == "dev"].reset_index(drop=True); H = P[P.split == "heldout"].reset_index(drop=True)
    res = dict(n_cells=int(P.tag.nunique()), n_held_passes=len(H), n_dev_passes=len(D), h_local_m=float(H.h_local_m.median()))
    # 1 decomposition
    for nm, S in (("all", P), ("held", H)):
        v = decompose(S); s = shares(v); res[f"shares_{nm}"] = s
        tags = S.tag.values; ci = {k: boot(tags, v, (lambda kk: lambda d: shares(d)[kk])(k), rng) for k in ("between", "drift", "white", "cross")}; res[f"shares_{nm}_ci"] = ci
        res[f"std_{nm}"] = {k: float(np.std(v[k])) for k in v}
    # 2 ceiling + kalman
    ap = pd.read_csv(pathlib.Path(root) / sub / pl / "adapt_params.csv").iloc[0]; ce = sj.ctrl_errors(H.assign(inw=True), ap.predictor)
    v = decompose(H); kp = fit_kalman(D); res["kalman_params"] = kp; Hm, Hs2, He, Ht = kalman(H, kp["r"], kp["q"], kp["m0"], kp["v0"])
    errs = dict(a_current=(ce.e.values, ce.tag.values), b_oracle_offset=(v["y"] - v["b"], H.tag.values), c_oracle_offset_drift=(v["z"], H.tag.values), d_kalman=(He, Ht))
    hl = res["h_local_m"]; ceil = {}
    for k, (e, tg) in errs.items():
        sp = spread(e); ceil[k] = dict(n=len(e), spread_ln=sp, spread_ci=boot(tg, dict(e=e), lambda d: spread(d["e"]), rng), std=float(np.std(e)), q_hi={str(p): float(np.quantile(e, 1 - p)) for p in EPS},
                                        q95=float(np.quantile(e, .95)), q05=float(np.quantile(e, .05)), Q_lo_at_015=float(0.15 * np.exp(np.quantile(e, .05))), Q_hi_at_015=float(0.15 * np.exp(np.quantile(e, .95))),
                                        corridor_ln_left_paper=LN5 - sp, altitude_shift_m=sp * hl)
    sa = ceil["a_current"]["spread_ln"]
    for k in ceil: ceil[k]["frac_of_current"] = ceil[k]["spread_ln"] / sa; ceil[k]["recovered_vs_current"] = 1 - ceil[k]["spread_ln"] / sa
    ceil["N2_irreducible_dominated"] = bool(ceil["c_oracle_offset_drift"]["spread_ln"] / sa >= 0.7); ceil["N4_kalman_recovers"] = bool(ceil["d_kalman"]["spread_ln"] / sa <= 0.9)
    res["ceiling"] = ceil
    # 3 chance-constrained rule
    _, Sd, Ed, _ = kalman(D, kp["r"], kp["q"], kp["m0"], kp["v0"]); zd = Ed / np.sqrt(Sd); rules = {}
    for kind in ("gaussian", "empirical"):
        rules[kind] = []
        for eps in EPS:
            tq = norm.ppf(1 - eps) if kind == "gaussian" else float(np.quantile(zd, 1 - eps))
            r = replay(Hm, Hs2, He, Ht, tq, rng, eps); r["t"] = float(tq); r["dev_expected_exceed_count"] = len(zd) * eps
            if kind == "empirical" and len(zd) * eps < 5: r["unresolvable_N3"] = True
            rules[kind].append(r)
    ecmp = ce.e.values; rules["ctrl_edge"] = dict(over_hi_freq=float((ecmp > ap.e_hi).mean()), e_hi=float(ap.e_hi), n=len(ecmp), ci_binom=cp(int((ecmp > ap.e_hi).sum()), len(ecmp)))
    rules["kalman_standardized_innov_lag1"] = float(sj.corr(*sj.lag1_pairs(H.assign(en=He / np.sqrt(Hs2)), "en")[:2]))
    res["rules"] = rules; res["flown_arms"] = flown(root, sub, pl)
    _, _, Ef, Tf = kalman(H, kp["r"], kp["q"], kp["m0"], kp["v0"]); res["fixed_arm_replay_check"] = dict(over_flown_inwindow=float((H.q_eff_w_cm2 > QMAX).mean()), hard_flown_inwindow=float((H.q_eff_w_cm2 > QHARD).mean()),
                                                                   low_flown_inwindow=float((H.q_eff_w_cm2 < QLOW).mean()), mean_Q_inwindow=float(H.q_eff_w_cm2.mean()))
    return res


def selftest():
    rng = np.random.default_rng(3); n = 60; rows = []
    for c in range(40):
        b = rng.normal(0, 0.5); y = b + rng.normal(0, 0.2, n)
        rows.append(pd.DataFrame(dict(tag=f"c{c}", split="dev" if c < 20 else "heldout", **{"pass": np.arange(n)}, y=y, z=y - b, ysys=b)))
    P = pd.concat(rows, ignore_index=True); s = shares(decompose(P)); assert s["between"] > 0.8 and abs(s["white"] - 0.138) < 0.05 and abs(s["cross"]) < 0.05 and s["drift"] < 1e-9, s
    D, H = P[P.split == "dev"].reset_index(drop=True), P[P.split == "heldout"].reset_index(drop=True); kp = fit_kalman(D); assert kp["rho"] <= 1e-2 and abs(kp["r"] - 0.04) < 0.01, kp
    Hm, Hs2, He, Ht = kalman(H, kp["r"], kp["q"], kp["m0"], kp["v0"]); late = np.concatenate([He[i * n + 10:(i + 1) * n] for i in range(20)])
    assert abs(np.std(late) - 0.2) < 0.03, np.std(late)                       # filter converges to white floor
    zz = He / np.sqrt(Hs2); rr = replay(Hm, Hs2, He, Ht, norm.ppf(0.95), rng, 0.05); assert rr["calibrated"] and abs(rr["over"] - 0.05) < 0.02, rr  # Gaussian rule is calibrated on Gaussian data
    assert cp(5, 100)[0] < 0.05 < cp(5, 100)[1]; print("selftest ok")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", default="run"); ap.add_argument("--f20"); ap.add_argument("--stagei"); ap.add_argument("--out"); a = ap.parse_args()
    if a.phase == "selftest": selftest(); sys.exit()
    rng = np.random.default_rng(SEED); R = {}
    for name, root, sub in (("I", a.stagei, "stagei"), ("F20", a.f20, "stagef20")): R[name] = {pl: planet(root, sub, pl, rng) for pl in PLANETS}
    if a.out: json.dump(R, open(a.out, "w"), indent=1, default=float)
