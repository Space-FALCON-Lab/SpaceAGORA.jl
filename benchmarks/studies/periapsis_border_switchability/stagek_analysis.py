"""Stage K analysis: natural (burn-free) periapsis drift as a free lever. Python only, offline, no flights.

Plan and pre-registered criteria K0..K6: tasks/2026-10-15_periapsis-region-border-targeting/stagek_plan.md (frozen before this ran on flight tables).

Usage:
  python3 -I stagek_analysis.py --phase selftest
  python3 -I stagek_analysis.py --phase maps    --f20 <stagef20 tree> --out <dir>     # numerical zonal maps per pair (cached CSV per planet)
  python3 -I stagek_analysis.py --phase analyze --f20 <stagef20 tree> --out <dir>     # criteria K0..K6, tables, json
--f20 is the stagef20 directory, which holds <planet>/cells_base.csv and <planet>/runs/open_mean_<tag>/passes.csv.
"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagek_lib as K

PLANETS = ["mars", "earth", "venus"]
E_MIN = 0.1                  # ASSUMED pair filter
SEED = 20261022              # ASSUMED bootstrap seed
NBOOT = 10000
DAY = 86400.0
# SOURCED textbook constants for the J2-dominance rule (not repo files)
K_B = {"mars": (2 * math.pi / (686.98 * DAY)) ** 2, "venus": (2 * math.pi / (224.70 * DAY)) ** 2,
       "earth": 4.9028e12 / 384400e3 ** 3 + 1.32712e20 / 1.496e11 ** 3}


def wrap(x): return (np.asarray(x) + math.pi) % (2 * math.pi) - math.pi


# ------------------------------------------------------------------ data
def load_pairs(f20, planet):
    root = pathlib.Path(f20) / planet; cells = pd.read_csv(root / "cells_base.csv"); rows = []
    for c in cells.itertuples():
        f = root / "runs" / f"open_mean_{c.tag}" / "passes.csv"
        if not f.exists(): continue
        d = pd.read_csv(f).sort_values("pass").reset_index(drop=True)
        for k in range(len(d) - 1):
            r0, r1 = d.iloc[k], d.iloc[k + 1]
            if r1["pass"] - r0["pass"] != 1 or r0.e < E_MIN or r1.e < E_MIN: continue
            rows.append(dict(planet=planet, tag=c.tag, split=c.split, epoch=int(c.epoch_idx), inc_deg=c.inc, raan0=c.raan, aop0=c.omega, k=int(r0["pass"]),
                             t0=r0.time_s, dt=r1.time_s - r0.time_s, a=r0.a_m, e=r0.e, i=r0.i_rad, raan=r0.raan_rad, w=r0.omega_rad,
                             a1=r1.a_m, e1=r1.e, w1=r1.omega_rad, r0=r0.radius_m, r1=r1.radius_m, dr=r1.radius_m - r0.radius_m, da=r1.a_m - r0.a_m,
                             dw_obs=float(wrap(r1.omega_rad - r0.omega_rad)), npass=len(d)))
    return pd.DataFrame(rows)


def add_maps(P, df):
    out = {k: [] for k in ("dr_Z", "da_Z", "dw_Z", "dr_Z2", "dr_A3", "dw_A")}
    for r in df.itertuples():
        o = K.zonal_map(P, r.a, r.e, r.i, r.raan, r.w, nmax=20)
        o2 = K.zonal_map(P, r.a, r.e, r.i, r.raan, r.w, nmax=2)
        T = 2 * math.pi * math.sqrt(r.a ** 3 / P["mu"]); _, wd, _ = K.zonal_secular(P, r.a, r.e, r.i, r.w)
        ed, _, _ = K.zonal_secular(P, r.a, r.e, r.i, r.w + 0.5 * wd * T, nlo=3)
        out["dr_Z"].append(o[0] * (1 - o[1]) - r.a * (1 - r.e)); out["da_Z"].append(o[0] - r.a); out["dw_Z"].append(float(wrap(o[4] - r.w)))
        out["dr_Z2"].append(o2[0] * (1 - o2[1]) - r.a * (1 - r.e)); out["dr_A3"].append(-r.a * ed * T); out["dw_A"].append(wd * T)
    for k, v in out.items(): df[k] = v
    df["dr_ZA"] = df.dr_Z2 + df.dr_A3
    return df


# ------------------------------------------------------------------ statistics
def ev(dr, model, centered=False):
    dr = np.asarray(dr); res = dr - np.asarray(model); base = dr - dr.mean() if centered else dr
    return 1 - (res ** 2).sum() / (base ** 2).sum() if (base ** 2).sum() > 0 else float("nan")


def cell_boot(df, fn, nboot=NBOOT, seed=SEED, strata=None):
    """Bootstrap over cells (resampling cell tags, within strata if given); fn(sub_df) -> float. Returns (point, lo, hi) 5/95 or 2.5/97.5 given by caller."""
    rng = np.random.default_rng(seed); groups = {t: g for t, g in df.groupby("tag")}; tags = np.array(list(groups))
    st = {t: strata(groups[t]) for t in tags} if strata else {t: 0 for t in tags}
    by = {}
    for t in tags: by.setdefault(st[t], []).append(t)
    vals = []
    for _ in range(nboot):
        pick = [t for s, ts in by.items() for t in rng.choice(ts, len(ts))]
        vals.append(fn(pd.concat([groups[t] for t in pick])))
    return np.array(vals)


def fit_beta(dev):
    x = dev.da - dev.da_Z; y = dev.dr - dev.dr_Z
    return float((x * y).sum() / (x * x).sum()) if (x * x).sum() > 0 else 0.0


def eps_3b(P, planet, df):
    """Tide-to-zonal drift ratio at the first pass of each cell: S_3B / |dr_Z|, S_3B = 2 pi a Pi_3B (15/4) e eta (the coefficient scale of the
    single-averaged quadrupole rate in shortperiod_analytic.jl; Pi_3B = K_b a^3 / mu), dr_Z the numerical zonal map. Ratio of the medians over cells."""
    f = df[df.k == df.groupby("tag").k.transform("min")]
    S = 2 * math.pi * f.a * (K_B[planet] * f.a ** 3 / P["mu"]) * 3.75 * f.e * np.sqrt(1 - f.e ** 2)
    return float(np.median(S) / max(np.median(np.abs(f.dr_Z)), 1e-3))


def lever(P, df): return np.array([K.lever_m_per_mps(P["mu"], a, e) for a, e in zip(df.a, df.e)])


def ang_window(P, a, e, inc, thr=0.25, step=5):
    ws = np.radians(np.arange(0, 360, step)); dr = np.array([(lambda o: o[0] * (1 - o[1]) - a * (1 - e))(K.zonal_map(P, a, e, inc, 0.0, w, nmax=20)) for w in ws])
    return np.degrees(ws), dr


def run_analysis(f20, out, planets=PLANETS):
    out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True); res = {"planets": {}}; lines = []
    say = lambda *a: (lines.append(" ".join(str(x) for x in a)), print(*a))
    for pl in planets:
        P = K.planet(pl.capitalize()); cache = out / f"pairs_model_{pl}.csv"
        df = pd.read_csv(cache)
        lev = lever(P, df); df["lev"] = lev; df["cancel_dv"] = np.abs(df.dr) / lev; df["pi_drift"] = K.pi_drift(df.a, df.e, df.dr)
        df["pi_drift_Z"] = K.pi_drift(df.a, df.e, df.dr_Z)
        dev, held = df[df.split == "dev"], df[df.split == "heldout"]
        beta = fit_beta(dev); df["dr_M"] = df.dr_Z + beta * (df.da - df.da_Z); dev, held = df[df.split == "dev"], df[df.split == "heldout"]
        R = dict(n_pairs=int(len(df)), n_cells=int(df.tag.nunique()), beta=beta, eps_3B=eps_3b(P, pl, df))
        R["zonal_dominated"] = bool(R["eps_3B"] <= 0.5)
        for nm, sub in (("all", df), ("dev", dev), ("heldout", held)):
            R[f"EV_{nm}"] = {m: ev(sub.dr, sub[c]) for m, c in (("Z2", "dr_Z2"), ("Z", "dr_Z"), ("ZA", "dr_ZA"), ("M", "dr_M"))}
            R[f"EVc_{nm}"] = {m: ev(sub.dr, sub[c], True) for m, c in (("Z2", "dr_Z2"), ("Z", "dr_Z"), ("M", "dr_M"))}
        b = cell_boot(held, lambda s: ev(s.dr, s.dr_M), nboot=NBOOT)
        R["EV_M_heldout_ci90"] = [float(np.percentile(b, 5)), float(np.percentile(b, 95))]
        # K0
        R["K0_median_rel"] = float(np.median(np.abs(df.dr_ZA - df.dr_Z) / np.maximum(np.abs(df.dr_Z), 5.0)))
        # K2 sign agreement
        R["K2"] = {}
        for sp, sub in (("heldout", held), ("dev", dev), ("all", df)):
            for db in (10, 25, 50):
                m = np.abs(sub.dr) >= db
                R["K2"][f"{sp}_{db}m"] = dict(n=int(m.sum()), frac=float((np.sign(sub.dr_M[m]) == np.sign(sub.dr[m])).mean()) if m.sum() else None)
        # K5 apsidal
        k5 = []
        for tag, g in df.groupby("tag"):
            pred = np.degrees(g.dw_A.sum()); obs = np.degrees(g.dw_obs.sum())
            if abs(pred) >= 5.0: k5.append(dict(tag=tag, pred_deg=pred, obs_deg=obs, rel=abs(obs / pred - 1)))
        R["K5"] = dict(n_cells=len(k5), median_rel=float(np.median([x["rel"] for x in k5])) if k5 else None, cells=k5)
        # K4 lever scaling (cells with >= 20 pairs)
        k4 = []
        for tag, g in df.groupby("tag"):
            if len(g) < 20: continue
            g = g.sort_values("k"); e5, l5 = g.iloc[:5], g.iloc[-5:]
            q = lambda s, c: float((np.abs(s[c]) / s.lev).mean())
            k4.append(dict(tag=tag, rho_l=float(l5.lev.mean() / e5.lev.mean()), rho_Q=q(l5, "dr") / q(e5, "dr") if q(e5, "dr") > 0 else None,
                           rho_Q_pred=q(l5, "dr_Z") / q(e5, "dr_Z") if q(e5, "dr_Z") > 0 else None))
        ok = [x for x in k4 if x["rho_Q"] and x["rho_Q_pred"]]
        R["K4"] = dict(n_cells=len(k4), median_rho_l=float(np.median([x["rho_l"] for x in k4])) if k4 else None,
                       median_abs_ln_ratio=float(np.median([abs(math.log(x["rho_Q_pred"] / x["rho_Q"])) for x in ok])) if ok else None,
                       median_rho_Q=float(np.median([x["rho_Q"] for x in ok])) if ok else None, median_rho_Q_pred=float(np.median([x["rho_Q_pred"] for x in ok])) if ok else None, cells=k4)
        # K3 AOP saving
        h30 = df[df.t0 <= 30 * DAY]; cellQ = []
        for tag, g in h30.groupby("tag"):
            span = (g.t0.max() + g.dt.iloc[-1] - g.t0.min()) / DAY
            cellQ.append(dict(tag=tag, inc_deg=g.inc_deg.iloc[0], aop0=g.aop0.iloc[0], Q=float(g.cancel_dv.sum() / span), Qnet=float(abs(g.dr.sum()) / g.lev.mean() / span)))
        cq = pd.DataFrame(cellQ); rank = {}
        for inc, g in df.groupby("inc_deg"):
            f1 = g[g.k == g.groupby("tag").k.transform("min")]; a_, e_ = float(f1.a.median()), float(f1.e.median())
            ths = {w: abs((lambda o: o[0] * (1 - o[1]) - a_ * (1 - e_))(K.zonal_map(P, a_, e_, math.radians(inc), 0.0, math.radians(w), nmax=20))) for w in (0, 90, 180, 270)}
            rank[inc] = dict(fav=min(ths, key=ths.get), unfav=max(ths, key=ths.get), pred_abs_dr={str(k): v for k, v in ths.items()}, a=a_, e=e_)
        cq["role"] = [("fav" if r.aop0 == rank[r.inc_deg]["fav"] else "unfav" if r.aop0 == rank[r.inc_deg]["unfav"] else "mid") for r in cq.itertuples()]

        def saving(sub, col):
            d = []
            for inc, g in sub.groupby("inc_deg"):
                f, u = g[g.role == "fav"][col], g[g.role == "unfav"][col]
                if len(f) and len(u): d.append(u.mean() - f.mean())
            return float(np.mean(d)) if d else float("nan")
        S = saving(cq, "Q"); Snet = saving(cq, "Qnet")
        rng = np.random.default_rng(SEED); bs = []; bn = []
        strat = {(r.inc_deg, r.role): [] for r in cq.itertuples()}
        for r in cq.itertuples(): strat[(r.inc_deg, r.role)].append(r.Index)
        for _ in range(NBOOT):
            idx = [i for ids in strat.values() for i in rng.choice(ids, len(ids))]; s = cq.loc[idx]
            bs.append(saving(s, "Q")); bn.append(saving(s, "Qnet"))
        bs, bn = np.array(bs), np.array(bn)
        unf = cq[cq.role == "unfav"].Q.mean(); fav = cq[cq.role == "fav"].Q.mean()
        R["K3"] = dict(S_mps_per_day=S, ci95=[float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))], pct_of_unfav=float(100 * S / unf) if unf else None,
                       mean_Q_unfav=float(unf), mean_Q_fav=float(fav), S_net=Snet, ci95_net=[float(np.nanpercentile(bn, 2.5)), float(np.nanpercentile(bn, 97.5))],
                       rank={str(k): v for k, v in rank.items()}, cells=cq.to_dict("records"))
        # K6 AOP curves (theory) and measured binned drift
        curves = {}
        for inc, rk in rank.items():
            w, drc = ang_window(P, rk["a"], rk["e"], math.radians(inc)); pic = K.pi_drift(rk["a"], rk["e"], drc); mx = np.abs(drc).max()
            curves[str(inc)] = dict(omega_deg=w.tolist(), dr_m=drc.tolist(), dv_mps=(drc / K.lever_m_per_mps(P["mu"], rk["a"], rk["e"])).tolist(), pi_drift=pic.tolist(),
                                    window_deg=float(5 * (np.abs(drc) < 0.25 * mx).sum()), max_abs_dr=float(mx),
                                    zero_crossings_deg=[float(w[i]) for i in range(len(w)) if drc[i] * drc[(i + 1) % len(w)] < 0])
        R["K6_curves"] = curves
        bins = np.floor(np.degrees(df.w % (2 * math.pi)) / 45).astype(int); df["wbin"] = bins
        tab = df.groupby(["inc_deg", "wbin"]).agg(n=("dr", "size"), mean_dr=("dr", "mean"), mean_dv=("cancel_dv", "mean"), mean_signed_dv=("dr", lambda s: 0.0)).reset_index()
        tab["mean_signed_dv"] = df.assign(sd=df.dr / df.lev).groupby(["inc_deg", "wbin"]).sd.mean().values
        tab.to_csv(out / f"aop_bins_{pl}.csv", index=False)
        # residual structure
        res_ = df.dr - df.dr_M; tot = ((res_ - res_.mean()) ** 2).sum()
        gr = lambda col: float(1 - ((res_ - res_.groupby(df[col]).transform("mean")) ** 2).sum() / tot) if tot > 0 else None
        R["K6_residual_variance_share"] = dict(raan=gr("raan0"), epoch=gr("epoch"), tag=gr("tag"))
        R["summary_dv"] = dict(median_abs_dr_m=float(np.median(np.abs(df.dr))), median_cancel_dv_mps_per_pass=float(np.median(df.cancel_dv)),
                               median_cancel_dv_mps_per_day=float(np.median(df.cancel_dv / df.dt * DAY)), median_abs_pi_drift=float(np.median(np.abs(df.pi_drift))))
        df.to_csv(out / f"pairs_full_{pl}.csv", index=False)
        res["planets"][pl] = R
        say(f"{pl}: pairs {R['n_pairs']} cells {R['n_cells']} beta {beta:.3e} eps_3B {R['eps_3B']:.3g} zonal-dominated {R['zonal_dominated']}")
        say(f"  EV heldout {R['EV_heldout']}  EV all {R['EV_all']}  M CI90 {R['EV_M_heldout_ci90']}")
    # verdicts
    m = res["planets"].get("mars")
    if m:
        v = {}
        v["K0"] = m["K0_median_rel"] <= 0.15
        v["K1"] = bool(m["EV_heldout"]["M"] >= 0.80)
        f = m["K2"]["heldout_10m"]; v["K2"] = bool(f["frac"] is not None and f["frac"] >= 0.90)
        v["K3"] = bool(m["K3"]["S_mps_per_day"] > 0 and m["K3"]["ci95"][0] > 0)
        v["K4"] = bool(m["K4"]["median_abs_ln_ratio"] is not None and m["K4"]["median_abs_ln_ratio"] <= math.log(1.5)) and bool(m["K4"]["median_rho_l"] is not None and m["K4"]["median_rho_l"] < 1)
        v["K5"] = bool(m["K5"]["median_rel"] is not None and m["K5"]["median_rel"] <= 0.10)
        res["verdicts_mars"] = v; say("Mars verdicts:", v)
    (out / "k_results.json").write_text(json.dumps(res, indent=1, default=float)); (out / "k_summary.txt").write_text("\n".join(lines) + "\n")
    return res


def run_maps(f20, out, planets=PLANETS):
    out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True)
    for pl in planets:
        P = K.planet(pl.capitalize()); df = load_pairs(f20, pl); df = add_maps(P, df); df.to_csv(out / f"pairs_model_{pl}.csv", index=False)
        print(pl, len(df), "pairs; zonal map done")


# ------------------------------------------------------------------ selftest on synthetic cells
def selftest():
    import tempfile
    P = K.planet("Mars"); tmp = pathlib.Path(tempfile.mkdtemp()); root = tmp / "stagef20" / "mars"; (root / "runs").mkdir(parents=True)
    rows = []; rng = np.random.default_rng(1)
    for n in range(12):
        inc = [25.2, 63.4][n % 2]; aop = [0, 90, 180, 270][n % 4]; tag = f"I{int(inc * 10):04d}_R000_W{aop:03d}_E{1 + n % 3}_{n}"
        rows.append((tag, inc, 0.0, aop, 1 + n % 3, "dev" if n % 3 else "heldout"))
        st = (1.6e7, 0.78, math.radians(inc), 0.3, math.radians(aop)); pas = []; t = 0.0
        for k in range(1, 31):
            pas.append(dict(**{"pass": k, "time_s": t, "radius_m": st[0] * (1 - st[1]), "a_m": st[0], "e": st[1], "i_rad": st[2], "raan_rad": st[3], "omega_rad": st[4]}))
            o = K.zonal_map(P, *st, nmax=20); t += o[5]; st = o[:5]
        d = pd.DataFrame(pas); (root / "runs" / f"open_mean_{tag}").mkdir(); d.to_csv(root / "runs" / f"open_mean_{tag}" / "passes.csv", index=False)
    pd.DataFrame(rows, columns=["tag", "inc", "raan", "omega", "epoch_idx", "split"]).to_csv(root / "cells_base.csv", index=False)
    df = add_maps(P, load_pairs(tmp / "stagef20", "mars")); assert len(df) == 12 * 29, len(df)
    assert np.abs(df.dr - df.dr_Z).max() < 1e-3, np.abs(df.dr - df.dr_Z).max()          # synthetic data IS the zonal map
    assert ev(df.dr, df.dr_Z) > 0.99999
    df.to_csv(tmp / "pairs_model_mars.csv", index=False)
    r = run_analysis(tmp / "stagef20", tmp, planets=["mars"])
    m = r["planets"]["mars"]; assert m["EV_heldout"]["M"] > 0.9999, m["EV_heldout"]; assert m["K2"]["heldout_10m"]["frac"] in (1.0, None)
    assert m["K0_median_rel"] < 0.1 and np.isfinite(m["eps_3B"])
    assert m["K5"]["n_cells"] > 0 and m["K5"]["median_rel"] < 0.05, m["K5"]
    print("SELFTEST OK", tmp)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True, choices=["selftest", "maps", "analyze"])
    ap.add_argument("--f20"); ap.add_argument("--out")
    a = ap.parse_args()
    if a.phase == "selftest": selftest()
    elif a.phase == "maps": run_maps(a.f20, a.out)
    else: run_analysis(a.f20, a.out)
