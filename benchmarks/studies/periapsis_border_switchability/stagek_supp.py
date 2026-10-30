"""Stage K supplements (UNFROZEN: added after the pre-registered criteria were evaluated; none changes a verdict). Python only, reads the tables written by
stagek_analysis.py and stagek_tess.py.

Usage: python3 -I stagek_supp.py --out <dir>
 S1  K4 re-evaluated with the zonal + tesseral + drag prediction (Mars)
 S2  measured cost by initial AOP (Q = pay-as-you-go cancel dV per day, Qnet = net-drift dV per day, first 30 d) with cell-bootstrap intervals, all planets
 S3  eta^2 (share of cell-to-cell variance) of AOP, inclination, RAAN, epoch for Q and for the signed 30-d net drift; Mars also for the 10-pass net
 S4  campaign-level arithmetic: dV over 30 d for pay-as-you-go and net, favorable-versus-unfavorable AOP by the measured medians, apsidal rotation rate by inclination (walk-in time)
 S5  sign statistics of Pi_drift (share positive, lag-1 autocorrelation of dr, ratio of mean drift to rms per-pass change)
"""
import argparse, json, math, pathlib
import numpy as np, pandas as pd

DAY = 86400.0; SEED = 20261022; NBOOT = 10000


def eta2(df, col, grp):
    tot = ((df[col] - df[col].mean()) ** 2).sum(); within = ((df[col] - df.groupby(grp)[col].transform("mean")) ** 2).sum()
    return float(1 - within / tot) if tot > 0 else None


def keeper_dv(f20, out):
    """S6: flown keeper (fixed 1000 m margin) total dV per cell from the Stage F20 F1 summaries (existing flights), by AOP role from the K3 ranking (flight-free)."""
    f20 = pathlib.Path(f20); k = json.load(open(pathlib.Path(out) / "k_results.json"))["planets"]; res = {}
    srcs = {"mars": "summary_f1_ws20.csv", "earth": "summary_f1_ws20.csv", "venus": "summary_f1_venusstagef.csv"}
    for pl in ("mars", "earth", "venus"):
        d = pd.read_csv(f20 / srcs[pl]); d = d[(d.planet == pl) & (d.ctrl == "fixed") & (d.status == "ok")]
        cells = pd.read_csv(f20 / pl / "cells_base.csv"); rk = k[pl]["K3"]["rank"]; R = {}
        for arm in ("mean", "pert"):
            x = d[d.arm == arm].merge(cells, on="tag")
            x["role"] = [("fav" if r.omega == rk[str(r.inc)]["fav"] else "unfav" if r.omega == rk[str(r.inc)]["unfav"] else "mid") for r in x.itertuples()]
            x["inc_key"] = x.inc
            def sav(sub):
                v = [g[g.role == "unfav"].total_dv_mps.mean() - g[g.role == "fav"].total_dv_mps.mean() for _, g in sub.groupby("inc_key") if (g.role == "fav").any() and (g.role == "unfav").any()]
                return float(np.mean(v)) if v else float("nan")
            rng = np.random.default_rng(SEED); strata = {kk: list(g.index) for kk, g in x.groupby(["inc_key", "role"])}
            b = np.array([sav(x.loc[[i for ids in strata.values() for i in rng.choice(ids, len(ids))]]) for _ in range(NBOOT)])
            R[arm] = dict(n=int(len(x)), mean_dv=float(x.total_dv_mps.mean()), sd_dv=float(x.total_dv_mps.std()), S_unfav_minus_fav=sav(x), ci95=[float(np.nanpercentile(b, 2.5)), float(np.nanpercentile(b, 97.5))],
                          by_aop={str(int(a)): dict(n=int(len(g)), mean=float(g.total_dv_mps.mean())) for a, g in x.groupby("omega")},
                          eta2_aop=eta2(x, "total_dv_mps", "omega"), eta2_inc=eta2(x, "total_dv_mps", "inc"), eta2_raan=eta2(x, "total_dv_mps", "raan"), eta2_epoch=eta2(x, "total_dv_mps", "epoch_idx"))
        res[pl] = R
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--f20"); a = ap.parse_args(); out = pathlib.Path(a.out); res = {}
    rng = np.random.default_rng(SEED)
    for pl in ("mars", "earth", "venus"):
        f = out / ("pairs_full_mars_tess.csv" if pl == "mars" else f"pairs_full_{pl}.csv"); df = pd.read_csv(f); R = {}
        h = df[df.t0 <= 30 * DAY]; rows = []
        for tag, g in h.groupby("tag"):
            g = g.sort_values("k"); span = (g.t0.max() + g.dt.iloc[-1] - g.t0.min()) / DAY
            rows.append(dict(tag=tag, inc=g.inc_deg.iloc[0], aop=g.aop0.iloc[0], raan=g.raan0.iloc[0], epoch=g.epoch.iloc[0], Q=g.cancel_dv.sum() / span,
                             Qnet=abs(g.dr.sum()) / g.lev.mean() / span, net_signed=g.dr.sum(), net_dv_signed=g.dr.sum() / g.lev.mean(), npairs=len(g), span=span,
                             wdot_deg_per_day=np.degrees(g.dw_A.sum()) / span, ZTD_net=(g["dr_ZTD"].sum() if "dr_ZTD" in g else np.nan)))
        c = pd.DataFrame(rows); R["cells_n"] = int(len(c))
        for col in ("Q", "Qnet"):
            R[f"eta2_{col}"] = {g: eta2(c, col, g) for g in ("aop", "inc", "raan", "epoch")}
        R["eta2_net_signed"] = {g: eta2(c, "net_dv_signed", g) for g in ("aop", "inc", "raan", "epoch")}
        by = {}
        for aop, g in c.groupby("aop"):
            b = [g.Q.sample(len(g), replace=True, random_state=int(rng.integers(1e9))).mean() for _ in range(2000)]
            bn = [g.Qnet.sample(len(g), replace=True, random_state=int(rng.integers(1e9))).mean() for _ in range(2000)]
            by[str(int(aop))] = dict(n=int(len(g)), Q_mean=float(g.Q.mean()), Q_ci95=[float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))],
                                     Qnet_mean=float(g.Qnet.mean()), Qnet_ci95=[float(np.percentile(bn, 2.5)), float(np.percentile(bn, 97.5))])
        R["by_aop"] = by
        qa = {k: v["Q_mean"] for k, v in by.items()}; qn = {k: v["Qnet_mean"] for k, v in by.items()}
        R["campaign_arithmetic"] = dict(span_days_median=float(c.span.median()), Q_mean_all=float(c.Q.mean()), Qnet_mean_all=float(c.Qnet.mean()),
                                        pay_as_you_go_30d_mps=float(30 * c.Q.mean()), net_30d_mps=float(30 * c.Qnet.mean()),
                                        best_aop_Q=min(qa, key=qa.get), worst_aop_Q=max(qa, key=qa.get), best_to_worst_Q_pct=float(100 * (1 - min(qa.values()) / max(qa.values()))),
                                        best_aop_Qnet=min(qn, key=qn.get), worst_aop_Qnet=max(qn, key=qn.get), best_to_worst_Qnet_pct=float(100 * (1 - min(qn.values()) / max(qn.values()))))
        R["wdot_deg_per_day_by_inc"] = {str(k): float(v) for k, v in c.groupby("inc").wdot_deg_per_day.median().items()}
        R["days_to_rotate_90deg_by_inc"] = {str(k): (float(90 / abs(v)) if abs(v) > 1e-6 else None) for k, v in c.groupby("inc").wdot_deg_per_day.median().items()}
        # S5 signs
        ac = [np.corrcoef(g.sort_values("k").dr.values[:-1], g.sort_values("k").dr.values[1:])[0, 1] for t, g in df.groupby("tag") if len(g) > 10]
        R["signs"] = dict(share_pi_drift_positive=float((df.pi_drift > 0).mean()), lag1_dr_median=float(np.median(ac)), mean_dr_m=float(df.dr.mean()), rms_dr_m=float(np.sqrt((df.dr ** 2).mean())),
                          mean_abs_over_rms=float(abs(df.dr.mean()) / np.sqrt((df.dr ** 2).mean())), per_cell_net_over_sum_abs_median=float(np.median([abs(g.dr.sum()) / g.dr.abs().sum() for t, g in df.groupby("tag")])))
        # S1 K4 with the full model (Mars)
        if pl == "mars":
            k4 = []
            for tag, g in df.groupby("tag"):
                if len(g) < 20: continue
                g = g.sort_values("k"); e5, l5 = g.iloc[:5], g.iloc[-5:]; q = lambda s, col: float((np.abs(s[col]) / s.lev).mean())
                k4.append((q(l5, "dr") / q(e5, "dr"), q(l5, "dr_ZTD") / q(e5, "dr_ZTD"), l5.lev.mean() / e5.lev.mean()))
            k4 = np.array(k4)
            R["K4_full_model"] = dict(n=int(len(k4)), median_abs_ln_ratio=float(np.median(np.abs(np.log(k4[:, 1] / k4[:, 0])))), median_rho_Q=float(np.median(k4[:, 0])),
                                      median_rho_Q_pred=float(np.median(k4[:, 1])), median_rho_l=float(np.median(k4[:, 2])), median_rho_Q_times_rho_l=float(np.median(k4[:, 0] * k4[:, 2])))
            # 10-pass nets: eta2 by AOP of the signed net
        res[pl] = R; print(pl, json.dumps({k: R[k] for k in ("eta2_Q", "eta2_Qnet", "eta2_net_signed", "campaign_arithmetic", "signs")}, indent=1, default=float))
    if a.f20: res["keeper_flown"] = keeper_dv(a.f20, out); print(json.dumps(res["keeper_flown"], indent=1, default=float))
    (out / "k_supp.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps({k: v for k, v in res["mars"].items() if k in ("K4_full_model", "by_aop", "wdot_deg_per_day_by_inc", "days_to_rotate_90deg_by_inc")}, indent=1, default=float))


if __name__ == "__main__":
    main()
