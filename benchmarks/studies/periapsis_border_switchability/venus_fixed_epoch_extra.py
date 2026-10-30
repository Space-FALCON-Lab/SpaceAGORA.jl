"""Extra offline checks for the Venus fixed-epoch study: (1) how fast the picked geometry degrades when the real epoch differs from the planned one (full ranking model, geometry held);
(2) the planner inside restricted inclination pools (where no near-zero-floor geometry exists). Usage: python3 -I venus_fixed_epoch_extra.py TRUTH_DIR"""
import json, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V
import venus_fixed_epoch as E


def main(td):
    an = pathlib.Path(td) / "analysis"; ctx, info = V.make_ctx(); res = {}
    pk = pd.read_csv(an / "planner_fine_pool_picks.csv"); rows = []
    for ep in E.HEAD:
        b = pk[(pk.epoch_off == ep) & (pk["rank"] == 0)].iloc[0]
        for d in (-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20):
            if ep + d < 0: continue
            c = dict(id="X", set="x", inc=float(b.inc), epoch_off=float(ep + d), omega=float(b.omega0), raan=float(b.raan0)); r, *_ = G.evaluate(ctx, c, None, None)
            rows.append(dict(planned_epoch=ep, epoch_error_d=d, geom="i%.0f/w%.0f/O%.0f" % (b.inc, b.omega0, b.raan0), floor_model=r["floor_mps"], R=r["R"]))
    S = pd.DataFrame(rows); S.to_csv(an / "extra_epoch_error.csv", index=False); print(S.pivot_table(index="planned_epoch", columns="epoch_error_d", values="floor_model").round(2).to_string())
    gr = pd.read_csv(an / "followup_grid_with_ev.csv"); rng = np.random.default_rng(0); pools = {"non-polar (i <= 70 or >= 110)": lambda d: (d.inc <= 70) | (d.inc >= 110), "critical only (63.4, 116.6)": lambda d: d.inc.isin([63.4, 116.6]),
                                                                                           "i = 39.2 only": lambda d: d.inc == 39.2, "i = 20 or 140 only": lambda d: d.inc.isin([20.0, 140.0])}
    out = []
    for nm, f in pools.items():
        for ep, g in gr.groupby("epoch_off"):
            g = g[f(g)]; best = g.floor.min(); pc = g.sort_values(["floor_cf", "R_cf"]).iloc[0]; pe = g.sort_values(["floor_ev", "R_ev"]).iloc[0]
            out.append(dict(pool=nm, epoch_off=ep, n=len(g), dp_best=best, plan_cf=pc.floor, plan_ev=pe.floor, random_mean=g.floor.mean(), random_best10=float(rng.choice(g.floor.values, (2000, 10)).min(1).mean()), vex_like_polar=np.nan))
    O = pd.DataFrame(out); O.to_csv(an / "extra_restricted_pools.csv", index=False)
    s = O.groupby("pool")[["dp_best", "plan_cf", "plan_ev", "random_mean", "random_best10"]].mean(); s["mean_regret_cf"] = (O.plan_cf - O.dp_best).groupby(O.pool).mean(); s["max_regret_cf"] = (O.plan_cf - O.dp_best).groupby(O.pool).max(); s["frac_pick_within_0.1_of_best"] = ((O.plan_cf - O.dp_best) <= 0.1).groupby(O.pool).mean()
    s["dp_best_min"] = O.groupby("pool").dp_best.min(); s["dp_best_max"] = O.groupby("pool").dp_best.max(); print(s.round(3).to_string()); s.to_csv(an / "extra_restricted_pools_summary.csv")


if __name__ == "__main__": main(sys.argv[1])
