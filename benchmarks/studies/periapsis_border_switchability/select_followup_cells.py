"""Follow-up strata, written BEFORE any follow-up run. python3 -I select_followup_cells.py [--results-root R]
Writes <R>/compare_baselines/closedloop_cells_followup.csv (new cells only) and closedloop_cells_all.csv (original 30 + new).
quiet_not_top: quiet-set cells ranked OUTSIDE the estimator-v2 top 10%; estimator_top_not_quiet: top-10% cells not in the quiet set;
quiet_not_top5 (fallback when quiet_not_top is empty): quiet-set cells outside the top 5%; quiet_in_top: quiet-set cells inside the top 5% (so inside the top 10%) and not already in the original selected_best. Disjoint from the original 30 and each other."""
import argparse, math, os, sys, tomllib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_root import add_results_root_argument, results_root, HERE
import compare_baselines as cb


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    root = results_root(a.results_root); base = os.path.join(root, "compare_baselines")
    cfg = tomllib.load(open(os.path.join(HERE, "baseline_closedloop_config.toml"), "rb"))["followup"]
    L, T, _ = cb.load(root); S, order = cb.build_sets(L, T)
    orig = pd.read_csv(os.path.join(base, "closedloop_cells.csv")); used = set(orig.tag)
    Li = L.set_index("tag"); Ti = T.set_index("tag")
    k = math.ceil(cfg["top_fraction"] * len(L)); rank = {t: i + 1 for i, t in enumerate(order.tag.values)}
    top = set(order.tag.values[:k]); quiet = set(S["B-quiet"])
    k5 = math.ceil(cfg["fallback_top_fraction"] * len(L)); top5 = set(order.tag.values[:k5])
    rng = np.random.default_rng(cfg["seed"]); rows = []; n = cfg["n_per_stratum"]

    def draw(name, pool):
        pool = sorted(t for t in pool if t not in used)
        take = list(rng.choice(pool, size=min(n, len(pool)), replace=False))
        print(f"{name}: {len(pool)} eligible, took {len(take)}" + ("  (fewer than requested: took all)" if len(pool) < n else ""))
        for t in take:
            used.add(t); r = Li.loc[t]
            calib = f"border/grid/{t[6:]}" if abs(r.inc_deg - 50.0) < 1e-6 else f"border2/grid/{t}"
            rows.append(dict(stratum=name, tag=t, inc_deg=r.inc_deg, raan_deg=r.raan_deg, omega_deg=r.omega_deg, epoch=int(r.epoch), calib_dir=calib,
                             passive_status=r.status, open_loop_dp_dv=r.dp_dv, open_loop_near_passive_dv=r.near_passive_dv, open_loop_near_passive_burns=r.near_passive_burns,
                             est_rank=rank[t], est_pred_lazy_dv=float(Ti.loc[t].v1_lazy_dv), est_pred_decayed=bool(Ti.loc[t].v1_decayed), in_top10=t in top, in_quiet=t in quiet))

    draw("quiet_not_top", quiet - top)   # empty: all 24 quiet cells are in the top 10% (reported, not hidden)
    draw("quiet_not_top5", quiet - top5)   # FALLBACK (named): quiet cells outside the estimator top 5%
    draw("estimator_top_not_quiet", top - quiet)
    draw("quiet_in_top", quiet & top5)   # inside the top 5% (hence inside the top 10%); contrasts rank-wise with quiet_not_top5
    D = pd.DataFrame(rows)
    for r in D.itertuples():
        assert os.path.isfile(os.path.join(root, r.calib_dir, "calibration.csv")), r.calib_dir
    D.to_csv(os.path.join(base, "closedloop_cells_followup.csv"), index=False)
    pd.concat([orig, D[orig.columns]], ignore_index=True).to_csv(os.path.join(base, "closedloop_cells_all.csv"), index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30); print(D.drop(columns=["calib_dir"]).to_string()); print("top10 size", k, "quiet size", len(quiet), "quiet in top", len(quiet & top))


if __name__ == "__main__":
    main()
