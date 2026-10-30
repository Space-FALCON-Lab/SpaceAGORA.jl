"""Write the stratified cell list for run_baseline_closedloop.jl BEFORE any truth run.
python3 -I select_closedloop_cells.py [--results-root R]  -> <R>/compare_baselines/closedloop_cells.csv
Cells come from the 336 truth cells (initial states = the stored calibration seeds). Seed and strata sizes: baseline_closedloop_config.toml."""
import argparse, os, sys, tomllib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_root import add_results_root_argument, results_root, HERE
import compare_baselines as cb


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    root = results_root(a.results_root); out = os.path.join(root, "compare_baselines"); os.makedirs(out, exist_ok=True)
    cfg = tomllib.load(open(os.path.join(HERE, "baseline_closedloop_config.toml"), "rb"))["selection"]
    L, T, _ = cb.load(root); S, order = cb.build_sets(L, T)
    Li = L.set_index("tag"); rng = np.random.default_rng(cfg["seed"]); used = set(); rows = []

    def take(stratum, tags):
        for t in tags:
            assert t not in used, t
            used.add(t); r = Li.loc[t]
            calib = f"border/grid/{t[6:]}" if abs(r.inc_deg - 50.0) < 1e-6 else f"border2/grid/{t}"   # layout of oracle_run.load_cells
            rows.append(dict(stratum=stratum, tag=t, inc_deg=r.inc_deg, raan_deg=r.raan_deg, omega_deg=r.omega_deg, epoch=int(r.epoch), calib_dir=calib,
                             passive_status=r.status, open_loop_dp_dv=r.dp_dv, open_loop_near_passive_dv=r.near_passive_dv, open_loop_near_passive_burns=r.near_passive_burns))

    def draw(pool, n):
        pool = sorted(t for t in pool if t not in used)
        return list(rng.choice(pool, size=n, replace=False))

    comp = L[L.status == "completed"]
    take("selected_best", [t for t in order.tag.values if Li.loc[t].status == "completed" and t not in used][:cfg["n_selected_best"]])
    take("worst_completed", list(comp.sort_values(["dp_dv", "tag"], ascending=[False, True]).tag.values[:cfg["n_worst_completed"]]))
    take("median_cost", draw(S["B-median"], cfg["n_median_cost"]))
    take("ref50", draw(S["B-ref50"], cfg["n_ref50"]))
    take("decayed_passive", draw(L.tag[L.status == "decayed"], cfg["n_decayed_passive"]))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(out, "closedloop_cells.csv"), index=False)
    for r in D.itertuples():   # initial states must exist
        assert os.path.isfile(os.path.join(root, r.calib_dir, "calibration.csv")), r.calib_dir
    print(D.to_string())


if __name__ == "__main__":
    main()
