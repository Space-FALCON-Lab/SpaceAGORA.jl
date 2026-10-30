"""Stage C3: closed-loop cell list for Mars, by the rule fixed in mars_generality_plan.md (written to disk before any closed-loop run).

Usage: python3 -I select_mars_closedloop_cells.py [--results-root DIR] [--table PATH] [--strata selected_best,typical,crit_not_top,odyssey_ref] [--out NAME]
--table: a selection_cell_table.csv whose "v2" columns hold the ranking predictor (default <root>/mars/selection_cell_table.csv; for estimator v3 use <root>/mars/v3/selection_cell_table.csv).
Reads <root>/mars/selection_cell_table.csv (estimator v2 predictions and true labels). Writes <root>/mars/closedloop_cells.csv
(columns tag, stratum, inc, raan, omega, epoch, rank_est, pred_dv, true_dv, true_class). Deterministic, no random draw.

Strata (disjoint, in this order):
  selected_best : the first n_selected_best cells by estimator-v2 rank (predicted lazy_s0 dV; predicted decays last; ties by predicted burns, then tag),
                  using ONLY estimator output ("completed" in the plan is read as predicted-complete: a selection may not use truth).
  typical       : n_typical truly completed cells whose true label is closest to the median true label of all truly completed cells (ties by tag).
  crit_not_top  : up to n_crit_not_top truly completed critical-inclination cells (63.4, 116.6) ranked outside the estimator's top ceil(top_fraction n), lowest tag first.
  odyssey_ref   : up to n_odyssey_ref truly completed cells at the Odyssey inclination closest to the median true label (ties by tag).
"""
import argparse, math, pathlib, sys, tomllib
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
MC = tomllib.load(open(HERE / "mars_generality_config.toml", "rb"))


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--table", default=None); ap.add_argument("--strata", default="selected_best,typical,crit_not_top,odyssey_ref"); ap.add_argument("--out", default="closedloop_cells.csv"); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    T = pd.read_csv(a.table or m / "selection_cell_table.csv"); strata = a.strata.split(","); cl = MC["closedloop"]
    T["pred_dv"] = np.where(T.v2_decayed, np.inf, T.v2_lazy_dv); T["pred_burns"] = np.where(T.v2_decayed, 10 ** 9, T.v2_lazy_burns)
    T = T.sort_values(["pred_dv", "pred_burns", "tag"], kind="stable").reset_index(drop=True); T["rank_est"] = np.arange(1, len(T) + 1)
    n_top = math.ceil(MC["analysis"]["top_fraction"] * len(T))
    comp = T[T.status == "completed"]; med = float(comp.true_cost_mps.median())
    used = set(); rows = []
    def take(df, stratum):
        for r in df.itertuples():
            used.add(r.tag)
            rows.append(dict(tag=r.tag, stratum=stratum, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch=int(r.epoch), rank_est=int(r.rank_est),
                             pred_dv=r.pred_dv if np.isfinite(r.pred_dv) else np.nan, true_dv=r.true_cost_mps, true_class=r.true_class))
    "selected_best" in strata and take(T.head(int(cl["n_selected_best"])), "selected_best")
    c = comp[~comp.tag.isin(used)].assign(d=lambda d: (d.true_cost_mps - med).abs()).sort_values(["d", "tag"], kind="stable")
    "typical" in strata and take(c.head(int(cl["n_typical"])), "typical")
    isc = comp.inc_deg.apply(lambda i: any(abs(i - ci) < 0.05 for ci in MC["analysis"]["critical_inclinations_deg"]))
    c = comp[isc & (comp.rank_est > n_top) & ~comp.tag.isin(used)].sort_values("tag", kind="stable")
    "crit_not_top" in strata and take(c.head(int(cl["n_crit_not_top"])), "crit_not_top")
    isod = comp.inc_deg.apply(lambda i: abs(i - MC["orbit"]["odyssey_inclination_deg"]) < 0.05)
    c = comp[isod & ~comp.tag.isin(used)].assign(d=lambda d: (d.true_cost_mps - med).abs()).sort_values(["d", "tag"], kind="stable")
    "odyssey_ref" in strata and take(c.head(int(cl["n_odyssey_ref"])), "odyssey_ref")
    out = pd.DataFrame(rows); out.to_csv(m / a.out, index=False)
    print(f"median true label (completed) {med:.4f} m/s; estimator top-10% cut: rank <= {n_top}; critical cells outside it: {int((isc & (comp.rank_est > n_top)).sum())}")
    print(out.to_string())


if __name__ == "__main__":
    main()
