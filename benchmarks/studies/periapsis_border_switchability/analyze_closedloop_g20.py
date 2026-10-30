"""Closed-loop check at Mars >= 20x20 (step 5): per-cell and per-stratum dV, burns, completion and band compliance of Arm B with the estimator v3 predictor
and with the frozen v2 predictor, beside the open-loop passive label of each cell.

Usage: python3 -I analyze_closedloop_g20.py --results-root BLIND_ROOT
Reads <root>/mars/{closedloop,closedloop_frozen?}/closedloop_rows.csv (run_mars_closedloop.jl; the default predictor writes closedloop/, --predictor=v3 writes closedloop_v3/),
<root>/mars/closedloop_cells.csv and <root>/mars/mars_cell_labels.csv. Writes <root>/mars/closedloop_g20_summary.csv.
Compliance: passes_in_band / passes (heat-rate band 0.14 to 0.32 W/cm2 at the pass), completed = engine return Success and no controller stop, exceeds_hard = any pass above 0.54.
"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); m = pathlib.Path(results_root(a.results_root)) / "mars"
    cells = pd.read_csv(m / "closedloop_cells.csv"); L = pd.read_csv(m / "mars_cell_labels.csv").set_index("tag"); rows = []
    for pred, d in (("v3", "closedloop_v3"), ("v2_frozen", "closedloop")):
        f = m / d / "closedloop_rows.csv"
        if not f.exists(): continue
        R = pd.read_csv(f).set_index("tag")
        for c in cells.itertuples():
            r = R.loc[c.tag]
            rows.append(dict(predictor=pred, tag=c.tag, stratum=c.stratum, inc=c.inc, pred_ol_dv=c.pred_dv, true_class_passive=L.loc[c.tag, "true_class"], passive_label_mps=L.loc[c.tag, "true_cost_label_mps"],
                             status=r.status, completed=bool(r.completed), t_end_days=r.t_end_days, total_dv_mps=r.total_dv_mps, burns=int(r.burns), passes=int(r.passes), passes_in_band=int(r.passes_in_band),
                             compliance=r.passes_in_band / max(1, r.passes), peak_q=r.peak_q_w_cm2, exceeds_hard=bool(r.exceeds_hard), controller_failures=str(r.controller_failures)[:80]))
    T = pd.DataFrame(rows); T.to_csv(m / "closedloop_g20_summary.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(T.drop(columns=["controller_failures"]).round(3).to_string(index=False))
    g = T.groupby(["predictor", "stratum"]).agg(n=("tag", "size"), completed=("completed", "sum"), dv_median=("total_dv_mps", "median"), dv_mean=("total_dv_mps", "mean"), burns_median=("burns", "median"),
                                                compliance_median=("compliance", "median"), compliance_min=("compliance", "min"), exceeds_hard=("exceeds_hard", "sum"))
    print(g.round(3).to_string()); print(T.groupby("predictor").agg(n=("tag", "size"), completed=("completed", "sum"), dv_median=("total_dv_mps", "median"), compliance_median=("compliance", "median"), exceeds_hard=("exceeds_hard", "sum")).round(3).to_string())


if __name__ == "__main__":
    main()
