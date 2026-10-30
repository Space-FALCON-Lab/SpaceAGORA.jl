"""Closed-loop table beside v3 and v4 (python3 -I analyze_v5_closedloop.py --results-root ROOT [--arms closedloop_v3,closedloop_v4,closedloop_v5]).

Per arm and stratum: runs completed of total, median (mean) total dV, median burns, median and minimum compliance (passes in the heat-rate band / passes), runs with a pass above the hard
limit. Criterion G of estimator_v5_plan.md: 12 of 12 complete, no hard breach, stratum median dV within +10% of v3, compliance median >= 0.94 over all runs."""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--arms", default="closedloop_v3,closedloop_v4,closedloop_v5")
a = ap.parse_args(); m = pathlib.Path(results_root(a.results_root)) / "mars"
rows = []; allrows = {}
for arm in a.arms.split(","):
    f = m / arm / "closedloop_rows.csv"
    if not f.exists(): continue
    R = pd.read_csv(f); R["compliance"] = R.passes_in_band / R.passes; allrows[arm] = R
    for st, g in list(R.groupby("stratum")) + [("all", R)]:
        rows.append(dict(arm=arm, stratum=st, completed=f"{int(g.completed.sum())}/{len(g)}", dv_median=g.total_dv_mps.median(), dv_mean=g.total_dv_mps.mean(), burns_median=g.burns.median(),
                         compliance_median=g.compliance.median(), compliance_min=g.compliance.min(), hard_breach_runs=int(g.exceeds_hard.sum()), peak_q_max=g.peak_q_w_cm2.max()))
T = pd.DataFrame(rows); pd.set_option("display.width", 250); print(T.round(3).to_string(index=False)); T.to_csv(m.parent / "v5" / "closedloop_compare.csv", index=False)
if "closedloop_v3" in allrows and "closedloop_v5" in allrows:
    for st in ("selected_best", "typical"):
        v3 = allrows["closedloop_v3"].query("stratum == @st").total_dv_mps.median(); v5 = allrows["closedloop_v5"].query("stratum == @st").total_dv_mps.median()
        print(f"stratum {st}: median dV v5 {v5:.3f} vs v3 {v3:.3f} (limit {1.1 * v3:.3f}): ratio {v5 / v3:.3f}")
    c = allrows["closedloop_v5"].compliance.median(); print(f"compliance median over all runs (v5): {c:.3f} (need >= 0.94)")
