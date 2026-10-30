"""UNFROZEN diagnostics of the D2/D4 failures (written after apoapsis_pi_check.py ran; never part of a verdict).
Reads <root>/apoapsis_pi/dev_pairs.csv; writes <root>/apoapsis_pi/diag_spread.csv.
Usage: python3 -I apoapsis_pi_diagnostics.py [--results-root DIR]
"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    out = pathlib.Path(results_root(a.results_root)) / "apoapsis_pi"
    P = pd.read_csv(out / "dev_pairs.csv"); P = P[P.keep].copy()
    P["resid_km"] = (P.obs_m - P.pred_m) / 1e3
    P["set"] = P.planet + np.where(P.file == "passes_deg3.csv", " deg3", " 20x20")
    rows = []
    for s, g in P.groupby("set"):
        cm = g.groupby("tag").ratio.transform("median")            # remove each cell's median: what is pass-to-pass scatter
        rows.append(dict(set=s, n=len(g), resid_km_p05=g.resid_km.quantile(.05), resid_km_p50=g.resid_km.median(), resid_km_p95=g.resid_km.quantile(.95),
                         resid_km_mad=(g.resid_km - g.resid_km.median()).abs().median(),
                         within5=np.mean(abs(g.ratio - 1) <= .05), within5_pred_ge_20km=np.mean(abs(g.ratio[abs(g.pred_m) >= 20e3] - 1) <= .05),
                         within5_after_cell_median=np.mean(abs(g.ratio / cm - 1) <= .05),
                         cell_median_p05=g.groupby("tag").ratio.median().quantile(.05), cell_median_p95=g.groupby("tag").ratio.median().quantile(.95)))
    # pair-level comparison of the same cells and passes, 20x20 minus deg3 residual (tesseral share of the scatter)
    m = P[P.planet == "Mars"].pivot_table(index=["tag", "pass_"], columns="file", values="resid_km").dropna()
    d = m["passes_3b.csv"] - m["passes_deg3.csv"]
    rows.append(dict(set="Mars 20x20 minus deg3, same tag and pass (km)", n=len(d), resid_km_p05=d.quantile(.05), resid_km_p50=d.median(), resid_km_p95=d.quantile(.95),
                     resid_km_mad=(d - d.median()).abs().median()))
    D = pd.DataFrame(rows); D.to_csv(out / "diag_spread.csv", index=False)
    pd.set_option("display.width", 250); print(D.to_string(index=False))
    for s, g in P.groupby("set"):    # residual against predicted size
        g = g.assign(b=pd.cut(abs(g.pred_m) / 1e3, [5, 20, 50, 100, 200, 1e4]))
        print(s); print(g.groupby("b", observed=True).agg(n=("ratio", "size"), within5=("ratio", lambda r: np.mean(abs(r - 1) <= .05)),
                                                         resid_mad_km=("resid_km", lambda r: (r - r.median()).abs().median())).to_string())


if __name__ == "__main__":
    main()
