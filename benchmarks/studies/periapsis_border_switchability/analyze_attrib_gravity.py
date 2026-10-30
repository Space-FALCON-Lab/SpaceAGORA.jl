"""Gravity-term attribution at Mars (step 2/3 of the >= 20x20 study; python3 -I analyze_attrib_gravity.py --results-root R).

Reads <root>/mars/{cells.csv, truth/<tag>/passes_{3b,zonal,deg3,j2only}.csv, attrib_gravity.csv} (run_mars_generality.jl --phase=diag --variants=zonal,deg3,j2only,
then run_estimator_v3.jl --phase=mars_attrib). Writes <root>/mars/{truth_variant_differences.csv, attrib_gravity_summary.csv}.
Part 1 (truth against truth, no estimator): RMS altitude difference over the common passes between gravity variants of the same cell:
  tesseral = full(20x20) - zonal(20x0);  zonal 4-20 = zonal(20x0) - deg3(3x0);  J3 = deg3(3x0) - j2only(2x0);  total vs the old 3x0 truth = full - deg3.
Part 2: RMS altitude error of each predictor against each truth variant (predictors see the matching truncation), by cell outcome.
"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); m = pathlib.Path(results_root(a.results_root)) / "mars"
    cells = pd.read_csv(m / "cells.csv"); rows = []
    for c in cells.itertuples():
        d = {v: pd.read_csv(m / "truth" / c.tag / f"passes_{v}.csv") for v in ("3b", "zonal", "deg3", "j2only")}

        def rms(x, y):
            n = min(len(d[x]), len(d[y])); e = d[x].alt_m.values[:n] - d[y].alt_m.values[:n]
            return float(np.sqrt(np.mean(e ** 2)) / 1e3), float(np.abs(e).max() / 1e3)
        r = dict(tag=c.tag, inc=c.inc, t_end_full=d["3b"].time_s.iloc[-1] / 86400, t_end_zonal=d["zonal"].time_s.iloc[-1] / 86400, t_end_deg3=d["deg3"].time_s.iloc[-1] / 86400, t_end_j2only=d["j2only"].time_s.iloc[-1] / 86400)
        for lab, (x, y) in {"tesseral": ("3b", "zonal"), "zonal4to20": ("zonal", "deg3"), "J3": ("deg3", "j2only"), "full_vs_3x0": ("3b", "deg3")}.items():
            r[lab + "_rms_km"], r[lab + "_max_km"] = rms(x, y)
        rows.append(r)
    T = pd.DataFrame(rows); T.to_csv(m / "truth_variant_differences.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    full_done = T[[f"t_end_{v}" for v in ("full", "zonal", "deg3", "j2only")]].min(axis=1) > 59.0
    cols = [c for c in T.columns if c.endswith("_rms_km")]
    print(f"Part 1: RMS altitude difference between truth variants (km); {int(full_done.sum())} cells completed in all four variants, {len(T)} cells in all")
    print(pd.concat({"completed in all": T[full_done][cols].describe().loc[["min", "50%", "max"]].T, "all cells": T[cols].describe().loc[["min", "50%", "max"]].T}, axis=1).round(3).to_string())
    print(T.groupby("inc")[cols].median().round(2).to_string())
    flips = {v: pd.crosstab(T.t_end_full > 59.5, T[f"t_end_{v}"] > 59.5).reindex(index=[False, True], columns=[False, True], fill_value=0).values.tolist() for v in ("zonal", "deg3", "j2only")}
    print("completion flips against the full truth ([[dec/dec, dec/comp], [comp/dec, comp/comp]] rows full, columns variant):", flips)
    A = pd.read_csv(m / "attrib_gravity.csv"); A["outcome"] = np.where(A.truth_t_end_days > 59.0, "completed", "decayed")
    S = A.groupby(["variant", "outcome", "predictor"]).rms_km.agg(n="count", median="median", mean="mean", p90=lambda x: x.quantile(0.9), max="max", n_above_1km=lambda x: int((x > 1).sum())).reset_index()
    S.to_csv(m / "attrib_gravity_summary.csv", index=False)
    print("Part 2: RMS altitude error (km) of each predictor against each truth variant"); print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
