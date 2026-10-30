"""Stage C2 diagnostics: choose the attribution subset (rule in mars_generality_plan.md, "Diagnostic budget").

Usage: python3 -I select_mars_diag_cells.py [--results-root DIR] [--n 12]
Rule: the n truly completed cells with the largest RMS altitude error of the frozen estimator v2 (truth seed) over the common passes; spread guard: at most
3 cells per inclination so that every geometry family is represented. Selection uses Mars truth and is for DIAGNOSIS only (it never changes a verdict).
Writes <root>/mars/diag_cells.csv (columns tag, inc, raan, omega, epoch_idx, rms_err_m).
"""
import argparse, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--n", type=int, default=12)
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    T = pd.read_csv(m / "selection_cell_table.csv"); rows = []
    for c in T[T.status == "completed"].itertuples():
        tr = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv"); pr = pd.read_csv(m / "pred" / f"{c.tag}.csv"); n = min(len(tr), len(pr))
        if n:
            rows.append(dict(tag=c.tag, inc=c.inc_deg, raan=c.raan_deg, omega=c.omega_deg, epoch_idx=int(c.epoch), rms_err_m=float(np.sqrt(np.mean((pr.alt_m.values[:n] - tr.alt_m.values[:n]) ** 2)))))
    D = pd.DataFrame(rows).sort_values("rms_err_m", ascending=False, kind="stable")
    keep = []; per = {}
    for r in D.itertuples():
        if per.get(r.inc, 0) < 3:
            keep.append(r.Index); per[r.inc] = per.get(r.inc, 0) + 1
        if len(keep) == a.n:
            break
    out = D.loc[keep]; out.to_csv(m / "diag_cells.csv", index=False); print(out.to_string())


if __name__ == "__main__":
    main()
