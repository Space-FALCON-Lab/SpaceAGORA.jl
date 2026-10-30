"""Revision 17 P3 in its registered form: prediction residual = measured ln q of pass j+1 minus the one-pass prediction z0 made before it (replan.csv: lnq_last_measured of the next row minus z0_pred_lnq of this row), median over passes 20..n per flight; plus the shift of y against the same seed at L0.
Usage: python3 -I planner_rev17_pred.py --planner-root R/planner --out F"""
import argparse, glob, json, pathlib
import numpy as np, pandas as pd


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); root = pathlib.Path(a.planner_root); rows = []
    for d in sorted(glob.glob(str(root / "rev17" / "ODY_PL_*"))):
        d = pathlib.Path(d); nm = d.name.split("_"); arm = nm[2]; lv = nm[-1]; seed = int(nm[5].lstrip("s")); R = pd.read_csv(d / "replan.csv"); P = pd.read_csv(d / "passes.csv"); M = pd.read_csv(d / "measured.csv")
        z = R.z0_pred_lnq.values.astype(float); lm = R.lnq_last_measured.values.astype(float); res = lm[1:] - z[:-1]; res = res[19:]; res = res[np.isfinite(res)]
        rows.append(dict(arm=arm, level=lv, seed=seed, resid_med=float(np.median(res)) if len(res) else float("nan"), resid_mean=float(np.mean(res)) if len(res) else float("nan"), y_med=float(np.nanmedian(M.y.values[4:])), n=len(R)))
    D = pd.DataFrame(rows); L0 = D[D.level == "L0"].set_index(["arm", "seed"]).y_med; D["dy_vs_L0"] = [r.y_med - L0.get((r.arm, r.seed), np.nan) for r in D.itertuples()]
    g = D.groupby(["arm", "level"]).agg(n=("seed", "size"), resid_med=("resid_med", "median"), resid_abs_med=("resid_med", lambda s: float(np.median(np.abs(s)))), dy_vs_L0=("dy_vs_L0", "median")).reset_index()
    pathlib.Path(a.out).write_text(g.to_string(index=False)); print(g.round(3).to_string(index=False)); D.to_csv(str(pathlib.Path(a.out).with_suffix(".csv")), index=False)


if __name__ == "__main__":
    main()
