"""Stage C2 diagnosis of the Mars failures (frozen estimator): which term, which cells.

Usage: python3 -I analyze_mars_diag.py [--results-root DIR]
Part 1 (all 128 cells, existing runs): RMS altitude error of the FROZEN estimator against the J3 exposure factor of the secular J3 rate of the eccentricity,
  g = sin(i) (1 - (5/4) sin^2 i) cos(omega)   (Kozai/Brouwer long-period J3 term, de/dt ~ -1.5 n J3 (R/p)^3 g; zero at the critical inclinations and at omega = 90, 270).
Part 2 (the 12 diagnostic cells, truth variants): estimator counterparts of the truth variants, RMS altitude error over the common passes
  full        frozen estimator vs the full truth (J2+J3, Sun, drag)
  j2only      frozen estimator vs truth with J3 removed (gravity degree 2): both lack J3, the estimator should match if J3 is the cause
  nosun       frozen estimator with the tide off vs truth without the Sun
  nodrag      frozen estimator with drag off vs truth without drag
  j3diag      UNFROZEN J3 estimator vs the full truth
Writes <root>/mars/{diag_table.csv, diag_summary.csv, j3_exposure.csv, fig_diag_j3.png}.
"""
import argparse, pathlib, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def rms_err(tr, pr):
    n = min(len(tr), len(pr))
    if n == 0:
        return np.nan, np.nan, 0
    e = pr.alt_m.values[:n] - tr.alt_m.values[:n]
    return float(np.sqrt(np.mean(e ** 2)) / 1e3), float(np.max(np.abs(e)) / 1e3), n


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    pd.set_option("display.width", 250)
    T = pd.read_csv(m / "selection_cell_table.csv"); rows = []
    for c in T.itertuples():
        tr = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv"); pr = pd.read_csv(m / "pred" / f"{c.tag}.csv"); pj = pd.read_csv(m / "pred_j3" / f"{c.tag}.csv")
        i = np.radians(c.inc_deg); w = np.radians(c.omega_deg); g = np.sin(i) * (1 - 1.25 * np.sin(i) ** 2) * np.cos(w)
        r, mx, n = rms_err(tr, pr); rj, mj, nj = rms_err(tr, pj)
        rows.append(dict(tag=c.tag, status=c.status, inc_deg=c.inc_deg, omega_deg=c.omega_deg, J3_exposure_g=g, abs_g=abs(g), rms_frozen_km=r, max_frozen_km=mx, n=n, rms_j3diag_km=rj, max_j3diag_km=mj))
    J = pd.DataFrame(rows); J.to_csv(m / "j3_exposure.csv", index=False)
    print("Spearman(|g|, frozen RMS error), all cells:", stats.spearmanr(J.abs_g, J.rms_frozen_km)[0], " completed:", stats.spearmanr(J[J.status == "completed"].abs_g, J[J.status == "completed"].rms_frozen_km)[0])
    J["g_class"] = pd.cut(J.abs_g, [-1e-9, 1e-6, 0.05, 0.15, 10], labels=["g=0 (critical inc. or omega=90/270)", "|g|<0.05", "0.05-0.15", ">0.15"])
    print(J.groupby("g_class", observed=True).agg(n=("tag", "size"), rms_frozen_km_median=("rms_frozen_km", "median"), rms_frozen_km_max=("rms_frozen_km", "max"),
                                                  rms_j3diag_km_median=("rms_j3diag_km", "median"), rms_j3diag_km_max=("rms_j3diag_km", "max")).round(3))
    # part 2
    D = pd.read_csv(m / "diag_cells.csv"); rows = []
    for c in D.itertuples():
        for variant, tfile, pdir in (("full", "passes_3b.csv", m / "pred_diag" / "full"), ("j2only", "passes_j2only.csv", m / "pred_diag" / "j2only"),
                                     ("nosun", "passes_nosun.csv", m / "pred_diag" / "nosun"), ("nodrag", "passes_nodrag.csv", m / "pred_diag" / "nodrag"),
                                     ("j3diag_vs_full", "passes_3b.csv", m / "pred_j3")):
            tr = pd.read_csv(m / "truth" / c.tag / tfile); pr = pd.read_csv(pdir / f"{c.tag}.csv"); r, mx, n = rms_err(tr, pr)
            ch_t = (tr.alt_m.values[min(len(tr), len(pr)) - 1] - tr.alt_m.values[0]) / 1e3 if min(len(tr), len(pr)) else np.nan
            rows.append(dict(tag=c.tag, inc_deg=c.inc, omega_deg=c.omega, variant=variant, rms_km=r, max_km=mx, n_common=n, n_truth=len(tr), n_pred=len(pr), truth_alt_change_km=ch_t))
    X = pd.DataFrame(rows); X.to_csv(m / "diag_table.csv", index=False)
    print(X.pivot(index="tag", columns="variant", values="rms_km").round(3).to_string())
    S = X.groupby("variant").rms_km.agg(["count", "median", "max"]).round(3); S.to_csv(m / "diag_summary.csv"); print(S)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for st, col in (("completed", "#0b6e99"), ("decayed", "#c4612b")):
        g = J[J.status == st]
        ax[0].scatter(g.abs_g + 1e-4, g.rms_frozen_km, s=14, c=col, alpha=0.7, label=f"frozen estimator ({st})"); ax[0].scatter(g.abs_g + 1e-4, g.rms_j3diag_km, s=14, c=col, marker="x", alpha=0.7, label=f"UNFROZEN +J3 ({st})")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("|J3 exposure g| (+1e-4)"); ax[0].set_ylabel("RMS altitude error (km)"); ax[0].legend(fontsize=6); ax[0].grid(alpha=0.2)
    piv = X.pivot(index="tag", columns="variant", values="rms_km")
    piv = piv.loc[D.tag]; xs = np.arange(len(piv)); w = 0.16
    for k, (v, col) in enumerate((("full", "#c4612b"), ("j2only", "#0b6e99"), ("nosun", "#5aa9c9"), ("nodrag", "#7f7f7f"), ("j3diag_vs_full", "#2a7f3f"))):
        ax[1].bar(xs + k * w, piv[v].values, w, label=v, color=col)
    ax[1].set_yscale("log"); ax[1].set_xticks(xs + 2 * w); ax[1].set_xticklabels(piv.index, rotation=90, fontsize=6); ax[1].set_ylabel("RMS altitude error (km)"); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.2, axis="y")
    ax[1].set_title("12 diagnostic cells: frozen estimator vs matching truth variant", fontsize=9)
    fig.tight_layout(); fig.savefig(m / "fig_diag_j3.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
