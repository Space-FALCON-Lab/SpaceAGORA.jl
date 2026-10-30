"""Stage C: Pi tables and the Pi-collapse view (Earth cells and Mars cells together).

Usage: python3 -I analyze_pi_collapse.py [--results-root DIR]
Reads <root>/mars/{pi_cells.csv, planet_constants.csv, body_distances.csv, mars_cell_labels.csv, selection_cell_table.csv}, the Earth truth
passes and labels (<root>/border*, <root>/oracle/oracle_cell_labels.csv) and the Earth v2 predictions (<root>/estimator_v2/pred/<tag>.csv).
Writes <root>/mars/{pi_table.csv, pi_table.md, pi_collapse_cells.csv, pi_collapse_fit.csv, fig_pi_collapse.png}.

Rule fixed in mars_generality_plan.md before any Mars truth run: fit log10(cost in band widths) on log10(Gamma) and log10(Pi_drag) with Earth cells
only (OLS), predict the Mars cells without refitting; "one collapses to the other" if the median |log10 prediction error| on Mars is <= 0.3.
"""
import argparse, importlib.util, pathlib, sys, tomllib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("oracle_lib_pi", HERE / "oracle_lib.py")
lib = importlib.util.module_from_spec(spec); spec.loader.exec_module(lib)
ORC = tomllib.load(open(HERE / "oracle_config.toml", "rb"))
MARGIN = 1000.0
EARTH_MU = lib.MU


def earth_dir(tag, inc, root):
    return root / "border" / "grid" / tag[6:] if inc == 50.0 else root / "border2" / "grid" / tag


def earth_cells(root):
    L = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv"); rows = []
    for c in L.itertuples():
        tr = pd.read_csv(earth_dir(c.tag, c.inc_deg, root) / "passes_3b.csv")
        if c.status == "completed":
            hq = tr.h_qhigh_m.values if "h_qhigh_m" in tr else None
        # shrunk band width: the label table stores the minimum; recompute the mean from the stored band columns where present
        if "h_qhigh_m" in tr:
            w = float(np.mean((tr.h_qlow_m - MARGIN) - (tr.h_qhigh_m + MARGIN)))
            bw = float(np.mean(tr.h_qlow_m - tr.h_qhigh_m))
        else:
            w = float(c.shrunk_width_km_min * 1000.0); bw = w + 2 * MARGIN
        Lk = lib.lever_kepler(tr.a_m.values, tr.e.values)
        cost = c.true_cost_label_mps if c.status == "completed" else np.nan
        pr = root / "estimator_v2" / "pred" / f"{c.tag}.csv"
        rms_err = np.nan
        if c.status == "completed" and pr.exists():
            p = pd.read_csv(pr); n = min(len(p), len(tr)); rms_err = float(np.sqrt(np.mean((p.alt_m.values[:n] - tr.alt_m.values[:n]) ** 2)))
        rows.append(dict(planet="Earth", tag=c.tag, status=c.status, cost_mps=cost, cost_band_widths=cost / (float(np.median(Lk)) * w) if np.isfinite(cost) else np.nan,
                         band_width_m=bw, rms_err_m=rms_err, rms_err_over_band=rms_err / bw if np.isfinite(rms_err) else np.nan))
    return pd.DataFrame(rows)


def mars_cells(root):
    m = root / "mars"; L = pd.read_csv(m / "mars_cell_labels.csv"); T = pd.read_csv(m / "selection_cell_table.csv").set_index("tag"); rows = []
    for c in L.itertuples():
        if c.status == "missing":
            continue
        rms_err = np.nan
        pr = m / "pred" / f"{c.tag}.csv"; tr = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv")
        if c.status == "completed" and pr.exists():
            p = pd.read_csv(pr); n = min(len(p), len(tr)); rms_err = float(np.sqrt(np.mean((p.alt_m.values[:n] - tr.alt_m.values[:n]) ** 2)))
        bw = float(np.mean(tr.h_qlow_m - tr.h_qhigh_m))
        pj = m / "pred_j3" / f"{c.tag}.csv"; rj = np.nan
        if c.status == "completed" and pj.exists():
            p2 = pd.read_csv(pj); n2 = min(len(p2), len(tr)); rj = float(np.sqrt(np.mean((p2.alt_m.values[:n2] - tr.alt_m.values[:n2]) ** 2)))
        rows.append(dict(planet="Mars", rms_err_j3diag_over_band=rj / bw if np.isfinite(rj) else np.nan, tag=c.tag, status=c.status, cost_mps=c.true_cost_label_mps, cost_band_widths=c.cost_band_widths, band_width_m=bw,
                         rms_err_m=rms_err, rms_err_over_band=rms_err / bw if np.isfinite(rms_err) else np.nan))
    return pd.DataFrame(rows)


def pi_tables(root):
    m = root / "mars"; P = pd.read_csv(m / "pi_cells.csv"); pc = pd.read_csv(m / "planet_constants.csv").set_index("planet")
    bd = pd.read_csv(m / "body_distances.csv")
    rows = []
    for planet, g in P.groupby("planet"):
        mu = float(pc.loc[planet, "mu_m3s2"])
        g = g.copy()
        for body in bd[bd.planet == planet].body.unique():
            d = bd[(bd.planet == planet) & (bd.body == body)].set_index("epoch")
            g[f"Pi_3B_{body}"] = [d.loc[int(ep), "mu_b_m3s2"] / mu * (a / d.loc[int(ep), "r_b_m"]) ** 3 for a, ep in zip(g.a_m, g.epoch)]
            g[f"a_over_rb_{body}"] = [a / d.loc[int(ep), "r_b_m"] for a, ep in zip(g.a_m, g.epoch)]
            g[f"nu_{body}"] = [d.loc[int(ep), "n_b_rad_s"] / n for n, ep in zip(g.n_rad_s, g.epoch)]
        rows.append(g)
    G = pd.concat(rows, ignore_index=True)
    G["J3"] = G.planet.map(lambda x: abs(float(pc.loc[x, "J3"]))); G["Re_m"] = G.planet.map(lambda x: float(pc.loc[x, "Re_m"]))
    G["Pi_J3"] = G.J3 * (G.Re_m / (G.a_m * (1 - G.e ** 2))) ** 3          # |J3| (R / p)^3 (framework document, section 3.1)
    G["Gamma_3B_J3"] = G.Pi_3B / G.Pi_J3
    return G


def fmt(x):
    return "n/a" if x is None or not np.isfinite(x) else (f"{x:.3g}")


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    G = pi_tables(root)
    cols = list(dict.fromkeys(["Pi_J2", "Pi_J3", "Pi_3B", "Gamma_3B_J2", "Gamma_3B_J3", "Pi_drag", "Pi_H", "a_over_rb_max"] + [c for c in G.columns if c.startswith(("Pi_3B_", "a_over_rb_", "nu_"))]))
    summ = G.groupby("planet")[cols].agg(["median", "min", "max"])
    summ.to_csv(m / "pi_table.csv")
    # compact markdown: medians with ranges, Earth cells vs Mars cells
    lines = ["| quantity | Earth (336 cells) median [min, max] | Mars (128 cells) median [min, max] |", "|---|---|---|"]
    for c in cols:
        e = G[G.planet == "Earth"][c].dropna(); mm = G[G.planet == "Mars"][c].dropna()
        f = lambda s: "n/a" if len(s) == 0 else f"{s.median():.3g} [{s.min():.3g}, {s.max():.3g}]"
        lines.append(f"| {c} | {f(e)} | {f(mm)} |")
    (m / "pi_table.md").write_text("\n".join(lines) + "\n"); print("\n".join(lines))
    # collapse view
    E = earth_cells(root); M = mars_cells(root)
    C = pd.concat([E, M], ignore_index=True).merge(G[["planet", "tag", "Gamma_3B_J2", "Pi_drag", "Pi_J2", "Pi_J3", "Pi_3B", "Gamma_3B_J3", "Pi_H", "a_over_rb_max"]], on=["planet", "tag"], how="left")
    C.to_csv(m / "pi_collapse_cells.csv", index=False)
    ok = C[(C.status == "completed") & np.isfinite(C.cost_band_widths) & (C.cost_band_widths > 0) & np.isfinite(C.Gamma_3B_J2)]
    fit_rows = []
    for ycol, lab in (("cost_band_widths", "cost (band widths of correction)"), ("rms_err_over_band", "RMS altitude error / band width")):
        d = ok[np.isfinite(ok[ycol]) & (ok[ycol] > 0)]
        e = d[d.planet == "Earth"]; mm = d[d.planet == "Mars"]
        if len(e) < 5 or len(mm) < 3:
            continue
        X = lambda df: np.column_stack([np.ones(len(df)), np.log10(df.Gamma_3B_J2), np.log10(df.Pi_drag)])
        beta, *_ = np.linalg.lstsq(X(e), np.log10(e[ycol]), rcond=None)
        pe = X(e) @ beta; pm = X(mm) @ beta; ym = np.log10(mm[ycol])
        r2e = 1 - np.sum((np.log10(e[ycol]) - pe) ** 2) / np.sum((np.log10(e[ycol]) - np.log10(e[ycol]).mean()) ** 2)
        r2m = 1 - np.sum((ym - pm) ** 2) / np.sum((ym - ym.mean()) ** 2)
        fit_rows.append(dict(quantity=lab, n_earth=len(e), n_mars=len(mm), b0=beta[0], b_logGamma=beta[1], b_logPidrag=beta[2], r2_earth_in_sample=r2e, r2_mars_out_of_sample=r2m,
                             median_abs_log10_error_mars=float(np.median(np.abs(ym - pm))), mean_log10_error_mars=float(np.mean(pm - ym)),
                             median_log10_mars=float(ym.median()), median_log10_earth=float(np.log10(e[ycol]).median()), collapse_rule_met=bool(np.median(np.abs(ym - pm)) <= 0.3)))
    F = pd.DataFrame(fit_rows); F.to_csv(m / "pi_collapse_fit.csv", index=False); print(F.round(3).to_string())
    # model-free summary: the dimensionless medians per planet and the Pi_J3 / Pi_3B view (supplementary, not under the pre-registered rule)
    sm = []
    for pl in ("Earth", "Mars"):
        d = ok[ok.planet == pl]
        sm.append(dict(planet=pl, n_completed=len(d), median_cost_band_widths=d.cost_band_widths.median(), median_rms_err_over_band=d.rms_err_over_band.median(),
                       median_rms_err_j3diag_over_band=d.get("rms_err_j3diag_over_band", pd.Series(dtype=float)).median() if pl == "Mars" else np.nan,
                       median_Pi_J3_over_Pi_3B=float((1 / d.Gamma_3B_J3).median()), median_Gamma_3B_J2=d.Gamma_3B_J2.median(), median_Pi_drag=d.Pi_drag.median()))
    SM = pd.DataFrame(sm); SM.to_csv(m / "pi_collapse_summary.csv", index=False); print(SM.round(4).to_string())
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8))
    for r, (ycol, lab) in enumerate((("cost_band_widths", "dV / (lever x band width)"), ("rms_err_over_band", "RMS altitude error / band width"))):
        for c_, (xcol, xl) in enumerate((("Gamma_3B_J2", "Gamma = Pi_3B / Pi_J2"), ("Pi_drag", "Pi_drag"))):
            ax = axes[r, c_]
            for pl, col in (("Earth", "#0b6e99"), ("Mars", "#c4612b")):
                d = ok[(ok.planet == pl) & np.isfinite(ok[ycol]) & (ok[ycol] > 0)]
                ax.scatter(d[xcol], d[ycol], s=12, c=col, alpha=0.6, label=f"{pl} (n={len(d)})")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(xl); ax.set_ylabel(lab); ax.grid(alpha=0.2); ax.legend(fontsize=8)
    fig.suptitle("Completed cells, Earth and Mars: dimensionless cost and error against Pi_3B/Pi_J2 and Pi_drag (two planets: a thin test)", fontsize=10)
    fig.tight_layout(); fig.savefig(m / "fig_pi_collapse.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
