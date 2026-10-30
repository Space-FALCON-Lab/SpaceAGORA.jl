"""Numerical predictor at Mars (option 3), Python half. Plan and thresholds: tasks/2026-10-15_periapsis-region-border-targeting/numerical_predictor_plan.md.

Usage (numpy, pandas; run after run_numerical_predictor.jl --phase=ladder, and after `analyze_mars_generality.py --variant NAME` for every NAME listed):
  python3 -I analyze_numerical_predictor.py --results-root BLIND_ROOT --phase hybrid --base num20x20 [--fractions 0.10,0.20]
        writes pred_numhybF_BASE/, pred_numhybF_BASE_selfcal/, pred_cells_numhybF_BASE.csv (F = percent; refined cells use N, the rest the v3 prediction) from pred_v3* and pred_<base>*;
        refined set = every cell whose v3 predicted dV (inf when predicted decayed) is <= the F-quantile of the cell scores (ties included; the actual fraction is printed)
  python3 -I analyze_numerical_predictor.py --results-root BLIND_ROOT --phase summary [--names num2x0,num3x0,...]
        consistency of N(20x20, truth tolerances) against the truth passes, per-bin error tables (pooled RMS = T5b, and the median over cells of the per-cell RMS),
        pass counts of the Stage C verdict table, selection metrics, the cost columns from num_timing.csv, and the share of the v3 failures N removes.
Reads <root>/mars/{cells.csv, truth/<tag>/passes_3b.csv, pred_<name>/, <name>/{verdict_table,lookahead_error}.csv, v3/verdict_table.csv, num_timing.csv}. Writes <root>/mars/numerical/*.csv.
"""
import argparse, pathlib, shutil, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

DAY = 86400.0
BINS = [0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]       # SOURCED: mars_generality_config_g20.toml [analysis] lookahead_bin_edges_days
N_SCAN = 4 * 22 * 24 * 3                                    # SOURCED: selection_skill_config.toml [scan]: 4 inclinations x 22 RAAN x 24 omega x 3 epochs = 6336
SCAN_BUDGET_CORE_H = 1.0                                    # ASSUMED (plan C3): one core-hour for the full candidate scan
V3_MS_PER_CANDIDATE = 8.9                                   # SOURCED: estimator_v3_results.md section 4 (median per 60-day prediction at Mars, single job)


def hybrid(m, base, fractions):
    T = pd.read_csv(m / "v3" / "selection_cell_table.csv")
    score = np.where(T.v2_decayed.values, np.inf, T.v2_lazy_dv.values.astype(float))
    pc3 = pd.read_csv(m / "pred_cells_v3.csv").set_index("tag"); pcn = pd.read_csv(m / f"pred_cells_{base}.csv").set_index("tag")
    cols = [c for c in pcn.columns if c in pc3.columns]
    for f in fractions:
        thr = np.quantile(score, f); refined = set(T.tag[score <= thr]); name = f"numhyb{int(round(100 * f))}_{base}"
        for suffix in ("", "_selfcal"):
            d = m / f"pred_{name}{suffix}"; d.mkdir(exist_ok=True)
            for tag in T.tag:
                src = m / (f"pred_{base}{suffix}" if tag in refined else f"pred_v3{suffix}") / f"{tag}.csv"
                if src.exists(): shutil.copy(src, d / f"{tag}.csv")
        rows = [(pcn if t in refined else pc3).loc[[t], cols] for t in T.tag if t in pc3.index and t in pcn.index]
        pd.concat(rows).reset_index().to_csv(m / f"pred_cells_{name}.csv", index=False)
        pd.DataFrame(dict(tag=sorted(refined))).to_csv(m / f"{name}_refined_cells.csv", index=False)
        print(f"{name}: refine {len(refined)} of {len(T)} cells ({100 * len(refined) / len(T):.1f}%), v3 score threshold {thr:.4f}")


def per_cell_bins(m, name, status):
    rows = []
    for tag, st in status.items():
        f = m / f"pred_{name}" / f"{tag}.csv"
        if not f.exists(): continue
        tr = pd.read_csv(m / "truth" / tag / "passes_3b.csv"); pr = pd.read_csv(f); n = min(len(tr), len(pr))
        if n < 1: continue
        e = pr.alt_m.values[:n] - tr.alt_m.values[:n]; td = tr.time_s.values[:n] / DAY
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            sel = (td > lo) & (td <= hi) if lo > 0 else (td >= 0) & (td <= hi)
            if sel.any(): rows.append(dict(tag=tag, status=st, bin=f"{lo:g}-{hi:g}", rms_km=float(np.sqrt(np.mean(e[sel] ** 2)) / 1e3), n=int(sel.sum())))
    return pd.DataFrame(rows)


def consistency(m, name, status):
    out = []
    for tag in status:
        f = m / f"pred_{name}" / f"{tag}.csv"
        if not f.exists(): continue
        tr = pd.read_csv(m / "truth" / tag / "passes_3b.csv"); pr = pd.read_csv(f)
        n = min(len(tr), len(pr)); d = pr.alt_m.values[:n] - tr.alt_m.values[:n]
        out.append(dict(tag=tag, n_truth=len(tr), n_pred=len(pr), n_common=n, max_abs_alt_m=float(np.max(np.abs(d))) if n else np.nan,
                        max_abs_time_s=float(np.max(np.abs(pr.t_s.values[:n] - tr.time_s.values[:n]))) if n else np.nan))
    return pd.DataFrame(out)


def closedloop(m, out):
    """Per-stratum closed-loop table: frozen v2 and v3 (run_mars_closedloop.jl) beside N variants (run_numerical_predictor.jl --phase=closedloop)."""
    sets = [("v2_frozen", "closedloop"), ("v3", "closedloop_v3"), ("N 20x20 T1 fixed band", "closedloop_num20x20"), ("N 20x20 T1 per-pass band", "closedloop_num20x20_perpass"),
            ("N 20x20 T2 fixed band", "closedloop_num20x20L")]
    cells = pd.read_csv(m / "closedloop_cells.csv"); rows = []
    for nm, d in sets:
        f = m / d / "closedloop_rows.csv"
        if not f.exists(): continue
        R = pd.read_csv(f).set_index("tag")
        for c in cells.itertuples():
            r = R.loc[c.tag]
            rows.append(dict(predictor=nm, tag=c.tag, stratum=c.stratum, completed=bool(r.completed), total_dv_mps=r.total_dv_mps, burns=int(r.burns), passes=int(r.passes),
                             compliance=r.passes_in_band / max(1, r.passes), peak_q=r.peak_q_w_cm2, exceeds_hard=bool(r.exceeds_hard),
                             predictor_evals=r.get("predictor_evals", np.nan), run_wall_s=r.get("run_wall_s", np.nan)))
    T = pd.DataFrame(rows); T.to_csv(out / "closedloop_cells_num.csv", index=False)
    agg = dict(n=("tag", "size"), completed=("completed", "sum"), dv_median=("total_dv_mps", "median"), dv_mean=("total_dv_mps", "mean"), burns_median=("burns", "median"),
               compliance_median=("compliance", "median"), compliance_min=("compliance", "min"), exceeds_hard=("exceeds_hard", "sum"), evals_median=("predictor_evals", "median"), run_wall_median_s=("run_wall_s", "median"))
    G = T.groupby(["predictor", "stratum"], sort=False).agg(**agg); A = T.groupby("predictor", sort=False).agg(**agg)
    G.to_csv(out / "closedloop_by_stratum_num.csv"); A.to_csv(out / "closedloop_all_num.csv")
    pd.set_option("display.width", 250); print(G.round(3).to_string()); print(A.round(3).to_string())
    P = T.pivot(index="tag", columns="predictor", values="total_dv_mps"); print(P.round(3).to_string())


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--phase", default="summary"); ap.add_argument("--base", default="num20x20")
    ap.add_argument("--fractions", default="0.10,0.20"); ap.add_argument("--names", default=None); a = ap.parse_args()
    root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    if a.phase == "hybrid":
        hybrid(m, a.base, [float(x) for x in a.fractions.split(",")]); return
    out = m / "numerical"; out.mkdir(exist_ok=True)
    if a.phase == "closedloop":
        closedloop(m, out); return
    cells = pd.read_csv(m / "cells.csv"); rs = pd.read_csv(m / "truth_run_summary.csv").set_index("tag")
    status = {t: ("completed" if rs.loc[t, "retcode"] == "Success" else "decayed") for t in cells.tag if t in rs.index and rs.loc[t, "status"] == "ok"}
    names = a.names.split(",") if a.names else sorted(p.name[len("pred_"):] for p in m.glob("pred_num*") if p.is_dir() and not p.name.endswith("_selfcal"))
    names = [n for n in names if (m / n / "verdict_table.csv").exists()]
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
    if "num20x20" in names or (m / "pred_num20x20").exists():
        C = consistency(m, "num20x20", status); C.to_csv(out / "consistency_20x20_T1.csv", index=False)
        print(f"consistency N(20x20,T1) vs truth: {len(C)} cells, max |alt diff| {C.max_abs_alt_m.max():.3g} m, max |time diff| {C.max_abs_time_s.max():.3g} s, "
              f"pass-count mismatches {(C.n_truth != C.n_pred).sum()}")
    V3 = pd.read_csv(m / "v3" / "verdict_table.csv"); v3fail = {(r.criterion, r.scope) for r in V3.itertuples() if not r.passed}
    refs = [("v3", V3)] + [(n, pd.read_csv(m / n / "verdict_table.csv")) for n in names]
    rows = []; bin_rows = []
    for n, VT in refs:
        pooled = VT[VT.scope.isin(["pooled", "completed"])].set_index("criterion"); get = lambda k: (pooled.loc[k, "value"] if k in pooled.index else np.nan)
        fails = {(r.criterion, r.scope) for r in VT.itertuples() if not r.passed}
        removed = len(v3fail - fails)
        T5b = VT[VT.criterion.str.startswith("T5b")].set_index("scope").value
        rows.append(dict(predictor=n, subcriteria=len(VT), passed=int(VT.passed.sum()), v3_failures_removed=removed, v3_failures=len(v3fail), newly_failing=len(fails - v3fail),
                         T1_spearman=get("T1 Spearman > max(B1, B2c)"), T2_recall=get("T2 top-10% recall > k/n"), T3_regret=get("T3 regret <= 0.76 m/s (predicted-complete pick)"),
                         T4_false_safe=get("T4 false-safe, controlled outcome (decay only), Q margin 10%"), T5a_pass1_km=get("T5a pass-1 offset RMS (completed, truth seed)"),
                         T6a_dv_ratio=get("T6a dV median predicted/true"), T6c_dv_spearman=get("T6c dV Spearman"), T7a_acc=get("T7a 3-class accuracy"), T7b_decay_recall=get("T7b decay recall"),
                         T7c_false_decay=get("T7c false-decay rate"), **{f"T5b_{k}": v for k, v in T5b.items()}))
        if n != "v3":
            B = per_cell_bins(m, n, status)
            if len(B):
                for (st, b), g in B.groupby(["status", "bin"]):
                    bin_rows.append(dict(predictor=n, status=st, bin=b, median_cell_rms_km=g.rms_km.median(), p90_cell_rms_km=g.rms_km.quantile(0.9), n_cells=len(g)))
    S = pd.DataFrame(rows); S.to_csv(out / "summary_verdicts.csv", index=False); print(S.round(3).to_string(index=False))
    if bin_rows:
        B = pd.DataFrame(bin_rows); B.to_csv(out / "median_cell_rms_by_bin.csv", index=False)
        for st in ("completed", "decayed"):
            P = B[B.status == st].pivot(index="predictor", columns="bin", values="median_cell_rms_km")[[f"{lo:g}-{hi:g}" for lo, hi in zip(BINS[:-1], BINS[1:])]]
            print(f"\nmedian over {st} cells of the per-cell RMS altitude error (km) by lookahead bin"); print(P.round(3).to_string())
    t = m / "num_timing.csv"
    if t.exists():
        T = pd.read_csv(t); g = T.groupby(["tol", "field"]).agg(median_s=("wall_s", "median"), q1=("wall_s", lambda v: v.quantile(0.25)), q3=("wall_s", lambda v: v.quantile(0.75)),
                                                                  n=("wall_s", "size"), load_start=("load_start", "first"), load_end=("load_end", "first")).reset_index()
        g["scan_core_h"] = g.median_s * N_SCAN / 3600.0; g["ratio_to_v3"] = g.median_s * 1e3 / V3_MS_PER_CANDIDATE
        g["within_budget"] = g.scan_core_h <= SCAN_BUDGET_CORE_H; g.to_csv(out / "cost_table.csv", index=False)
        print(f"\ncost (single job, serial; scan = {N_SCAN} candidates; v3 {V3_MS_PER_CANDIDATE} ms; budget {SCAN_BUDGET_CORE_H} core-h ASSUMED)"); print(g.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
