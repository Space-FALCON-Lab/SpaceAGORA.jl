"""Estimator v4 scoring beside v3 and the numerical predictor (python3 -I analyze_estimator_v4.py --results-root ROOT [--variants v3,v4,...] [--tag NAME]).

Reads, for each variant V, <root>/mars/V/verdict_table.csv (written by `analyze_mars_generality.py --variant V`, the Stage C T-list) and the prediction folder
<root>/mars/pred_V (truth seed). Writes <root>/v4/compare_<tag>_*.csv and prints:
  1  pass counts of the Stage C sub-criteria (pooled + folds) per variant, and the sub-criteria v4 gains / loses against v3
  2  the pooled selection metrics (T1-T4, T6-T7) and T5a/T5b/T5c per lookahead bin
  3  median and 90th percentile over completed cells of the per-cell RMS altitude error per lookahead bin (truth seed)
  4  pass-time error: all passes per bin, and at commensurability crossings of the TRUTH history (m nu crosses an integer, m = 2, 3, 4)
  5  pass-1 offset, decayed-cell statistics
Spin rate for nu: 7.08823596e-5 rad/s (SOURCED: planets.jl Mars struct, as tesseral_param.py).
"""
import argparse, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

SPIN = 7.08823596e-5
BINS = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60)]
DAY = 86400.0


def counts(vt):
    ok = vt.passed.astype(bool)
    return int(ok.sum()), len(vt)


def percell(root, var, labels):
    m = root / "mars"; rows = []
    for c in labels.itertuples():
        tf = m / "truth" / c.tag / "passes_3b.csv"; pf = m / f"pred_{var}" / f"{c.tag}.csv"
        if not (tf.exists() and pf.exists()): continue
        tr = pd.read_csv(tf); p = pd.read_csv(pf); n = min(len(tr), len(p))
        if n == 0: continue
        e = p.alt_m.values[:n] - tr.alt_m.values[:n]; dt = p.t_s.values[:n] - tr.time_s.values[:n]; t = tr.time_s.values[:n] / DAY
        nu = SPIN * np.diff(tr.time_s.values) / (2 * math.pi)
        row = dict(tag=c.tag, status=c.status, n_truth=len(tr), n_pred=len(p), pass1_err_m=float(e[0]), t_end_d=float(tr.time_s.values[-1] / DAY))
        for lo, hi in BINS:
            s = (t > lo) & (t <= hi) if lo > 0 else (t <= hi)
            row[f"rms_{lo}_{hi}"] = float(np.sqrt(np.mean(e[s] ** 2)) / 1e3) if s.any() else np.nan
            row[f"dt_rms_{lo}_{hi}"] = float(np.sqrt(np.mean(dt[s] ** 2))) if s.any() else np.nan
        # crossings of the truth commensurability: m nu_k passes an integer between pass k and k+1
        cross = []
        for m_ in (2, 3, 4):
            x = m_ * nu
            for k in range(len(x) - 1):
                if math.floor(x[k]) != math.floor(x[k + 1]): cross.append((k + 1, m_, int(max(math.floor(x[k]), math.floor(x[k + 1])))))   # pass index (1-based) after which nu has crossed
        row["n_cross"] = len(cross); errs = []
        for kc, m_, j in cross:
            lo_i, hi_i = kc, min(kc + 5, n)       # passes kc+1 .. kc+5 (0-based kc .. kc+4)
            if hi_i - lo_i >= 3: errs.append(float(np.sqrt(np.mean(dt[lo_i:hi_i] ** 2))))
        row["cross_dt_rms_s"] = float(np.max(errs)) if errs else np.nan      # worst crossing of the cell
        row["cross_dt_rms_mean_s"] = float(np.mean(errs)) if errs else np.nan
        row["cross_dt_before_s"] = float(np.mean([abs(dt[max(kc - 2, 0)]) for kc, _, _ in cross if kc - 2 < n])) if cross else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--variants", default="v3,v4,num20x20L"); ap.add_argument("--tag", default="dev"); ap.add_argument("--ref", default="v3")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); (root / "v4").mkdir(exist_ok=True)
    vs = a.variants.split(","); labels = pd.read_csv(root / "mars" / "mars_cell_labels.csv")
    VT = {v: pd.read_csv(root / "mars" / v / "verdict_table.csv") for v in vs if (root / "mars" / v / "verdict_table.csv").exists()}
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    # 1
    rows = []
    for v, t in VT.items():
        n, tot = counts(t); pooled = t[t.scope == "pooled"]
        rows.append(dict(variant=v, passed=n, of=tot, pooled_passed=int(pooled.passed.sum()), pooled_of=len(pooled)))
    C = pd.DataFrame(rows); print("== 1. Stage C sub-criteria passed =="); print(C.to_string(index=False)); C.to_csv(root / "v4" / f"compare_{a.tag}_counts.csv", index=False)
    ref = VT.get(a.ref)
    if ref is not None:
        key = lambda t: t.criterion + " | " + t.scope
        rk = set(key(ref)[ref.passed.astype(bool)])
        for v, t in VT.items():
            if v == a.ref: continue
            tk = set(key(t)[t.passed.astype(bool)]); print(f"-- {v} vs {a.ref}: gained {len(tk - rk)}, lost {len(rk - tk)}")
            for x in sorted(rk - tk): print("   LOST  ", x)
    # 2
    out = []
    for v, t in VT.items():
        g = lambda crit, scope="pooled": (lambda r: float(r.value.iloc[0]) if len(r) else np.nan)(t[(t.criterion.str.startswith(crit)) & (t.scope == scope)])
        r = dict(variant=v, T1=g("T1"), T2=g("T2"), T3=g("T3"), T4=g("T4 "), T5a=g("T5a"), T6a=g("T6a", "completed"), T6b=g("T6b", "completed"), T6c=g("T6c", "completed"), T7a=g("T7a"), T7b=g("T7b"), T7c=g("T7c"))
        for lo, hi in BINS:
            sc = f"({lo if lo else -0.001}, {float(hi)}]" if lo == 0 else f"({float(lo)}, {float(hi)}]"
            r[f"T5b_{lo}_{hi}"] = g("T5b", sc)
        out.append(r)
    M = pd.DataFrame(out); print("== 2. pooled metrics and T5b (pooled RMS altitude error, km, completed cells) =="); print(M.round(3).T.to_string()); M.to_csv(root / "v4" / f"compare_{a.tag}_metrics.csv", index=False)
    # 3, 4
    PCs = {v: percell(root, v, labels) for v in vs if (root / "mars" / f"pred_{v}").exists()}
    rows = []
    for v, P in PCs.items():
        P = P[P.status == "completed"]
        r = dict(variant=v, n_completed=len(P))
        for lo, hi in BINS:
            r[f"med_{lo}_{hi}"] = P[f"rms_{lo}_{hi}"].median(); r[f"p90_{lo}_{hi}"] = P[f"rms_{lo}_{hi}"].quantile(0.9)
        rows.append(r)
    Q = pd.DataFrame(rows); print("== 3. per-cell RMS altitude error (km), completed cells: median / 90th percentile =="); print(Q.round(3).T.to_string()); Q.to_csv(root / "v4" / f"compare_{a.tag}_percell_rms.csv", index=False)
    rows = []
    for v, P in PCs.items():
        P = P[P.status == "completed"]; x = P[P.n_cross > 0]; y = P[P.n_cross == 0]
        r = dict(variant=v, n_cells_cross=len(x), n_cells_nocross=len(y), cross_dt_rms_median_s=x.cross_dt_rms_s.median(), cross_dt_rms_p90_s=x.cross_dt_rms_s.quantile(0.9),
                 frac_cross_cells_le_1000s=float((x.cross_dt_rms_s <= 1000).mean()) if len(x) else np.nan, cross_dt_before_median_s=x.cross_dt_before_s.median())
        for lo, hi in BINS: r[f"dt_{lo}_{hi}"] = float(np.sqrt(np.nanmean(P[f"dt_rms_{lo}_{hi}"] ** 2)))
        rows.append(r)
    T = pd.DataFrame(rows); print("== 4. pass-time error (s): rms over completed cells per bin, and the worst-crossing RMS over the 5 passes after a truth commensurability crossing =="); print(T.round(0).T.to_string()); T.to_csv(root / "v4" / f"compare_{a.tag}_passtime.csv", index=False)
    for v, P in PCs.items(): P.to_csv(root / "v4" / f"percell_{a.tag}_{v}.csv", index=False)
    # 5
    rows = []
    for v, P in PCs.items():
        c = P[P.status == "completed"]; d = P[P.status != "completed"]
        rows.append(dict(variant=v, pass1_rms_m_completed=float(np.sqrt(np.mean(c.pass1_err_m ** 2))), n_decayed=len(d),
                         decayed_rms_0_20=float(np.sqrt(np.nanmean(np.r_[d.rms_0_5, d.rms_5_10, d.rms_10_20] ** 2))) if len(d) else np.nan,
                         passcount_diff_median=float((P.n_pred - P.n_truth).median())))
    print("== 5 =="); print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
