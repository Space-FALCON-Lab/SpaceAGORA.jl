"""Earth altitude and pass-time error of estimator versions against the truth (python3 -I analyze_earth_v5.py --results-root ROOT [--variants estimator_v3,estimator_v5] [--seed pred|pred_selfcal]).

Reads the Earth cell list <root>/oracle/oracle_cell_labels.csv, the truth pass tables (<root>/border/grid/<tag minus 'I0500_'> for the 50 deg cells, <root>/border2/grid/<tag> otherwise)
and, for each variant folder <root>/<variant>/<seed>/<tag>.csv, the predicted pass table. Prints, per lookahead bin (0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60 d) and for completed cells and all cells:
pooled RMS altitude error (km), pooled RMS pass-time error (s), and the per-cell medians; the pass-count difference; and the cells where a variant changes the altitude history of another
(byte-identity of the prediction files is checked separately with cmp). Writes <root>/v5/earth_compare_<tag>.csv."""
import argparse, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

BINS = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60)]
DAY = 86400.0


def truth_dir(root, tag, inc):
    return root / "border" / "grid" / tag[7:] if inc == 50.0 else root / "border2" / "grid" / tag


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--variants", default="estimator_v3,estimator_v5"); ap.add_argument("--seed", default="pred"); ap.add_argument("--tag", default="dev")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); (root / "v5").mkdir(exist_ok=True)
    labels = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv"); vs = a.variants.split(",")
    rows = []
    for c in labels.itertuples():
        tf = truth_dir(root, c.tag, c.inc_deg) / "passes_3b.csv"
        if not tf.exists(): continue
        tr = pd.read_csv(tf)
        for v in vs:
            pf = root / v / a.seed / f"{c.tag}.csv"
            if not pf.exists(): continue
            p = pd.read_csv(pf); n = min(len(tr), len(p))
            if n == 0: continue
            e = p.alt_m.values[:n] - tr.alt_m.values[:n]; dt = p.t_s.values[:n] - tr.time_s.values[:n]; t = tr.time_s.values[:n] / DAY
            for lo, hi in BINS:
                s = (t > lo) & (t <= hi) if lo > 0 else (t <= hi)
                if s.any():
                    rows.append(dict(tag=c.tag, status=c.status, variant=v, bin=f"{lo}-{hi}", n=int(s.sum()), se2=float((e[s] ** 2).sum()), sdt2=float((dt[s] ** 2).sum()),
                                     rms_km=float(np.sqrt(np.mean(e[s] ** 2)) / 1e3), dt_rms_s=float(np.sqrt(np.mean(dt[s] ** 2)))))
            rows.append(dict(tag=c.tag, status=c.status, variant=v, bin="count", n=len(p) - len(tr), se2=np.nan, sdt2=np.nan, rms_km=np.nan, dt_rms_s=np.nan))
    D = pd.DataFrame(rows); D.to_csv(root / "v5" / f"earth_compare_{a.tag}.csv", index=False)
    pd.set_option("display.width", 250)
    for scope in ("completed", "all"):
        d = D if scope == "all" else D[D.status == "completed"]
        print(f"== {scope} cells: pooled RMS altitude error (km) / pooled RMS pass-time error (s) / per-cell median altitude RMS (km) / per-cell median |pass-time RMS| (s) by lookahead bin ==")
        out = []
        for v in vs:
            for b in [f"{lo}-{hi}" for lo, hi in BINS]:
                g = d[(d.variant == v) & (d.bin == b)]
                if len(g) == 0: continue
                out.append(dict(variant=v, bin=b, cells=len(g), pooled_alt_km=math.sqrt(g.se2.sum() / g.n.sum()) / 1e3, pooled_dt_s=math.sqrt(g.sdt2.sum() / g.n.sum()), med_alt_km=g.rms_km.median(), med_dt_s=g.dt_rms_s.median()))
        print(pd.DataFrame(out).round(3).to_string(index=False))
        for v in vs:
            g = d[(d.variant == v) & (d.bin == "count")]
            print(f"  {v}: pass-count difference (pred - truth): median {g.n.median()}, mean {g.n.mean():.2f}, cells {len(g)}")
    print("\n== cells where the first variant is worse than the second in the 40-60 d bins (per-cell altitude RMS) ==")
    for b in ("40-50", "50-60"):
        x = D[(D.bin == b) & (D.status == "completed")].pivot(index="tag", columns="variant", values="rms_km").dropna()
        if len(vs) >= 2 and len(x):
            r = x[vs[1]] / x[vs[0]]
            print(f"  {b}: cells {len(x)}; ratio {vs[1]}/{vs[0]} median {r.median():.3f}, share improved {float((r < 1).mean()):.2f}, share worse by >10% {float((r > 1.1).mean()):.2f}")


if __name__ == "__main__":
    main()
