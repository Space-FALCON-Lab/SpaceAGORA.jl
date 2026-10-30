"""Pass/fail table of the pre-registered criteria of estimator_v5_plan.md for one Mars results root (python3 -I analyze_v5_criteria.py --results-root ROOT --tag TAG [--set dev|blind|seen]).

Needs, per root, the outputs of `analyze_mars_generality.py --variant V` for V in v3, v4, v5, num20x20L (<root>/mars/V/verdict_table.csv) and of
`analyze_estimator_v5.py --tag TAG` (<root>/v5/compare_TAG_{counts,metrics,passtime}.csv). Criteria (thresholds from the plan, written before the blind runs):
A  T5a <= 0.25 km and T5b <= 1.0 km in all seven bins (pooled, completed cells)
B  sub-criteria passed >= 44 of 56 (blind sets) / >= 36 of 46 (development)
C  no sub-criterion that v3 passes is failed by v5 (against v4 also reported)
D  every cell with a truth commensurability crossing has an RMS pass-time error <= 1,000 s over the 5 passes after the crossing (fraction reported)
TT1 pooled pass-time RMS <= 1,000 s in the bins 0-5, 5-10, 10-20 d
TT2 pooled pass-time RMS of v5 <= 0.5 x that of v4 in each of the seven bins"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

BINS = ["0_5", "5_10", "10_20", "20_30", "30_40", "40_50", "50_60"]


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--tag", required=True); ap.add_argument("--set", default="blind", choices=["dev", "blind", "seen"])
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"; v5d = root / "v5"
    VT = {v: pd.read_csv(m / v / "verdict_table.csv") for v in ("v3", "v4", "v5", "num20x20L") if (m / v / "verdict_table.csv").exists()}
    key = lambda t: t.criterion + " | " + t.scope
    ok = {v: set(key(t)[t.passed.astype(bool)]) for v, t in VT.items()}
    n_tot = len(VT["v5"]); need = 36 if a.set == "dev" else 44
    M = pd.read_csv(v5d / f"compare_{a.tag}_metrics.csv").set_index("variant"); P = pd.read_csv(v5d / f"compare_{a.tag}_passtime.csv").set_index("variant")
    rows = []
    t5a = M.loc["v5", "T5a"]; t5b = [M.loc["v5", f"T5b_{b}"] for b in BINS]
    rows.append(("A accuracy", f"T5a {t5a:.3f} km; T5b " + " / ".join(f"{x:.2f}" for x in t5b) + " km", bool(t5a <= 0.25 and all(x <= 1.0 for x in t5b))))
    rows.append(("B selection", f"{len(ok['v5'])} of {n_tot} (need {need}); v3 {len(ok['v3'])}, v4 {len(ok['v4'])}, N20 {len(ok['num20x20L'])}", len(ok["v5"]) >= need))
    lost3 = sorted(ok["v3"] - ok["v5"]); lost4 = sorted(ok["v4"] - ok["v5"])
    rows.append(("C no loss vs v3", f"lost {len(lost3)}: " + "; ".join(lost3) + f"  | vs v4 lost {len(lost4)}", len(lost3) == 0))
    fr = P.loc["v5", "frac_cross_cells_le_1000s"]
    rows.append(("D crossings <= 1000 s", f"fraction of cells {fr:.2f}; median worst-crossing RMS {P.loc['v5', 'cross_dt_rms_median_s']:.0f} s (v4 {P.loc['v4', 'cross_dt_rms_median_s']:.0f})", bool(fr >= 1.0)))
    dt5 = [P.loc["v5", f"dt_{b}"] for b in BINS]; dt4 = [P.loc["v4", f"dt_{b}"] for b in BINS]
    rows.append(("TT1 pass time <= 1000 s (0-20 d)", "v5 " + " / ".join(f"{x:.0f}" for x in dt5[:3]) + " s (v4 " + " / ".join(f"{x:.0f}" for x in dt4[:3]) + ")", all(x <= 1000 for x in dt5[:3])))
    ratios = [x / y for x, y in zip(dt5, dt4)]
    rows.append(("TT2 v5 <= 0.5 x v4 in every bin", "ratios " + " / ".join(f"{r:.2f}" for r in ratios), all(r <= 0.5 for r in ratios)))
    T = pd.DataFrame(rows, columns=["criterion", "value", "pass"]); pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 200)
    print(f"== {a.tag} ({a.set}) =="); print(T.to_string(index=False)); T.to_csv(v5d / f"criteria_{a.tag}.csv", index=False)
    print("\npass-time RMS by bin (s):  v3 / v4 / v5 / N20x20 / oracle-time variants where present")
    for v in P.index: print(f"  {v:16s}", " ".join(f"{P.loc[v, f'dt_{b}']:>9.0f}" for b in BINS))
    print("\nT5b by bin (km):")
    for v in M.index: print(f"  {v:16s}", " ".join(f"{M.loc[v, f'T5b_{b}']:>9.3f}" for b in BINS), f"| T5a {M.loc[v, 'T5a']:.3f}")


if __name__ == "__main__":
    main()
