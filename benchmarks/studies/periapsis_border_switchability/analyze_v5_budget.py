"""Error-budget tables of the pass-timing study (python3 -I analyze_v5_budget.py --dir DIR [--planet mars|earth]).

Reads <dir>/{planet}_inject_<stepper>.csv and {planet}_chain_<stepper>.csv written by diag_v5_timing.jl / diag_v5_earth.jl and prints, per force subset (variant) and stepper:
  injection  median |E_full| (s per fresh periapsis-to-periapsis period, fresh start from the reference apoapsis state), median |E_half|, median |dr2| (periapsis radius error of the second pass, m),
             by pass-index bin; and the relative error of the drag-induced period change
  chain      median and 90th percentile of the pass-time error (s) and of the periapsis radius error (m) after j = 5, 10, 20, 40 (Earth 30) passes
Tables go to <dir>/budget_<planet>_*.csv. Diagnostic only; never a verdict."""
import argparse, pathlib
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--planet", default="mars"); a = ap.parse_args()
d = pathlib.Path(a.dir); P = a.planet; bins = [1, 6, 16, 41, 100000]; labels = ["1-5", "6-15", "16-40", "41+"]
nodrag_of = {"full": "nodrag", "drag_notess": "nodrag_notess", "j2drag": "j2only", "drag2b": "point", "zonaldrag": "zonalonly", "j2sundrag": "j2sun_nodrag"}
steppers = [s for s in ("v3", "v4", "v5nokp", "v5", "v5zsp") if (d / f"{P}_inject_{s}.csv").exists()]
pd.set_option("display.width", 300)
rows = []
for s in steppers:
    D = pd.read_csv(d / f"{P}_inject_{s}.csv"); D["bin"] = pd.cut(D.k, bins, right=False, labels=labels)
    for v, g in D.groupby("variant"):
        for b, gb in g.groupby("bin", observed=True):
            rows.append(dict(stepper=s, variant=v, bin=b, n=len(gb), E_full_med_abs=np.median(np.abs(gb.E_full)), E_full_med_signed=np.median(gb.E_full), E_full_p90_abs=np.percentile(np.abs(gb.E_full), 90),
                             E_half_med_abs=np.median(np.abs(gb.E_half)), dr1_med_abs=np.median(np.abs(gb.dr1_m)), dr2_med_abs=np.median(np.abs(gb.dr2_m)), dr2_p90_abs=np.percentile(np.abs(gb.dr2_m), 90)))
    key = D.set_index(["tag", "k", "variant"])
    for (tag, k, v), r in key.iterrows():
        nv = nodrag_of.get(v)
        if nv is not None and (tag, k, nv) in key.index:
            dT = r.T_C - key.loc[(tag, k, nv)].T_C
            if abs(dT) > 1.0: rows.append(dict(stepper=s, variant=v, bin=pd.cut([k], bins, right=False, labels=labels)[0], n=-1, E_full_med_abs=np.nan, E_full_med_signed=100 * r.E_full / dT))
I = pd.DataFrame(rows); I.to_csv(d / f"budget_{P}_inject.csv", index=False)
J = I[I.n > 0]
print(f"== {P}: injection, per fresh period: median |E_full| (s) [signed median] by pass-index bin ==")
for v in ("drag2b", "j2drag", "zonaldrag", "j2sundrag", "tessdrag", "drag_notess", "full", "nodrag", "j2sun_nodrag", "zonalonly", "j2only", "point"):
    g = J[J.variant == v]
    if g.empty: continue
    print(f"  {v:13s}")
    for s in steppers:
        gs = g[g.stepper == s]
        if gs.empty: continue
        print(f"     {s:7s}", " | ".join(f"{r.bin}: {r.E_full_med_abs:7.3f} [{r.E_full_med_signed:7.3f}]" for r in gs.itertuples()))
print(f"\n== {P}: relative error of the drag-induced period change (percent of the change; median signed) ==")
R = I[I.n == -1]
for v in ("drag2b", "j2drag", "zonaldrag", "j2sundrag", "drag_notess", "full"):
    g = R[R.variant == v]
    if g.empty: continue
    for s in steppers:
        gs = g[g.stepper == s]
        if gs.empty: continue
        print(f"  {v:13s} {s:7s}", " | ".join(f"{b}: {gb.E_full_med_signed.median():7.3f}" for b, gb in gs.groupby("bin", observed=True)))
print(f"\n== {P}: injection, periapsis radius error of the second pass |dr2| (m): median [p90] ==")
for v in ("nodrag", "full", "tessdrag", "zonaldrag", "j2drag"):
    g = J[J.variant == v]
    if g.empty: continue
    for s in steppers:
        gs = g[g.stepper == s]
        if gs.empty: continue
        print(f"  {v:10s} {s:7s}", " | ".join(f"{r.bin}: {r.dr2_med_abs:6.1f} [{r.dr2_p90_abs:6.1f}]" for r in gs.itertuples()))
rows = []
for s in steppers:
    f = d / f"{P}_chain_{s}.csv"
    if not f.exists(): continue
    C = pd.read_csv(f)
    for v, g in C.groupby("variant"):
        for j in (5, 10, 20, 30, 40):
            x = g[g.j == j]
            if len(x): rows.append(dict(stepper=s, variant=v, j=j, n=len(x), dt_med_abs=np.median(np.abs(x.dt)), dt_p90_abs=np.percentile(np.abs(x.dt), 90), drp_med_abs=np.median(np.abs(x.drp)), drp_p90_abs=np.percentile(np.abs(x.drp), 90)))
Ch = pd.DataFrame(rows); Ch.to_csv(d / f"budget_{P}_chain.csv", index=False)
print(f"\n== {P}: chain of passes from the reference apoapsis state after truth pass 1: median [p90] |pass-time error| (s) and |periapsis radius error| (m) after j passes ==")
for v in ("drag2b", "j2drag", "zonaldrag", "j2sundrag", "tessdrag", "full", "nodrag", "j2sun_nodrag"):
    g = Ch[Ch.variant == v]
    if g.empty: continue
    print(f"  {v}")
    for s in steppers:
        gs = g[g.stepper == s]
        if gs.empty: continue
        print(f"     {s:7s}", " | ".join(f"j={r.j}: dt {r.dt_med_abs:8.1f} [{r.dt_p90_abs:8.1f}] drp {r.drp_med_abs:6.1f} [{r.drp_p90_abs:7.1f}]" for r in gs.itertuples()))
