"""Venus with GRAM truth, part (b): the apoapsis-decay checks D1V and D2V (venus_gram_plan.md sections 2c and 6.3; criteria of the parent plan), with the definitions of
apoapsis_pi_check.py imported, not copied (observable, prediction formula (1), pair filter |pred| >= 5 km). H = the GRAM local scale height of the pass (`h_local_m`, PRIMARY); companion H = the cell-start single-fit H. No altitude-layer filter.
Usage: python3 -I venus_apoapsis_pi_check.py [--results-root DIR] -> <root>/venus/apoapsis_pi/{dev_pairs.csv, criteria.csv, diagnostics.csv}
"""
import argparse, glob, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import apoapsis_pi_check as ck

VC = tomllib.load(open(HERE / "venus_generality_config.toml", "rb"))


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); v = root / "venus"
    out = v / "apoapsis_pi"; out.mkdir(exist_ok=True)
    pc = pd.read_csv(v / "planet_constants.csv").iloc[0]; veh = VC["vehicle"]
    area = veh["bus_dims_m"][0] * veh["bus_dims_m"][2] + 2 * veh["panel_dims_m"][1] * veh["panel_dims_m"][2]
    mu, omega = float(pc.mu_m3s2), abs(float(pc.spin_rad_s)); Hfit = pd.read_csv(v / "cells.csv").set_index("tag").fit_H_m.to_dict()
    rs = pd.read_csv(v / "truth_run_summary.csv").set_index("tag"); rows = []
    for p in sorted(glob.glob(str(v / "truth" / "*" / "passes_3b.csv"))):
        tag = pathlib.Path(p).parent.name; T = pd.read_csv(p); hl = dict(zip(T.alt_m, T.h_local_m))
        r = ck.pairs(p, "Venus", veh["cd"], area, mu, omega, lambda h: hl[h]); rf = ck.pairs(p, "Venus", veh["cd"], area, mu, omega, lambda h: Hfit[tag])
        for x, y in zip(r, rf): x["status"] = "completed" if rs.loc[x["tag"], "retcode"] == "Success" else "decayed"; x["pred_fit_m"] = y["pred_m"]; x["H_fit_m"] = y["H_m"]
        rows += r
    P = pd.DataFrame(rows); P["ratio"] = P.obs_m / P.pred_m; P["ratio_raw"] = P.obs_raw_m / P.pred_m; P["ratio_fit"] = P.obs_m / P.pred_fit_m; P.to_csv(out / "dev_pairs.csv", index=False)
    K = P[P.keep]; crit, diag = [], []
    crit.append(dict(id="D1V", statement="Venus 20x20: median obs/pred", value=float(K.ratio.median()), lo=0.97, hi=1.03))
    crit.append(dict(id="D2V", statement="Venus 20x20: fraction |ratio-1|<=0.05", value=ck.frac_within(K.ratio, 0.05), lo=0.90, hi=1.0))
    Kf = P[P.pred_fit_m.abs() >= ck.MIN_PRED_M]
    diag.append(dict(item="COMPANION (cell-start single-fit H), pairs kept: median ratio, fraction within 5%", value=f"{len(Kf)}, {Kf.ratio_fit.median():.4f}, {ck.frac_within(Kf.ratio_fit, 0.05):.4f}"))
    diag.append(dict(item="H_local of kept pairs (km) min/median/max; H_fit of cells min/median/max", value=" / ".join("%.2f" % (x / 1e3) for x in (K.H_m.min(), K.H_m.median(), K.H_m.max())) + "; " + " / ".join("%.2f" % (x / 1e3) for x in (K.H_fit_m.min(), K.H_fit_m.median(), K.H_fit_m.max()))))
    for c in crit: c["passed"] = bool(c["lo"] <= c["value"] <= c["hi"])
    diag.append(dict(item="pairs total / kept / cells", value=f"{len(P)} / {len(K)} / {K.tag.nunique()}"))
    diag.append(dict(item="ratio p05/p25/p50/p75/p95", value=" / ".join("%.4f" % x for x in K.ratio.quantile([.05, .25, .5, .75, .95]))))
    diag.append(dict(item="raw-apoapsis ratio p05/p50/p95 (secondary)", value=" / ".join("%.4f" % x for x in K.ratio_raw.quantile([.05, .5, .95]))))
    diag.append(dict(item="Laplace c = a e / H min/median/max", value=" / ".join("%.0f" % x for x in (K.c.min(), K.c.median(), K.c.max()))))
    diag.append(dict(item="|pred| (km) min/median/max of kept", value=" / ".join("%.1f" % (x / 1e3) for x in (K.pred_m.abs().min(), K.pred_m.abs().median(), K.pred_m.abs().max()))))
    for st in ("completed", "decayed"):
        g = K[K.status == st]
        if len(g): diag.append(dict(item=f"UNFROZEN {st} cells: n, median ratio, fraction within 5%", value=f"{len(g)}, {g.ratio.median():.4f}, {ck.frac_within(g.ratio, 0.05):.4f}"))
    for lo, hi in ((0.0, 0.84), (0.84, 0.842), (0.842, 0.8435), (0.8435, 1.0)):
        s = K[(K.e >= lo) & (K.e < hi)]
        if len(s): diag.append(dict(item=f"UNFROZEN e in [{lo},{hi}): n, median ratio, within 5%", value=f"{len(s)}, {s.ratio.median():.4f}, {ck.frac_within(s.ratio, 0.05):.4f}"))
    # spread of the pair residual in km (UNFROZEN): obs - pred
    res = (K.obs_m - K.pred_m) / 1e3
    diag.append(dict(item="UNFROZEN residual obs - pred (km): p05/p50/p95, and its rms", value=" / ".join("%.3f" % x for x in res.quantile([.05, .5, .95])) + f"; rms {np.sqrt((res**2).mean()):.3f}"))
    pd.DataFrame(crit).to_csv(out / "criteria.csv", index=False); pd.DataFrame(diag).to_csv(out / "diagnostics.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 140)
    print(pd.DataFrame(crit).to_string(index=False)); print(pd.DataFrame(diag).to_string(index=False))


if __name__ == "__main__":
    main()
