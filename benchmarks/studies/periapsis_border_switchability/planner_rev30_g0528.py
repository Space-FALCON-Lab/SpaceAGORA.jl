"""Revision 30.4: G0528 position, rev23geo (OLRT3-W3, EI 2000) against rev26val2 (OLRT3-W3L, EI 150), seed by seed. Analysis only; Python, read-only on flight data.
Decomposition (a) timing, (b) predictor, (c) history (planner_rev21_position.flight_position), arrival, dV, guard decisions (decisions.csv: burns, reasons, first differing pass), arm configuration (summary columns), Fisher exact test on
'within 1 deg with L1' and the systematic shift of the history term (EI 150 minus EI 2000, paired by seed, bootstrap CI of the median, 10,000 resamples, seed 20261030). Usage: python3 -I planner_rev30_g0528.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev21_position as PS
import planner_rev28_t as T28
import planner_rev20_score as SC


def fisher(a, b, c, d):
    """Two-sided Fisher exact p for [[a, b], [c, d]]."""
    n = a + b + c + d; r1, c1 = a + b, a + c
    pr = lambda x: math.comb(r1, x) * math.comb(n - r1, c1 - x) / math.comb(n, c1)
    p0 = pr(a); return float(sum(pr(x) for x in range(max(0, c1 - (n - r1)), min(r1, c1) + 1) if pr(x) <= p0 + 1e-12))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    F1 = V.load("rev23geo", None); F2 = V.load("rev26val2", None); S1 = pd.read_csv(next(V.ROOT.glob("summary_rev23geo_*.csv"))); S2 = pd.read_csv(next(V.ROOT.glob("summary_rev26val2_*.csv")))
    r1 = S1[S1.geometry == "G0528"].iloc[0]; r2 = S2[S2.geometry == "G0528"].iloc[0]; skip = {"run_id", "tag", "wall_s", "run_wall_s", "passes", "t_end_days", "total_dv_mps", "burns", "nu_cum", "n_decisions", "n_replan_rows", "planner_dp_total_s", "planner_dp_fraction", "host", "gram_core_epoch", "safety_alt0_m", "safety_h_local_m", "collapsed_decisions", "error", "controller_failures"}
    cfg_diff = {c: (str(r1[c]), str(r2[c])) for c in S1.columns if c in S2.columns and c not in skip and str(r1[c]) != str(r2[c])}; only2 = [c for c in S2.columns if c not in S1.columns]
    rows = []
    for s in (13, 14, 15, 16):
        f1 = F1[("OLRT3W3", "G0528", "L0", s)]; f2 = F2[("OLRT3W3L", "G0528", "L0", s)]; row = dict(seed=s)
        for tag, f in (("a", f1), ("b", f2)):
            q = PS.flight_position(f["dir"], "G0528", f["tstar"]); D = pd.read_csv(f["dir"] / "decisions.csv"); b = D[D.delta_v_mps.abs() > 0]
            row[tag] = dict(raw=f["pos"]["raw"], l1=f["pos"]["l1"], tau=f["pos"]["tau"], a=q["a_timing_deg"], b=q["b_pred_deg"], c=q["c_hist_deg"], t_arr=f["t_arr_d"], dt=q["dt_d"], dv=f["dv_sum"], burns=len(b), passes=f["n"], err=f["err"], reasons=b.reason.value_counts().to_dict(), rank_ok=bool(f["ok"]), win=bool(f["window"]),
                            dv_vec=None)
        # first decision pass that differs (target_pass, executed dv)
        D1 = pd.read_csv(F1[("OLRT3W3", "G0528", "L0", s)]["dir"] / "decisions.csv"); D2 = pd.read_csv(F2[("OLRT3W3L", "G0528", "L0", s)]["dir"] / "decisions.csv"); m = min(len(D1), len(D2)); dd = np.where(~np.isclose(D1.delta_v_mps.values[:m], D2.delta_v_mps.values[:m], atol=1e-12))[0]
        row["first_dv_diff_pass"] = int(D1.target_pass.values[dd[0]]) if len(dd) else None; rows.append(row)
    ca = np.array([r["a"]["c"] for r in rows]); cb = np.array([r["b"]["c"] for r in rows]); dc = cb - ca
    sh = T28.paired(dc); da = np.array([r["b"]["a"] - r["a"]["a"] for r in rows]); dbb = np.array([r["b"]["b"] - r["a"]["b"] for r in rows]); dl1 = np.array([r["b"]["l1"] - r["a"]["l1"] for r in rows])
    n1 = sum(r["a"]["l1"] < 1 for r in rows); n2 = sum(r["b"]["l1"] < 1 for r in rows); r1_ = sum(r["a"]["raw"] < 1 for r in rows); r2_ = sum(r["b"]["raw"] < 1 for r in rows)
    out_ = dict(rows=rows, cfg_diff=cfg_diff, only_in_rev26=only2, history_shift=sh, fisher_l1=fisher(n1, 4 - n1, n2, 4 - n2), fisher_raw=fisher(r1_, 4 - r1_, r2_, 4 - r2_), counts=dict(l1=(n1, n2), raw=(r1_, r2_)),
                med_shift=dict(c=float(np.median(dc)), a=float(np.median(da)), b=float(np.median(dbb)), l1=float(np.median(dl1)), c_mean=float(np.mean(dc)), c_sd=float(np.std(dc, ddof=1))),
                pooled_8=dict(l1=n1 + n2, raw=r1_ + r2_), drag_expect_shift_deg=float(0.001 * 63.6 * 0.6))
    (out / "rev30_g0528.json").write_text(json.dumps(SC.nan_none(out_), indent=1, default=float))
    print(cfg_diff, only2); print(sh); print("fisher l1", out_["fisher_l1"], "raw", out_["fisher_raw"], out_["counts"], out_["med_shift"])
    for r in rows: print(r["seed"], {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r["a"].items()}, "\n   ", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r["b"].items()}, r["first_dv_diff_pass"])


if __name__ == "__main__":
    main()
