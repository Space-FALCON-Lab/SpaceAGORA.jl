"""Revision 31 extension (rev31ext): seeds 17-20 at G0486 (W4L5, W3E5) and G1038 (W4L5), pooled with rev31pos s13-16 (n = 8); G1038 v3 comparator rev29aim s13-20. Analysis only; Python, read-only on flight data.
Per cell at n = 8: within 1 deg with L1 (Wilson), ALL (success + hard 0 + window + L1 < 1 deg; corridor reported), paired dV v4 minus v3. Usage: python3 -I planner_rev31_ext.py --out-dir DIR"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev31_pos as P31
import planner_rev20_score as SC


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev31ext", "rev31pos", "rev29aim")}
    def g(arm, geom, blks):
        d = {}
        for b in blks: d.update({k[3]: f for k, f in B[b].items() if k[0] == arm and k[1] == geom and k[2] == "L0"})
        return [d[s] for s in sorted(d)]
    cells = {"G0486": dict(v4=g("OLRT3W4L5", "G0486", ("rev31pos", "rev31ext")), v3=g("OLRT3W3E5", "G0486", ("rev31pos", "rev31ext"))), "G1038": dict(v4=g("OLRT3W4L5", "G1038", ("rev31pos", "rev31ext")), v3=g("OLRT3W3E5", "G1038", ("rev29aim",)))}
    R = {}
    for geom, c in cells.items():
        s4 = P31.stats(c["v4"]); s3 = P31.stats(c["v3"]); n = min(len(c["v4"]), len(c["v3"])); d = [x["dv_sum"] - y["dv_sum"] for x, y in zip(c["v4"][:n], c["v3"][:n])]
        R[geom] = dict(v4=s4, v3=s3, dv=T28.paired(d), fisher_l1=None)
    cc = SC.compare_passes(V.ROOT / "rev31ext" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=cc["equal"], ncols=cc["ncols"], rows=(cc["rows_a"], cc["rows_b"]), first=cc["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev31ext_*.csv"))); rr = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(rr.ei_km); reg["stop"] = str(rr.stop_reason)
    (out / "rev31ext.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1, default=float))
    for g_, r in R.items():
        for t in ("v4", "v3"):
            s = r[t]; print(g_, t, f"n {s['n']} ok {s['ok']} hard {s['hard']} win {s['win']} corr {s['corr']} l1med {s['l1_med']:.2f} l11 {s['l11']} W=({s['wilson'][0]:.2f},{s['wilson'][1]:.2f}) ALL {s['all']} W=({s['wilson_all'][0]:.2f},{s['wilson_all'][1]:.2f}) dV {s['dv_med']:.2f} cross {s['cross_med']:.2f} share {s['cross_share']:.2f}")
        p = r["dv"]; print("  dV", np.round(p["d"], 2).tolist(), round(p["med"], 2), np.round(p["ci"], 2).tolist(), p["verdict"], p["p"])
    print(reg)


if __name__ == "__main__":
    main()
