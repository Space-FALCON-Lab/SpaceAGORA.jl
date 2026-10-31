"""Revision 31b (rev31ext2): OLRT3-W4L5 at G3257_AODY (80 deg) and G0140_AODY (20 deg), seeds 17-20, pooled with rev31cls s13-16 (n = 8). Analysis only; Python, read-only on flight data.
Usage: python3 -I planner_rev31_ext2.py --out-dir DIR"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev31_pos as P31
import planner_rev20_score as SC


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev31cls", "rev31ext2")}; R = {}
    for geom in ("G3257_AODY", "G0140_AODY"):
        d = {}
        for b in ("rev31cls", "rev31ext2"): d.update({k[3]: f for k, f in B[b].items() if k[0] == "OLRT3W4L5" and k[1] == geom and k[2] == "L0"})
        fs = [d[s] for s in sorted(d)]; s = P31.stats(fs); s.update(seeds=sorted(d), err=[float(f["err"]) for f in fs], dt=[float(f["dt_nom_w"]) for f in fs]); R[geom] = s
    cc = SC.compare_passes(V.ROOT / "rev31ext2" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=cc["equal"], ncols=cc["ncols"], rows=(cc["rows_a"], cc["rows_b"]), first=cc["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev31ext2_*.csv"))); rr = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(rr.ei_km); reg["stop"] = str(rr.stop_reason)
    (out / "rev31ext2.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1, default=float))
    for g, s in R.items(): print(g, f"n {s['n']} ok {s['ok']} hard {s['hard']} win {s['win']} corr {s['corr']} raw1 {s['raw1']} l11 {s['l11']} W=({s['wilson'][0]:.2f},{s['wilson'][1]:.2f}) ALL {s['all']} W=({s['wilson_all'][0]:.2f},{s['wilson_all'][1]:.2f}) dV {s['dv_med']:.2f} l1 {np.round(s['l1'],2).tolist()} win {s['win']}")
    print(reg)


if __name__ == "__main__":
    main()
