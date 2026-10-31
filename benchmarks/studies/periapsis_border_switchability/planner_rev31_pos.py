"""Revision 31.1 (rev31pos): OLRT3-W4L5 (v4 + L1) and OLRT3-W3E5 (v3) at G1, G0975M, G0063M, G0890, G0486 (and W4L5 at G1038), L0, seeds 13-16, 5 km aim, EI 150. Analysis only; Python, read-only on flight data.
Per class: success, hard, window, corridor (reported), position raw / L1, within 1 deg (Wilson), cross-track versus along-track share of the remaining error (squared, after L1: cross^2 / L1^2), paired dV v4 minus v3, ALL (success + hard 0 + window + L1 < 1 deg; corridor reported).
Usage: python3 -I planner_rev31_pos.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev30_seed as S30
import planner_rev20_score as SC


def stats(fs):
    n = len(fs); l1 = [f["pos"]["l1"] for f in fs]; raw = [f["pos"]["raw"] for f in fs]; cr = np.array([f["pos"]["cross"] for f in fs]); al = np.array([f["pos"]["along"] for f in fs]); l1a = np.array(l1)
    al_after = np.sqrt(np.maximum(l1a ** 2 - cr ** 2, 0.0)); allc = sum(bool(f["ok"] and not f["hard"] and f["window"] and f["pos"]["l1"] < 1.0) for f in fs)
    return dict(n=n, ok=int(sum(f["ok"] for f in fs)), hard=int(sum(f["hard"] for f in fs)), win=int(sum(f["window"] for f in fs)), corr=int(sum(f["corr"] for f in fs)), raw_med=float(np.median(raw)), l1_med=float(np.median(l1)), raw1=int(sum(x < 1 for x in raw)), l11=int(sum(x < 1 for x in l1)),
                wilson=S30.wilson(int(sum(x < 1 for x in l1)), n), all=int(allc), wilson_all=S30.wilson(int(allc), n), cross_med=float(np.median(cr)), along_raw_med=float(np.median(np.abs(al))), along_after_med=float(np.median(al_after)),
                cross_share=float(np.sum(cr ** 2) / max(np.sum(l1a ** 2), 1e-12)), dv_med=float(np.median([f["dv_sum"] for f in fs])), dv=[float(f["dv_sum"]) for f in fs], l1=[float(x) for x in l1], raw=[float(x) for x in raw], cross=[float(x) for x in cr], along_after=[float(x) for x in al_after])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev31pos", "rev29aim")}
    sel = lambda blk, arm, geom, seeds=(13, 14, 15, 16): [f for k, f in sorted(B[blk].items(), key=lambda t: t[0][3]) if k[0] == arm and k[1] == geom and k[2] == "L0" and k[3] in seeds]
    R = {}
    for geom in ("G1", "G0975M", "G0063M", "G0890", "G0486", "G1038"):
        f4 = sel("rev31pos", "OLRT3W4L5", geom); f3 = sel("rev31pos", "OLRT3W3E5", geom) if geom != "G1038" else sel("rev29aim", "OLRT3W3E5", geom)
        d = [x["dv_sum"] - y["dv_sum"] for x, y in zip(f4, f3)]
        R[geom] = dict(v4=stats(f4), v3=stats(f3), dv_pair=T28.paired(d), v3_src="rev31pos" if geom != "G1038" else "rev29aim", seeds_23=bool(stats(f4)["l11"] in (2, 3)))
    c = SC.compare_passes(V.ROOT / "rev31pos" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev31pos_*.csv"))); rr = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(rr.ei_km); reg["stop"] = str(rr.stop_reason)
    (out / "rev31pos.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1, default=float))
    for g, r in R.items():
        for tag in ("v4", "v3"):
            s = r[tag]; print(g, tag, f"ok {s['ok']} hard {s['hard']} win {s['win']} corr {s['corr']} raw {s['raw_med']:.2f} l1 {s['l1_med']:.2f} raw1 {s['raw1']} l11 {s['l11']} ALL {s['all']} cross {s['cross_med']:.2f} alongafter {s['along_after_med']:.2f} share {s['cross_share']:.2f} dV {s['dv_med']:.2f}")
        p = r["dv_pair"]; print("   dV v4-v3", np.round(p["d"], 2).tolist(), round(p["med"], 2), np.round(p["ci"], 2).tolist(), p["verdict"])
    print(reg)


if __name__ == "__main__":
    main()
