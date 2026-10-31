"""Revision 31.3 (rev31cls): OLRT3-W4L5 (v4 + L1, 5 km aim, EI 150) at four unflown cells (20 and 80 deg classes, Odyssey radius and 2,000 km), seeds 13-16, plus the regression flight. Analysis only; Python, read-only on flight data.
Tables planner_ground_tables_rev31.json and IC plans headline_rev31.json are found by planner_rev21_position (table_entry / headline_ic globs). Cross-check against the logged pos_err_deg / l1_err_deg.
Usage: python3 -I planner_rev31_cls.py --out-dir DIR"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev31_pos as P31
import planner_rev20_score as SC


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    F = V.load("rev31cls", None); S = pd.read_csv(next(V.ROOT.glob("summary_rev31cls_*.csv"))).set_index("run_id"); R = {}; dmax = [0.0, 0.0]
    for geom in ("G0140_AODY", "G0042_A2000", "G3257_AODY", "G3161_A2000"):
        fs = [f for k, f in sorted(F.items(), key=lambda t: t[0][3]) if k[1] == geom and k[0] == "OLRT3W4L5"]; s = P31.stats(fs)
        lg = [(float(S.loc[f["dir"].name].pos_err_deg), float(S.loc[f["dir"].name].l1_err_deg)) for f in fs]
        for f, (lr, ll) in zip(fs, lg): dmax[0] = max(dmax[0], abs(f["pos"]["raw"] - lr)); dmax[1] = max(dmax[1], abs(f["pos"]["l1"] - ll))
        s.update(logged_l1=[x[1] for x in lg], logged_raw=[x[0] for x in lg], tstar=float(fs[0]["tstar"]), dt=[float(f["dt_nom_w"]) for f in fs], err=[float(f["err"]) for f in fs], q=[float(f["q_last"]) for f in fs], passes=[int(f["n"]) for f in fs]); R[geom] = s
    c = SC.compare_passes(V.ROOT / "rev31cls" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"], ei=float(S[S.arm == "OLRT3"].iloc[0].ei_km), stop=str(S[S.arm == "OLRT3"].iloc[0].stop_reason))
    (out / "rev31cls.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg, crosscheck_max_dev=dmax)), indent=1, default=float))
    print("cross-check max |raw|, |L1| vs logged:", dmax)
    for g, s in R.items(): print(g, f"ok {s['ok']} hard {s['hard']} win {s['win']} corr {s['corr']} raw {s['raw_med']:.2f} l1 {s['l1_med']:.2f} raw1 {s['raw1']} l11 {s['l11']} ALL {s['all']} cross {s['cross_med']:.2f} alongafter {s['along_after_med']:.2f} share {s['cross_share']:.2f} dV {s['dv_med']:.2f}", np.round(s['l1'], 2).tolist(), np.round(s['dt'], 2).tolist(), np.round(s['err'], 1).tolist(), np.round(s['q'], 3).tolist())
    print(reg)


if __name__ == "__main__":
    main()
