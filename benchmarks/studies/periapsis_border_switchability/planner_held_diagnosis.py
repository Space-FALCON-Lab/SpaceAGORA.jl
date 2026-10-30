"""POST HOC diagnosis of the held ladder at Odyssey (L0, L1, L2, L3): why L0 and L1 cost less dV than L2 and L3. Not pre-registered. From replan.csv/passes.csv of the held flights, with the CALIBRATED tube each arm flew
(q_cap 0.2503; width = ground-table W for L1 to L3, l0_width_ln for L0). 'forced' burn: the one-pass prediction z0 is outside the full tube; 'anticipatory': inside (horizon, reserve-shrunk tube, feed-forward, early weights).
Usage: python3 -I planner_held_diagnosis.py --planner-root R/planner [--out F]"""
import argparse, glob, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB


def one(d, lo, hi):
    R = pd.read_csv(d / "replan.csv"); P = pd.read_csv(d / "passes.csv"); ex = R.executed_dv_mps.values; z0 = R.z0_pred_lnq.values; burn = np.abs(ex) >= 1e-3
    forced = burn & ((z0 > hi) | (z0 < lo)); anti = burn & ~forced; sg = np.sign(ex[burn])
    return dict(burns=int(burn.sum()), dv=float(np.abs(ex).sum()), forced_burns=int(forced.sum()), forced_dv=float(np.abs(ex[forced]).sum()), antic_burns=int(anti.sum()), antic_dv=float(np.abs(ex[anti]).sum()),
                down_dv=float(-ex[ex < 0].sum()), up_dv=float(ex[ex > 0].sum()), reversals=int((np.diff(sg) != 0).sum()) if burn.sum() > 1 else 0, q_mean=float(P.q_w_cm2.mean()), frac_above_cap=float((np.log(P.q_w_cm2) > hi).mean()), passes=len(P))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--out"); a = ap.parse_args(); cfg = PLB.load_cfg(); root = pathlib.Path(a.planner_root) / "runs"
    td = PLB.task_dir(); W = json.loads((td / cfg["planner"]["ground_table"]["path"]).read_text())["geometries"]["ODY"]["W"]; hi = math.log(cfg["planner"]["q_cap_mean_w_cm2"]); out = {}
    for arm in ("L0", "L1", "L2", "L3"):
        lo = hi - (cfg["planner"]["l0_width_ln"] if arm == "L0" else W); v = {int(str(d).split("_s")[-1]): one(pathlib.Path(d), lo, hi) for d in sorted(glob.glob(str(root / f"ODY_PL_{arm}_pert_ODY_s*")))}
        df = pd.DataFrame(v).T; out[arm] = dict(n=len(df), sums=df.sum().to_dict(), medians=df.median().to_dict())
    txt = json.dumps(out, indent=1, default=float); print(txt)
    if a.out: pathlib.Path(a.out).write_text(txt)


if __name__ == "__main__":
    main()
