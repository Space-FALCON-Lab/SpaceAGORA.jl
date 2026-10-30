"""POST HOC diagnosis of the dev result that L3 at Odyssey costs more dV than L1 (Revision 12 request; not a pre-registered analysis). Reads dev replan.csv / passes.csv and the frozen dev tube
(q_cap 0.234, W 0.7186 from the pre-calibration ground table, the config the dev flights flew): which L3 burns are 'forced' (the predicted next-pass ln q z0 lies outside the full tube, so L1 would burn as well)
and which are 'anticipatory' (inside the tube: caused by the horizon, the drift feed-forward and the early-timing weights), burn-size and sign statistics, and the mean flown heat level per arm.
Usage: python3 -I planner_dev_diagnosis.py --results-root R [--out F]"""
import argparse, glob, json, math, pathlib
import numpy as np, pandas as pd

Q_CAP_DEV = 0.234; W_DEV = 0.7185552648324176      # dev config (planner_config.toml before the calibration freeze; ODY ground table W at d_ref 64.2)


def one(d):
    R = pd.read_csv(d / "replan.csv"); P = pd.read_csv(d / "passes.csv")
    hi = math.log(Q_CAP_DEV); lo = hi - W_DEV
    ex = R.executed_dv_mps.values; z0 = R.z0_pred_lnq.values; burn = np.abs(ex) >= 1e-3
    forced = burn & ((z0 > hi) | (z0 < lo)); antic = burn & ~forced
    sgn = np.sign(ex[burn]); rev = int((np.diff(sgn) != 0).sum()) if burn.sum() > 1 else 0
    return dict(orbits=len(R), burns=int(burn.sum()), dv=float(np.abs(ex).sum()), forced_burns=int(forced.sum()), forced_dv=float(np.abs(ex[forced]).sum()), antic_burns=int(antic.sum()), antic_dv=float(np.abs(ex[antic]).sum()),
                antic_down_dv=float(np.abs(ex[antic & (ex < 0)]).sum()), up_dv=float(ex[ex > 0].sum()), down_dv=float(-ex[ex < 0].sum()), sign_reversals=rev, median_burn=float(np.median(np.abs(ex[burn]))) if burn.any() else 0.0,
                mean_lnq=float(np.log(P.q_w_cm2).mean()), q_mean=float(P.q_w_cm2.mean()), frac_above_hi=float((np.log(P.q_w_cm2) > hi).mean()), z0_mean_offset_from_mid=float(np.nanmean(z0) - 0.5 * (lo + hi)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--out"); a = ap.parse_args()
    root = pathlib.Path(a.results_root) / "planner" / "dev"; out = {}
    for arm in ("L1", "L3"):
        for d in sorted(glob.glob(str(root / f"ODY_PL_{arm}_pert_ODY_s*"))):
            out.setdefault(arm, {})[int(d.split("_s")[-1])] = one(pathlib.Path(d))
    summ = {}
    for arm, v in out.items():
        df = pd.DataFrame(v).T; summ[arm] = df.median().to_dict(); summ[arm + "_sum"] = df.sum().to_dict()
    res = dict(per_flight=out, medians={k: v for k, v in summ.items() if not k.endswith("_sum")}, totals={k: v for k, v in summ.items() if k.endswith("_sum")})
    txt = json.dumps(res, indent=1, default=float); print(txt)
    if a.out: pathlib.Path(a.out).write_text(txt)


if __name__ == "__main__":
    main()
