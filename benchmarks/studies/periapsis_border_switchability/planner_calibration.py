"""Calibration from the dev flights (plan U13, U14, Revision 12): read-only. It does NOT edit any config.
  sigma_p (U14): sd of the one-pass prediction error e_p of ln q_mean, e_p = ln q_meas(pass j) - (z0_pred(j) - dv_j / w0_j) from replan.csv and passes.csv of the dev flights (the planner's own prediction of the
  next pass after its burn, against the measured mean-field heat rate); compared with the ASSUMED 0.10 behind q_cap_mean 0.234.
  D_ref (U13): per dev seed min(days B-ODY, days F-O) at the Odyssey geometry (B-ODY from this dev phase, F-O from the Stage O dev1 flights, seeds 1 to 8); compared with the ASSUMED 64.2 d (DERIVED, R1.2).
Usage: python3 -I planner_calibration.py --results-root R [--stageo-dev DIR] [--out F]"""
import argparse, glob, json, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB


def ep_of(run_dir):
    R = pd.read_csv(run_dir / "replan.csv"); P = pd.read_csv(run_dir / "passes.csv").set_index("pass")
    out = []
    for _, r in R.iterrows():
        j = int(r["pass"])
        if j in P.index and np.isfinite(r["z0_pred_lnq"]):
            out.append(float(np.log(P.loc[j, "q_w_cm2"]) - (r["z0_pred_lnq"] - r["executed_dv_mps"] / r["w0"])))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--stageo-dev"); ap.add_argument("--out"); a = ap.parse_args()
    cfg = PLB.load_cfg(); root = pathlib.Path(a.results_root) / "planner"; res = dict(assumed=dict(sigma_p=cfg["planner"]["sigma_p"], q_cap_mean=cfg["planner"]["q_cap_mean_w_cm2"], d_ref_days=cfg["planner"]["d_ref_days"]))
    S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(root / "summary_dev_mars_*.csv")))]); S = S[S.status == "ok"]
    eps = {}; days = {}
    for r in S.itertuples():
        d = root / "dev" / r.run_id; key = f"{r.arm}@{'ODY' if r.geometry == 'ODY' else 'G1'}"
        if (d / "replan.csv").exists():
            e = ep_of(d); eps.setdefault(key, []).append(e)
        if r.arm == "BODY":
            days[int(str(r.run_id).split("_pert_")[1].split("_")[1].lstrip("s"))] = float(r.t_end_days)
    res["e_p"] = {k: dict(n=int(sum(len(x) for x in v)), mean=float(np.mean(np.concatenate(v))), sd=float(np.std(np.concatenate(v), ddof=1)), n_flights=len(v)) for k, v in eps.items()}
    allv = np.concatenate([np.concatenate(v) for v in eps.values()]) if eps else np.array([]); res["sigma_p_measured_pooled"] = float(np.std(allv, ddof=1)) if len(allv) > 2 else None
    sd = pathlib.Path(a.stageo_dev) if a.stageo_dev else PLB.task_dir() / "stageo_data" / "dev"
    fo = {}
    for f in sorted(glob.glob(str(sd / "summary_dev1_*.csv"))):
        T = pd.read_csv(f)
        for r in T[(T.status == "ok") & (T.arm.astype(str).str.replace("-", "") == "FO")].itertuples():
            fo[int(str(r.run_id).split("_pert_")[1].split("_")[1].lstrip("s"))] = float(r.t_end_days)
    seeds = sorted(set(days) & set(fo)); dref = {k: min(days[k], fo[k]) for k in seeds}
    res["d_ref_per_seed"] = dref; res["body_days"] = {k: days[k] for k in sorted(days)}; res["fo_days_dev"] = {k: fo[k] for k in seeds}
    res["d_ref_median"] = float(np.median(list(dref.values()))) if dref else None; res["body_days_median"] = float(np.median(list(days.values()))) if days else None
    res["verdict_vs_assumed"] = dict(sigma_p_ratio=(res["sigma_p_measured_pooled"] / cfg["planner"]["sigma_p"]) if res["sigma_p_measured_pooled"] else None,
                                     d_ref_ratio=(res["d_ref_median"] / cfg["planner"]["d_ref_days"]) if res["d_ref_median"] else None)
    txt = json.dumps(res, indent=1, default=float); print(txt)
    if a.out:
        pathlib.Path(a.out).write_text(txt)


if __name__ == "__main__":
    main()
