"""P25-5 measurement (analysis only; Python, read-only on flight data): flown apoapsis drop per unit q over the model's drop per unit q at the same apoapsis radius
(Revision 19 B state-dependent method of planner_rev18_calibrate: model drop = d0(ra_prev) * q / q_req, d0 interpolated on the table's ra; kappa = 1 is the model), on the rev25low flights with the mean-field rev25mf tables.
Per flight: kappa = sum(flown drops) / sum(model drops) over the selected passes; per pass: flown/model. Selection: ra before the pass < 7,408 km; and additionally e before the pass < 0.36 (e from the previous apoapsis and this pass's periapsis radius).
Usage: python3 -I planner_rev25_droplaw.py --out DIR"""
import argparse, json, math, os, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
TD = pathlib.Path(PLB.task_dir()); ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
VE = 220.0 * 9.80665; RA0 = 7.408e6


def flight(d, tab):
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); q = P.q_w_cm2.values.astype(float); rp = P.rp_m.values.astype(float)
    d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); ra_prev = np.concatenate(([float(tab["ra"][0])], ra[:-1]))      # start apoapsis of pass 1 = table row 0 (as planner_rev18_calibrate)
    model = np.interp(ra_prev, tab["ra"][::-1], d0[::-1]) * q / tab["q_req"]; flown = ra_prev - ra; e_pre = (ra_prev - rp) / (ra_prev + rp)
    ok = np.isfinite(flown) & (model > 0) & np.isfinite(e_pre)
    return dict(P=P, ra_prev=ra_prev, flown=flown, model=model, e_pre=e_pre, ok=ok)


def summ(x):
    x = np.asarray([v for v in x if np.isfinite(v)]); return dict(n=len(x), med=float(np.median(x)), q25=float(np.percentile(x, 25)), q75=float(np.percentile(x, 75)), mn=float(x.min()), mx=float(x.max())) if len(x) else dict(n=0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--block", default="rev25low"); ap.add_argument("--arm", default="OLRT3-W3"); ap.add_argument("--tag", default="rev25low"); a = ap.parse_args(); out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tabs = json.load(open(TD / "planner_ground_tables_rev25mf.json"))["geometries"]; S = pd.read_csv(next(ROOT.glob(f"summary_{a.block}_*.csv"))).set_index("run_id"); R = {}
    for tok in ("G0486", "G1230", "G1038"):
        g = tabs[tok]; tab = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in g.items()}; fl = []
        for dd in sorted((ROOT / a.block).glob(f"{tok}_PL_{a.arm}_pert_ODY_s*_L0")):
            F = flight(dd, tab); P = F["P"]; sr = S.loc[dd.name]; n = len(P); m = P.mass_kg.values.astype(float)
            r = dict(run=dd.name, passes=n, days=float(P.time_s.iloc[-1] / 86400), stop=str(sr.stop_reason), dv_mass=VE * math.log(m[0] / m[-1]), dv_summary=float(sr.total_dv_mps), ra_last_km=float(P.ra_m.dropna().iloc[-1] / 1e3), alt_last_km=float((P.ra_m.dropna().iloc[-1] - 3396190.0) / 1e3), n_nan_ra=int(P.ra_m.isna().sum()), e_last=float(F["e_pre"][np.isfinite(F["e_pre"])][-1]), e_first_sel=float(F["e_pre"][F["ok"] & (F["ra_prev"] < RA0)][0]))
            for tag, sel in (("ra<7408", F["ok"] & (F["ra_prev"] < RA0)), ("e<0.36", F["ok"] & (F["e_pre"] < 0.36) & (F["ra_prev"] < RA0)), ("e 0.30-0.36", F["ok"] & (F["e_pre"] < 0.36) & (F["e_pre"] >= 0.30)), ("e 0.25-0.30", F["ok"] & (F["e_pre"] < 0.30) & (F["e_pre"] >= 0.25)), ("e 0.20-0.25", F["ok"] & (F["e_pre"] < 0.25) & (F["e_pre"] >= 0.20)), ("e 0.17-0.20", F["ok"] & (F["e_pre"] < 0.20) & (F["e_pre"] >= 0.17)), ("e 0.10-0.17", F["ok"] & (F["e_pre"] < 0.17) & (F["e_pre"] >= 0.10)), ("e<0.10", F["ok"] & (F["e_pre"] < 0.10))):
                r[tag] = dict(n=int(sel.sum()), kappa=float(F["flown"][sel].sum() / F["model"][sel].sum()) if sel.any() else float("nan"), pass_ratio=(F["flown"][sel] / F["model"][sel]).tolist(), e_min=float(F["e_pre"][sel].min()) if sel.any() else float("nan"), e_max=float(F["e_pre"][sel].max()) if sel.any() else float("nan"))
            fl.append(r)
        cell = dict(tok=tok, flights=fl)
        for tag in ("ra<7408", "e<0.36", "e 0.30-0.36", "e 0.25-0.30", "e 0.20-0.25", "e 0.17-0.20", "e 0.10-0.17", "e<0.10"):
            allr = [x for r in fl for x in r[tag]["pass_ratio"]]
            cell[tag] = dict(flight_kappa=summ([r[tag]["kappa"] for r in fl if r[tag]["n"]]), pass_ratio=summ(allr), npass=len(allr), e_min=min([r[tag]["e_min"] for r in fl if r[tag]["n"]] or [float("nan")]), e_max=max(r[tag]["e_max"] for r in fl if r[tag]["n"]) if any(r[tag]["n"] for r in fl) else float("nan"))
        R[tok] = cell
    (out / f"droplaw_{a.tag}.json").write_text(json.dumps(R, indent=1)); 
    for t, c in R.items():
        for r in c["flights"]: print(t, r["run"][-6:], r["passes"], round(r["days"], 2), r["stop"], round(r["dv_mass"], 2), r["dv_summary"], round(r["alt_last_km"]), r["n_nan_ra"], round(r["e_last"], 3), r["ra<7408"]["n"], round(r["ra<7408"]["kappa"], 3), r["e<0.36"]["n"])


if __name__ == "__main__":
    main()
