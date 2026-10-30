"""Rebuild summary rows for flights that completed (passes.csv present) but have no row in any summary_*.csv, e.g. when a campaign was stopped at a flight boundary (summaries are written only when a campaign
finishes). NEW file, not a frozen file. A rebuilt row has the columns planner_analysis.collect reads (run_id, arm, geometry, status, total_dv_mps, burns, t_end_days, stop_reason, run_wall_s) and is MARKED
(column rebuilt = 1): dV = sum of |decided burns| (equals the flown dV to 1e-11 m/s in the median of 80 real rows; the last decision may not have been executed when the run ended at the next apoapsis, so column dv_uncertainty_mps = |last decision| bounds the error; maximum observed 0.34 m/s), burns = nonzero decisions, t_end_days = last periapsis time + half the osculating period of the last pass (the
engine's end time is the closing apoapsis; error well below 0.1 d, checked against real rows by --check), stop_reason = end_radius when the last logged apoapsis is at or below the target 7,408 km, else 'interrupted'.
run_wall_s is NaN. Partial flights (no passes.csv) are not touched. Usage:
  python3 -I planner_rebuild_summary.py --planner-root R/planner --phase held --out summary_held_mars_ws_h3_rebuilt.csv [--check]     (writes only rows for run folders absent from all summaries)
--check compares the reconstruction with the existing real rows and prints the differences instead of writing."""
import argparse, glob, math, pathlib
import numpy as np, pandas as pd

MU = 42828.375e9; ISP = 220.0; G0 = 9.80665; RA_T = 7.408e6; M0 = 461.0


def rebuild(d):
    P = pd.read_csv(d / "passes.csv"); Dc = d / "decisions.csv"
    D = pd.read_csv(Dc) if Dc.exists() and Dc.stat().st_size > 2 else pd.DataFrame(dict(delta_v_mps=[]))
    last = P.iloc[-1]; a = float(last.a_m); T = 2 * math.pi * math.sqrt(a ** 3 / MU)
    ra = P.ra_m.dropna().values; ended = len(ra) > 0 and ra[-1] <= RA_T
    rid = d.name; arm = rid.split("_PL_")[1].split("_pert_")[0]; gid = rid.split("_PL_")[0]
    return dict(run_id=rid, arm=arm, geometry=gid, status="ok", passes=len(P), t_end_days=(float(last.time_s) + 0.5 * T) / 86400.0, total_dv_mps=float(D.delta_v_mps.abs().sum()), dv_uncertainty_mps=float(abs(D.delta_v_mps.iloc[-1])) if len(D) else 0.0,
                burns=int((D.delta_v_mps.abs() >= 1e-3).sum()), stop_reason="end_radius" if ended else "interrupted", run_wall_s=float("nan"), rebuilt=1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--phase", default="held"); ap.add_argument("--out"); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    root = pathlib.Path(a.planner_root); S = [pd.read_csv(f) for f in sorted(glob.glob(str(root / f"summary_{a.phase}_*.csv")))]
    have = pd.concat(S) if S else pd.DataFrame(columns=["run_id"]); have = have[have.status == "ok"] if "status" in have else have
    runs = sorted(glob.glob(str(root / "runs" / "*" / "passes.csv")))
    if a.check:
        worst = dict(t=0.0, dv=0.0, burns=0, stop=0); n = 0
        for r in have.itertuples():
            d = root / "runs" / r.run_id
            if not (d / "passes.csv").exists():
                continue
            b = rebuild(d); n += 1
            worst["t"] = max(worst["t"], abs(b["t_end_days"] - r.t_end_days)); worst["dv"] = max(worst["dv"], abs(b["total_dv_mps"] - r.total_dv_mps)); worst["burns"] = max(worst["burns"], abs(b["burns"] - r.burns)); worst["stop"] += b["stop_reason"] != r.stop_reason
        print(f"checked {n} real rows: max |t_end error| {worst['t']:.4f} d, max |dV error| {worst['dv']:.5f} m/s, max |burn count error| {worst['burns']}, stop_reason mismatches {worst['stop']}"); return
    known = set(have.run_id); rows = [rebuild(pathlib.Path(p).parent) for p in runs if pathlib.Path(p).parent.name not in known]
    print(f"{len(runs)} completed run folders, {len(known)} in summaries, {len(rows)} rebuilt"); 
    if rows and a.out:
        pd.DataFrame(rows).to_csv(root / a.out, index=False); print("wrote", root / a.out)


if __name__ == "__main__":
    main()
