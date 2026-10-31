"""Revision 27 B27-1 pairing (analysis only; Python, read-only on flight data): OLRT3-W3E at ODY L0 seeds 13-16 with ei_km = 150 (rev27ei) against OLRT3-W3 at the same seeds with ei_km = 2000 (rev21v).
Per pair: success, hard, window (against the logged nominal t* = aim minus the 0.15 d aim offset), t_arr, aerobraking dV, final apoapsis error, passes, first-pass q, burn comparison pass by pass on decisions.csv (executed dv and request time),
per-pass apoapsis drop; then the Revision 27 smallness rule. Usage: python3 -I planner_rev27_ei.py --out-dir DIR"""
import argparse, json, math, os, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev20_score as SC
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
OFF = 0.15      # ASSUMED aim offset (days) at Odyssey


def fl(block, arm, seed, pat):
    S = SC.load_summary(ROOT, pat); d = ROOT / block / f"ODY_PL_{arm}_pert_ODY_s{seed}_L0"; f = SC.flight(d, S.get(d.name), dv_pre_fallback=0.0); sr = S.get(d.name) or {}
    nom = f["tt_aim_d"] - OFF; f["nom"] = nom; f["dt_nom_w"] = (f["t_arr_d"] - nom) / f["wmax_d"]; f["window"] = bool(abs(f["dt_nom_w"]) <= 1.0); f["dir"] = d; f["wall_s"] = float(sr.get("run_wall_s", float("nan"))); f["ei_km"] = float(sr.get("ei_km", 2000.0) if sr.get("ei_km") == sr.get("ei_km") else 2000.0)
    return f


def burns(d):
    D = pd.read_csv(d / "decisions.csv"); b = D[D.delta_v_mps.abs() > 0]; return {int(r.target_pass): (float(r.delta_v_mps), float(r.request_time_s)) for r in b.itertuples()}, D


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True); rows = []
    for s in (13, 14, 15, 16):
        a_ = fl("rev21v", "OLRT3-W3", s, "summary_rev21v_*.csv"); b_ = fl("rev27ei", "OLRT3-W3E", s, "summary_rev27ei_*.csv")
        Pa = pd.read_csv(a_["dir"] / "passes.csv").sort_values("pass").reset_index(drop=True); Pb = pd.read_csv(b_["dir"] / "passes.csv").sort_values("pass").reset_index(drop=True); n = min(len(Pa), len(Pb))
        ra_a = np.concatenate(([28.533e6], Pa.ra_m.values[:-1])) if False else None
        da = -np.diff(np.concatenate(([np.nan], Pa.ra_m.values.astype(float)))); db = -np.diff(np.concatenate(([np.nan], Pb.ra_m.values.astype(float))))
        ba, Da = burns(a_["dir"]); bb, Db = burns(b_["dir"])
        # first pass where a burn differs (set, dv, or request time)
        allp = sorted(set(ba) | set(bb)); first = None
        for p in allp:
            x, y = ba.get(p), bb.get(p)
            if x is None or y is None or abs(x[0] - y[0]) > 1e-12 or abs(x[1] - y[1]) > 1e-6: first = (p, x, y); break
        # first pass where decisions differ in dv (including zeros) and the request-time offset between the runs on the shared prefix
        m = min(len(Da), len(Db)); dv_a = Da.delta_v_mps.values[:m]; dv_b = Db.delta_v_mps.values[:m]; idx = np.where(~np.isclose(dv_a, dv_b, atol=1e-12))[0]; fdv = int(Da.target_pass.values[idx[0]]) if len(idx) else None
        pref = fdv - 1 if fdv else m          # passes whose executed dv are identical (shared prefix)
        dt_req = (Da.request_time_s.values[:m] - Db.request_time_s.values[:m])
        ddrop = (db[1:n] - da[1:n]) / 1e3
        q1 = (float(Pa.q_w_cm2.iloc[0]), float(Pb.q_w_cm2.iloc[0]))
        rows.append(dict(seed=s, ok=(a_["ok"], b_["ok"]), hard=(a_["hard"], b_["hard"]), window=(a_["window"], b_["window"]), dt_nom_w=(a_["dt_nom_w"], b_["dt_nom_w"]), t_arr=(a_["t_arr_d"], b_["t_arr_d"]), d_tarr=b_["t_arr_d"] - a_["t_arr_d"], dv=(a_["dv_sum"], b_["dv_sum"]), d_dv=b_["dv_sum"] - a_["dv_sum"],
                         err=(a_["err"], b_["err"]), d_err=b_["err"] - a_["err"], passes=(a_["n"], b_["n"]), burns=(a_["burns"], b_["burns"]), q1=q1, q1_rel=q1[1] / q1[0] - 1, wall=(a_["wall_s"], b_["wall_s"]), stop=(a_["stop"], b_["stop"]),
                         first_burn_diff=None if first is None else dict(pass_=first[0], a=first[1], b=first[2]), first_dv_diff_pass=fdv, shared_prefix_passes=pref, dreq_prefix_s=dt_req[:max(pref, 1)].tolist(),
                         drop1=(da[1] / 1e3 if False else None), ddrop_first10_km=ddrop[:10].tolist(), drop_a_first10_km=(da[1:11] / 1e3).tolist(), drop_b_first10_km=(db[1:11] / 1e3).tolist(), ddrop_pref_med_abs_km=float(np.median(np.abs(ddrop[:max(pref - 1, 1)]))), ddrop_pref_max_abs_km=float(np.max(np.abs(ddrop[:max(pref - 1, 1)]))),
                         ra1=(float(Pa.ra_m.iloc[0] / 1e3), float(Pb.ra_m.iloc[0] / 1e3)), ddrop_all_med_km=float(np.median(ddrop)), n_burn_a=len(ba), n_burn_b=len(bb), rdv=sum(1 for p in ba if p in bb)))
    med = lambda k: float(np.median([r[k] for r in rows]))
    smallness = dict(outcome_changes=int(sum(1 for r in rows if r["ok"][0] != r["ok"][1] or r["hard"][0] != r["hard"][1] or r["window"][0] != r["window"][1])), med_abs_dt_d=float(np.median([abs(r["d_tarr"]) for r in rows])), med_abs_ddv=float(np.median([abs(r["d_dv"]) for r in rows])), med_abs_derr=float(np.median([abs(r["d_err"]) for r in rows])))
    smallness["small"] = bool(smallness["outcome_changes"] == 0 and smallness["med_abs_dt_d"] <= 0.1 and smallness["med_abs_ddv"] <= 0.3 and smallness["med_abs_derr"] <= 5.0)
    (out / "rev27ei_pairs.json").write_text(json.dumps(dict(rows=rows, smallness=smallness, med=dict(d_tarr=med("d_tarr"), d_dv=med("d_dv"), d_err=med("d_err"))), indent=1, default=float)); print(json.dumps(smallness)); print({k: med(k) for k in ("d_tarr", "d_dv", "d_err")})
    for r in rows: print(r["seed"], r["ok"], r["hard"], r["window"], "dt %.3f dV %.2f err %.1f->%.1f passes %s q1 %.5f first_burn_diff %s fdv %s pref %s ddrop first10 %s wall %s" % (r["d_tarr"], r["d_dv"], r["err"][0], r["err"][1], r["passes"], r["q1_rel"], r["first_burn_diff"], r["first_dv_diff_pass"], r["shared_prefix_passes"], np.round(r["ddrop_first10_km"], 1).tolist(), np.round(r["wall"], 0).tolist()))


if __name__ == "__main__":
    main()
