"""Revision 26 validation (rev26val2) and the final acceptance table inputs (analysis only; Python, read-only on flight data).
Per flight of any block: success ([-25, 0] km of the geometry's target radius), hard, window (first passage, 3 final orbits, against the logged nominal t*), corridor (AC5), aerobraking dV to the target pass (ve ln m0/m), the apoapsis position error raw and with L1
(Revision 26: the rendezvous time t*' is chosen in [t* - W, t* + W] and within W of the arrival, error |E - rate tau s_hat|), mirroring planner_l1 / planner_pos_final with planner_rev21_position (J2Track nominal, target direction at t*, flight direction at arrival).
Cross-check: rev26val2 flights log pos_err_deg, l1_err_deg, pos_t_arr_d, pos_tstar_d in the summary. Mission success = success and window and L1 position < 1.0 deg (the Revision 26 user requirement).
Usage: python3 -I planner_rev26_val2.py --out-dir DIR"""
import argparse, json, math, os, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev20_score as SC
import planner_rev21_position as PS
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
OFF = {"ODY": 0.15}; MU = PS.MU; TOLP = 1.0      # aim offset ASSUMED (days); position tolerance 1.0 deg (user requirement, Revision 26)


def per(a, rp): return 2 * math.pi * math.sqrt(((a + rp) / 2) ** 3 / MU)


def position(d, tok, tstar_d, rt, ic_x0=None, nom=None):
    """Raw and L1 position error of a flight directory (mirror of planner_pos_final). Returns dict or None when the target radius is not reached."""
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= rt))[0]
    if not len(hit): return None
    n = int(hit[0]); row = P.iloc[n]; a, e, inc = float(row.a_m), float(row.e), float(row.i_rad); dom, dOm = PS.j2_orbit(a, e, inc)
    u_fl = PS.apo_hat(float(row.omega_rad) + 0.5 * dom, float(row.raan_rad) + 0.5 * dOm, inc); nom = nom or PS.nominal(tok, ic_x0); nt = nom["track"]; ts = tstar_d * 86400.0
    u_tgt = nt.apo(ts); E = PS.logmap(u_tgt, u_fl); rv = PS.tproj(u_tgt, nt.rate_vec(ts)); s = rv / np.linalg.norm(rv); rate_dd = math.degrees(float(np.linalg.norm(nt.rate_vec(ts)))) * 86400.0
    t_arr = float(row.time_s) + 0.5 * 2 * math.pi * math.sqrt(a ** 3 / MU); dt = t_arr / 86400.0 - tstar_d
    ra_f = float(row.ra_m) if np.isfinite(row.ra_m) else a * (1 + e); rp_f = float(row.rp_m) if np.isfinite(row.rp_m) else a * (1 - e); W = 3.0 * per(ra_f, rp_f) / 86400.0
    al = math.degrees(float(np.dot(E, s))); cr = math.degrees(float(np.linalg.norm(E - np.dot(E, s) * s))); lo = max(-W, dt - W); hi = min(W, dt + W)
    if lo > hi: l1, tau = math.degrees(float(np.linalg.norm(E))), 0.0
    else:
        tau = float(np.clip(al / rate_dd, lo, hi)) if rate_dd > 0 else 0.0; l1 = math.hypot(al - rate_dd * tau, cr)
    return dict(raw=math.degrees(PS.ang(u_fl, u_tgt)), l1=l1, tau=tau, t_arr_d=t_arr / 86400.0, along=al, cross=cr)


def load(block, glob_pat, arms=None):
    S = SC.load_summary(ROOT, f"summary_{block}_*.csv"); out = {}
    for k, d in SC.scan(ROOT, block).items():
        if not SC.finished(d): continue
        tok = PS.GEOM_ID.get(SC.GID.get(k[1], k[1]), k[1]); rt = PS.target_radius(tok); sr = S.get(d.name) or {}
        f = SC.flight(d, sr or None, dv_pre_fallback=0.0, target=rt); f["dir"] = d; f["tok"] = tok; f["rt"] = rt; f["block"] = block; f["summ"] = sr
        aim = f["tt_aim_d"]; icx = float(sr.get("ic_x0")) if sr.get("ic_x0") is not None and np.isfinite(float(sr.get("ic_x0"))) else None
        if np.isfinite(aim): nomt = float(sr.get("pos_tstar_d")) if sr.get("pos_tstar_d") is not None and np.isfinite(float(sr.get("pos_tstar_d"))) else aim - OFF.get(k[1], 0.0)
        else:
            try: nomt = PS.nominal(tok)["days"]
            except Exception: nomt = float("nan")
        f["tstar"] = nomt
        if f["reached"] and np.isfinite(nomt):
            f["dt_nom_w"] = (f["t_arr_d"] - nomt) / f["wmax_d"]; f["window"] = bool(abs(f["dt_nom_w"]) <= 1.0)
            try: f["pos"] = position(d, tok, nomt, rt, icx if str(sr.get("ic_scenario")) == "S-free" else None)
            except Exception as ex: f["pos"] = None; f["pos_err"] = str(ex)
        else: f["dt_nom_w"] = float("nan"); f["window"] = False; f["pos"] = None
        f["corr"] = bool(f["reached"] and f["corridor"] and f["qs_last_ok"]); f["l1ok"] = bool(f["pos"] and f["pos"]["l1"] < TOLP); f["rawok"] = bool(f["pos"] and f["pos"]["raw"] < TOLP)
        f["all"] = bool(f["ok"] and not f["hard"] and f["window"] and f["corr"] and f["l1ok"]); f["all_raw"] = bool(f["ok"] and not f["hard"] and f["window"] and f["corr"] and f["rawok"]); out[k] = f
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    blocks = {b: load(b, None) for b in ("rev26val2", "rev27ei", "rev21v", "rev20val", "rev23geo", "rev27low", "rev21c", "rev17")}
    rows = []
    for b, F in blocks.items():
        for k, f in F.items():
            sm = f["summ"]
            rows.append(dict(block=b, arm=k[0], geom=k[1], level=k[2], seed=k[3], reached=f["reached"], ok=f["ok"], hard=f["hard"], window=f["window"], corr=f["corr"], dv=f["dv_sum"], dv_total=f["dv_sum"] + (f["dv_pre"] if np.isfinite(f["dv_pre"]) else 0.0), burns=f["burns"], err=f["err"], dt_nom_w=f["dt_nom_w"], t_arr=f["t_arr_d"], tstar=f["tstar"],
                             pos_raw=None if not f["pos"] else f["pos"]["raw"], pos_l1=None if not f["pos"] else f["pos"]["l1"], tau=None if not f["pos"] else f["pos"]["tau"], all=f["all"], all_raw=f["all_raw"], log_pos=float(sm.get("pos_err_deg", float("nan"))) if sm else float("nan"), log_l1=float(sm.get("l1_err_deg", float("nan"))) if sm else float("nan"),
                             log_t_arr=float(sm.get("pos_t_arr_d", float("nan"))) if sm else float("nan"), log_tstar=float(sm.get("pos_tstar_d", float("nan"))) if sm else float("nan"), stop=f["stop"], ei=float(sm.get("ei_km", 2000.0)) if sm else 2000.0))
    (out / "val2_rows.json").write_text(json.dumps(SC.nan_none(rows), indent=1)); df = pd.DataFrame(rows)
    v = df[df.block == "rev26val2"]; v = v[v.log_pos.notna()]
    print("cross-check rev26val2 (mine vs logged): max |raw|", float((v.pos_raw - v.log_pos).abs().max()), "max |L1|", float((v.pos_l1 - v.log_l1).abs().max()), "max |t_arr|", float((v.t_arr - v.log_t_arr).abs().max()), "max |tstar|", float((v.tstar - v.log_tstar).abs().max()))
    print(v[["arm", "geom", "level", "seed", "pos_raw", "log_pos", "pos_l1", "log_l1", "t_arr", "log_t_arr"]].round(3).to_string())


if __name__ == "__main__":
    main()
