"""Revision 30.2: position of the default-t* arm (rev28t OLRT3-W3T) under an ASSUMED convention: the target apoapsis direction is the planner's own predicted apoapsis direction at the chosen t*.
Construction (DERIVED, exact): the model-world open-loop plan of the flights (planner arm OLR, table of the geometry, flown kappa, pinned IC level x0 of the S-fixed plan, burn penalty pen) is re-solved with the Revision 18 time price: a linear cost -lam * w_k on the level
nu_k of pass k (w_k = apoapsis drop weight of pass k relative to the arrival pass on the lam = 0 path), lam (m/s per orbit of arrival time) bisected until the model arrival (planner_rev21_position.arrival_history, kappa) equals t* to within one final orbit.
The plan's (a, e) history drives the J2 track (planner_rev21_position.J2Track, i = i0) and the target direction at t*. Variants (to test construction dependence): pen 0.01 (flight) / 0; kappa flight / 1; the plan arriving just before / just after t*; the constant-shift rebuild of results_rev28t.md.
Usage: python3 -I planner_rev30_conv.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev21_position as PS
import planner_rev20_score as SC


def plan_for_tstar(tok, tstar_d, pen=None, kappa=None, side="nearest"):
    RP = PS._rp(); sw, P = RP.arm_def("OLR"); tab = PS.get_tab(tok); x0 = math.log(RP.QFIRST_CAL.get(tok, PS.O.Q_FIRST)) + RP.ic_shift(PS.headline_ic(tok)); rt = PS.target_radius(tok)
    pen = P["dp_burn_penalty_mps"] if pen is None else pen; kap = PS.kappa_for(tok) if kappa is None else kappa
    inp = PLB.ol_inputs(sw, tab, P, x0); lvl = np.array(inp["level"]); lo = np.array(inp["lo"]); hi = np.array(inp["hi"]); w = np.array(inp["w"]); tw = PLB.tie_weights(w, P["dp_tie_eps"]); A = lo - lvl; B = hi - lvl
    arr = lambda x: PS.arrival_history(tab, x, x0, kap, rt)
    nu0 = RP.dpf(A, B, tw, pen=pen, lin=None); x_ref = np.concatenate(([x0], lvl + nu0)); h0 = arr(x_ref)
    d0 = tab["d0"]; n = tab["n"]; dk = np.array([d0[min(k, n) - 1] * math.exp(x_ref[min(k, len(x_ref)) - 1]) / tab["q_req"] for k in range(1, len(x_ref) + 1)]); wT = dk[1:1 + len(lvl)] / dk[min(len(dk) - 1, h0["N"] - 1)]
    def at(lam):
        nu = RP.dpf(A, B, tw, pen=pen, lin=-lam * wT)
        if nu is None: return None
        x = np.concatenate(([x0], lvl + nu)); h = arr(x); return dict(lam=lam, x=x, h=h, days=h["days"])
    p0 = dict(lam=0.0, x=x_ref, h=h0, days=h0["days"]); need_early = p0["days"] > tstar_d; sg = 1.0 if need_early else -1.0
    lo_l, hi_l = 0.0, 1.0e-3; ph = at(sg * hi_l)
    while ph is not None and ((ph["days"] > tstar_d) if need_early else (ph["days"] < tstar_d)) and hi_l < 200.0:
        lo_l = hi_l; hi_l *= 2.0; ph = at(sg * hi_l)
    if ph is None: return None
    plans = {}                                    # bisect the discontinuous arrival(lam) to the two sides of t*
    pa, pb = (p0 if lo_l == 0.0 else at(sg * lo_l)), ph
    for _ in range(50):
        mid = 0.5 * (lo_l + hi_l); pm = at(sg * mid)
        if pm is None: break
        if ((pm["days"] > tstar_d) if need_early else (pm["days"] < tstar_d)): lo_l, pa = mid, pm
        else: hi_l, pb = mid, pm
    before = pa if (pa["days"] <= tstar_d) else pb; after = pb if (pb["days"] >= tstar_d) else pa
    cand = {"before": before, "after": after}; best = min(cand.values(), key=lambda p: abs(p["days"] - tstar_d)); sel = best if side == "nearest" else cand[side]
    i0, om0, Om0 = PS.geometry_elements(tok); h = sel["h"]; tr = PS.J2Track(h["t"], h["a"], h["e"], np.full(len(h["t"]), i0), om0, Om0, i0)
    return dict(track=tr, days=h["days"], hist=h, tok=tok, N=h["N"], lam=sel["lam"], days_before=before["days"], days_after=after["days"], pen=pen, kappa=kap)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    F = V.load("rev28t", None); rows = []; var = {}
    for (arm, geom, lv, seed), f in sorted(F.items()):
        if arm != "OLRT3W3T": continue
        key = (geom, round(f["tstar"], 6))
        if key not in var:
            tok = f["tok"]; ts = f["tstar"]; v = {}
            v["price"] = plan_for_tstar(tok, ts)                                      # primary: flight penalty, flight kappa
            v["pen0"] = plan_for_tstar(tok, ts, pen=0.0)
            v["kappa1"] = plan_for_tstar(tok, ts, kappa=1.0)
            v["before"] = plan_for_tstar(tok, ts, side="before"); v["after"] = plan_for_tstar(tok, ts, side="after")
            v["shift"] = T28.nominal_for_tstar(tok, ts)                               # earlier constant-shift construction
            v["unshifted"] = PS.nominal(tok)                                          # logged convention: nominal plan at its own arrival
            var[key] = v
        v = var[key]; r = dict(geom=geom, seed=seed, ok=f["ok"], window=f["window"], corr=f["corr"], hard=f["hard"], dv=f["dv_sum"], err=f["err"], tstar=f["tstar"], t_arr=f["t_arr_d"])
        for name, nm in v.items():
            if nm is None: r[name] = None; continue
            q = V.position(f["dir"], f["tok"], f["tstar"], f["rt"], None, nm); r[name] = dict(raw=q["raw"], l1=q["l1"])
        rows.append(r)
    # target-direction spread between variants at t* (deg)
    spread = {}
    for (geom, ts), v in var.items():
        u = {k: nm["track"].apo(ts * 86400.0) for k, nm in v.items() if nm is not None}; base = u["price"]
        spread[geom] = {k: float(math.degrees(PS.ang(base, x))) for k, x in u.items()}
        spread[geom]["plan_days"] = {k: float(nm["days"]) for k, nm in v.items() if nm is not None}; spread[geom]["lam"] = float(v["price"]["lam"]); spread[geom]["days_before_after"] = [float(v["price"]["days_before"]), float(v["price"]["days_after"])]
    (out / "rev30_conv.json").write_text(json.dumps(SC.nan_none(dict(rows=rows, spread=spread)), indent=1, default=float))
    df = pd.DataFrame(rows); print(json.dumps(spread, indent=1))
    for g in ("ODY", "G1230", "G1038", "G0975M"):
        x = df[df.geom == g]
        print(g, {k: (round(float(np.median([r[k]["raw"] for r in x.to_dict("records")])), 2), round(float(np.median([r[k]["l1"] for r in x.to_dict("records")])), 2), int(sum(r[k]["l1"] < 1 for r in x.to_dict("records")))) for k in ("unshifted", "price", "pen0", "kappa1", "before", "after", "shift")})


if __name__ == "__main__":
    main()
