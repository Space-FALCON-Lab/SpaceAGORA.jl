"""Revision 28C (analysis only: Python, no Julia, no flights, nothing in src/): insertion-plane / node matching.

Idea: with a target plane (i_t, Omega_t at t*) the mission inserts at inclination i_t and at the node Omega_0 = Omega_t - dOmega_J2, where dOmega_J2 is the J2 node drift over the campaign computed from the nominal plan's (a, e)
history (per pass: -3 pi J2 (R/p)^2 cos i, p = a (1 - e^2); same law as geometry_rank.schedule). The aerobraking geometry is then the DP over the remaining free elements (omega0, epoch, x0).

Parts: (a) required dOmega and Omega_0 per class and target (25 targets of target_orbits.json); (b) residual node error at arrival from the flown spread (rev21v, rev23geo, rev25low, rev26val2, rev27low), split into the J2 history term,
the timing term and the model bias; (c) residual plane-change cost (exact apoapsis-node law, rev22 tools); (d) DP aerobraking dV over (omega0, epoch) with i and Omega_0 fixed, as proxies from the sampled Omega_0 (nearest and
interpolated) against the free-geometry best; (e) totals against the in-orbit plane change (Revision 26a).
Approach-side constraints are NAMED ASSUMED parameters with NO value (None): the achievable insertion inclination from the arrival-asymptote declination and the approach cost of a node choice. Results are conditional on them.
Usage: python3 -I planner_rev28_plane.py selftest | run"""
import argparse, glob, hashlib, json, math, os, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev23_geometry as X
import planner_rev20_ic as R
import planner_rev22_rendezvous as R22
import planner_select_geometry as SG

TD = pathlib.Path(PLB.task_dir()); OUT = TD / "rev28_data"
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
# ---- named ASSUMED parameters (no invented values: None = not provided, result conditional)
APPROACH_MIN_INC_DEG = None          # ASSUMED: lowest/highest insertion inclination achievable from the arrival-asymptote declination (not sourced; None = unconstrained)
APPROACH_NODE_COST_MPS = None        # ASSUMED: approach-side cost of choosing the insertion node (not sourced; None = not charged: results are conditional on it)
TARGET_NODES_DEG = tuple(range(0, 360, 45))   # ASSUMED sampling of the target node Omega_t (not given by any file)
EXAMPLE_NODE_DEG = 33.67             # ASSUMED example target node (Odyssey model arrival node, Revision 23); used only to print an Omega_0 example
TOL_FREEDOM = X.TOL_FREEDOM          # 0.25 m/s ASSUMED (Revision 23)
NODE_TOL_DEG = 1.0                   # SOURCED P28-6 / Revision 26 (position tolerance 1 deg)
IN_ORBIT_PLANE_MPS = (17.0, 25.0)    # SOURCED Revision 26a: in-orbit plane change at G1 (first order, rev22 tools)
IN_ORBIT_PHASING_MPS = (15.0, 78.0)  # SOURCED Revision 26a: phasing at G1
KAPPA = {"ODY": 0.9967520763650508, "G5863": 0.879308559055656}    # MEASURED planner_config kappa; 1.0 (ASSUMED) elsewhere
TOK_ALT_KM = {"G0486": 2000.0, "G1230": 1000.0, "G1038": 500.0}   # Revision 25 block cells; every other token ends at 7,408 km
TABLE_FILES = [("planner_ground_tables.json", ("ODY", "G5863")), ("planner_ground_tables_generic.json", ("G0975", "G0063", "G0390", "G1216")), ("planner_ground_tables_rev23.json", ("G0528", "G0890", "G0975M", "G0063M")),
               ("planner_ground_tables_rev25mf.json", ("G0486", "G1230", "G1038"))]
BLOCKS = ("rev21v", "rev23geo", "rev25low", "rev26val2", "rev27low")
CLASS_TAGS = {"ODY": 93.522, "HIGH": 80.0, "CRIT63": 63.4, "CRIT117": 116.6, "EQ": 20.0}


def target_ra(tok):
    return SG.radius_of(TOK_ALT_KM[tok]) if tok in TOK_ALT_KM else float(SG.AN.RA_TARGET)


def load_tabs():
    out = {}
    for f, keys in TABLE_FILES:
        p = TD / f
        if not p.exists(): continue
        g = json.load(open(p))["geometries"]
        for k in keys:
            if k in g: out[k] = {a: (np.array(b, float) if isinstance(b, list) else b) for a, b in g[k].items()}
    return out


# ----------------------------------------------------------------------------------------------------- the plan's (a, e) history and the J2 node drift

def j2_node_step(ctx, a, e, inc_rad):
    p = a * (1 - e * e); return -3 * math.pi * ctx.J2 * (ctx.Rg / p) ** 2 * math.cos(inc_rad)


def plan_for(tab, ra_t, kappa=1.0, scenario="S-fixed"):
    """S-fixed DP plan on the table to closing radius ra_t (arrival horizon, pen 0). Returns dict(x, N, days) or None."""
    with X.patched(tab, ra_t=ra_t, horizon_d=600.0):
        try:
            c = R.Case("X", scenario, 250e3, False, kappa, "arrival", 0.0, 0.0); p = c.nominal
        except Exception:
            return None
    return None if p is None else dict(x=np.asarray(p["x"], float), N=int(p["N"]), days=float(p["days"]), s=float(p["s"]), dv=float(p["dv_campaign"] + p["dv_pre"]))


def history(tab, x, kappa, ra_t, q_first_ln=None):
    """Replicates planner_rev18_offline.arrival and records the (ra, rp, a, e, t) of every pass up to the closing pass."""
    q1 = math.log(R.Q_FIRST) if q_first_ln is None else q_first_ln; O = R.O; n = tab["n"]; q_req = tab["q_req"]; ra = float(tab["ra"][0]); t = float(tab["t"][0]); path = tab["path"]
    d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); ra_t_tab = tab["ra"][::-1]; d0_r = d0[::-1]; rows = []
    for k in range(1, int(n * 4) + 1):
        kk = min(k, n); xk = float(x[min(k, len(x)) - 1]); nu = xk - (q1 + float(path[kk - 1] - path[0])); rp = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nu; a = 0.5 * (ra + rp)
        drop = kappa * float(np.interp(ra, ra_t_tab, d0_r)) * math.exp(xk) / q_req; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2 * math.pi * math.sqrt(amid ** 3 / O.MU)
        rows.append((k, ra, rp, a, (ra - rp) / (ra + rp), t))
        if ra_next <= ra_t: return dict(rows=np.array(rows), N=k, days=(t + 0.5 * Tk) / 86400.0, T_final=Tk)
        ra = ra_next; t += Tk
    return None


def node_drift_deg(ctx, hist, inc_deg):
    """J2 node drift from the first pass to the closing apoapsis (deg): sum of the per-pass steps over passes 1..N-1 plus half of the last one."""
    r = hist["rows"]; st = np.array([j2_node_step(ctx, r[k, 3], r[k, 4], math.radians(inc_deg)) for k in range(len(r))])
    return math.degrees(st[:-1].sum() + 0.5 * st[-1]), math.degrees(st[-1]) / (hist["T_final"] / 86400.0)       # deg, final-orbit rate deg per day


# ----------------------------------------------------------------------------------------------------- (b) flown residuals

def finished(d):
    return all((d / f).exists() for f in ("passes.csv", "decisions.csv", "measured.csv", "replan.csv"))


def flown_residuals(ctx, tabs):
    rows = []; plans = {}
    for blk in BLOCKS:
        for d in sorted((ROOT / blk).glob("*")) if (ROOT / blk).exists() else []:
            nm = d.name; gid = nm.split("_PL_")[0]; tok = {"G5863": "G5863"}.get(gid, gid)
            if tok not in tabs or not finished(d) or "_OLRT3" not in nm: continue
            P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra_t = target_ra(tok); hit = np.where(P.ra_m.values <= ra_t)[0]
            if not len(hit): continue
            n = int(hit[0]) + 1; inc = P.i_rad.values; raan = np.unwrap(P.raan_rad.values); a = P.a_m.values; e = P.e.values; kap = KAPPA.get(tok, 1.0)
            if tok not in plans:
                pl = plan_for(tabs[tok], ra_t, kap); plans[tok] = None if pl is None else dict(pl, hist=history(tabs[tok], pl["x"], kap, ra_t))
            pl = plans[tok]
            if pl is None: continue
            inc_deg = math.degrees(inc[0]); nom, rate = node_drift_deg(ctx, pl["hist"], inc_deg)
            st = np.array([j2_node_step(ctx, a[k], e[k], inc[k]) for k in range(n)]); j2_flown_hist = math.degrees(st[:-1].sum() + 0.5 * st[-1]); flown = math.degrees(raan[n - 1] - raan[0]) + 0.5 * math.degrees(st[-1])
            t_arr = (P.time_s.values[n - 1] + 0.5 * 2 * math.pi * math.sqrt(a[n - 1] ** 3 / R.O.MU)) / 86400.0; dt = t_arr - pl["days"]
            res = flown - nom; timing = rate * dt; hist_term = (j2_flown_hist - nom) - timing; bias = flown - j2_flown_hist
            di_node = 2 * math.degrees(math.asin(math.sin(math.radians(inc_deg)) * math.sin(math.radians(abs(res)) / 2)))
            di_inc = math.degrees(inc[n - 1] - inc[0]); di = math.hypot(di_node, di_inc); rp = float(P.rp_m.values[n - 1]); ra_c = float(P.ra_m.values[n - 1])
            rows.append(dict(block=blk, run=nm, token=tok, inc_deg=inc_deg, n=n, N_nom=pl["N"], t_arr_d=t_arr, t_nom_d=pl["days"], dt_d=dt, drift_nom_deg=nom, drift_flown_deg=flown, residual_deg=res, timing_deg=timing, history_deg=hist_term, bias_deg=bias,
                             rate_final_deg_d=rate, di_deg=di, di_node_deg=di_node, di_inc_deg=di_inc, plane_cost_mps=SG.exact_plane_increment(di, ra_c, rp), ra_close_m=ra_c, rp_close_m=rp))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------- (a) targets and (d) proxies

def target_table(ctx):
    """Required node offset per target (best geometry of the target from target_orbits.json: S-fixed plan, plan history, J2 drift)."""
    T = json.load(open(TD / "target_orbits.json"))["targets"]; G = X.init(); out = []; cache = {}
    for e in T:
        cid = e["S-fixed"]["best_id"]; ra_t = e["radius_m"]; inc = CLASS_TAGS[e["class_tag"]]
        c = next((c for c in G["grid"] + G["fine"] if c["id"] == cid), None)
        if c is None: continue
        key = (cid, ra_t)
        if key not in cache:
            ody = abs(ra_t - float(SG.AN.RA_TARGET)) < 1.0
            if ody:
                tab, info = SG.candidate_table(c, None, preview=False)
            else:
                tab, info = SG.candidate_table(c, ra_t, preview=True)
            pl = plan_for(tab, ra_t, 1.0); cache[key] = None if pl is None else (pl, history(tab, pl["x"], 1.0, ra_t), info)
        v = cache[key]
        if v is None: out.append(dict(target_id=e["target_id"], feasible=False)); continue
        pl, h, info = v; dO, rate = node_drift_deg(ctx, h, inc)
        out.append(dict(target_id=e["target_id"], class_tag=e["class_tag"], inc_deg=inc, alt_km=e["alt_km"], radius_m=ra_t, best_id=cid, N=pl["N"], t_arr_d=pl["days"], node_drift_deg=dO, rate_final_deg_d=rate, field_mode=info["field_mode"],
                        omega0_example_deg=(EXAMPLE_NODE_DEG - dO) % 360.0, plan_aero_mps=pl["dv"]))
    return out


def fixed_node_proxy(targets):
    """(d) DP aerobraking dV with i and Omega_0 fixed, from the sampled Omega_0 of each (omega0, epoch) cell: nearest sample and linear interpolation between the two bracketing samples (periodic), for the 8 target nodes; the
    free best is the class best over all geometries. S-fixed aero from the Revision 23 (7,408 km, field) and Revision 24 (ladder, preview) raw runs."""
    r23 = json.load(open(TD / "rev23_data" / "raw_cut256.json")) + json.load(open(TD / "rev23_data" / "raw_grid1536_rest.json")); r24 = json.load(open(TD / "rev24_data" / "raw_rev24.json")); rows = []
    A = {}
    for r in r23:
        v = r["dp"]["S-fixed"]; A.setdefault(("ODY", r["id"]), dict(inc=r["inc"], omega0=r["omega0"], raan0=r["raan0"], epoch=r["epoch_off"], aero=v["aero"] if v.get("feasible") else float("nan"), src="rev23 field, Omega_0 step 45 deg (90 deg at 80/20)"))
    for r in r24:
        for lv, d in r["levels"].items():
            p = d.get("field") if lv == "ODY" and "field" in d else d["preview"]; v = p["S-fixed"]
            if lv == "ODY" and ("ODY", r["id"]) in A: continue
            A[(lv, r["id"])] = dict(inc=r["inc"], omega0=r["omega0"], raan0=r["raan0"], epoch=r["epoch"], aero=v["aero"] if v.get("feasible") else float("nan"), src="rev24 %s, Omega_0 step 90 deg" % ("field" if lv == "ODY" else "preview"))
    df = pd.DataFrame([dict(level=k[0], id=k[1], **v) for k, v in A.items()]); df = df[df.id != "ODY"]
    for t in targets:
        if not t.get("class_tag"): continue
        lv = "ODY" if abs(t["alt_km"] - 4011.81) < 1 else "%d" % int(t["alt_km"]); g = df[(df.level == lv) & ((df.inc - t["inc_deg"]).abs() < 0.1)]
        if not len(g): continue
        free = float(g.aero.min()); res = dict(target_id=t["target_id"], level=lv, n_cand=int(len(g)), free_best=free, node_steps=sorted(set(g.raan0)), src=g.src.iloc[0], per_node={})
        nst = sorted(set(g.raan0)); cells = g.groupby(["omega0", "epoch"])
        for Ot in TARGET_NODES_DEG:
            req = (Ot - t["node_drift_deg"]) % 360.0; best_n = []; best_i = []
            for _, c in cells:
                c = c.sort_values("raan0"); xs = c.raan0.values; ys = c.aero.values
                if not np.isfinite(ys).any(): continue
                dist = np.abs(((xs - req + 180) % 360) - 180); k = int(np.argmin(dist)); best_n.append(ys[k])
                xs2 = np.concatenate((xs, [xs[0] + 360])); ys2 = np.concatenate((ys, [ys[0]])); rq = req if req >= xs[0] else req + 360
                best_i.append(float(np.interp(rq, xs2, ys2)))
            bn = np.array(best_n); bi = np.array(best_i)
            res["per_node"][str(Ot)] = dict(omega0_req=req, best_nearest=float(np.nanmin(bn)), best_interp=float(np.nanmin(bi)), share_cells_nearest_le_tol=float(np.mean(bn <= free + TOL_FREEDOM)), share_cells_interp_le_tol=float(np.mean(bi <= free + TOL_FREEDOM)))
        pn = res["per_node"].values(); res["median_best_nearest"] = float(np.median([p["best_nearest"] for p in pn])); res["worst_best_nearest"] = float(np.max([p["best_nearest"] for p in pn]))
        res["median_best_interp"] = float(np.median([p["best_interp"] for p in pn])); res["worst_best_interp"] = float(np.max([p["best_interp"] for p in pn]))
        res["median_share_good_nearest"] = float(np.median([p["share_cells_nearest_le_tol"] for p in pn])); res["median_share_good_interp"] = float(np.median([p["share_cells_interp_le_tol"] for p in pn]))
        rows.append(res)
    return rows


def run():
    X.init(); ctx = X.G["ctx"]; OUT.mkdir(parents=True, exist_ok=True); tabs = load_tabs()
    fr = flown_residuals(ctx, tabs); fr.to_csv(OUT / "flown_residuals.csv", index=False)
    cls = fr.assign(cls=lambda d: d.token.map(lambda t: {"ODY": "ODY 93.5", "G0890": "93.5 (ω0 90)", "G0975": "93.5 (G0975)", "G1038": "93.5 (500 km)", "G0528": "63.4", "G0486": "63.4 (2000 km)", "G5863": "116.6 (G1)", "G1230": "116.6 (1000 km)", "G0063": "25.2 (G0063)", "G0063M": "25.2 (G0063M)"}.get(t, t)))
    agg = {}
    for tok, g in fr.groupby("token"):
        agg[tok] = dict(n=int(len(g)), inc=float(g.inc_deg.iloc[0]), drift_nom=float(g.drift_nom_deg.median()), res_med=float(g.residual_deg.median()), res_abs_med=float(g.residual_deg.abs().median()), res_abs_max=float(g.residual_deg.abs().max()), res_sd=float(g.residual_deg.std(ddof=1)) if len(g) > 1 else None,
                        timing_abs_med=float(g.timing_deg.abs().median()), history_med=float(g.history_deg.median()), history_abs_med=float(g.history_deg.abs().median()), bias_med=float(g.bias_deg.median()), bias_sd=float(g.bias_deg.std(ddof=1)) if len(g) > 1 else None,
                        dt_sd_d=float(g.dt_d.std(ddof=1)) if len(g) > 1 else None, di_med=float(g.di_deg.median()), di_max=float(g.di_deg.max()), cost_med=float(g.plane_cost_mps.median()), cost_max=float(g.plane_cost_mps.max()),
                        within_1deg=int((g.residual_deg.abs() < NODE_TOL_DEG).sum()), res_after_bias_abs_med=float((g.residual_deg - g.bias_deg.median()).abs().median()))
    tt = target_table(ctx); px = fixed_node_proxy(tt)
    json.dump(dict(flown=agg, targets=tt, proxy=px, params=dict(APPROACH_MIN_INC_DEG=APPROACH_MIN_INC_DEG, APPROACH_NODE_COST_MPS=APPROACH_NODE_COST_MPS, TARGET_NODES_DEG=list(TARGET_NODES_DEG), TOL_FREEDOM=TOL_FREEDOM)), open(OUT / "rev28_plane.json", "w"), indent=1, default=float)
    print(json.dumps(agg, indent=0, default=float)[:3000]); print(len(tt), len(px))


def selftest():
    G = X.init(); ctx = G["ctx"]; tabs = load_tabs(); tab = tabs["ODY"]; pl = plan_for(tab, float(SG.AN.RA_TARGET), KAPPA["ODY"]); h = history(tab, pl["x"], KAPPA["ODY"], float(SG.AN.RA_TARGET))
    ref = R.select_ic("ODY", "S-fixed"); assert h["N"] == ref["N"] == pl["N"] and abs(h["days"] - ref["t_arr_d"]) < 1e-9, (h["N"], ref["N"], h["days"], ref["t_arr_d"])
    d, rate = node_drift_deg(ctx, h, 93.522); assert 4.5 < d < 6.5 and abs(rate) > 0.0, (d, rate)       # Odyssey: about 5.4 deg over the campaign (Revision 23)
    assert node_drift_deg(ctx, h, 63.4)[0] < 0 < node_drift_deg(ctx, h, 116.6)[0]                       # prograde regresses, retrograde advances
    assert abs(SG.exact_plane_increment(0.0, 7.408e6, 3.504e6)) < 1e-12
    OUT.mkdir(parents=True, exist_ok=True); (OUT / "selftest28.txt").write_text("planner_rev28_plane selftest ok\n"); print("planner_rev28_plane selftest ok")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "run"]); a = ap.parse_args()
    selftest() if a.cmd == "selftest" else run()
