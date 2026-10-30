"""Revision 24 (WP24; analysis only: Python, no Julia, no flights): a set of circular target orbits (Odyssey's plane, high inclination, near-critical 63.4 / 116.6 deg, near-equatorial) at apoapsis altitudes from 500 km to the Odyssey value.

Each target is a circle whose radius equals the closing apoapsis radius of the campaign (the apoapsis boundary condition becomes the circle radius). Per target and per candidate aerobraking geometry:
  aerobraking dV = the Revision 20 DP (planner_rev20_ic.Case via planner_rev23_geometry.dp, S-fixed and S-free, free arrival) to that final apoapsis from the same initial apoapsis (28,533 km), q_req and tube width W of the
  7,408 km campaign kept ("same corridor"); plane term = exact apoapsis-node burn increment at the target's circle speed (inclination difference, node free); circle match = v_circ - v_a at r_a = circle radius.
Extension beyond the existing tables (DERIVED, MODEL EXTRAPOLATION): geometry_rank.schedule / channels run with orbit.ra_target_m lowered; the mean field is NOT available along the longer track, so every ladder level uses the
preview table (fallback scale height, no field term) and the Odyssey level is also given with the field (Revision 23 table) to measure the gap.
Imports are read-only except planner_rev23_geometry (additive keyword arguments ra_t, horizon_d, key rp_final).
Usage: python3 -I planner_rev24_targets.py selftest | run [--jobs 4] | summarize"""
import argparse, contextlib, hashlib, json, math, os, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev23_geometry as X
import planner_rev20_ic as R
import planner_rev22_rendezvous as R22

TD = pathlib.Path(PLB.task_dir())
OUT = TD / "rev24_data"
MU = R22.MU                                   # SOURCED planner config mu (Mars)
R_EQ = 3396190.0                              # SOURCED geometry_rank_config.toml [mars] r_eq_m (pck00011)
RA_ODY = 7408000.0                            # SOURCED stageo_config.toml [stageo] ra_target_m; a RADIUS (see note)
RP_FINAL = 3.504e6                            # ASSUMED (Revision 22): final periapsis radius of the circle-match ellipse; sensitivity ladder RP_SENS_M
RP_SENS_M = (3.5e6, 3.504e6, 3.6e6, 3.646e6)  # ASSUMED sensitivity of the circle match to the final periapsis radius (3.646e6 = 250 km altitude, the Revision 20 h_des)
H_MIN_ALT = 200e3                             # ASSUMED (Revision 22) lowest circle altitude outside the drag pass
ALT_LADDER_KM = (500.0, 1000.0, 2000.0, 3000.0, None)    # ASSUMED ladder (5 levels, ~ geometric then Odyssey); None = the Odyssey value r = 7,408 km
HIGH_INC_DEG = 80.0                           # ASSUMED high inclination distinct from Odyssey's (93.522 deg is already high); present in the fine mean-field grid
NEAR_EQ_DEG = 20.0                            # ASSUMED near-equatorial inclination (lowest of the fine grid; the 1,536 grid minimum is 25.2); present in the fine mean-field grid
CRIT_DEG = (63.4, 116.6)                      # SOURCED geometry_rank_config.toml inclination classes (the Revision 23 zero-cost classes)
NODE_POLICY = "free"                          # ASSUMED: target node unconstrained (node-fixed needs Omega0-matched samples, missing)
MAX_PASSES_EXT = 5000                         # ASSUMED guard of the schedule loop for the extended campaigns
HORIZON_EXT_D = 600.0                         # ASSUMED arrival horizon (days) of the extended campaigns (planner horizon is 130 d)
TOL_FREEDOM = X.TOL_FREEDOM                   # 0.25 m/s ASSUMED, as Revision 23
CHEAP_MPS = 0.5                               # ASSUMED: a target is "cheap" if the best net aerobraking + plane dV is at most this (Revision 23 falsifier level)
INC_TARGETS = [("ODY", 93.522, "SOURCED stageo_config.toml [stageo] inclination_deg (odyssey_tolson.toml reference elements)"),
               ("HIGH", HIGH_INC_DEG, "ASSUMED named parameter HIGH_INC_DEG (Odyssey's 93.522 deg is itself high; 80 deg is the distinct high-inclination class, sampled in meanfield_fine)"),
               ("CRIT63", 63.4, "SOURCED geometry_rank_config.toml inclinations_deg (near-critical, Revision 23 zero-cost class)"),
               ("CRIT117", 116.6, "SOURCED geometry_rank_config.toml inclinations_deg (near-critical retrograde, Revision 23 zero-cost class)"),
               ("EQ", NEAR_EQ_DEG, "ASSUMED named parameter NEAR_EQ_DEG (lowest inclination of the fine mean-field grid)")]


def alt_of(a_km):
    return RA_ODY - R_EQ if a_km is None else a_km * 1e3


def circle_costs(alt_m, rp_m=RP_FINAL):
    """Circle radius r = R_EQ + alt = apoapsis radius; v_circ, v_a of the final ellipse (r_p given), V0 = v_circ - v_a at the apoapsis. DERIVED."""
    r = R_EQ + alt_m; vc = math.sqrt(MU / r); va = math.sqrt(MU * (2 / r - 2 / (r + rp_m)))
    return dict(r_m=r, v_circ=vc, v_a=va, V0=vc - va, e_final=(r - rp_m) / (r + rp_m), reachable=bool(rp_m < r and alt_m >= H_MIN_ALT))


def plane_exact(di_deg, vc, va):
    """Increment over the coplanar match of an apoapsis-node burn with the target velocity rotated by di: sqrt((vc - va)^2 + 4 vc va sin^2(di/2)) - (vc - va). DERIVED."""
    return math.sqrt((vc - va) ** 2 + 4 * vc * va * math.sin(math.radians(di_deg) / 2) ** 2) - (vc - va)


@contextlib.contextmanager
def orbit(ctx, ra_t):
    o = ctx.cfg["orbit"]; old = (o["ra_target_m"], o["max_passes"]); o["ra_target_m"] = float(ra_t); o["max_passes"] = MAX_PASSES_EXT
    try:
        yield
    finally:
        o["ra_target_m"], o["max_passes"] = old


def build_ext(cand, field, ra_t, preview=True):
    """Channels of the campaign to apoapsis radius ra_t at the q_req (heat level) of the 7,408 km campaign (same corridor); preview = no field, fallback scale height."""
    G = X.init(); ctx = G["ctx"]; GR = G["GR"]; b = ctx.cfg["band"]; H = None if (preview or field is None) else field["H"]
    q_req = GR.solve_q_req(ctx, cand, b["d_ref_days"], H)
    with orbit(ctx, ra_t):
        S = GR.schedule(ctx, cand, q_req, H); d = GR.channels(ctx, cand, S, None if preview else field)
    return dict(q_req=q_req, S=S, d=d)


def level_row(cand, field, a_km, preview):
    ra_t = R_EQ + alt_of(a_km); bd = build_ext(cand, field, ra_t, preview); tab = X.table_from(bd); G = X.G
    row = dict(n=int(tab["n"]), days_q_req=float(G["GR"].campaign_days(bd["S"])), q_req=float(bd["q_req"]), W=float(tab["W"]), e_end=float(bd["S"]["e"][-1]), rp_end=float(bd["S"]["rp"][-1]))
    for sc in ("S-fixed", "S-free"):
        r = X.dp(tab, sc, ra_t=ra_t, horizon_d=HORIZON_EXT_D)
        row[sc] = {k: r.get(k) for k in ("feasible", "why", "aero", "campaign", "pre", "N", "days", "burns", "rp_final")}
    return row


def analyze(task):
    G = X.init(); cand = task["cand"]; src = task["src"]; mf = G["mf"][src] if src != "fine" else X.load_fine(); field = mf.get(task["fid"]) if task["fid"] else None
    t0 = time.time(); out = dict(id=task["id"], fid=task["fid"], src=src, inc=cand["inc"], omega0=cand["omega"], raan0=cand["raan"], epoch=cand["epoch_off"], role=task.get("role"), levels={})
    for a_km in ALT_LADDER_KM:
        key = "ODY" if a_km is None else "%d" % a_km
        out["levels"][key] = dict(preview=level_row(cand, field, a_km, True))
        if a_km is None and field is not None:
            out["levels"][key]["field"] = level_row(cand, field, a_km, False)
    out["wall_s"] = time.time() - t0; return out


def tasks():
    G = X.init(); GR = G["GR"]; T = X.cut256_tasks(); seen = {t["id"] for t in T}
    for inc in (HIGH_INC_DEG, NEAR_EQ_DEG):                 # cut-style subsets (omega0, Omega0 every 90 deg, 4 epochs) of the fine grid at the added inclinations
        for c in G["fine"]:
            if c["inc"] == inc and abs(c["omega"] / 90 - round(c["omega"] / 90)) < 1e-9 and abs(c["raan"] / 90 - round(c["raan"] / 90)) < 1e-9:
                T.append(dict(id=c["id"], fid=c["id"], cand=c, src="fine", level="full"))
    for t in X.anchor_tasks():
        if t["id"] not in seen: T.append(t)
    return T


def cmd_run(jobs):
    X.init(); X.load_fine(); T = tasks(); print("candidates", len(T), flush=True)
    rows = X.pmap(analyze, T, jobs, "rev24"); sha = X.write_json(OUT / "raw_rev24.json", rows); print("sha", sha[:16])


def selftest():
    G = X.init(); ctx = G["ctx"]
    c = circle_costs(RA_ODY - R_EQ)
    Fae = R22.ell(RP_FINAL, RA_ODY); assert abs(R22.apoapsis_cost(Fae, (RA_ODY, 0.0)) - c["V0"]) < 0.05      # same number as Revision 22 (read-only)
    assert abs(c["V0"] - 477.54) < 0.05, c                                           # Revision 22b: v_circ 2,404.4 - v_a 1,926.9
    assert abs(plane_exact(1.0, c["v_circ"], c["v_a"]) + c["V0"] - 479.0) < 0.1
    v = [circle_costs(a_ * 1e3 if a_ else RA_ODY - R_EQ)["V0"] for a_ in ALT_LADDER_KM]; assert all(x < y for x, y in zip(v, v[1:])), v
    # extension reproduces the standard table at 7,408 km and the DP state is restored
    cand = G["GR"].anchor_cand(ctx.cfg); fld = G["mf"]["cut"]["ODY"]
    t0 = X.table_from(X.build(cand, fld)); t1 = X.table_from(build_ext(cand, fld, RA_ODY, preview=False))
    assert np.array_equal(t0["path"], t1["path"]) and t0["n"] == t1["n"]
    a0 = X.dp(t0, "S-fixed"); a1 = X.dp(t1, "S-fixed", ra_t=RA_ODY, horizon_d=HORIZON_EXT_D); assert abs(a0["aero"] - a1["aero"]) < 1e-12
    assert R.RA_T == RA_ODY and R.HORIZON_D == 130.0 and ctx.cfg["orbit"]["ra_target_m"] == RA_ODY and ctx.cfg["orbit"]["max_passes"] == 400
    # a lower target: more passes, finite, arrival at the lower radius
    t2 = X.table_from(build_ext(cand, fld, 3.896e6, preview=True)); assert t2["n"] > t0["n"] and t2["ra"][-1] < 3.97e6
    OUT.mkdir(parents=True, exist_ok=True); (OUT / "selftest24.txt").write_text("planner_rev24_targets selftest ok\n"); print("planner_rev24_targets selftest ok")


# ----------------------------------------------------------------------------------------------------- summaries

def _v(d):
    return d["aero"] if d and d.get("feasible") else float("nan")


def long_table(rows):
    out = []
    for r in rows:
        for key, lv in r["levels"].items():
            for tv, p in lv.items():
                for sc in ("S-fixed", "S-free"):
                    d = p[sc]; out.append(dict(id=r["id"], role=r.get("role"), inc=r["inc"], omega0=r["omega0"], raan0=r["raan0"], epoch=r["epoch"], level=key, table=tv, sc=sc, aero=_v(d), feasible=bool(d.get("feasible")), why=d.get("why"),
                                               N=d.get("N"), days=d.get("days"), rp_final=d.get("rp_final"), n=p["n"], days_q_req=p["days_q_req"], e_end=p["e_end"]))
    return pd.DataFrame(out)


def spearman(x, y):
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[m], y[m])[0]) if m.sum() > 3 else float("nan")


def summarize():
    rows = json.load(open(OUT / "raw_rev24.json")); L = long_table(rows); S = {}; ody_alt = RA_ODY - R_EQ
    levels = [("ODY" if a is None else "%d" % a, alt_of(a)) for a in ALT_LADDER_KM]
    # --- circle costs and the dependence on the final periapsis (P24-1)
    S["circle"] = {k: dict(alt_km=alt / 1e3, **circle_costs(alt), V0_by_rp={("%.0f" % (rp / 1e3)): circle_costs(alt, rp)["V0"] for rp in RP_SENS_M}) for k, alt in levels}
    v0 = [S["circle"][k]["V0"] for k, _ in levels]; S["P24_1"] = dict(V0_500=v0[0], monotonic=bool(all(x < y for x, y in zip(v0, v0[1:]))), V0_ody=v0[-1])
    # --- existing tables (P24-2)
    mins = {}
    for f in ("planner_ground_tables.json", "planner_ground_tables_generic.json", "planner_ground_tables_rev23.json"):
        for k, g in json.load(open(TD / f))["geometries"].items(): mins["%s:%s" % (f, k)] = float(g["ra"][-1])
    S["P24_2"] = dict(min_final_ra_km=min(mins.values()) / 1e3, n_tables=len(mins), any_below_7000=bool(min(mins.values()) < 7.0e6), per_table_km={k: v / 1e3 for k, v in mins.items()})
    # --- campaign length (P24-3): same corridor (q_req of the 7,408 km campaign), preview tables
    pv = L[(L.table == "preview") & (L.sc == "S-fixed")].drop_duplicates(["id", "level"]); w = pv.pivot(index="id", columns="level", values=["n", "days_q_req"])
    S["P24_3"] = {k: dict(passes_ratio_median=float((w["n"][k] / w["n"]["ODY"]).median()), days_ratio_median=float((w["days_q_req"][k] / w["days_q_req"]["ODY"]).median()), days_ratio_min=float((w["days_q_req"][k] / w["days_q_req"]["ODY"]).min()),
                       days_ratio_max=float((w["days_q_req"][k] / w["days_q_req"]["ODY"]).max()), n_median=float(w["n"][k].median()), days_median=float(w["days_q_req"][k].median())) for k, _ in levels}
    ao = L[(L.id == "ODY") & (L.table == "preview") & (L.sc == "S-fixed")].set_index("level")
    S["P24_3"]["odyssey_anchor"] = {k: dict(n=int(ao.loc[k, "n"]), days_q_req=float(ao.loc[k, "days_q_req"]), N_dp=ao.loc[k, "N"], days_dp=ao.loc[k, "days"]) for k, _ in levels}
    # --- field vs preview at the Odyssey level
    fo = L[(L.level == "ODY") & (L.table == "field")]; po = L[(L.level == "ODY") & (L.table == "preview")]; m = fo.merge(po, on=["id", "sc"], suffixes=("_f", "_p"))
    S["field_vs_preview_ody"] = {sc: dict(n=int(len(g)), spearman=spearman(g.aero_f, g.aero_p), median_shift_preview_minus_field=float((g.aero_p - g.aero_f).median()), max_abs_shift=float((g.aero_p - g.aero_f).abs().max())) for sc, g in m.groupby("sc")}
    # --- P24-4: best per inclination class at 7,408 km with the field, node free
    best = {}
    for sc in ("S-fixed", "S-free"):
        for inc, g in fo[fo.sc == sc].groupby("inc"):
            best.setdefault(str(inc), {})[sc] = dict(best=float(g.aero.min()), best_id=str(g.loc[g.aero.idxmin(), "id"]), n=int(g.aero.notna().sum()), median=float(g.aero.median()))
    S["P24_4"] = best
    # --- per target
    T = []
    for tag, inc_t, prov in INC_TARGETS:
        for k, alt in levels:
            tc = circle_costs(alt); ent = dict(target_id="%s_%s" % (tag, k if k == "ODY" else "%dkm" % int(alt / 1e3)), class_tag=tag, inc_deg=inc_t, inc_provenance=prov, alt_km=alt / 1e3, circle_r_km=tc["r_m"] / 1e3,
                                               alt_provenance=("DERIVED: ra_target_m 7,408 km (a radius, stageo_config.toml) minus R_eq 3,396.19 km" if k == "ODY" else "ASSUMED ladder ALT_LADDER_KM"), node_deg=None, node_policy=NODE_POLICY,
                                               reachable=tc["reachable"], e_final=tc["e_final"], V0_circle_match=tc["V0"], v_circ=tc["v_circ"], v_a=tc["v_a"], V0_by_final_rp=S["circle"][k]["V0_by_rp"])
            tbl = "field" if k == "ODY" else "preview"
            for sc in ("S-fixed", "S-free"):
                g = L[(L.level == k) & (L.table == tbl) & (L.sc == sc)].copy(); g["di"] = (g.inc - inc_t).abs(); g["plane"] = [plane_exact(d, tc["v_circ"], tc["v_a"]) for d in g.di]; g["net"] = g.aero + g.plane
                ok = g[np.isfinite(g.net)]; b = ok.loc[ok.net.idxmin()]; cls = ok[ok.di < 0.1]; ody = g[g.id == "ODY"].iloc[0]; bc = cls.loc[cls.aero.idxmin()] if len(cls) else None
                ent[sc] = dict(n_feasible=int(len(ok)), n_candidates=int(len(g)), best_id=str(b.id), best_inc=float(b.inc), best_aero=float(b.aero), best_plane=float(b.plane), best_net=float(b.net), best_days=(None if b.days is None else float(b.days)),
                               class_n=int(len(cls)), class_best_id=(None if bc is None else str(bc.id)), class_best_aero=(None if bc is None else float(bc.aero)), class_median=(None if bc is None else float(cls.aero.median())),
                               class_share_within_tol=(None if bc is None else float((cls.aero <= bc.aero + TOL_FREEDOM).mean())), grid_share_within_tol=float((ok.net <= b.net + TOL_FREEDOM).mean()),
                               ody_geometry_aero=float(ody.aero), ody_geometry_plane=float(ody.plane), ody_geometry_net=float(ody.net), saving_vs_ody_geometry=float(ody.net - b.net), n_beat_ody=int((ok.net < ody.net).sum()),
                               total_best=float(b.net + tc["V0"]), total_ody_geometry=float(ody.net + tc["V0"]), aero_only_best=float(b.aero))
            cheap = ent["S-fixed"]["best_net"] <= CHEAP_MPS
            ent["classification"] = "not reachable within the corridor" if (not tc["reachable"] or ent["S-fixed"]["n_feasible"] == 0) else ("cheap (aerobraking + plane <= %.1f m/s)" % CHEAP_MPS if cheap else "reachable, aerobraking-costly (best net > %.1f m/s)" % CHEAP_MPS)
            ent["plane_change_if_node_fixed"] = "node must be matched by Omega0 (continuous); the sampled grid pays >= 257 m/s at 93.5 deg (Revision 23); matched-node values need new mean-field samples"
            T.append(ent)
    S["targets"] = T
    # --- missing samples and cost (derived from the extended path lengths)
    n500 = pv[pv.level == "500"].set_index("id").n; nod = pv[pv.level == "ODY"].set_index("id").n; extra = (n500 - nod).clip(lower=0)
    S["missing_samples"] = dict(candidates=int(len(extra)), extra_passes_median=float(extra.median()), extra_passes_max=int(extra.max()), calls_per_pass=3, extra_calls_total=int(3 * extra.sum()), est_call_seconds=float(3 * extra.sum() * 0.157e-3),
                                note="3 altitudes per pass (SOURCED sampling alt_offsets_m), 0.157 ms per call (SOURCED probe); the tracks of the shallower targets are prefixes of the 500 km track")
    X.write_json(OUT / "summary24.json", S)
    ent_json = [{k: v for k, v in e.items()} for e in T]
    X.write_json(TD / "target_orbits.json", dict(schema="target_orbits/1", note="Revision 24; circular targets, radius = closing apoapsis radius; costs DERIVED (model world), preview tables below the Odyssey altitude", RA_ody_m=RA_ODY, R_eq_m=R_EQ, RP_FINAL_m=RP_FINAL, targets=ent_json))
    md = ["| target | i (deg) | alt (km) | r (km) | e final | circle match V0 (m/s) | best geometry (S-fixed) | aero (m/s) | plane (m/s) | total = aero + plane + V0 | Odyssey geometry total | class |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for e in T:
        f = e["S-fixed"]; md.append("| %s | %.3f | %.0f | %.0f | %.3f | %.1f | %s (i %.1f) | %.3f | %.3f | %.1f | %.1f | %s |" % (e["target_id"], e["inc_deg"], e["alt_km"], e["circle_r_km"], e["e_final"], e["V0_circle_match"], f["best_id"], f["best_inc"], f["best_aero"], f["best_plane"], f["total_best"], f["total_ody_geometry"], e["classification"]))
    (TD / "target_orbits.md").write_text("# Revision 24 target orbits (circular; radius = closing apoapsis radius)\n\nAerobraking and plane costs DERIVED from the Revision 20 DP (S-fixed, free arrival, node free); preview tables below the Odyssey altitude; V0 at r_p = 3,504 km (ASSUMED). Provenance per entry in target_orbits.json.\n\n" + "\n".join(md) + "\n")
    print("summary written"); return S


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "run", "summarize"]); ap.add_argument("--jobs", type=int, default=4); a = ap.parse_args()
    if a.cmd == "selftest": selftest()
    elif a.cmd == "run": cmd_run(a.jobs)
    else: summarize()
