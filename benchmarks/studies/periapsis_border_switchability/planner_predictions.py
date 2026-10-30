"""WP-P1 / F9 for the planner arms: the ground tables the flight code reads (planner_ground_tables.json) and the pre-registered predictions per arm (planner_predictions.json), both hashed.
Python only; no Julia, no flights. Plan: tasks/2026-10-15_periapsis-region-border-targeting/analytic_framework_and_transfer_plan.md R4.5, R4.10 (F9), R7.3, R8.1, R8.2.

  --dry-run            the flight list (274) against the plan (planner_lib.dry_run)
  --tables             write planner_ground_tables.json (+ .sha256) for the geometries the planner flies (ODY and the G1 geometry of geometry_gate_table.json)
  --predict            write planner_predictions.json (+ .sha256): days, passes, dV, burns, floors per arm and geometry, the OL schedules, the attribution expectations, the time-law feasibility check
  --selftest           tables / predictions invariants on the Odyssey geometry (about 15 s)
Where the numbers come from. The drift path is geometry_rank.py's burn-free ln-q path (J2 apse and node rates, zonal de/dt, the T4a tesseral map, the Sun third-body channel, drag descent Lambda H, velocity loss,
oblate surface, GRAM mean-field term where sampled) at the Q-hold schedule of q_req (the heat rate at which the campaign takes d_ref_days), as stageg3_analysis.model_path. The dV of each arm is the planner's own DP replayed in the
NOISE-FREE model world (planner_lib.simulate): a model-consistent prediction, not a noise-inclusive one. The floors are geometry_rank's tube floors (deterministic and the 12-member phase-jitter ensemble).
Labels: SOURCED (file named), DERIVED (computed here from sourced inputs), ASSUMED (named choice, flagged).
"""
import argparse, json, math, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB

TABLES_NAME = "planner_ground_tables.json"
PRED_NAME = "planner_predictions.json"


def _imports():
    import geometry_rank as GR
    import stageg3_analysis as G3A
    return GR, G3A


def make_ctx_cal(cfg3, cfg):
    """stageg3_analysis.make_ctx with the calibrated planner values (q_cap_mean_w_cm2, d_ref_days from planner_config.toml, the dev-measured ones) applied to the ranking config in memory BEFORE the anchor context is built;
    the frozen geometry_rank_config.toml is not edited. Everything else is make_ctx."""
    GR, G3A = _imports(); gcfg = GR.load_cfg(); gcfg["band"]["q_cap_mean_w_cm2"] = cfg["planner"]["q_cap_mean_w_cm2"]; gcfg["band"]["d_ref_days"] = cfg["planner"]["d_ref_days"]
    ctx = GR.Ctx(gcfg, float(cfg3["stageg3"]["burn_response_ratio"])); GR.anchor_ctx(ctx); return ctx


def geometries(cfg, cfg3):
    """token -> dict(id, cand, g): ODY (anchor), G1, G2, Gworst from the hash-verified gate table."""
    GR, G3A = _imports()
    out = {"ODY": dict(id="ODY", g=dict(cfg["planner"]["geometry"]["ODY"]), cand=GR.anchor_cand(GR.load_cfg()))}
    for g in cfg3["stageg3"]["geometry"]:
        out[g["role"]] = dict(id=g["id"], g=g, cand=G3A.cand_of(g))
    return out


def ground_table(ctx, cand, field, cfg):
    """Ground table of one geometry (lists, JSON-able) and the ranking row. The tube width W uses the ranking's slack_band_factor and q_cap (asserted equal to the planner config's)."""
    GR, G3A = _imports()
    b = ctx.cfg["band"]; P = cfg["planner"]
    assert b["q_cap_mean_w_cm2"] == P["q_cap_mean_w_cm2"] and b["slack_band_factor"] == P["slack_band_factor"] and b["d_ref_days"] == P["d_ref_days"], "ranking context was not calibrated from planner_config.toml"
    assert abs(GR.fast_width(ctx.cfg, ctx.ratio) - P["fast_reserve_ln"]) < 1e-9 and ctx.ratio == P["burn_response_ratio"]
    m = G3A.model_path(ctx, cand, field); S = m["S"]
    W = b["slack_band_factor"] * math.log(b["q_cap_mean_w_cm2"] / m["q_req"])
    path = np.asarray(m["path"], float); slow = np.asarray(m["slow"], float)
    tab = dict(n=int(m["n"]), q_req=float(m["q_req"]), W=float(W), path=path.tolist(), slow=slow.tolist(), d_slow=np.diff(slow).tolist(), fast=(path - slow).tolist(), w=np.asarray(m["w"], float).tolist(),
               H=np.asarray(S["H"], float).tolist(), a=np.asarray(S["a"], float).tolist(), e=np.asarray(S["e"], float).tolist(), t=np.asarray(S["t"], float).tolist(),
               rp=np.asarray(S["rp"], float).tolist(), ra=np.asarray(S["ra"], float).tolist(), vrel=np.asarray(S["V"], float).tolist(), frot=np.asarray(S["frot"], float).tolist())
    return tab, m


def build_tables(cfg, cfg3, ctx, geoms):
    G3A = _imports()[1]
    out = {}
    for tok in ("ODY", "G1"):
        g = geoms[tok]; field, src = G3A.meanfield_for(cfg3, ctx, None, g["g"]); tab, _ = ground_table(ctx, g["cand"], field, cfg)
        tab.update(token=tok, field_mode="meanfield" if field is not None else "preview_no_field", field_source=src, inc_deg=g["cand"]["inc"], omega_deg=g["cand"]["omega"], raan_deg=g["cand"]["raan"], epoch_off_d=g["cand"]["epoch_off"])
        out[g["id"]] = tab
    return out


def write_tables(cfg, cfg3, ctx, geoms, outdir):
    T = build_tables(cfg, cfg3, ctx, geoms)
    doc = dict(schema="planner_ground_tables/1", note="Burn-free ln-q path of the geometry_rank.py model at the Q-hold schedule of q_req (d_ref_days); read by planner_common.jl (planner_table). Index = pass - 1.",
               labels="DERIVED from geometry_rank.py (validated terms, analytic_results.md); constants in planner_config.toml and geometry_rank_config.toml", geometries=T)
    p = pathlib.Path(outdir) / TABLES_NAME; p.write_text(json.dumps(doc, indent=0, default=float)); h = PLB.sha256_file(p)
    (pathlib.Path(outdir) / (TABLES_NAME + ".sha256")).write_text(h + "  " + TABLES_NAME + "\n"); print(f"wrote {p} sha256 {h}"); return h


def _sw(cfg, arm):
    return PLB.arm_switches(cfg, arm)


def predictions(cfg, cfg3, ctx, geoms, tables, with_ens=True):
    GR, G3A = _imports()
    P = PLB.flat_P(cfg); B = cfg["planner"]["body"]; q1 = math.log(cfg["stageo"]["bands"]["odyssey"]["q_target_w_cm2"]); arms = {}; sched_out = {}; rows = {}
    for tok, g in geoms.items():
        field, src = G3A.meanfield_for(cfg3, ctx, None, g["g"])
        rows[tok] = G3A.pre_flight_row(ctx, g["g"], field) if with_ens else dict(floor_det_mps=None, floor_ens_mps=None, floor_ens_sd=None, floor_driftfree_mps=None, regime=None, days=None, n_passes=None, frot_mean=None, vrel_start_mps=None, vrel_end_mps=None, W_eff=None, R_ln=None)
    fo, n2 = cfg["planner"]["arms"]["FO"], cfg["planner"]["arms"]["N2I"]
    for tok in ("ODY", "G1"):
        g = geoms[tok]; tab = tables[g["id"]]; r = rows[tok]
        for arm in ("BODY", "L0", "L1", "L2", "L3", "OL", "OLR", "L3Z", "OLRZ"):
            sw = _sw(cfg, arm) if arm not in ("L3Z", "OLRZ") else _sw(cfg, arm); Pa = PLB.flat_P(cfg, arm); sim = PLB.simulate(sw, tab, Pa, B, q1)
            e = dict(arm=arm, burn_penalty_mps=Pa["dp_burn_penalty_mps"], burn_penalty_label="ASSUMED (user 2026-10-27; one model replay)" if Pa["dp_burn_penalty_mps"] > 0 else "ASSUMED sensitivity arm (none)", geometry=tok, gid=g["id"], mode=sw["mode"], days=float(r["days"]) if r["days"] is not None else tab["t"][-1] / 86400.0, passes=int(tab["n"]),
                     dv_model_world_mps=sim["dv"], burns_model_world=sim["burns"], burn_passes_model_world=sim["burn_passes"], model_violations_above_cap=sim["viol_hi"], model_violations_below_floor=sim["viol_lo"],
                     floor_det_mps=r["floor_det_mps"], floor_ens_mps=r["floor_ens_mps"], floor_ens_sd_mps=r["floor_ens_sd"], floor_driftfree_mps=r["floor_driftfree_mps"], regime=r["regime"], q_req_w_cm2=tab["q_req"], tube_W_ln=tab["W"],
                     basis="DERIVED: planner DP (burn penalty as burn_penalty_mps) replayed in the noise-free model world (planner_lib.simulate) on the ground table; floors from geometry_rank.py")
            if arm == "BODY":
                e.update(days=cfg["planner"]["d_ref_days"], dv_reference_mps=2.924, burns_reference=13, burns_upper_bound=int(cfg["planner"]["d_ref_days"]),
                         basis="B-ODY defines D_ref (64.2 d, DERIVED plan R1.2); passes from the model table (R1.2 says about 200); dV has NO analytic prediction: its burns react to noise that the model world lacks "
                               "(dv_model_world is only an illustration). Reference: the flown Odyssey main phase 2.924 m/s in 13 burns (SOURCED, stageo_config.toml [stageo.flown]), an ASSUMED anchor; burns are bounded by the one-per-24-h cadence")
            if arm == "OL":
                s, _ = PLB.ol_schedule(sw, tab, Pa, q1, Pa["dp_tie_eps"]); sched_out[g["id"]] = [[int(k), float(v)] for k, v in sorted(s.items())]
                e.update(burns_schedule=len(s), dv_schedule_mps=float(sum(abs(v) for v in s.values())))
            e["c4_window_mps"] = [None, 2.0 * e["dv_model_world_mps"]]
            arms[f"{arm}@{tok}"] = e
    for arm, src in (("FO", "F-O"), ("N2I", "N2-I")):
        arms[f"{arm}@ODY"] = dict(arm=arm, geometry="ODY", mode="baseline", days={"FO": 84.9, "N2I": 70.9}[arm], passes={"FO": 267, "N2I": None}[arm], dv_mps={"FO": 3.86, "N2I": 4.21}[arm], burns=None,
                                  basis="SOURCED: Stage O dev arm means (plan section 7 and R1.1: F-O 267 passes, 84.9 d, 3.86 m/s; N2-I 70.9 d, 4.21 m/s); passes of N2-I and burn counts not in the plan text (null)")
    d = lambda a, g: arms[f"{a}@{g}"]["dv_model_world_mps"]
    gL1 = d("L1", "ODY") - d("L1", "G1"); gL3 = d("L3", "ODY") - d("L3", "G1")
    attrib = dict(basis="DERIVED: model-world dV differences (planner_lib.simulate); noise-driven parts (C, N) have no analytic prediction",
                  P_ff=d("L1", "ODY") - d("L2", "ODY"), P_t=d("L2", "ODY") - d("L3", "ODY"), N=d("L0", "ODY") - d("L1", "ODY"), geometry_effect_L1=gL1, geometry_effect_L3=gL3,
                  geometry_main_effect=0.5 * (gL1 + gL3), geometry_interaction=gL3 - gL1, P_g_total=gL3, floor_ens_ODY=rows["ODY"]["floor_ens_mps"], floor_ens_G1=rows["G1"]["floor_ens_mps"])
    sens = {}
    for tok in ("ODY", "G1"):
        for pen in (0.0, 0.01, 0.03):
            for arm in ("L1", "L3", "OL", "OLR"):
                r_ = PLB.simulate(_sw(cfg, arm), tables[geoms[tok]["id"]], dict(P, dp_burn_penalty_mps=pen), B, q1, fast=True)
                sens[f"{arm}@{tok}@pen{pen}"] = dict(dv_mps=r_["dv"], burns=r_["burns"], model_violations_above_cap=r_["viol_hi"])
    feas = {}
    for tok in ("ODY", "G1"):
        r = rows[tok]; feas[tok] = dict(days_model=r["days"], d_ref_days=cfg["planner"]["d_ref_days"], days_at_qcap=None, frot_mean=r["frot_mean"], vrel_start_mps=r["vrel_start_mps"], vrel_end_mps=r["vrel_end_mps"],
                                       note="time law re-checked with f_rot(i) and V_rel: geometry_rank.schedule applies vrel_parts (V_rel, f_rot) per pass at the candidate's inclination (R4.4); campaign time is the model's, not the B-ODY flown time")
    gate_rows = {tok: dict(floor_det_mps=rows[tok]["floor_det_mps"], floor_ens_mps=rows[tok]["floor_ens_mps"], regime=rows[tok]["regime"], gid=geoms[tok]["id"]) for tok in geoms}
    body = dict(schema="planner_predictions/1", labels="SOURCED / DERIVED / ASSUMED per entry (basis field); NOT frozen until Freeze B", arms=arms, ol_schedules=sched_out, attribution_expectations=attrib,
                feasibility=feas, geometry_floors=gate_rows,
                burn_penalty_sensitivity=dict(basis="ASSUMED illustrative penalties (m/s per burn) on the noise-free replay with the vectorized DP (oracle_lib.dp_schedule); the config default is 0.0 (none)", rows=sens),
                constants=dict(q_cap_mean_w_cm2=dict(value=cfg["planner"]["q_cap_mean_w_cm2"], label="DERIVED (R1.2, sigma_p 0.10 ASSUMED U14)"), q_never_exceed_qs_w_cm2=dict(value=cfg["planner"]["q_never_exceed_qs_w_cm2"], label="SOURCED value (Hanna), ASSUMED use as never-exceed (U12)"),
                               d_ref_days=dict(value=cfg["planner"]["d_ref_days"], label="DERIVED (R1.2)"), q_first_w_cm2=dict(value=cfg["stageo"]["bands"]["odyssey"]["q_target_w_cm2"], label="DERIVED geometric mean of the Odyssey band (Stage F rule)"),
                               burn_response_ratio=dict(value=cfg["planner"]["burn_response_ratio"], label="DERIVED (halflever_twin_results.md)")))
    return body


def write_predictions(cfg, cfg3, ctx, geoms, tables, outdir):
    body = predictions(cfg, cfg3, ctx, geoms, tables)
    doc = dict(body); doc["body_sha256"] = PLB.sha256_bytes(PLB.canon(body).encode())
    p = pathlib.Path(outdir) / PRED_NAME; p.write_text(json.dumps(doc, indent=1, default=float)); h = PLB.sha256_file(p)
    (pathlib.Path(outdir) / (PRED_NAME + ".sha256")).write_text(h + "  " + PRED_NAME + "\n"); print(f"wrote {p} sha256 {h}"); return doc, h


def selftest():
    cfg = PLB.load_cfg(); GR, G3A = _imports(); cfg3 = G3A.load3(); ctx = make_ctx_cal(cfg3, cfg); geoms = geometries(cfg, cfg3)
    assert geoms["G1"]["id"] == "G5863" and geoms["G2"]["id"] == "G2348" and geoms["Gworst"]["id"] == "G3503" and geoms["ODY"]["cand"]["inc"] == 93.522
    g = geoms["ODY"]; field, _ = G3A.meanfield_for(cfg3, ctx, None, g["g"]); tab, m = ground_table(ctx, g["cand"], field, cfg)
    n = tab["n"]
    assert all(len(tab[k]) == n for k in ("path", "slow", "fast", "w", "H", "a", "e", "t")) and len(tab["d_slow"]) == n - 1 and abs(tab["path"][0]) < 1e-12 and tab["w"][0] == 0.0 and all(x > 0 for x in tab["w"][1:])
    assert np.allclose(np.asarray(tab["path"]) - np.asarray(tab["slow"]), tab["fast"]) and 190 <= n <= 260 and 0.12 < tab["q_req"] < 0.20 and 0.40 < tab["W"] < 1.20
    # lever weight of the table equals planner_lib.burn_weight with the table's own apsides
    mu = ctx.mu
    for k in (1, 50, n - 1):
        assert abs(tab["w"][k] - PLB.burn_weight(tab["H"][k], ctx.ratio, mu, tab["a"][k - 1], tab["e"][k - 1])) < 1e-9 * tab["w"][k]
    # replay: the OL schedule keeps the model-world truth inside the tube, and the dV of the DP on the slow path is close to the ranking floor's order of magnitude
    P = PLB.flat_P(cfg); B = cfg["planner"]["body"]; q1 = math.log(cfg["stageo"]["bands"]["odyssey"]["q_target_w_cm2"])
    for arm in ("L1", "L2", "L3", "OLR"):
        s = PLB.simulate(PLB.arm_switches(cfg, arm), tab, P, B, q1); assert s["viol_hi"] == 0, (arm, s["viol_hi"], s["dv"])
    ol = PLB.simulate(PLB.arm_switches(cfg, "OL"), tab, P, B, q1)           # blind: the plan keeps the SLOW path inside the reserve-shrunk tube, so the fast part (5-95% width 0.257) crosses the cap on some passes
    assert ol["viol_hi"] <= 0.25 * n and ol["viol_lo"] <= 0.25 * n, (ol["viol_hi"], ol["viol_lo"])
    d3 = PLB.simulate(PLB.arm_switches(cfg, "L3"), tab, P, B, q1)["dv"]; d2 = PLB.simulate(PLB.arm_switches(cfg, "L2"), tab, P, B, q1)["dv"]
    assert d3 <= d2 * 1.05 + 1e-9, (d3, d2)                                             # prepayment does not cost more than lazy timing (5% slack for the horizon)
    # JSON round trip and body hash
    body = dict(a=1.5, b=[1, 2]); assert PLB.sha256_bytes(PLB.canon(body).encode()) == PLB.sha256_bytes(PLB.canon(json.loads(json.dumps(body))).encode())
    print(f"planner_predictions selftest ok: ODY table n={n}, q_req={tab['q_req']:.4f}, W={tab['W']:.4f}; lever weights equal planner_lib.burn_weight; L1/L2/L3/OL-R stay inside the tube in the model world; "
          f"dV(L3)={d3:.3f} <= dV(L2)={d2:.3f} m/s; OL blind model-world passes above the cap {ol['viol_hi']}/{n}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dry-run", action="store_true"); ap.add_argument("--tables", action="store_true"); ap.add_argument("--predict", action="store_true")
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--out-dir"); a = ap.parse_args()
    if a.dry_run:
        ok, _ = PLB.dry_run(PLB.load_cfg()); sys.exit(0 if ok else 1)
    if a.selftest:
        return selftest()
    if not (a.tables or a.predict):
        ap.error("--dry-run | --tables | --predict | --selftest")
    cfg = PLB.load_cfg(); GR, G3A = _imports(); cfg3 = G3A.load3(); ctx = make_ctx_cal(cfg3, cfg); geoms = geometries(cfg, cfg3)
    outdir = pathlib.Path(a.out_dir) if a.out_dir else PLB.task_dir()
    tables = build_tables(cfg, cfg3, ctx, geoms)
    if a.tables:
        write_tables(cfg, cfg3, ctx, geoms, outdir)
    if a.predict:
        write_predictions(cfg, cfg3, ctx, geoms, tables, outdir)


if __name__ == "__main__":
    main()
