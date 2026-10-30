"""Stage G3 analysis (WP-G3, gate V-G): are the free-drift predictions of geometry_rank.py right at other geometries? Plan: tasks/2026-10-15_periapsis-region-border-targeting/
analytic_framework_and_transfer_plan.md, R4.3 (design, criteria, outcomes), R4.4 and R5.1 (U15: heat rate on the atmosphere-relative speed), R4.1 (exact lever, burn removal). Python only, no Julia.
Imports geometry_rank (model, tube DP, mean field loader; never edited) and stageo_analysis (inertial Q_s form qs_of, Odyssey tiers).
  selftest                                  closed forms, the U15 rotation identity, criteria and outcomes on synthetic cases, burn removal against geometry_rank.dev_natural, epoch feeding identity
  predict   [--meanfield DIR] [--out F]     pre-flight prediction table per geometry (deterministic and jitter-ensemble floors, drift range R, days, passes, regime): no flights needed
  gate      --results-root R [--meanfield DIR]   criteria (a) to (c) per flight with the model fed the flown pass epochs, geometry verdicts, the R4.3 outcome
  integrity --results-root R                route, failures, epoch, NaN, burns, stop reasons, scoring table (Q_s on V_rel and, for comparison, inertial)
How each U15 piece is implemented (R5.1: relative speed in BOTH the keeper's decision quantity and the scoring):
  decision: stageg3_common.jl StageGPolicy (V_rel at the osculating periapsis from periapsis_vrel, planet rotation vector); scoring: qs_rel() below on the per-pass vrel_mps column (the exact relative speed at the
  periapsis event, the one q_w_cm2 uses); the inertial qs_of of Stage O is reported next to it. The rotation correction equals omega r_p cos(i) along track at periapsis (checked in both selftests).
Usage: python3 -I stageg3_analysis.py selftest | predict | gate --results-root R | integrity --results-root R"""
import argparse, contextlib, glob, io, json, math, pathlib, sys, tempfile, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import geometry_rank as GR
import stageo_analysis as SOA

CFG3_PATH = HERE / "stageg3_config.toml"


def table_path(t):
    """Location of the geometry table: env SPACEAGORA_BORDER_TASK_DIR, else the task directory found by walking up from this file."""
    import os
    d = os.environ.get("SPACEAGORA_BORDER_TASK_DIR"); d = pathlib.Path(d) if d else SOA.task_dir()
    if d is None: raise FileNotFoundError("task directory not found (set SPACEAGORA_BORDER_TASK_DIR)")
    return d / t["path"]


def load3(path=CFG3_PATH):
    """Config with the geometry list injected from the final table (geometry_gate_table.json). The table's sha256 is verified against [stageg3.table].sha256 and a mismatch is an ERROR.
    cfg["stageg3"]["geometry"] = default-set geometries (plus G3 only when include_g3), each with id, role, inc_deg, omega_deg, raan_deg, epoch_off_d, epoch_utc, field_id (= id)."""
    import hashlib
    cfg = tomllib.loads(pathlib.Path(path).read_text()); t = cfg["stageg3"]["table"]; fp = table_path(t); raw = fp.read_bytes()
    h = hashlib.sha256(raw).hexdigest()
    if h != t["sha256"]: raise ValueError(f"geometry table sha256 mismatch: {fp} is {h}, config has {t['sha256']}")
    T = json.loads(raw); G = []
    for g in T["geometries"]:
        if not (g["in_default_set"] or (t["include_g3"] and g["role"] == "G3")): continue
        G.append(dict(id=g["id"], role=g["role"], inc_deg=g["inclination_deg"], omega_deg=g["omega0_deg"], raan_deg=g["raan0_deg"], epoch_off_d=g["epoch_offset_days"], epoch_utc=g["epoch_utc"], field_id=g["id"]))
    cfg["stageg3"]["geometry"] = G; cfg["stageg3"]["table_start"] = T["start"]; cfg["stageg3"]["table_hash"] = h; return cfg


# ---- heat rate on the relative speed (U15) ------------------------------------------------------------------------------------------------------------

def qs_rel(rho_peak, vrel):
    """Scoring heat rate on the atmosphere-relative speed, W/cm2: 0.5 rho_peak V_rel^3 / 1e4 (SI inputs); vrel is the exact speed at the periapsis event (passes.csv vrel_mps). Inertial form: SOA.qs_of."""
    return 0.5 * np.asarray(rho_peak, float) * np.asarray(vrel, float) ** 3 / 1e4


def rotation_identity(inc_deg, raan_deg, argp_deg, rp, ra, mu, w):
    """(along-track rotation correction, omega r_p cos i, V_rel exact, V_rel closed form) for the periapsis of an orbit with the given elements, planet rotation about +z at rate w. Independent of the Julia code."""
    i, O, om = map(math.radians, (inc_deg, raan_deg, argp_deg)); a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
    Rz = lambda x: np.array([[math.cos(x), -math.sin(x), 0], [math.sin(x), math.cos(x), 0], [0, 0, 1]]); Rx = np.array([[1, 0, 0], [0, math.cos(i), -math.sin(i)], [0, math.sin(i), math.cos(i)]])
    R = Rz(O) @ Rx @ Rz(om); vp = math.sqrt(mu * (2 / rp - 1 / a))
    rvec = R @ np.array([rp, 0, 0]); vvec = R @ np.array([0, vp, 0]); omega = np.array([0, 0, w]); rot = np.cross(omega, rvec)
    along = float(rot @ (vvec / np.linalg.norm(vvec))); vrel = float(np.linalg.norm(vvec - rot)); sphi = math.sin(i) * math.sin(om)
    closed = math.sqrt(vp ** 2 - 2 * w * rp * vp * math.cos(i) + (w * rp) ** 2 * (1 - sphi ** 2))
    return along, w * rp * math.cos(i), vrel, closed


# ---- natural path of a flown run (burn removal with the exact lever, R4.1) --------------------------------------------------------------------------

def flown_natural(ctx, run, ratio):
    """Natural (burn-removed) ln q_mean path of one flight and the burn-cost weights (m/s per ln unit), as geometry_rank.dev_natural but tolerant of the empty decisions.csv of a zero-burn flight.
    Burn k acts between passes k-1 and k; the burn's effect on ln q is -dv / w with w = H / (ratio x lever) (the exact vis-viva lever for ratio 1, R4.1)."""
    run = pathlib.Path(run); P = pd.read_csv(run / "passes.csv").dropna(subset=["ra_m"]).reset_index(drop=True)
    try: D = pd.read_csv(run / "decisions.csv")
    except (FileNotFoundError, pd.errors.EmptyDataError): D = pd.DataFrame(columns=["target_pass", "delta_v_mps"])
    n = len(P); lnq = np.log(P.q_w_cm2.values); H = P.h_local_m.values; rp = P.rp_m.values; ra = P.ra_m.values; burn = np.zeros(n); w = np.zeros(n)
    for i in range(1, n):
        a = 0.5 * (ra[i - 1] + rp[i - 1]); e = (ra[i - 1] - rp[i - 1]) / (ra[i - 1] + rp[i - 1]); w[i] = H[i] / (ratio * GR.lever(ctx, a, e))
    for r in D.itertuples():
        i = int(r.target_pass) - 1
        if r.delta_v_mps != 0 and 1 <= i < n: burn[i] = -r.delta_v_mps / w[i]
    return P, lnq - np.cumsum(burn), w


# ---- the model fed the flown pass epochs (R4.3) ------------------------------------------------------------------------------------------------------

def _trunc(S, n):
    out = {}
    for k, v in S.items(): out[k] = v[:n] if isinstance(v, np.ndarray) else v
    out["n"] = n; return out


def model_path(ctx, cand, field=None, t_flown=None):
    """Burn-free ln q path of the model for `cand`. With t_flown (flown periapsis times, s from the start) the pass epochs are replaced by the flown ones for the first n = min(model, flown) passes, which removes the timing
    error the ranking's jitter ensemble covered; everything else (apsides, scale heights, mean field by pass index) stays the model's. Returns dict(path, slow, w, S, q_req, n)."""
    cfg = ctx.cfg; b = cfg["band"]; H_of_k = None if field is None else field["H"]
    q_req = GR.solve_q_req(ctx, cand, b["d_ref_days"], H_of_k); S = GR.schedule(ctx, cand, q_req, H_of_k)
    if t_flown is not None:
        n = min(S["n"], len(t_flown)); S = _trunc(S, n); S["t"] = np.asarray(t_flown, float)[:n]
    d = GR.channels(ctx, cand, S, field); path = GR.path_of(d); slow = GR.split_slow_fast(path, b["fast_window_passes"])[0]
    n = S["n"]; w = np.zeros(n); w[1:] = S["H"][1:] / (ctx.ratio * np.array([GR.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, n)]))
    return dict(path=path, slow=slow, w=w, S=S, q_req=q_req, n=n)


def cand_of(g):
    """geometry_rank candidate dict from a [[stageg3.geometry]] row (the id is a label only; the elements are what is used)."""
    return dict(id=str(g["id"]), inc=float(g["inc_deg"]), epoch_off=float(g["epoch_off_d"]), omega=float(g["omega_deg"]), raan=float(g["raan_deg"]))


def pre_flight_row(ctx, g, field):
    """Pre-flight prediction (ranking evaluate, deterministic and jitter ensemble) of one geometry."""
    r, n, S, d = GR.evaluate(ctx, cand_of(g), field, ctx.s_t_ody); fl = GR.jitter_floors(ctx, cand_of(g), field)
    return dict(id=g["id"], role=g["role"], inc_deg=g["inc_deg"], omega0_deg=g["omega_deg"], raan0_deg=g["raan_deg"], epoch_off_d=g["epoch_off_d"], n_passes=r["n_passes"], days=r["days"], q_req=r["q_req"],
                W=r["W"], fast_w=r["fast_w"], W_eff=r["W_eff"], floor_det_mps=r["floor_mps"], floor_ens_mps=float(np.mean(fl)), floor_ens_sd=float(np.std(fl)), floor_driftfree_mps=r["floor_driftfree_mps"],
                R_ln=r["R"], N_ln=r["N"], regime=r["regime"], field_mode=r["field_mode"], frot_mean=r["frot_mean"], vrel_start_mps=r["vrel_start_mps"], vrel_end_mps=r["vrel_end_mps"])


# ---- criteria (R4.3, U19) ----------------------------------------------------------------------------------------------------------------------------

def criteria(pred_slow, fl_slow, F_pred, F_fl, F_loss, W_eff, cr, regime_tol):
    """R4.3 pass criteria on one flight. (a) Pearson correlation of the smoothed paths >= min_correlation; (b) drift range within range_factor of the flown, or within +-range_abs_ln when both are below W_eff;
    (c) the floor on the flown path within +-max(floor_abs, floor_rel x predicted) of the predicted AND on the same side of the neutral / interferes boundary (F > F_loss + tol). Returns a dict of values and booleans."""
    a = np.asarray(pred_slow, float); b = np.asarray(fl_slow, float)
    corr = float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else float("nan")
    Rp, Rf = float(np.ptp(a)), float(np.ptp(b)); lo, hi = min(Rp, Rf), max(Rp, Rf)
    ratio = hi / lo if lo > 0 else float("inf")
    b_factor = bool(ratio <= cr["range_factor"]); b_abs = bool(abs(Rp - Rf) <= cr["range_abs_ln"] and hi < W_eff)
    tol_f = max(cr["floor_abs_mps"], cr["floor_rel"] * abs(F_pred)); c_floor = bool(abs(F_fl - F_pred) <= tol_f)
    side_p, side_f = bool(F_pred > F_loss + regime_tol), bool(F_fl > F_loss + regime_tol)
    out = dict(corr=corr, a_pass=bool(np.isfinite(corr) and corr >= cr["min_correlation"]), R_pred=Rp, R_flown=Rf, R_ratio=ratio, b_factor=b_factor, b_abs=b_abs, b_pass=b_factor or b_abs,
               F_pred=float(F_pred), F_flown=float(F_fl), F_tol=tol_f, c_floor=c_floor, regime_pred="interferes" if side_p else "neutral_or_helps", regime_flown="interferes" if side_f else "neutral_or_helps",
               c_side=side_p == side_f)
    out["c_pass"] = bool(c_floor and out["c_side"]); out["pass"] = bool(out["a_pass"] and out["b_pass"] and out["c_pass"]); return out


def verdict(per_geom, rule="all"):
    """Geometry verdict from its flights' criteria dicts (None = insufficient, counts as not passed). rule 'all': every flight passes."""
    if not per_geom: return None
    ok = [(c is not None and c["pass"]) for c in per_geom]
    return bool(all(ok)) if rule == "all" else bool(sum(ok) > len(ok) / 2)


def outcome(g1, g2):
    """R4.3 outcomes, with the gap in R4.3 closed by a USER DECISION dated 2026-10-26: G1 passes and G2 fails -> PROCEED WITH G1 ONLY, G2 reported as diagnostic.
    G1 and G2 pass: geometry arms proceed. G1 fails, G2 passes: G2 becomes the headline geometry. Both fail: geometry arms dropped. Gworst informs D4 only and never enters."""
    if g1 is None or g2 is None: return "INCOMPLETE"
    if g1 and g2: return "PROCEED_G1_AND_G2"
    if g1 and not g2: return "PROCEED_G1_ONLY (user decision 2026-10-26; G2 diagnostic)"
    if (not g1) and g2: return "G2_HEADLINE"
    return "GEOMETRY_ARMS_DROPPED"


def score_flight(ctx, cfg3, g, field, run, ratio):
    """Criteria for one flown run directory at geometry row g, with the model fed the flown epochs. Returns dict (criteria + context) or an INSUFFICIENT record."""
    cr = cfg3["stageg3"]["criteria"]; win = ctx.cfg["band"]["fast_window_passes"]
    P, nat, wfl = flown_natural(ctx, run, ratio); rec = dict(run=str(pathlib.Path(run).name), n_flown=len(P))
    if len(P) < cr["min_passes_scored"]: rec.update(insufficient=True, reason="fewer than min_passes_scored flown passes"); return rec
    m = model_path(ctx, cand_of(g), field, P.time_s.values); n = min(m["n"], len(P)); rec["n_scored"] = n
    if n < cr["min_passes_scored"]: rec.update(insufficient=True, reason="model schedule shorter than min_passes_scored"); return rec
    nat = nat[:n]; slow_f, fast_f = GR.split_slow_fast(nat, win); pred = m["slow"][:n]
    # prediction row for W, the model fast width and the drift-free floor (the ranking's own evaluate, unfed); W_eff of each side uses its own fast width, as the pre-check compares the model with the oracle
    row = GR.evaluate(ctx, cand_of(g), field, ctx.s_t_ody)[0]; levels = ctx.cfg["floor"]["n_levels"]
    W_eff_pred = row["W_eff"]; F_pred = GR.tube_floor(pred, m["w"][:n], W_eff_pred, levels)
    W_eff_fl = row["W"] - GR.width_5_95(fast_f); F_fl = GR.tube_floor(slow_f, wfl[:n], W_eff_fl, levels) if W_eff_fl > 0 else float("inf")
    c = criteria(pred, slow_f, F_pred, F_fl, row["floor_driftfree_mps"], W_eff_pred, cr, ctx.cfg["floor"]["regime_tol_mps"])
    rec.update(c); rec.update(insufficient=False, W=row["W"], W_eff_pred=W_eff_pred, W_eff_flown=W_eff_fl, fast_w_flown=GR.width_5_95(fast_f), F_driftfree=row["floor_driftfree_mps"], model_days=row["days"], model_passes=row["n_passes"]); return rec


# ---- results layout ---------------------------------------------------------------------------------------------------------------------------------

def sd(root): return pathlib.Path(root) / "stageg3"
def adir(root): d = sd(root) / "analysis"; d.mkdir(parents=True, exist_ok=True); return d
def dump(root, name, res): (adir(root) / name).write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))


def summaries(root, phase="gate"):
    fs = sorted(glob.glob(str(sd(root) / f"summary_{phase}_mars_*.csv"))); return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def meanfield_for(cfg3, ctx, meanfield, g):
    """Mean field of geometry g (key field_id) from --meanfield, else the config directories (relative to the task directory); None (and a note) when not found."""
    td = SOA.task_dir(); cr = cfg3["stageg3"]["criteria"]; dirs = [meanfield] if meanfield else [td / cr["meanfield_dir"], td / cr["meanfield_dir_fallback"]] if td else []
    for d in dirs:
        try: F = GR.load_meanfield(d, ctx.cfg)
        except Exception: continue
        if g["field_id"] in F: return F[g["field_id"]], str(d)
    return None, None


def make_ctx(cfg3):
    cfg = GR.load_cfg(); ratio = float(cfg3["stageg3"]["burn_response_ratio"])
    ctx = GR.Ctx(cfg, ratio); GR.anchor_ctx(ctx); return ctx


def cmd_predict(cfg3, meanfield, out):
    ctx = make_ctx(cfg3); rows = []
    for g in cfg3["stageg3"]["geometry"]:
        field, src = meanfield_for(cfg3, ctx, meanfield, g); r = pre_flight_row(ctx, g, field); r["meanfield_source"] = src; rows.append(r); print(json.dumps(r, default=float))
    if out: pathlib.Path(out).write_text(json.dumps(rows, indent=1, default=float))
    return rows


def cmd_gate(cfg3, root, meanfield):
    ctx = make_ctx(cfg3); S = summaries(root, "gate"); G3 = cfg3["stageg3"]; ok = S[S.status == "ok"] if len(S) else S; per = {}; recs = []
    for g in G3["geometry"]:
        field, src = meanfield_for(cfg3, ctx, meanfield, g); rows = ok[ok.geometry == g["id"]] if len(ok) else ok; per[g["id"]] = []
        for r in rows.itertuples():
            rec = score_flight(ctx, cfg3, g, field, sd(root) / "runs" / r.run_id, float(G3["burn_response_ratio"]))
            rec.update(geometry=g["id"], role=g["role"], keeper=r.keeper, tag=r.tag, stop_reason=r.stop_reason, meanfield=src, field_mode="meanfield" if field is not None else "NO_FIELD (prediction without the mean-field term)")
            recs.append(rec); per[g["id"]].append(None if rec.get("insufficient") else rec)
    ver = {g["id"]: verdict(per[g["id"]], G3["criteria"]["flights_rule"]) for g in G3["geometry"]}; role = {g["role"]: g["id"] for g in G3["geometry"]}
    res = dict(flights=recs, verdict={gid: dict(passed=v, role=next(g["role"] for g in G3["geometry"] if g["id"] == gid)) for gid, v in ver.items()},
               outcome=outcome(ver[role["G1"]], ver[role["G2"]]), gworst=dict(id=role["Gworst"], passed=ver[role["Gworst"]], note="informs D4 only; cannot block"),
               flights_rule=G3["criteria"]["flights_rule"], criteria=G3["criteria"], burn_response_ratio=G3["burn_response_ratio"])
    dump(root, "gate_verdict.json", res); return res


def cmd_integrity(cfg3, root):
    S = summaries(root, "gate"); G3 = cfg3["stageg3"]; out = dict(n_rows=len(S), expected=G3["plan_counts"]["gate"], problems=[]); P = out["problems"]
    if not len(S): P.append("no gate summary rows"); dump(root, "integrity.json", out); return out
    out["rows_ok"] = int((S.status == "ok").sum()); (S.status != "ok").any() and P.append("rows not ok: " + ", ".join(S[S.status != "ok"].run_id))
    len(S) != G3["plan_counts"]["gate"] and P.append(f"row count {len(S)} != {G3['plan_counts']['gate']}")
    for f in sorted(glob.glob(str(sd(root) / "mc_gate_mars_*.toml"))):
        t = tomllib.loads(pathlib.Path(f).read_text()); t.get("mc_route") != "process" and P.append(f"route {t.get('mc_route')} in {pathlib.Path(f).name}, not process"); t.get("mc_failed", 0) and P.append(f"mc_failed={t['mc_failed']}")
    td = {g["id"]: g for g in G3["geometry"]}; tab = []
    for r in S[S.status == "ok"].itertuples():
        g = td[r.geometry]; ep = pd.Timestamp(GR_epoch(G3, g)); core = f"{ep.year}-{ep.month}-{ep.day} {ep.hour}:{ep.minute}:{ep.second}"
        not str(r.gram_core_epoch).startswith(core) and P.append(f"{r.run_id}: GRAM core epoch {r.gram_core_epoch} != {core}")
        r.keeper == "zero_burn" and r.burns != 0 and P.append(f"{r.run_id}: zero-burn flight has {r.burns} burns")
        r.keeper == "fixed" and r.stop_reason not in ("end_radius", "horizon") and P.append(f"{r.run_id}: keeper stop_reason {r.stop_reason}")
        r.keeper == "zero_burn" and r.stop_reason in ("controller_stopped",) and P.append(f"{r.run_id}: zero-burn stop_reason {r.stop_reason}")
        D = pd.read_csv(sd(root) / "runs" / r.run_id / "passes.csv"); d = D.dropna(subset=["ra_m"])
        d[["rho_peak_kg_m3", "q_w_cm2", "vrel_mps", "rp_m"]].isna().any().any() and P.append(f"{r.run_id}: NaN in rho_peak / q / vrel / rp")
        qr = qs_rel(d.rho_peak_kg_m3, d.vrel_mps); qi = SOA.qs_of(d.rho_peak_kg_m3, d.rp_m, d.ra_m)
        tab.append(dict(run_id=r.run_id, geometry=r.geometry, role=r.role, keeper=r.keeper, passes=len(D), days=r.t_end_days, dv_mps=r.total_dv_mps, burns=r.burns, stop_reason=r.stop_reason, safety_reason=r.safety_reason,
                        qs_rel_max=float(qr.max()), qs_inertial_max=float(qi.max()), qs_rel_over_032=float((qr > 0.32).mean()), qs_inertial_over_032=float((qi > 0.32).mean()), ratio_rel_over_inertial_median=float(np.median(qr / qi)),
                        any_rel_over_065=bool((qr > 0.65).any()), any_rel_over_054=bool((qr > 0.54).any())))
    out["scoring"] = tab; dump(root, "integrity.json", out); return out


def GR_epoch(G3, g):
    return pd.Timestamp(G3["base_epoch_utc"]) + pd.Timedelta(days=float(g["epoch_off_d"]))


# ---- selftest ----------------------------------------------------------------------------------------------------------------------------------------

def _fake_run(root, ctx, cand, S, burns, noise, tamper=None):
    """Write a synthetic flight (passes.csv, decisions.csv) from a model schedule S: q follows `noise` (a ln q path) plus the cumulative effect of `burns` {pass index: dv} with the exact lever weights."""
    n = S["n"]; w = np.zeros(n); w[1:] = S["H"][1:] / (ctx.ratio * np.array([GR.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, n)]))
    eff = np.zeros(n)
    for i, dv in burns.items(): eff[i] = -dv / w[i]
    lnq = np.log(0.2) + noise + np.cumsum(eff)
    P = pd.DataFrame(dict(pass_=np.arange(1, n + 1), time_s=S["t"], q_w_cm2=np.exp(lnq), h_local_m=S["H"], rp_m=S["rp"], ra_m=S["ra"], rho_peak_kg_m3=1e-8, vrel_mps=S["V"])).rename(columns={"pass_": "pass"})
    D = pd.DataFrame(dict(target_pass=[i + 1 for i in burns], delta_v_mps=list(burns.values())))
    run = pathlib.Path(root); run.mkdir(parents=True, exist_ok=True); P.to_csv(run / "passes.csv", index=False)
    (D.to_csv(run / "decisions.csv", index=False) if len(D) else (run / "decisions.csv").write_text("")); return w


def selftest():
    cfg3 = load3(); G3 = cfg3["stageg3"]; cr = G3["criteria"]; gs = G3["geometry"]
    # config
    assert G3["burn_response_ratio"] == 1.0 and [g["role"] for g in gs] == ["G1", "G2", "Gworst"] and G3["plan_counts"]["gate"] == 10 and G3["plan_counts"]["keeper"] + G3["plan_counts"]["zero_burn"] == 10
    assert (cr["min_correlation"], cr["range_factor"], cr["range_abs_ln"], cr["floor_abs_mps"], cr["floor_rel"], cr["smooth_passes"]) == (0.85, 1.5, 0.15, 0.3, 0.30, 15)
    assert G3["safety"]["q_stop_w_cm2"] == 0.65 and G3["safety"]["n_hloc"] == 3.0
    n_keeper = len(gs) * G3["keeper_seeds_per_geometry"]; n_zb = sum(g["role"] in G3["zero_burn_roles"] for g in gs) * G3["zero_burn_seeds"]; assert (n_keeper, n_zb) == (6, 4)
    assert [g["id"] for g in gs] == ["G5863", "G2348", "G3503"] and [g["epoch_utc"] for g in gs] == ["2002-02-04T19:00:32", "2001-11-06T19:00:32", "2001-11-06T19:00:32"] and [g["inc_deg"] for g in gs] == [116.6, 63.4, 88.6]
    for g in gs: assert pd.Timestamp(g["epoch_utc"]) == pd.Timestamp(G3["base_epoch_utc"]) + pd.Timedelta(days=g["epoch_off_d"])
    assert cfg3["stageg3"]["table_start"]["ra0_m"] == 28533000.0 and cfg3["stageg3"]["table_start"]["ra_target_m"] == 7408000.0 and cfg3["stageg3"]["table_start"]["true_anomaly_deg"] == 180.0
    import tempfile as _tf
    with _tf.TemporaryDirectory() as _d:
        pth = pathlib.Path(_d) / "bad.toml"; pth.write_text(CFG3_PATH.read_text().replace(G3["table"]["sha256"], "0" * 64))
        try: load3(pth); raise SystemExit("hash mismatch was not rejected")
        except ValueError as e: assert "mismatch" in str(e)
    with _tf.TemporaryDirectory() as _d:
        pth = pathlib.Path(_d) / "g3.toml"; pth.write_text(CFG3_PATH.read_text().replace("include_g3 = false", "include_g3 = true")); assert [g["role"] for g in load3(pth)["stageg3"]["geometry"]] == ["G1", "G2", "Gworst", "G3"]
    cfgr = GR.load_cfg(); assert cfgr["floor"]["burn_response_ratio"] == G3["burn_response_ratio"], "geometry_rank_config burn_response_ratio differs from the gate config"
    assert cfgr["band"]["fast_window_passes"] == cr["smooth_passes"]
    # U15: rotation identity along track and V_rel closed form (independent numpy implementation), R4.4 table
    mu, w = G3["rotation"]["mu_m3s2"], G3["rotation"]["w_rot_rad_s"]; worst = worst_v = 0.0; cnt = 0
    for inc in (25.2, 63.4, 93.5, 93.522, 116.6, 150.0):
        for raan in (0.0, 90.0, 270.0):
            for argp in (0.0, 45.0, 90.0, 200.0):
                al, wr, vr, cl = rotation_identity(inc, raan, argp, 3.4838e6, 28533e3, mu, w); worst = max(worst, abs(al - wr)); worst_v = max(worst_v, abs(vr - cl)); cnt += 1
    assert worst < 1e-6 and worst_v < 1e-6, (worst, worst_v)
    rp0, ra0 = 3.4838e6, 28533e3; vp0 = math.sqrt(mu * (2 / rp0 - 2 / (ra0 + rp0)))
    for inc, want in ((63.4, 0.976), (93.5, 1.003), (116.6, 1.024)): assert abs((vp0 - w * rp0 * math.cos(math.radians(inc))) / vp0 - want) < 0.002
    assert abs(float(qs_rel(2 * 0.65e4 / 4700 ** 3, 4700.0)) - 0.65) < 1e-12 and abs(float(SOA.qs_of(1e-8, 3.5e6, 2.8e7)) - 0.5e-8 * (mu * (2 / 3.5e6 - 2 / 3.15e7)) ** 1.5 / 1e4) < 1e-12
    # criteria on synthetic paths
    t = np.linspace(0, 1, 200); a = 0.3 * np.sin(2 * np.pi * 1.5 * t); good = criteria(a, 1.1 * a, 0.5, 0.6, 0.1, 0.46, cr, 0.05); assert good["pass"] and good["corr"] > 0.999
    assert not criteria(a, -a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["a_pass"]                                    # anticorrelated
    assert not criteria(a, 3 * a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["b_pass"]                                 # range ratio 3 > 1.5, difference 0.9 > 0.15
    assert criteria(0.05 * a, 0.08 * a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["b_pass"]                            # factor 1.6 fails but both below W_eff and within 0.15 ln
    assert criteria(a, a, 0.5, 0.79, 0.1, 0.46, cr, 0.05)["c_floor"] and not criteria(a, a, 0.5, 0.81, 0.1, 0.46, cr, 0.05)["c_floor"]   # larger of 0.3 and 30% of 0.5 is 0.3
    assert criteria(a, a, 2.0, 2.59, 0.1, 0.46, cr, 0.05)["c_floor"] and not criteria(a, a, 2.0, 2.61, 0.1, 0.46, cr, 0.05)["c_floor"]  # 30% of 2.0 = 0.6
    assert criteria(a, a, 0.12, 0.30, 0.1, 0.46, cr, 0.05)["c_floor"] and not criteria(a, a, 0.12, 0.30, 0.1, 0.46, cr, 0.05)["c_side"]   # within tolerance but across the neutral / interferes boundary
    assert (outcome(True, True), outcome(False, True), outcome(False, False)) == ("PROCEED_G1_AND_G2", "G2_HEADLINE", "GEOMETRY_ARMS_DROPPED") and outcome(True, False).startswith("PROCEED_G1_ONLY") and outcome(None, True) == "INCOMPLETE"
    assert verdict([good, good]) and not verdict([good, None]) and not verdict([good, dict(good, **{"pass": False})]) and verdict([good, dict(good, **{"pass": False})], "majority") is False
    # model feeding and burn removal on the Odyssey anchor (no field): feeding the model's own epochs reproduces its path; a 3% longer clock changes the tesseral path; burn removal inverts the lever weights
    ctx = GR.Ctx(cfgr, G3["burn_response_ratio"]); anc = GR.anchor_cand(cfgr)
    full = model_path(ctx, anc); n0 = min(full["n"], 60); own = model_path(ctx, anc, None, full["S"]["t"][:n0]); assert own["n"] == n0 and np.allclose(own["path"], full["path"][:n0], atol=1e-12)
    skew = model_path(ctx, anc, None, 1.03 * full["S"]["t"][:n0]); assert not np.allclose(skew["path"], own["path"], atol=1e-6)
    with tempfile.TemporaryDirectory() as td:
        rng = np.random.default_rng(7); Sf = _trunc(full["S"], n0); noise = np.concatenate([[0.0], np.cumsum(np.diff(full["path"][:n0]))]) + 0.0
        burns = {5: 0.07, 6: -0.04, 20: 0.1}; w_ = _fake_run(pathlib.Path(td) / "k", ctx, anc, Sf, burns, noise)
        P, nat, wf = flown_natural(ctx, pathlib.Path(td) / "k", G3["burn_response_ratio"]); assert np.allclose(nat - nat[0], noise - noise[0], atol=1e-9) and np.allclose(wf, w_)
        P2, nat2, w2 = GR.dev_natural(ctx, str(pathlib.Path(td) / "k"), G3["burn_response_ratio"]); assert np.allclose(nat2, nat) and np.allclose(w2, wf)          # same as the ranking's own reconstruction
        _fake_run(pathlib.Path(td) / "z", ctx, anc, Sf, {}, noise); Pz, natz, _ = flown_natural(ctx, pathlib.Path(td) / "z", 1.0); assert np.allclose(natz - natz[0], noise - noise[0], atol=1e-9)   # empty decisions.csv
        # end to end on the synthetic zero-burn flight: the flight IS the model path, so it passes (a), (b); the floor criterion compares like with like
        gA = dict(id="ODY", role="G1", inc_deg=anc["inc"], omega_deg=anc["omega"], raan_deg=anc["raan"], epoch_off_d=0.0, field_id="ODY")
        cfg_small = dict(cfg3); cfg_small["stageg3"] = dict(G3, criteria=dict(cr, min_passes_scored=30))
        rec = score_flight(ctx, cfg_small, gA, None, pathlib.Path(td) / "z", 1.0); assert not rec["insufficient"] and rec["a_pass"] and rec["b_pass"] and rec["corr"] > 0.999, rec
        _fake_run(pathlib.Path(td) / "bad", ctx, anc, Sf, {}, -noise); recb = score_flight(ctx, cfg_small, gA, None, pathlib.Path(td) / "bad", 1.0); assert not recb["a_pass"] and not recb["pass"]
        _fake_run(pathlib.Path(td) / "short", ctx, anc, _trunc(full["S"], 12), {}, noise[:12]); assert score_flight(ctx, cfg_small, gA, None, pathlib.Path(td) / "short", 1.0)["insufficient"]
    print(f"stageg3 selftest ok: config (10 flights = 6 keeper + 4 zero-burn), U15 rotation correction = omega r_p cos i along track (max error {worst:.2e} m/s over {cnt} orbits), V_rel closed form (max error {worst_v:.2e} m/s), "
          f"R4.4 table, criteria (a) to (c), outcomes, burn removal = geometry_rank.dev_natural, epoch feeding identity")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "predict", "gate", "integrity"]); ap.add_argument("--results-root"); ap.add_argument("--meanfield"); ap.add_argument("--out"); ap.add_argument("--config", default=str(CFG3_PATH))
    a = ap.parse_args(); cfg3 = load3(a.config)
    if a.cmd == "selftest": return selftest()
    if a.cmd == "predict": return cmd_predict(cfg3, a.meanfield, a.out)
    a.results_root or ap.error("--results-root required")
    return cmd_gate(cfg3, a.results_root, a.meanfield) if a.cmd == "gate" else cmd_integrity(cfg3, a.results_root)


if __name__ == "__main__": main()
