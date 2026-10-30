"""Revision 23 (WP23; analysis only: Python, no Julia, no flights): the initial geometry as a planner decision (passive corridor maintenance).

What this file does. It connects the geometry ranking (geometry_rank.py, GRAM mean-field samples in analytic_data/) to the Revision 20 DP (planner_rev20_ic.Case / select_ic) for an ARBITRARY sampled
geometry: the ground table is rebuilt from the geometry's candidate dict and its mean field exactly as planner_predictions.ground_table does (bit-identical for ODY and G5863, selftest), then planner_rev20_ic.Case
runs the exact DP on it. Everything else (ranking, robustness, drivers, rendezvous coupling, validation) is built on that one function. All other modules are imported READ-ONLY; two module globals of
planner_rev20_ic (load_tab, X_REF, P) are swapped inside a context manager (`patched`) and restored.

DP metric ("aerobraking dV", Revision 20d: no final periapsis maneuver): campaign dV through the arrival pass + the pre-burn where it is charged. S-fixed = capture periapsis given at the reference start
(a change pays w_pre |nu_pre|); S-free = the start level is free at no campaign cost (select_ic scenarios). horizon "arrival", pen 0 (the Revision 18/19/20 model-world convention).

Usage:  python3 -I planner_rev23_geometry.py selftest
        python3 -I planner_rev23_geometry.py run  [--jobs 3] [--part cut|anchors|grid]    (writes rev23_data/raw_*.json)
        python3 -I planner_rev23_geometry.py summarize                                    (reads raw_*.json and flights; writes the summary JSONs)
Env: SPACEAGORA_BORDER_TASK_DIR (task dir, planner_lib.task_dir), SPACEAGORA_PLANNER_ROOT (flight root; fallback ~/Documents/SpaceAGORA-planner-data/planner).
Labels in the output: SOURCED (file), DERIVED (computed here), MEASURED (flights), ASSUMED (named constant below)."""
import argparse, contextlib, hashlib, json, math, os, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev18_offline as O
import planner_rev20_ic as R

TD = pathlib.Path(PLB.task_dir())
OUT = TD / "rev23_data"
FLIGHT_ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))

OFFSET_LN = 0.2231            # SOURCED: Revision 17a, measured mean ln(rho/rho_base) of the +/- dust levels L2/L3
L1_SCALE = 1.5                # SOURCED: Revision 17a decision (1): L1 noise scale 1.5 (GRAM perturbation scale)
RESERVE0 = 0.257              # SOURCED: planner_config.toml fast_reserve_ln (5-95% width of the fast part at scale 1)
KAPPA_DEFAULT = 1.0           # ASSUMED: planner_config.toml kappa_default (no flight calibration); ODY/G1 MEASURED values are used only where stated
KAPPA_G1 = 0.879308559055656  # MEASURED: planner_config.toml [planner.terminal_periapsis.kappa] G5863
JIT_FRAC, JIT_N = 0.05, 12    # SOURCED: geometry_rank_config.toml phase_jitter_frac / phase_jitter_n (pass-epoch time-scale jitter)
JIT_NT = 6                    # ASSUMED: members of the jitter ensemble in the (slower) common-arrival mode; subsampled from 12 for run time
V_A, V_C = 1926.9, 2404.4     # DERIVED, Revision 22b: v_a of the final ellipse (r_p 3,504 km ASSUMED) and v_circ at r_a = 7,408 km
V0 = 477.5                    # DERIVED, Revision 22b: coplanar circular match at r_a
GOOD_FRAC = 0.5               # ASSUMED: a geometry is "good" if its S-fixed DP dV is at most this fraction of the Odyssey anchor's
TIE_MPS = 0.005              # ASSUMED: DP dV values closer than this are ties (reported; ranks use dV, then the jitter-ensemble mean, then the id)
ABL = ("tesseral", "third_body", "zonal", "field")      # channels of geometry_rank.channels removed one at a time (drivers)
G = {}                        # worker globals (filled by init, inherited by fork)


# ----------------------------------------------------------------------------------------------------- setup

def init():
    if G:
        return G
    import geometry_rank as GR, planner_predictions as PP
    cfg = PLB.load_cfg(); _, G3A = PP._imports(); cfg3 = G3A.load3(); ctx = PP.make_ctx_cal(cfg3, cfg)
    gcfg = ctx.cfg
    cut = GR.load_meanfield(TD / "analytic_data" / "meanfield", gcfg)
    G.update(GR=GR, PP=PP, G3A=G3A, cfg=cfg, cfg3=cfg3, ctx=ctx, gcfg=gcfg, mf=dict(cut=cut), grid=GR.make_grid(gcfg), fine=GR.make_grid(gcfg, fine=True))
    G["fine_key"] = {key(c): c["id"] for c in G["fine"]}
    return G


def key(c):
    return (round(float(c["inc"]), 3), round(float(c["epoch_off"]), 3), round(float(c["omega"]) % 360.0, 3), round(float(c["raan"]) % 360.0, 3))


def load_fine():
    if "fine" not in G["mf"]:
        G["mf"]["fine"] = G["GR"].load_meanfield(TD / "analytic_data" / "meanfield_fine", G["gcfg"])
    return G["mf"]["fine"]


# ----------------------------------------------------------------------------------------------------- the ground table of an arbitrary geometry

def build(cand, field, jitter=0.0, q_req=None, preview=False):
    """Burn-free channels of one geometry (geometry_rank model at the calibrated planner context). preview: no mean-field term and the fallback scale height (what the generic tables were)."""
    ctx = G["ctx"]; GR = G["GR"]; b = ctx.cfg["band"]; H = None if (field is None or preview) else field["H"]
    if q_req is None:
        q_req = GR.solve_q_req(ctx, cand, b["d_ref_days"], H)
    S = GR.schedule(ctx, cand, q_req * (1.0 + jitter), H)
    d = GR.channels(ctx, cand, S, None if preview else field)
    return dict(q_req=q_req, S=S, d=d)


@contextlib.contextmanager
def orbit_target(ctx, ra_t, max_passes=None):
    """Revision 25 (additive): run geometry_rank.schedule with the campaign stopping at apoapsis radius ra_t (and, if given, max_passes); restores the config. ra_t None = unchanged."""
    o = ctx.cfg["orbit"]; old = (o["ra_target_m"], o["max_passes"])
    if ra_t is not None: o["ra_target_m"] = float(ra_t)
    if max_passes is not None: o["max_passes"] = int(max_passes)
    try:
        yield
    finally:
        o["ra_target_m"], o["max_passes"] = old


def table_from(bd, keys=None):
    """Ground-table dict (planner_predictions.ground_table layout) from channels; keys = subset of channels for the path (None = all, the real table)."""
    ctx = G["ctx"]; GR = G["GR"]; b = ctx.cfg["band"]; S = bd["S"]; n = S["n"]
    path = GR.path_of(bd["d"], keys); slow = GR.split_slow_fast(path, b["fast_window_passes"])[0]
    w = np.zeros(n); w[1:] = S["H"][1:] / (ctx.ratio * np.array([GR.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, n)]))
    W = b["slack_band_factor"] * math.log(b["q_cap_mean_w_cm2"] / bd["q_req"])
    return dict(n=int(n), q_req=float(bd["q_req"]), W=float(W), path=path, slow=slow, d_slow=np.diff(slow), fast=path - slow, w=w, H=np.asarray(S["H"], float), a=np.asarray(S["a"], float),
                e=np.asarray(S["e"], float), t=np.asarray(S["t"], float), rp=np.asarray(S["rp"], float), ra=np.asarray(S["ra"], float), vrel=np.asarray(S["V"], float), frot=np.asarray(S["frot"], float))


# ----------------------------------------------------------------------------------------------------- the DP on an arbitrary table

@contextlib.contextmanager
def patched(tab, delta=0.0, reserve=None, harden=(0.0, 0.0), ra_t=None, horizon_d=None):
    """planner_rev20_ic.Case reads load_tab, X_REF and P as module globals. Inside: the table is `tab`; the pass-1 reference level is ln Q_FIRST + delta (an INFORMED uniform density offset delta in ln: the realized
    ln q of every pass is higher by delta at the same periapsis; pre-burn weights and the periapsis-shift bookkeeping of the arrival are unchanged); the tube reserve is `reserve` (ln); harden = (dlo, dhi): the
    planned tube [lo, hi] becomes [lo + dlo, hi - dhi] (a HARDENED tube: the plan keeps an unknown offset of up to dhi (+) or dlo (-) inside the original tube)."""
    old = (R.load_tab, R.X_REF, R.P, PLB.ol_inputs, R.RA_T, O.RA_T, R.HORIZON_D)      # Revision 24 (additive): ra_t = closing apoapsis radius of the campaign, horizon_d = arrival horizon
    if ra_t is not None: R.RA_T = O.RA_T = float(ra_t)
    if horizon_d is not None: R.HORIZON_D = float(horizon_d)
    t = dict(tab); t["token"] = "X"
    R.load_tab = lambda tok: dict(t)
    R.X_REF = math.log(R.Q_FIRST) + float(delta)
    if reserve is not None:
        Pm = dict(R.P); Pm["fast_reserve_ln"] = float(reserve); R.P = Pm
    if harden != (0.0, 0.0):
        orig = PLB.ol_inputs
        def hardened(sw, tb, P_, q1):
            inp = orig(sw, tb, P_, q1); inp["lo"] = [x + harden[0] for x in inp["lo"]]; inp["hi"] = [x - harden[1] for x in inp["hi"]]; return inp
        PLB.ol_inputs = hardened
    try:
        yield
    finally:
        R.load_tab, R.X_REF, R.P, PLB.ol_inputs, R.RA_T, O.RA_T, R.HORIZON_D = old


def dp(tab, scenario, kappa=KAPPA_DEFAULT, pen=0.0, delta=0.0, reserve=None, horizon="arrival", harden=(0.0, 0.0), tstar=None, ra_t=None, horizon_d=None):
    """Aerobraking-dV plan of one scenario ('ref' | 'S-fixed' | 'S-free') on a table. tstar (days): arrival window |t - tstar| <= W_max (3 final periods, ASSUMED Revision 18) by the Revision 18 price
    (Case.window_plan); None = free arrival (Revision 20 as written). Returns a flat dict; feasible False (with a reason) if no plan exists."""
    if harden[0] + harden[1] >= tab["W"] - (RESERVE0 if reserve is None else reserve):
        return dict(feasible=False, why="tube_closed")
    with patched(tab, delta, reserve, harden, ra_t, horizon_d):
        try:
            c = R.Case("X", scenario, 250e3, False, kappa, horizon, 0.0, pen); p = c.nominal; win = None
            if tstar is not None:
                w = c.window_plan(float(tstar), R.W_MAX_PERIODS * p["T_final_s"] / 86400.0); p = w["plan"]; win = bool(w["feasible"])
        except Exception as ex:                                    # no DP solution (no self-consistent arrival horizon, start outside the tube): infeasible
            return dict(feasible=False, why="dp:" + type(ex).__name__)
    if p is None or not math.isfinite(p["days"]):
        return dict(feasible=False, why="no_arrival")
    ok = (not p["beyond_table"]) and p["days"] <= (R.HORIZON_D if horizon_d is None else horizon_d) and (win is None or win)
    return dict(feasible=bool(ok), why="" if ok else ("beyond_table" if p["beyond_table"] else ("window" if win is False else "after_horizon")), aero=float(p["dv_campaign"] + p["dv_pre"]), campaign=float(p["dv_campaign"]),
                pre=float(p["dv_pre"]), pre_equiv=float(p["dv_pre_equiv"]), N=int(p["N"]), days=float(p["days"]), burns=int(p["burns"]), x0=float(p["x0"]), nu_pre=float(p["s"]), margin_x0=float(p["margin_hi_x0"]),
                path_margin_min=float(p["path_margin_min"]), fb_ignored=float(p["fb"]), rp_final=float(p["rp_final"]))


# ----------------------------------------------------------------------------------------------------- one candidate

def chan_stats(bd):
    GR = G["GR"]; win = G["ctx"].cfg["band"]["fast_window_passes"]; out = {}
    for k, v in bd["d"].items():
        p = GR.path_of({k: v}); s = GR.split_slow_fast(p, win)[0]
        out[k] = dict(net=float(p[-1]), range_slow=float(s.max() - s.min()), rms_fast=float(np.std(p - s)))
    return out


def analytic_params(cand, bd):
    """Analytic perturbation parameters of the geometry (DERIVED from the schedule S): J2 apse rotation and node drift over the campaign, critical-inclination proximity, latitude of periapsis."""
    S = bd["S"]; inc = math.radians(cand["inc"])
    return dict(apse_rot_total_deg=float(math.degrees(S["omega"][-1] - S["omega"][0])), apse_rot_per_orbit_deg=float(math.degrees((S["omega"][-1] - S["omega"][0]) / max(S["n"] - 1, 1))),
                node_drift_total_deg=float(math.degrees(S["raan"][-1] - S["raan"][0])), raan_arrival_deg=float(math.degrees(S["raan"][-1]) % 360.0), omega_arrival_deg=float(math.degrees(S["omega"][-1]) % 360.0),
                crit_incl_proximity=float(abs(math.sin(inc) ** 2 - 0.8)), lat_pk_range_deg=float(math.degrees(S["lat"].max() - S["lat"].min())), lat_pk_abs_max_deg=float(math.degrees(np.abs(S["lat"]).max())),
                frot_mean=float(np.mean(S["frot"])), vrel_start=float(S["V"][0]))


def feas_flags(tab, d_s):
    """Feasibility terms. tube: W_eff = W - fast reserve > 0; q_cap: q_req < q_cap; never65: q_cap * exp(0.856) < 0.65 (same constants for every geometry, ASSUMED transfer); target/arrival from the DP plan."""
    b = G["ctx"].cfg["band"]; Weff = tab["W"] - RESERVE0
    return dict(tube=bool(Weff > 0), W_eff=float(Weff), q_cap=bool(tab["q_req"] < b["q_cap_mean_w_cm2"]), never65=bool(b["q_cap_mean_w_cm2"] * math.exp(b["ln_qs_over_qmean_max"]) < b["q_never_exceed_qs_w_cm2"]),
                target_7408=bool(tab["ra"][-1] > 0 and tab["n"] < G["gcfg"]["orbit"]["max_passes"]), arrival=bool(d_s.get("feasible", False)))


def run_base(tab, kappa=KAPPA_DEFAULT, tstar=None):
    return {sc: dp(tab, sc, kappa, tstar=tstar) for sc in ("ref", "S-fixed", "S-free")}


def both(tab, scs=("S-fixed", "S-free"), **kw):
    """Free-arrival and common-arrival (t* = D_ref) DP of the same variant."""
    T = G["cfg"]["planner"]["d_ref_days"]
    return dict(free={sc: dp(tab, sc, **kw) for sc in scs}, tstar={sc: dp(tab, sc, tstar=T, **kw) for sc in scs})


def aerov(r):
    return r["aero"] if r.get("feasible") else float("nan")


def analyze(task):
    """task = dict(id, fid, cand, src ('cut'|'fine'), level ('base'|'full')). Returns a flat JSON-able row. Every DP number is the aerobraking dV (m/s) of its scenario; 'free' = arrival chosen by the plan
    (Revision 20), 'tstar' = common arrival t* = D_ref = planner_config d_ref_days within W_max."""
    G_ = init(); cand = task["cand"]; src = task["src"]
    mf = G_["mf"][src] if src != "fine" else load_fine()
    fid = task["fid"]; field = mf.get(fid) if fid else None
    t0 = time.time(); T = G_["cfg"]["planner"]["d_ref_days"]
    bd = build(cand, field); tab = table_from(bd)
    row = dict(id=task["id"], fid=fid, src=src, inc=cand["inc"], omega0=cand["omega"], raan0=cand["raan"], epoch_off=cand["epoch_off"], field_found=field is not None, n=tab["n"], q_req=tab["q_req"], W=tab["W"],
               days_model=float(G_["GR"].campaign_days(bd["S"])), R_slow=float(tab["slow"].max() - tab["slow"].min()), N_net=float(tab["path"][-1]), role=task.get("role"))
    row.update(analytic_params(cand, bd)); row["chan"] = chan_stats(bd)
    base = run_base(tab); row["dp"] = base; row["dpt"] = {sc: dp(tab, sc, tstar=T) for sc in ("S-fixed", "S-free")}
    row["feas"] = feas_flags(tab, base["S-fixed"])
    row["feasible"] = bool(all(row["feas"][k] for k in ("tube", "q_cap", "never65", "target_7408")) and base["S-free"].get("feasible", False) and base["S-fixed"].get("feasible", False))
    row["feasible_tstar"] = bool(all(row["feas"][k] for k in ("tube", "q_cap", "never65", "target_7408")) and all(row["dpt"][sc].get("feasible", False) for sc in row["dpt"]))
    if task["level"] == "full":
        row["dp_pen01"] = {sc: dp(tab, sc, pen=0.01) for sc in ("S-fixed", "S-free")}
        row["informed"] = {("%+.4f" % s): {sc: dp(tab, sc, delta=s) for sc in ("S-fixed", "S-free")} for s in (+OFFSET_LN, -OFFSET_LN)}
        row["harden"] = {"%+.4f" % (+OFFSET_LN): both(tab, harden=(0.0, OFFSET_LN)), "%+.4f" % (-OFFSET_LN): both(tab, harden=(OFFSET_LN, 0.0))}
        row["reserveL1"] = both(tab, reserve=RESERVE0 * L1_SCALE)
        row["dp_kappa_g1"] = {sc: dp(tab, sc, kappa=KAPPA_G1) for sc in ("S-fixed", "S-free")}
        row["dp_ablate"] = {c: {sc: dp(table_from(bd, [k for k in bd["d"] if k != c]), sc) for sc in ("S-fixed", "S-free")} for c in ABL}
        row["dp_ablate"]["only_drag_speed_geometry"] = {sc: dp(table_from(bd, ["drag", "speed", "geometry"]), sc) for sc in ("S-fixed", "S-free")}
        bp = build(cand, field, preview=True); tp = table_from(bp); row["dp_preview"] = {sc: dp(tp, sc) for sc in ("S-fixed", "S-free")}; row["preview_q_req"] = tp["q_req"]
        ens = {"free": {"S-fixed": [], "S-free": []}, "tstar": {"S-fixed": [], "S-free": []}}
        for j in np.linspace(-JIT_FRAC, JIT_FRAC, JIT_N):
            bj = build(cand, field, jitter=float(j), q_req=bd["q_req"]); tj = table_from(bj)
            for sc in ("S-fixed", "S-free"):
                ens["free"][sc].append(aerov(dp(tj, sc)))
        for j in np.linspace(-JIT_FRAC, JIT_FRAC, JIT_NT):
            bj = build(cand, field, jitter=float(j), q_req=bd["q_req"]); tj = table_from(bj)
            for sc in ("S-fixed", "S-free"):
                ens["tstar"][sc].append(aerov(dp(tj, sc, tstar=T)))
        row["jitter"] = {m: {sc: dict(mean=float(np.nanmean(v)) if np.isfinite(v).any() else float("nan"), sd=float(np.nanstd(v)) if np.isfinite(v).any() else float("nan"), n_ok=int(np.isfinite(v).sum()), vals=[float(x) for x in v]) for sc, v in d.items()} for m, d in ens.items()}
        # geometry_rank floor in the SAME calibrated context (comparison only; the frozen ranking used q_cap 0.234 / D_ref 64.2)
        r_, *_ = G_["GR"].evaluate(G_["ctx"], cand, field, G_["ctx"].s_t_ody); row["floor_cal_mps"] = float(r_["floor_mps"]); row["floor_driftfree_cal_mps"] = float(r_["floor_driftfree_mps"]); row["regime_cal"] = r_["regime"]
    row["wall_s"] = time.time() - t0
    return row


# ----------------------------------------------------------------------------------------------------- candidate lists

def cut256_tasks():
    G_ = init(); GR = G_["GR"]; out = []
    for c in G_["grid"]:
        if GR.in_cut(c, G_["gcfg"], "256"):
            out.append(dict(id=c["id"], fid=c["id"], cand=c, src="cut", level="full"))
    out.append(dict(id="ODY", fid="ODY", cand=GR.anchor_cand(G_["gcfg"]), src="cut", level="full"))
    return out


def anchor_tasks():
    """Flown/anchor geometries with the mean field of the fine sampling (found by their elements, not by the id of the generic tables, which has no field)."""
    G_ = init(); cfg = G_["cfg"]; out = []
    for gid, o in cfg["planner"]["geometry"].items():
        if gid == "ODY":
            continue
        c = dict(id=gid, inc=float(o["inc_deg"]), epoch_off=float(o["epoch_off_d"]), omega=float(o["omega_deg"]), raan=float(o["raan_deg"]))
        out.append(dict(id=gid, fid=G_["fine_key"].get(key(c)), cand=c, src="fine", level="full", role=o["role"]))
    for g in G_["cfg3"]["stageg3"]["geometry"]:          # G1 (G5863), G2 (G2348), Gworst (G3503) of the gate table
        c = G_["G3A"].cand_of(g)
        out.append(dict(id=g["id"], fid=g["id"], cand=c, src="fine", level="full", role=g["role"]))
    return out


def grid_tasks(skip):
    """The 1,536 grid (R2.3) minus the ids already in `skip`: base DP only, fine mean field (the same geometry under its fine id)."""
    G_ = init(); out = []
    for c in G_["grid"]:
        if c["id"] in skip:
            continue
        out.append(dict(id=c["id"], fid=G_["fine_key"].get(key(c)), cand=c, src="fine", level="base"))
    return out


def pmap(fn, items, jobs, label=""):
    import multiprocessing as mp
    res = []; t0 = time.time()
    if jobs <= 1:
        for i, x in enumerate(items):
            res.append(fn(x))
            if (i + 1) % 10 == 0: print("%s %d/%d %.0f s" % (label, i + 1, len(items), time.time() - t0), flush=True)
        return res
    with mp.get_context("fork").Pool(jobs) as pool:
        for i, r in enumerate(pool.imap(fn, items, chunksize=1)):
            res.append(r)
            if (i + 1) % 10 == 0: print("%s %d/%d %.0f s" % (label, i + 1, len(items), time.time() - t0), flush=True)
    return res


def write_json(p, obj):
    p = pathlib.Path(p); p.parent.mkdir(parents=True, exist_ok=True); p.write_text(json.dumps(obj, indent=0, default=float)); return hashlib.sha256(p.read_bytes()).hexdigest()


def cmd_run(part, jobs):
    init()
    if part in ("anchors", "grid"):
        load_fine()
    if part == "cut":
        rows = pmap(analyze, cut256_tasks(), jobs, "cut256"); print("sha", write_json(OUT / "raw_cut256.json", rows)[:16])
    elif part == "anchors":
        rows = pmap(analyze, anchor_tasks(), jobs, "anchors"); print("sha", write_json(OUT / "raw_anchors.json", rows)[:16])
    elif part == "grid":
        skip = {t["id"] for t in cut256_tasks()}; rows = pmap(analyze, grid_tasks(skip), jobs, "grid1536"); print("sha", write_json(OUT / "raw_grid1536_rest.json", rows)[:16])


# ----------------------------------------------------------------------------------------------------- rendezvous plane coupling (pure functions)

def rel_inc_deg(i1, O1, i2, O2):
    """Angle between two orbit planes (deg) from inclinations and nodes (deg): cos di = cos i1 cos i2 + sin i1 sin i2 cos(O1 - O2)."""
    a, b, d = map(math.radians, (i1, i2, O1 - O2)); c = math.cos(a) * math.cos(b) + math.sin(a) * math.sin(b) * math.cos(d)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def plane_cost_quadrature(di_deg):
    """Increment over the coplanar match at r_a, as the contract states it (Revision 22b): sqrt(V0^2 + (2 v_a sin(di/2))^2) - V0, V0 = 477.5, 2 v_a sin(0.5 deg) = 33.63 m/s per deg. DERIVED."""
    return math.sqrt(V0 ** 2 + (2 * V_A * math.sin(math.radians(di_deg) / 2)) ** 2) - V0


def plane_cost_exact(di_deg):
    """Same increment from the exact vis-viva burn at the apoapsis node, |v_circ rotated by di - v_a| = sqrt((v_c - v_a)^2 + 4 v_c v_a sin^2(di/2)) (DERIVED; agrees with the Revision 22b sweep 479.02 m/s at 1 deg)."""
    return math.sqrt((V_C - V_A) ** 2 + 4 * V_C * V_A * math.sin(math.radians(di_deg) / 2) ** 2) - (V_C - V_A)


def breakeven_di_deg(saving_mps):
    """Relative inclination at which the quadrature plane cost equals a given aerobraking saving: (2 v_a sin(di/2))^2 = 2 V0 S + S^2. DERIVED."""
    return math.degrees(2 * math.asin(math.sqrt(2 * V0 * saving_mps + saving_mps ** 2) / (2 * V_A)))


# ----------------------------------------------------------------------------------------------------- selftest

def selftest():
    G_ = init(); GR = G_["GR"]; PP = G_["PP"]; cfg = G_["cfg"]; geoms = PP.geometries(cfg, G_["cfg3"]); mfs = dict(ody=G_["mf"]["cut"]["ODY"]); lines = []
    head = json.load(open(TD / "rev20_data" / "ic" / "headline.json"))
    kap = R.load_kappa(); worst = 0.0
    for tok in ("ODY", "G1"):
        g = geoms[tok]; field, src = G_["G3A"].meanfield_for(G_["cfg3"], G_["ctx"], None, g["g"])
        # (1) the rebuilt table is the stored table, bit for bit
        tab = table_from(build(g["cand"], field)); old = O.load_tab(tok)
        for k in ("path", "w", "H", "ra", "rp", "t", "a", "e", "slow", "d_slow"):
            assert np.array_equal(tab[k], np.asarray(old[k], float)), (tok, k)
        assert tab["n"] == old["n"] and tab["q_req"] == old["q_req"] and tab["W"] == old["W"]
        # (2) the DP on the rebuilt table equals select_ic on the stored table, and the Revision 20 headline file, within 0.01 m/s
        for sc in ("ref", "S-fixed", "S-free"):
            ic = R.select_ic(tok, sc); mine = dp(tab, sc, kap[tok]); ref = ic["campaign_dv_mps"] + ic["dv_pre_mps"]
            hp = [p for p in head["plans"] if p["geom"] == tok and p["scenario"] == sc and p["fb_term"] is False and p["horizon"] == "arrival" and abs(p["delta_ln"]) < 1e-12 and abs(p["h_des_m"] - 250000.0) < 1e-6 and p["t_star_d"] is None]
            assert len(hp) == 1, (tok, sc, len(hp)); hv = hp[0]["campaign_dv_mps"] + hp[0]["dv_pre_mps"]
            d1, d2 = abs(mine["aero"] - ref), abs(mine["aero"] - hv); worst = max(worst, d1, d2)
            assert d1 <= 0.01 and d2 <= 0.01 and mine["N"] == ic["N"], (tok, sc, mine["aero"], ref, hv)
            lines.append("%s %-8s DP %.6f  select_ic %.6f  headline.json %.6f  N %d  (diff %.1e / %.1e)" % (tok, sc, mine["aero"], ref, hv, mine["N"], d1, d2))
    lines.append("anchor check ok: worst difference %.2e m/s (threshold 0.01)" % worst)
    # (3) module state is restored by `patched`; offset 0 is the nominal; reserve patch changes only the reserve
    assert R.X_REF == math.log(R.Q_FIRST) and R.load_tab.__name__ == "load_tab" and R.P["fast_reserve_ln"] == RESERVE0
    t_ody = table_from(build(geoms["ODY"]["cand"], G_["mf"]["cut"]["ODY"]))
    assert dp(t_ody, "S-free", delta=0.0)["aero"] == dp(t_ody, "S-free")["aero"]
    assert dp(t_ody, "S-fixed", reserve=RESERVE0 * 1.5)["aero"] >= dp(t_ody, "S-fixed")["aero"] - 1e-12      # a larger reserve shrinks the tube: never cheaper
    assert R.X_REF == math.log(R.Q_FIRST) and R.P["fast_reserve_ln"] == RESERVE0
    # (4) plane geometry: equal-inclination planes di = 2 asin(sin i sin(dOmega/2)); coplanar 0; quadrature break-even for 3.6 m/s is about 1.75 deg (hand, Revision 23); exact law agrees at 1 deg with Rev 22b
    for inc, dO in ((93.5, 10.0), (63.4, 3.0), (116.6, 40.0)):
        assert abs(rel_inc_deg(inc, 0.0, inc, dO) - math.degrees(2 * math.asin(math.sin(math.radians(inc)) * math.sin(math.radians(dO) / 2)))) < 1e-9
    assert rel_inc_deg(60.0, 100.0, 60.0, 100.0) == 0.0 and abs(rel_inc_deg(60.0, 0.0, 70.0, 33.0) - rel_inc_deg(70.0, 33.0, 60.0, 0.0)) < 1e-12
    assert abs(breakeven_di_deg(3.6) - 1.747) < 0.01 and abs(plane_cost_quadrature(breakeven_di_deg(3.6)) - 3.6) < 1e-6
    assert abs(2 * V_A * math.sin(math.radians(0.5)) - 33.63) < 0.02 and abs((V0 + plane_cost_exact(1.0)) - 479.02) < 0.1
    import planner_rev22_rendezvous as R22
    F = dict(i=math.radians(93.5), Om=math.radians(10.0), om=0.0); T = dict(i=math.radians(93.5), Om=math.radians(12.0), om=0.0)
    assert abs(math.degrees(R22.relative_plane(F, T)[0]) - rel_inc_deg(93.5, 10.0, 93.5, 12.0)) < 1e-9
    lines.append("plane formulas ok (quadrature break-even for 3.6 m/s: %.3f deg; exact-law 1 deg total %.2f m/s vs Rev 22b 479.02)" % (breakeven_di_deg(3.6), V0 + plane_cost_exact(1.0)))
    # (5) a cut-256 field equals the same geometry's fine-field entry (the two sampling jobs agree), one candidate
    c = next(c for c in G_["grid"] if G_["GR"].in_cut(c, G_["gcfg"], "256") and c["inc"] == 63.4); fid = G_["fine_key"][key(c)]; mfn = load_fine()
    for k in ("H", "lnrho0"):
        assert np.allclose(G_["mf"]["cut"][c["id"]][k], mfn[fid][k], rtol=0, atol=1e-9), (c["id"], fid, k)
    lines.append("cut-256 sample %s equals fine sample %s (H, ln rho): ok" % (c["id"], fid))
    OUT.mkdir(parents=True, exist_ok=True); (OUT / "selftest.txt").write_text("\n".join(lines) + "\nplanner_rev23_geometry selftest ok\n"); print("\n".join(lines)); print("planner_rev23_geometry selftest ok")


# ----------------------------------------------------------------------------------------------------- summaries (read raw_*.json; write the summary JSONs and tables)

MODES = ("free", "tstar"); SCS = ("S-fixed", "S-free"); SC_SHORT = {"S-fixed": "fx", "S-free": "fr"}
TOL_FREEDOM = 0.25            # ASSUMED: a geometry is "good" in a case if its net dV is within this of the best net dV of the case's reference (roughly the jitter sd of one realization)
NAN = float("nan")


def _a(d, sc):
    x = (d or {}).get(sc) or {}
    return x["aero"] if x.get("feasible") else NAN


def flat(rows):
    """One row per candidate with the aerobraking dV (m/s) of every variant: columns <variant>_<mode>_<fx|fr>; NaN where the DP is infeasible."""
    out = []
    for r in rows:
        f = dict(id=r["id"], fid=r["fid"], role=r.get("role"), inc=r["inc"], omega0=r["omega0"], raan0=r["raan0"], epoch=r["epoch_off"], n=r["n"], q_req=r["q_req"], W=r["W"], days_model=r["days_model"],
                 R_slow=r["R_slow"], N_net=r["N_net"], feasible=r["feasible"], feasible_tstar=r["feasible_tstar"], raan_arr=r["raan_arrival_deg"], omega_arr=r["omega_arrival_deg"],
                 node_drift=r["node_drift_total_deg"], apse_rot=r["apse_rot_total_deg"], crit=r["crit_incl_proximity"], lat_rng=r["lat_pk_range_deg"], field_found=r["field_found"])
        for k in ("tube", "q_cap", "never65", "target_7408", "arrival"):
            f["feas_" + k] = r["feas"][k]
        for sc in SCS:
            s = SC_SHORT[sc]; f["base_free_" + s] = _a(r["dp"], sc); f["base_tstar_" + s] = _a(r["dpt"], sc)
            f["days_free_" + s] = (r["dp"][sc].get("days", NAN)); f["N_free_" + s] = (r["dp"][sc].get("N", NAN)); f["burns_free_" + s] = r["dp"][sc].get("burns", NAN)
            f["nu_pre_free_" + s] = r["dp"][sc].get("nu_pre", NAN)
        f["ref_free"] = _a(r["dp"], "ref")
        if "harden" in r:
            for m in MODES:
                for sg, nm in (("+0.2231", "hp"), ("-0.2231", "hm")):
                    for sc in SCS:
                        f["%s_%s_%s" % (nm, m, SC_SHORT[sc])] = _a(r["harden"][sg][m], sc)
                for sc in SCS:
                    f["rl1_%s_%s" % (m, SC_SHORT[sc])] = _a(r["reserveL1"][m], sc); f["jm_%s_%s" % (m, SC_SHORT[sc])] = r["jitter"][m][sc]["mean"]; f["js_%s_%s" % (m, SC_SHORT[sc])] = r["jitter"][m][sc]["sd"]
                    f["jn_%s_%s" % (m, SC_SHORT[sc])] = r["jitter"][m][sc]["n_ok"]
            for sc in SCS:
                s = SC_SHORT[sc]
                f["ip_free_" + s] = _a(r["informed"]["+0.2231"], sc); f["im_free_" + s] = _a(r["informed"]["-0.2231"], sc); f["pen_free_" + s] = _a(r["dp_pen01"], sc)
                f["kg1_free_" + s] = _a(r["dp_kappa_g1"], sc); f["prev_free_" + s] = _a(r["dp_preview"], sc)
                for c, v in r["dp_ablate"].items():
                    f["abl_%s_free_%s" % (c, s)] = _a(v, sc)
            f["floor_cal"] = r["floor_cal_mps"]; f["floor_driftfree_cal"] = r["floor_driftfree_cal_mps"]; f["regime"] = r["regime_cal"]
            for c, v in r["chan"].items():
                f["rng_" + c] = v["range_slow"]; f["net_" + c] = v["net"]
        out.append(f)
    return pd.DataFrame(out)


def spearman(x, y):
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4: return NAN
    return float(spearmanr(x[m], y[m])[0])


def order(df, col, tie="id"):
    """Ranking order of a column (ascending dV; NaN last): ties are broken by the jitter-ensemble mean of the same variant when present, then the id (so the set is deterministic)."""
    d = df.copy(); d["_k"] = (d[col] / TIE_MPS).round().fillna(1e9)      # values within TIE_MPS are ties
    keys = ["_k"] + ([tie] if tie in d.columns else []) + ["id"]
    return d.sort_values(keys)


def spread(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x): return {}
    q = np.percentile(x, [0, 10, 25, 50, 75, 90, 100]); return dict(n=int(len(x)), min=q[0], p10=q[1], p25=q[2], median=q[3], p75=q[4], p90=q[5], max=q[6], n_le_0p01=int((x <= 0.01).sum()), n_le_0p1=int((x <= 0.1).sum()),
                                                                     n_le_0p5=int((x <= 0.5).sum()), n_ge_1=int((x >= 1).sum()), n_ge_3=int((x >= 3).sum()), n_tie_best=int((x <= x.min() + TIE_MPS).sum()))


def eta2(df, col, by):
    d = df[[col, by]].dropna(); g = d.groupby(by)[col]; sst = ((d[col] - d[col].mean()) ** 2).sum()
    return float((g.size() * (g.mean() - d[col].mean()) ** 2).sum() / sst) if sst > 0 else NAN


def top_overlap(df, col_nom, col_var, k=10, tie_nom=None):
    """Top-k sets of the nominal and of the variant ranking; ties (within TIE_MPS) are broken by the SAME key in both, the nominal jitter-ensemble mean tie_nom, then the id. Also: how many nominal top-k members have no
    feasible plan in the variant (NaN)."""
    a = list(order(df, col_nom, tie_nom).head(k).id); b = list(order(df, col_var, tie_nom).head(k).id)
    return dict(overlap=int(len(set(a) & set(b))), nominal=a, variant=b, n_nan_variant=int(df[col_var].isna().sum()), top_nominal_nan_in_variant=int(df.set_index("id").loc[a, col_var].isna().sum()))


def band_overlap(df, col_nom, col_var, k=10, band=0.1):
    """Of the nominal top-k, how many are still within `band` m/s of the variant's best value (an ASSUMED tolerance, robust to ties)."""
    a = order(df, col_nom).head(k); vmin = df[col_var].min()
    return int((df.set_index("id").loc[a.id, col_var] <= vmin + band).sum())


def _fb_one(t):
    G_ = init(); T = G_["cfg"]["planner"]["d_ref_days"]; r = t["row"]; cand = dict(id=r["id"], inc=r["inc"], epoch_off=r["epoch_off"], omega=r["omega0"], raan=r["raan0"])
    field = (G_["mf"]["cut"] if r["src"] == "cut" else load_fine()).get(r["fid"])
    return dp(table_from(build(cand, field)), t["sc"], horizon="full", tstar=T)


def fallback(rows, df, jobs=4):
    """Where the free-arrival plan of a scenario does not exist (the cost-free optimum arrives beyond the table + tail, Revision 20 caveat 9), the ranking value is the deadline plan: horizon 'full', arrival window
    t* = D_ref = 70.4 d (planner_config d_ref_days), W_max. base_free_* then holds the completed value; strict_free_* keeps the strict one (NaN where no free plan exists); comp_flag_* marks the completed rows."""
    byid = {r["id"]: r for r in rows}; tasks = []; where = []
    for sc in SCS:
        s = SC_SHORT[sc]; df["strict_free_" + s] = df["base_free_" + s]; df["comp_flag_" + s] = False
        for i in df.index[df["base_free_" + s].isna()]:
            tasks.append(dict(row=byid[df.at[i, "id"]], sc=sc)); where.append((i, s))
    load_fine(); res = pmap(_fb_one, tasks, jobs if len(tasks) > 8 else 1, "fallback")
    for (i, s), v in zip(where, res):
        if v.get("feasible"):
            df.at[i, "base_free_" + s] = v["aero"]; df.at[i, "comp_flag_" + s] = True
    return df


def summarize():
    init(); rows = json.load(open(OUT / "raw_cut256.json")); df = fallback(rows, flat(rows)); S = {}
    ody = df[df.id == "ODY"].iloc[0]; cut = df[df.id != "ODY"].reset_index(drop=True); S["n_cut"] = int(len(cut))
    S["labels"] = "DERIVED unless stated; model world (noise free, exact burns, superposition); kappa 1.0 ASSUMED for every row; ODY row at kappa 1.0"
    # --- feasibility
    S["feasibility"] = dict(n=int(len(cut)), tube=int(cut.feas_tube.sum()), q_cap=int(cut.feas_q_cap.sum()), never65=int(cut.feas_never65.sum()), target_7408=int(cut.feas_target_7408.sum()), feasible_free=int(cut.feasible.sum()),
                            feasible_tstar=int(cut.feasible_tstar.sum()), infeasible_free_ids=list(cut[~cut.feasible].id), infeasible_tstar_ids=list(cut[~cut.feasible_tstar].id), W_eff_min=float((cut.W - RESERVE0).min()))
    # --- spread per mode and scenario (feasible only)
    S["spread"] = {m + "_" + SC_SHORT[sc]: spread(cut["base_%s_%s" % (m, SC_SHORT[sc])]) for m in MODES for sc in SCS}
    S["odyssey_anchor"] = {m + "_" + SC_SHORT[sc]: float(ody["base_%s_%s" % (m, SC_SHORT[sc])]) for m in MODES for sc in SCS}
    S["odyssey_anchor"]["ref_free"] = float(ody["ref_free"])
    S["odyssey_rank"] = {m + "_" + SC_SHORT[sc]: int((cut["base_%s_%s" % (m, SC_SHORT[sc])] < ody["base_%s_%s" % (m, SC_SHORT[sc])]).sum()) for m in MODES for sc in SCS}
    # --- by inclination class
    S["by_inc"] = {}
    for inc, g in cut.groupby("inc"):
        S["by_inc"][str(inc)] = {m + "_" + SC_SHORT[sc]: dict(median=float(g["base_%s_%s" % (m, SC_SHORT[sc])].median()), min=float(g["base_%s_%s" % (m, SC_SHORT[sc])].min()), max=float(g["base_%s_%s" % (m, SC_SHORT[sc])].max()))
                                 for m in MODES for sc in SCS}
    # --- variance shares of the geometry axes (eta^2)
    S["eta2"] = {m + "_" + SC_SHORT[sc]: {by: eta2(cut, "base_%s_%s" % (m, SC_SHORT[sc]), by) for by in ("inc", "omega0", "raan0", "epoch")} for m in MODES for sc in SCS}
    # --- correlations with the geometry_rank floor and with the drift range R (P23-1, P23-3)
    fine = pd.read_csv(TD / "analytic_data" / "rank_fine" / "geometry_rank.csv").set_index("id"); G_ = init()
    cutc = {c["id"]: c for c in G_["grid"]}
    cut["fine_id"] = [G_["fine_key"][key(cutc[i])] for i in cut.id]
    cut["floor_fine_det"] = [float(fine.loc[i, "floor_mps"]) for i in cut.fine_id]; cut["floor_fine_ens"] = [float(fine.loc[i, "floor_ens_mps"]) if np.isfinite(fine.loc[i, "floor_ens_mps"]) else NAN for i in cut.fine_id]
    cut["R_fine"] = [float(fine.loc[i, "R"]) for i in cut.fine_id]
    S["corr"] = {}
    for m in MODES:
        for sc in SCS:
            c_ = "base_%s_%s" % (m, SC_SHORT[sc]); S["corr"][m + "_" + SC_SHORT[sc]] = dict(floor_fine_det=spearman(cut[c_], cut.floor_fine_det), floor_cal=spearman(cut[c_], cut.floor_cal), n_ens=int(cut.floor_fine_ens.notna().sum()),
                                                                                           floor_fine_ens=spearman(cut[c_], cut.floor_fine_ens), R_slow=spearman(cut[c_], cut.R_slow), R_fine=spearman(cut[c_], cut.R_fine),
                                                                                           node_drift=spearman(cut[c_], cut.node_drift.abs()), apse_rot_abs=spearman(cut[c_], cut.apse_rot.abs()), crit=spearman(cut[c_], cut.crit),
                                                                                           lat_range=spearman(cut[c_], cut.lat_rng), jitter_mean=spearman(cut[c_], cut["jm_%s_%s" % (m, SC_SHORT[sc])]))
    S["corr"]["free_vs_tstar"] = {SC_SHORT[sc]: spearman(cut["base_free_" + SC_SHORT[sc]], cut["base_tstar_" + SC_SHORT[sc]]) for sc in SCS}
    S["corr"]["fixed_vs_free"] = {m: spearman(cut["base_%s_fx" % m], cut["base_%s_fr" % m]) for m in MODES}
    # --- drivers: ablation (channel removed) and per-channel slow range
    chans = [c for c in ABL] + ["only_drag_speed_geometry"]; S["ablation"] = {}
    for c in chans:
        for sc in SCS:
            s = SC_SHORT[sc]; base = cut["base_free_" + s]; ab = cut["abl_%s_free_%s" % (c, s)]; dlt = ab - base
            S["ablation"][c + "_" + s] = dict(n_ok=int(dlt.notna().sum()), median_change=float(dlt.median()), mean_change=float(dlt.mean()), median_dV_without=float(ab.median()), share_improved_by_gt_0p1=float((dlt < -0.1).mean()),
                                                share_worsened_by_gt_0p1=float((dlt > 0.1).mean()), spearman_without_vs_base=spearman(base, ab))
    S["channel_range_corr"] = {c + "_" + SC_SHORT[sc]: spearman(cut["base_free_" + SC_SHORT[sc]], cut["rng_" + c]) for c in ("zonal", "tesseral", "third_body", "field", "drag", "speed", "geometry") for sc in SCS}
    S["channel_range_median_ln"] = {c: float(cut["rng_" + c].median()) for c in ("zonal", "tesseral", "third_body", "field", "drag", "speed", "geometry")}
    # --- robustness of the top 10 (P23-2) and of the whole ranking
    S["robust"] = {}
    for m in MODES:
        for sc in SCS:
            s = SC_SHORT[sc]; nom = "base_%s_%s" % (m, s); tie = "jm_%s_%s" % (m, s)
            R_ = {}
            for nm, col in (("harden_plus", "hp_%s_%s" % (m, s)), ("harden_minus", "hm_%s_%s" % (m, s)), ("reserve_L1", "rl1_%s_%s" % (m, s)), ("jitter_mean", "jm_%s_%s" % (m, s))):
                R_[nm] = dict(**top_overlap(cut, nom, col, 10, tie), band_overlap_0p1=band_overlap(cut, nom, col), spearman=spearman(cut[nom], cut[col]), n_feasible=int(cut[col].notna().sum()),
                              median_shift=float((cut[col] - cut[nom]).median()), median_shift_all_feasible=float((cut[col] - cut[nom]).dropna().median()))
            if m == "free":
                for nm, col in (("informed_plus", "ip_free_" + s), ("informed_minus", "im_free_" + s)):
                    R_[nm] = dict(**top_overlap(cut, nom, col, 10, tie), band_overlap_0p1=band_overlap(cut, nom, col), spearman=spearman(cut[nom], cut[col]), n_feasible=int(cut[col].notna().sum()), median_shift=float((cut[col] - cut[nom]).median()))
            S["robust"][m + "_" + s] = R_
    S["sens"] = {}
    for sc in SCS:
        s = SC_SHORT[sc]
        for nm in ("pen", "kg1", "prev"):
            col = "%s_free_%s" % (nm, s); S["sens"]["%s_%s" % (nm, s)] = dict(spearman=spearman(cut["base_free_" + s], cut[col]), median_shift=float((cut[col] - cut["base_free_" + s]).median()), max_abs_shift=float((cut[col] - cut["base_free_" + s]).abs().max()),
                                                                                  n_ok=int(cut[col].notna().sum()))
    # --- ranking tables
    tabs = {}
    for m in MODES:
        for sc in SCS:
            s = SC_SHORT[sc]; col = "base_%s_%s" % (m, s); o = order(cut, col, "jm_%s_%s" % (m, s)); keep = ["id", "inc", "omega0", "raan0", "epoch", col, "jm_%s_%s" % (m, s), "js_%s_%s" % (m, s), "hp_%s_%s" % (m, s), "hm_%s_%s" % (m, s), "rl1_%s_%s" % (m, s), "R_slow", "floor_cal", "days_free_" + s, "comp_flag_" + s] if m == "free" else ["id", "inc", "omega0", "raan0", "epoch", col, "jm_%s_%s" % (m, s), "js_%s_%s" % (m, s), "hp_%s_%s" % (m, s), "hm_%s_%s" % (m, s), "rl1_%s_%s" % (m, s), "R_slow", "floor_cal", "days_free_" + s]
            tabs[m + "_" + s] = dict(top10=o.head(10)[keep].to_dict("records"), bottom10=o[o[col].notna()].tail(10)[keep].to_dict("records"), n_best_ties=int((cut[col] <= cut[col].min() + TIE_MPS).sum()))
    S["tables"] = tabs
    # --- break-even (P23-4) from the largest saving of each ranking
    S["breakeven"] = {}
    for m in MODES:
        for sc in SCS:
            s = SC_SHORT[sc]; col = "base_%s_%s" % (m, s); x = cut[col].dropna(); sv = float(x.max() - x.min()); sv_o = float(ody[col] - x.min()); sv_med = float(ody[col] - x.median())
            S["breakeven"][m + "_" + s] = dict(largest_saving=sv, di_deg=breakeven_di_deg(sv), saving_vs_ody_best=sv_o, di_deg_vs_ody=breakeven_di_deg(max(sv_o, 0.0)), exact_law_increment_at_di=plane_cost_exact(breakeven_di_deg(sv)),
                                               saving_ody_vs_median=sv_med)
    write_json(OUT / "cut256_table.json", cut.to_dict("records")); cut.to_csv(OUT / "cut256_table.csv", index=False)
    S["files"] = {}
    # --- anchors, validation, coupling
    anchors(S, df, cut)
    coupling(S, cut)
    sha = write_json(OUT / "summary.json", S); print("summary sha", sha[:16]); return S


def _anchor_rows():
    p = OUT / "raw_anchors.json"
    return json.load(open(p)) if p.exists() else []


def anchors(S, df, cut):
    rows = _anchor_rows(); A = {}
    if not rows:
        return
    ad = flat(rows)
    for _, r in ad.iterrows():
        d = {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in r.items() if k in ("id", "role", "inc", "omega0", "raan0", "epoch", "base_free_fx", "base_free_fr", "base_tstar_fx", "base_tstar_fr", "ref_free", "hp_free_fx", "hp_free_fr", "hm_free_fx",
                                                                                                          "hm_free_fr", "rl1_free_fx", "rl1_free_fr", "jm_free_fx", "jm_free_fr", "js_free_fx", "js_free_fr", "prev_free_fx", "prev_free_fr", "R_slow", "floor_cal", "feasible", "field_found", "fid", "days_free_fx", "days_free_fr")}
        for sc in SCS:
            d["rank_free_" + SC_SHORT[sc]] = int((cut["base_free_" + SC_SHORT[sc]] < r["base_free_" + SC_SHORT[sc]]).sum()) + 1
        A[r["id"]] = d
    S["anchors"] = A
    write_json(OUT / "anchors.json", A)
    validation(S, ad)


def flown_table():
    """Flown aerobraking dV (m/s, dV to the last counted pass, Revision 17a) per (arm, geometry) cell from FINISHED flight directories (all four csv present and not touched for 15 min). Success = the Revision 13/14 rule of planner_rev20_score."""
    import glob
    import planner_rev20_score as SC
    rows = []
    for blk in ("rev17", "rev20val", "rev21c", "gen1", "devR", "devT", "dev", "held"):
        p = FLIGHT_ROOT / blk
        if not p.exists(): continue
        summ = {}
        for f in glob.glob(str(FLIGHT_ROOT / ("summary_%s_mars_*.csv" % blk))):
            for r in pd.read_csv(f).itertuples(): summ[r.run_id] = r._asdict()
        for d in sorted(p.glob("*")):
            k = SC.parse_id(d.name)
            if not k or not SC.finished(d): continue
            if time.time() - max((d / x).stat().st_mtime for x in SC.FILES) < 900: continue          # not a running flight
            try:
                f = SC.flight(d, summ.get(d.name))
            except Exception:
                continue
            rows.append(dict(block=blk, run=d.name, geom=k["geom"], arm=k["arm"], level=k["level"], seed=k["seed"], ok=bool(f["ok"]), dv=float(f["dv_sum"]), dv_pre=(float(f["dv_pre"]) if np.isfinite(f["dv_pre"]) else 0.0), days=float(f["days"])))
    # held block: no directories under the results root; the corrected dV of the dv_correction table (Revision 17a), success = stop_reason end_radius of the summary row
    fj = json.load(open(TD / "dv_correction" / "flights.json")); sm = {}
    for r in fj:
        if r["phase"] != "held": continue
        if r["summary"] not in sm:
            f_ = FLIGHT_ROOT / r["summary"]; sm[r["summary"]] = pd.read_csv(f_).set_index("run_id") if f_.exists() else None
        if sm[r["summary"]] is None: continue                 # summary file of another host not present here: success unknown, flight skipped
        sr = sm[r["summary"]].loc[r["run_id"]] if r["run_id"] in sm[r["summary"]].index else None
        rows.append(dict(block="held", run=r["run_id"], geom={"G5863": "G1"}.get(r["geometry"], r["geometry"]), arm=r["arm"], level="", seed=0, ok=bool(sr is not None and str(sr.get("stop_reason")) == "end_radius"), dv=float(r["dv_corrected"]), dv_pre=0.0, days=float("nan")))
    return pd.DataFrame(rows)


def validation(S, ad):
    FT = flown_table(); V = {}; pred = {}
    for tok in ("ODY", "G1", "G0975", "G0063", "G0390", "G1216"):
        pred[tok] = {}
        for sc in ("ref", "S-fixed", "S-free"):
            ic = R.select_ic(tok, sc); pred[tok][sc] = dict(aero=ic["campaign_dv_mps"] + ic["dv_pre_mps"], N=ic["N"], days=ic["t_arr_d"], note="stored ground table (G1 = G5863 mean field; ODY mean field; generic: preview, no field), Revision 20 convention")
    V["predicted_ol_convention"] = {tok: dict(aero=R.select_ic(tok, "ref", horizon="full", pen=0.01)["campaign_dv_mps"], burns=R.select_ic(tok, "ref", horizon="full", pen=0.01)["burns"], note="horizon full, burn penalty 0.01 (the held OL arm's convention), reference start")
                                    for tok in ("ODY", "G1")}
    field_pred = {}
    for _, r in ad.iterrows():
        field_pred[r["id"]] = dict(ref=r["ref_free"], S_fixed=r["base_free_fx"], S_free=r["base_free_fr"], tstar_S_fixed=r["base_tstar_fx"], tstar_S_free=r["base_tstar_fr"], kappa="1.0 (ASSUMED)", field_found=bool(r["field_found"]))
    V["predicted_stored_tables"] = pred; V["predicted_field_tables"] = field_pred
    cells = [("OLRT3", "L0", "ODY", "rev17"), ("OLRT3", "L0", "G1", "rev20val"), ("OLRT3", "", "G0975", "gen1"), ("OLRT3", "L0", "G0975", "rev21c"), ("OLRT35", "", "ODY", "devR"), ("OLRT35", "", "G1", "devR"), ("OLRT35", "", "G0975", "gen1"),
             ("OLRT3W", "L0", "ODY", "rev20val"), ("OLRT3W", "L0", "G1", "rev20val"), ("L3", "", "ODY", "dev"), ("L3", "", "G1", "dev"), ("L3", "", "ODY", "held"), ("L3", "", "G1", "held"), ("L1", "", "ODY", "held"), ("L1", "", "G1", "held"),
             ("OLR", "", "ODY", "held"), ("OLR", "", "G1", "held"), ("OL", "", "ODY", "held"), ("OL", "", "G1", "held")]
    out = []
    for arm, lv, geom, blk in cells:
        d = FT[(FT.arm == arm) & (FT.geom == geom) & (FT.block == blk) & ((FT.level == lv) if lv else True)] if len(FT) else FT
        if not len(d): out.append(dict(arm=arm, level=lv, geom=geom, block=blk, n=0)); continue
        ok = d[d.ok]
        out.append(dict(arm=arm, level=lv, geom=geom, block=blk, n=int(len(d)), n_ok=int(len(ok)), dv_median_ok=float(ok.dv.median()) if len(ok) else None, dv_pre_median=float(ok.dv_pre.median()) if len(ok) else None,
                        dv_plus_pre_median_ok=float((ok.dv + ok.dv_pre).median()) if len(ok) else None, dv_median_all=float(d.dv.median()), days_median_ok=float(ok.days.median()) if len(ok) else None))
    V["flown"] = out; V["flown_blocks_available"] = {b: int((FT.block == b).sum()) for b in (FT.block.unique() if len(FT) else [])}
    S["validation"] = V; write_json(OUT / "validation.json", V)


def coupling(S, cut):
    """Item 3: plane-mismatch trade and the share of the geometry freedom that survives, on the 1,536 grid (cut-256 + the rest) for named ASSUMED target planes."""
    rest_p = OUT / "raw_grid1536_rest.json"
    if not rest_p.exists():
        S["coupling"] = dict(note="raw_grid1536_rest.json missing"); return
    rows = json.load(open(OUT / "raw_cut256.json")) + json.load(open(rest_p)); g = fallback(rows, flat(rows)); ody_row = g[g.id == "ODY"].iloc[0]; g = g[g.id != "ODY"].reset_index(drop=True)
    C = {"n_grid": int(len(g))}; INC_T = {"25.2": 25.2, "63.4": 63.4, "ODY 93.522": 93.522, "116.6": 116.6}
    for sc in SCS:
        s = SC_SHORT[sc]
        for m in MODES:
            col = "base_%s_%s" % (m, s); key_ = "%s_%s" % (m, s); aero = g[col]; bA = float(aero.min()); good = aero <= bA + TOL_FREEDOM; res = dict(best_A=bA, n_good_A=int(good.sum()), frac_good_A=float(good.mean()), ody=float(ody_row[col]))
            # case B: target inclination fixed, node free (the target node is unconstrained, so the node can be matched; the cost is the inclination difference)
            B = {}
            for nm, it in INC_T.items():
                di = (g.inc - it).abs(); J = aero + di.apply(plane_cost_quadrature); onc = g.inc.sub(it).abs() < 0.1
                bX = float(J.min()); surv = int(((J <= bA + TOL_FREEDOM) & good).sum())
                B[nm] = dict(n_on_class=int(onc.sum()), frac_grid_on_class=float(onc.mean()), best_in_class=float(aero[onc].min()), median_in_class=float(aero[onc].median()), best_net=bX, best_net_id=str(g.loc[J.idxmin(), "id"]), best_net_inc=float(g.loc[J.idxmin(), "inc"]),
                             best_net_plane_cost=float(plane_cost_quadrature(float(di[J.idxmin()]))), survive_good=surv, survive_frac_of_good=float(surv / max(good.sum(), 1)), n_good_in_class=int((good & onc).sum()),
                             n_net_le_ody=int((J <= ody_row[col]).sum()), frac_net_le_ody=float((J <= ody_row[col]).mean()))
            res["B"] = B
            # case C: inclination AND node fixed: Odyssey's model plane (i 93.522, node at arrival of the Odyssey anchor) as an ASSUMED example; plane mismatch from each candidate's arrival node
            it = 93.522; Ot = float(ody_row["raan_arr"]); di = np.array([rel_inc_deg(i, o, it, Ot) for i, o in zip(g.inc, g.raan_arr)]); J = aero + np.array([plane_cost_quadrature(x) for x in di])
            nodes_fixed = dict(target_inc=it, target_node_arrival_deg=Ot)
            onc = (g.inc.sub(it).abs() < 0.1).values
            res["C_grid"] = dict(**nodes_fixed, best_net=float(np.nanmin(J)), best_net_id=str(g.loc[int(np.nanargmin(J)), "id"]), n_net_le_ody=int((J <= ody_row[col]).sum()), frac_net_le_ody=float((J <= ody_row[col]).mean()),
                                 survive_good=int(((J <= bA + TOL_FREEDOM) & good).sum()), survive_frac_of_good=float(((J <= bA + TOL_FREEDOM) & good).sum() / max(good.sum(), 1)), median_di_on_class_deg=float(np.median(di[onc])) if onc.any() else None,
                                 min_di_on_class_deg=float(di[onc].min()) if onc.any() else None, di_p10_all=float(np.percentile(di, 10)))
            # case C with the node matched by choosing the insertion node: per (inclination class, omega0, epoch) cell the best (min over the sampled Omega0) and the worst aero; required Omega0 = Ot - node drift
            cell = g[onc].copy(); cell["cell"] = cell.omega0.astype(str) + "|" + cell.epoch.astype(str); gg = cell.groupby("cell")[col]
            drift = float(cell.node_drift.median()); need = (Ot - drift) % 360.0
            nearest = cell.assign(d=(((cell.raan0 - need + 180) % 360) - 180).abs()).sort_values("d").groupby("cell").head(2)       # two nearest sampled Omega0 per cell
            res["C_cells"] = dict(n_cells=int(gg.ngroups), omega0_needed_deg=need, node_drift_deg=drift, within_cell_range_median=float((gg.max() - gg.min()).median()), within_cell_range_max=float((gg.max() - gg.min()).max()),
                                  cells_good_best_omega0=int((gg.min() <= bA + TOL_FREEDOM).sum()), cells_good_nearest2=int((nearest.groupby("cell")[col].min() <= bA + TOL_FREEDOM).sum()),
                                  cells_good_all_omega0=int((gg.max() <= bA + TOL_FREEDOM).sum()), cells_le_ody_nearest2=int((nearest.groupby("cell")[col].min() <= ody_row[col]).sum()),
                                  nearest_grid_omega0_error_deg=float(nearest.d.min()), best_per_cell_median=float(gg.min().median()), nearest2_per_cell_median=float(nearest.groupby("cell")[col].min().median()))
            cg = g.groupby(["inc", "omega0", "epoch"])[col]; rg = (cg.max() - cg.min()).dropna()
            res["raan_profile"] = dict(n_cells=int(len(rg)), median_range=float(rg.median()), p90_range=float(np.percentile(rg, 90)), max_range=float(rg.max()), frac_cells_range_le_tol=float((rg <= TOL_FREEDOM).mean()),
                                       by_inc={str(i_): float(rg[rg.index.get_level_values(0) == i_].median()) for i_ in sorted(g.inc.unique())})
            C[key_] = res
    # break-even di for the largest saving of the cut-256 (P23-4) is in S["breakeven"]; plane term table
    C["plane_cost_table"] = {str(d): dict(quadrature=plane_cost_quadrature(d), exact=plane_cost_exact(d)) for d in (0.1, 0.25, 0.5, 1.0, 1.75, 2.0, 3.0, 5.0, 10.0)}
    # node timing: the rate of the node (deg per day) of the cut classes, from the campaign node drift over the campaign days
    C["node_rate_deg_per_day"] = {str(i): float((gi.node_drift / gi.days_model).median()) for i, gi in g.groupby("inc")}
    # Omega0 sensitivity of the DP inside a cell (all 1,536): the surviving-freedom proxy rests on it
    S["coupling"] = C; write_json(OUT / "coupling.json", C)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "run", "summarize"]); ap.add_argument("--jobs", type=int, default=3); ap.add_argument("--part", default="cut", choices=["cut", "anchors", "grid"]); a = ap.parse_args()
    if a.cmd == "selftest":
        selftest()
    elif a.cmd == "run":
        cmd_run(a.part, a.jobs)
    else:
        summarize()
