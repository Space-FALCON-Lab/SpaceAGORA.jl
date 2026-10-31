"""Revision 20, WP2 (analysis only, Python, no Julia, no flights): closed-loop counterfactual replay of the OLRT3 planner, with and without the Revision 18 time-to-go guard, on the realized
natural (burn-free) paths of the existing flights.

COUNTERFACTUAL ASSUMPTION (stated, ASSUMED, tested in `validate`). The GRAM noise a flight meets is independent of the policy that flies it. Per flight and per pass p the extraction keeps
  xnat_p  natural ln q = flown ln q_p minus the cumulative executed offset nu_p (replan.csv nu_cum);   residual rx_p = xnat_p - [xnat_1 + table path at the flown state];
  e_p     one-pass prediction error of the planner (flown ln q_p - predicted post-burn level, as part2_derivations_check.flight);
  eps_p   drop noise = flown apoapsis drop / (kappa d0(ra_{p-1}) q_p / q_req)  (kappa: Revision 19 B, MEASURED per geometry);
  s_p     ln(rho_peak / rho) so that Q_s,p = exp(x_p + s_p) (U15, V_rel basis: the flight's own pass-maximum ratio);
  rr_p    periapsis residual about the table at the flown state; rho_p flown/model period ratio (about 1.000).
Two indexings of the natural level are implemented: noise='pass' (xnat_p used as flown at pass index p, whatever the state) and noise='state' (the table level at the CURRENT replayed state, k(ra_{p-1}),
plus the flown residual). The choice is tested in `validate` (own path: exact by construction; cross-arm: arm A replayed on arm B's natural path of the same seed against the flown arm A). Beyond the
last flown pass the residual series are block-bootstrapped (ASSUMED, fixed seed) from the last flown block; `ext` counts such passes.
State dependent world (Revision 19 B): drop_p = eps_p kappa d0(ra_{p-1}) exp(x_p) / q_req, rp_p = rpnat_p - H_p u_p, period from Kepler on the osculating orbit, anchored by rho_p.

Policy (mirror of planner_common.jl planner_decide for OLRT3): a priori schedule (ol_schedule) -> per-apoapsis decision with the terminal aim (planner_terminal_aim, c from the last 6 drop/q), else the
tube DP (plan_inputs + exact DP) with the drift update bias, the replan-deviation threshold dev_thr_mps = 0.15 (snap to the a priori burn unless outside the full tube). GUARD (Revision 18 design (b)):
at each apoapsis outside the terminal aim, predict the closing-apoapsis time from the current state, the measured drop/q ratio and the DP path (level forecast anchored on the mean residual of the last
anchor_k passes, drift bias persisted bias_steps steps); if it is outside [t* - s, t* + s] re-solve the DP with a larger "price" and bisect it, at most 12 solves in all, until the prediction is inside
[t* - aim_frac s, t* + aim_frac s] or the tube binds (flag time_infeasible, best-effort plan at the largest price). Actuator: "price" = Revision 18's linear term -lambda wT_k nu_k (bang-bang, chatters);
"tube" (the tuned headline, a DEVIATION from the Revision 18 text) = the same bisection variable tightens the tube edge on the needed side (never raising the upper edge). The tube, the terminal aim
and Q_s 0.65 are never relaxed. All guard parameters are ASSUMED and named in GUARD_DEFAULT / GUARD_TUNED (chosen on DEV paths, seeds 1-8, only); GUARD_ROUND2 is a post hoc aim-point shift.

Usage: python3 -I planner_rev20_replay.py selftest | noisecorr | validate | evaluate --set dev|held [--scenarios ref,S-fixed,S-free] [--noise pass|state|time] [--round 1|2] [--tag T] | curve | rev17 | gen1
(data to tasks/.../rev20_data/replay/; dev tuning: tune_stage() with a list of guard-parameter dicts)
"""
import argparse, glob, hashlib, json, math, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev18_offline as O
import planner_rev18_flights as FL
import planner_analysis as AN

INF = math.inf
MU = O.MU
RA_T = AN.RA_TARGET                                    # SOURCED 7,408 km
CFG = PLB.load_cfg()
TD = pathlib.Path(PLB.task_dir())
ROOT = pathlib.Path(__import__("os").environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
OUT = TD / "rev20_data" / "replay"
Q_FIRST = O.Q_FIRST                                    # MEASURED calibration target 0.212
QFIRST_CAL = {"ODY": 0.2119115774128761, "G1": 0.2119998026887692, "G0975": 0.211971}      # SOURCED: q_first_w_cm2 of summary_cal_mars_ws_cal.csv / summary_calG_mars_ws_calG.csv (the planner's w.q_first)
GUARD_W_RP_M = 3.48e6                                  # ASSUMED: periapsis radius of the nominal closing orbit, used only for the guard's slack basis W_max
QS_LIMIT = CFG["planner"]["q_never_exceed_qs_w_cm2"]   # SOURCED 0.65
TOL_KM = CFG["planner"]["terminal"]["tolerance_km"]    # ASSUMED 25 km
HORIZON_D = CFG["planner"]["horizon_days"]             # ASSUMED 130 d
W_MAX_PERIODS = 3.0                                    # ASSUMED (Revision 18)
H_DES_M = 250e3                                        # ASSUMED (Revision 18)
TIME_SERIES = ("e", "eps", "s")          # noise series that are indexed by time when noise="time"
EXT_BLOCK = 10                                         # ASSUMED bootstrap block length (passes) for noise beyond the flown path
EXT_FROM = 30                                          # ASSUMED: blocks drawn from the last 30 flown passes
P_MAX_EXTRA = 220                                      # NUMERICAL cap on replayed passes beyond the flown path
GUARD_DEFAULT = dict(s_frac=0.5,       # ASSUMED: guard slack s = s_frac x W_max + rel_slack x (time to go), W_max = 3 final-orbit periods of the nominal closing orbit
                     rel_slack=0.02,   # ASSUMED: the slack grows with the predicted time to go (the one-pass-ahead prediction of a far arrival is uncertain; a tight band chatters, see results note)
                     edge_margin_ln=0.0,  # ASSUMED: guard solves use the tube shrunk by this margin at both edges (0 = the planner's own tube); tube edges are never relaxed
                     aim_frac=0.5,     # ASSUMED: the bisection stops when the prediction is within t* +/- aim_frac x s
                     lam_max=1.0,      # ASSUMED m/s per orbit (the Revision 18 grid ends at 0.5; 1.0 is past where the tube binds)
                     lam_min=1.0e-4,   # ASSUMED lower end of the log bisection (Revision 18 grid starts at 3.2e-4)
                     max_solves=12,    # SOURCED Revision 18 design (b)
                     k_window=6,       # ASSUMED: median of the last 6 drop ratios (= term_c_window)
                     bias_steps=15,    # ASSUMED: the observed drift bias is persisted this many steps in the guard's forecast (15 = the planner's own cap)
                     infeasible_mode="pin",  # ASSUMED: when the tube binds ("time-infeasible") fly the largest-price plan ("pin") or the planner's own plan ("plain")
                     anchor_k=15,      # ASSUMED: the guard anchors its level forecast on the mean natural-level residual of the last anchor_k passes (0 = the planner's one-pass anchor z0 - fast_j)
                     actuator="price", # ASSUMED: "price" = Revision 18's linear price lambda in the DP; "tube" = tighten the tube edge on the needed side (never raising the upper edge), a lazy actuator
                     eps_max=1.0,      # ASSUMED: tube actuator: fraction of the tube width the edge may move at the largest price
                     t_offset_d=0.0,   # ASSUMED (round 2 only): the guard aims at t* + offset; value = minus the median arrival error of the round-1 guard on the DEV paths, per geometry (float or dict by geometry)
                     term_model=False, # ASSUMED switch: the guard's arrival prediction emulates the terminal aim (equal drops of at most term_step_m inside term_activate_steps steps of the target)
                     kappa_stat="median",  # ASSUMED: statistic of the measured drop ratio window (the planner's terminal aim uses the median; the mean is unbiased for the summed drop)
                     ema_beta=0.0,     # ASSUMED: exponential smoothing factor of the predicted arrival over decisions (0 = off)
                     min_burn_g=0.0,   # ASSUMED (v2): guard burns smaller than this (m/s) are replaced by the planner's own burn
                     rev_dead_g=0.0,   # ASSUMED (v2): reversal deadband (m/s): a guard burn smaller than this and opposite to the previous guard burn within rev_win_g passes is not flown
                     rev_win_g=6,      # ASSUMED (v2): passes of the reversal window
                     dev_thr_g=0.0,    # ASSUMED: the guard's burn replaces the planner's own burn only when it differs by more than this (m/s); 0 = always
                     pen_g=None,       # ASSUMED: fixed cost per burn (m/s) in the guard's DP solves (None = the planner's 0.01)
                     start_pass=2)     # ASSUMED: guard active from the first decision


# ------------------------------------------------------------------------------------------------ tables and calibration

def load_tab(token):
    if token in ("ODY", "G1"):
        tab = O.load_tab(token)
    elif token == "G0975":
        g = json.load(open(TD / "planner_ground_tables_generic.json"))["geometries"]["G0975"]
        tab = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in g.items()}
    else:
        raise ValueError(token)
    tab["token"] = token
    n = tab["n"]; ra = tab["ra"]
    tab["d0"] = np.concatenate((np.diff(-ra), [float(ra[-2] - ra[-1])]))
    tab["_ra_inc"] = ra[::-1]; tab["_d0_inc"] = tab["d0"][::-1]; tab["_k_inc"] = np.arange(n, 0, -1, dtype=float)
    tab["_kk"] = np.arange(1, n + 1, dtype=float)
    return tab


def load_kappa():
    """MEASURED (rev19_data/rev18_calibrate.json, Revision 19 B): drop-per-q calibration per geometry; G0975 ASSUMED 1 (no flight calibration of the generic table)."""
    d = json.load(open(TD / "rev19_data" / "rev18_calibrate.json"))["geometry"]
    k = {g["geom"]: float(g["kappa_med"]) for g in d}
    k.setdefault("G0975", 1.0)
    return k


KAPPA = load_kappa()
_TABS = {}


def get_tab(token):
    if token not in _TABS:
        _TABS[token] = load_tab(token)
    return _TABS[token]


def d0_at(tab, ra):
    return float(np.interp(ra, tab["_ra_inc"], tab["_d0_inc"]))


def kidx(tab, ra):
    """continuous table pass index (1..n) of an apoapsis radius"""
    return float(np.interp(ra, tab["_ra_inc"], tab["_k_inc"]))


def tfun(tab, key, k):
    return float(np.interp(k, tab["_kk"], tab[key]))


def Tk(ra, rp):
    a = 0.5 * (ra + rp)
    return 2.0 * math.pi * math.sqrt(a ** 3 / MU)


# ------------------------------------------------------------------------------------------------ fast exact DP (vectorized; identical objective to planner_lib.dp_exact, plus a linear term)

def _run(c):
    m = np.minimum.accumulate(c)
    idx = np.arange(len(c))
    ia = np.maximum.accumulate(np.where(c <= m, idx, 0))
    return m, ia


def dpf(A, B, L, pen=0.0, lin=None, tol=1.0e-9):
    """min sum_k L_k |nu_k - nu_{k-1}| + pen 1(move) + lin_k nu_k, nu_{-1} = 0, A_k <= nu_k <= B_k, on the vertex set {0} U {A} U {B}. Returns nu or None. Same relaxation / penalty / mask order as
    planner_lib.dp_exact and planner_rev18_offline.dp_lin (checked in selftest)."""
    A = np.asarray(A, float); B = np.asarray(B, float); L = np.asarray(L, float); n = len(A)
    val = np.unique(np.concatenate(([0.0], A, B))); ns = len(val); i0 = int(np.argmin(np.abs(val)))
    prev = np.full(ns, INF); prev[i0] = 0.0
    parent = np.empty((n, ns), dtype=np.int32); ar = np.arange(ns)
    for k in range(n):
        Lk = L[k]
        cm, ia = _run(prev - Lk * val)
        fwd = Lk * val + cm
        cm2r, ia2r = _run((prev + Lk * val)[::-1])
        bwd = (-Lk * val + cm2r[::-1]); ab = (ns - 1 - ia2r)[::-1]
        useb = bwd < fwd
        cur = np.where(useb, bwd, fwd); arg = np.where(useb, ab, ia)
        keep = prev <= cur
        cur = np.where(keep, prev, cur); arg = np.where(keep, ar, arg)
        if pen > 0.0:
            stay = prev <= pen + cur
            cur = np.where(stay, prev, pen + cur); arg = np.where(stay, ar, arg)
        if lin is not None:
            cur = cur + lin[k] * val
        cur = np.where((val < A[k] - tol) | (val > B[k] + tol), INF, cur)
        parent[k] = arg; prev = cur
    j = int(np.argmin(prev))
    if not np.isfinite(prev[j]):
        return None
    nu = np.zeros(n)
    for k in range(n - 1, -1, -1):
        nu[k] = val[j]; j = parent[k][j]
    return nu


# ------------------------------------------------------------------------------------------------ flight extraction (natural path and noise series)

def run_dir(run_id):
    for ph in ("dev", "devT", "devR", "gen1", "sens", "rev17", "runs", "rev20val", "rev21c", "rev21v", "rev20olt", "../trx_pull4/planner/devR"):
        c = ROOT / ph / run_id
        if (c / "passes.csv").exists():
            return c
    return None


def find_flight(run_id):
    d = run_dir(run_id)
    return d


def extract(run_id, tab, kappa, Hloc=None):
    """Natural path and noise series of one flown flight (all arrays indexed 1..n at [p-1]). Returns dict."""
    d = run_dir(run_id)
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True)
    R = pd.read_csv(d / "replan.csv")
    ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= RA_T))[0]
    reached = len(hit) > 0
    n = int(hit[0] + 1) if reached else len(P)
    q = P.q_w_cm2.values.astype(float)[:n]; x = np.log(q); rp = P.rp_m.values.astype(float)[:n]; ra = ra[:n]; t = P.time_s.values.astype(float)[:n]
    rowp = {int(p): i for i, p in enumerate(R["pass"].values)}
    nu = np.zeros(n); cur = 0.0; e = np.zeros(n); dvx = np.zeros(n); w0 = np.full(n, np.nan); z0p = np.full(n, np.nan)
    for p in range(2, n + 1):
        i = rowp.get(p)
        if i is not None:
            cur = float(R.nu_cum.values[i]); dv = float(R.executed_dv_mps.values[i]); w = float(R.w0.values[i]); z = float(R.z0_pred_lnq.values[i])
            nu[p - 1] = cur; dvx[p - 1] = dv; w0[p - 1] = w; z0p[p - 1] = z
            e[p - 1] = x[p - 1] + dv / w - z
        else:
            nu[p - 1] = cur
    e_sur = False
    if not np.all(np.isfinite(e)):                                         # open-loop (OL) flights never predict: surrogate prediction errors N(0, sigma_p) (sigma_p SOURCED planner_config, ASSUMED to transfer), fixed seed
        rng_e = np.random.default_rng(int(hashlib.sha256(run_id.encode()).hexdigest()[:8], 16)); bad = ~np.isfinite(e); e = np.where(bad, rng_e.normal(0.0, CFG["planner"]["sigma_p"], n), e); e[0] = 0.0; e_sur = True
    xnat = x - nu
    H_pass = P.h_local_m.values.astype(float)[:n]
    rpnat = rp + H_pass * nu
    ra_prev = np.concatenate(([float(tab["ra"][0])], ra[:-1]))
    drop = ra_prev - ra
    d0 = np.array([d0_at(tab, r) for r in ra_prev])
    model = kappa * d0 * np.exp(x) / tab["q_req"]
    eps = np.where(np.isfinite(drop) & (model > 0), drop / model, 1.0)
    rho = P.rho.values.astype(float)[:n]; rpk = P.rho_peak_kg_m3.values.astype(float)[:n]
    s = np.where((rho > 0) & (rpk > 0), np.log(rpk / np.where(rho > 0, rho, 1.0)), 0.0)
    # state residuals
    kb = np.array([kidx(tab, r) for r in ra_prev])
    qb = xnat[0] + np.array([tfun(tab, "path", k) for k in kb]) - float(tab["path"][0])
    rx = xnat - qb
    rpb = np.array([tfun(tab, "rp", k) for k in kb])
    rr = rpnat - rpb
    # periods: model T over the orbit leaving pass p (half before/after the burn at the next apoapsis)
    rho_T = np.ones(n)
    rpn = np.concatenate((rp[1:], [rp[-1]]))
    for p in range(1, n):
        Tm = 0.5 * Tk(ra[p - 1], rp[p - 1]) + 0.5 * Tk(ra[p - 1], rpn[p - 1])
        rho_T[p - 1] = (t[p] - t[p - 1]) / Tm
    a_m = float(P.a_m.values[n - 1]); t_arr_ref = t[n - 1] + 0.5 * 2 * math.pi * math.sqrt(a_m ** 3 / MU)      # Revision 18 definition (planner_rev18_flights.flight)
    rho_T[n - 1] = (t_arr_ref - t[n - 1]) / (0.5 * Tk(ra[n - 1], rp[n - 1]))
    if Hloc is None:
        lev = np.array([PLB.lever_m_per_mps(MU, a, ee) for a, ee in zip(R.a_m.values, R.e.values)])
        Hloc = float(np.nanmedian(R.w0.values * lev)) if len(R) else float("nan")
    dra = np.zeros(n)                                                  # osculating apoapsis at the decision (a (1+e), replan.csv) minus the post-pass apoapsis of passes.csv: tens of metres, anchors the terminal aim
    lw = np.ones(n)                                                    # flown / model lever (osculating a, e at the apoapsis against the table-free post-pass orbit): anchors w0
    for p in range(2, n + 1):
        i = rowp.get(p)
        if i is not None:
            am = 0.5 * (ra[p - 2] + rp[p - 2]); em = (ra[p - 2] - rp[p - 2]) / (ra[p - 2] + rp[p - 2])
            dra[p - 1] = float(R.a_m.values[i]) * (1 + float(R.e.values[i])) - ra[p - 2]
            lw[p - 1] = PLB.lever_m_per_mps(MU, float(R.a_m.values[i]), float(R.e.values[i])) / PLB.lever_m_per_mps(MU, am, em)
    xslow = np.convolve(np.pad(xnat, (7, 7), mode="edge"), np.ones(15) / 15, mode="valid")[:n]
    qs = np.exp(x + s)
    return dict(e_surrogate=e_sur, xslow=xslow, hp=xnat - xslow, lw=lw, dra=dra, run_id=run_id, dir=str(d), n=n, reached=reached, x=x, xnat=xnat, nu=nu, rp=rp, ra=ra, t=t, e=e, eps=eps, s=s, rx=rx, rr=rr, rpnat=rpnat, H=H_pass, rhoT=rho_T,
                dvx=dvx, w0=w0, Hloc=Hloc, t_arr_d=t_arr_ref / 86400.0, alt=P.alt_m.values.astype(float)[:n], qs=qs, ra_prev=ra_prev, mass0=float(P.mass_kg.values[0]),
                dv_flown=float(np.abs(dvx[1:n]).sum()), final_err_m=float(ra[n - 1] - RA_T), qs_max=float(np.max(qs)), q1=float(q[0]))


def extend(F, pmax, seed):
    """Fill the residual series beyond the flown path by a block bootstrap of the last EXT_FROM flown passes (ASSUMED). Returns F with arrays of length pmax and F['ext_from']=n."""
    n = F["n"]
    if pmax <= n:
        return F
    rng = np.random.default_rng(seed); lo = max(2, n - EXT_FROM); src = np.arange(lo, n)           # 0-based indices of flown passes usable as noise (not pass 1)
    need = pmax - n; idx = []
    while len(idx) < need:
        b = int(rng.integers(0, len(src))); idx.extend(src[b:b + EXT_BLOCK].tolist())
    idx = np.array(idx[:need])
    G = dict(F)
    for k in ("e", "eps", "s", "rx", "rr", "rhoT", "H", "lw", "dra", "hp"):
        G[k] = np.concatenate((F[k], F[k][idx]))
    for k in ("xnat", "rpnat"):
        G[k] = np.concatenate((F[k], np.full(need, np.nan)))
    G["ext_n"] = need
    return G


# ------------------------------------------------------------------------------------------------ the arm: switches and the a priori schedule

_SCHED = {}


def arm_def(arm):
    a = CFG["planner"]["arms"][arm]
    sw = PLB.arm_switches(CFG, arm)
    sw.update(terminal=bool(a.get("terminal", False)), dev_thr=float(a.get("dev_thr_mps", 0.0)), lower_margin=float(a.get("lower_margin_ln", 0.0)),
              rev_window=int(a.get("rev_window_orbits", 0)), rev_deadband=float(a.get("rev_deadband_ln", 0.0)), min_burn_arm=float(a.get("min_burn_arm_mps", 0.0)), smooth_beta=float(a.get("smooth_beta", 0.0)))
    return sw, PLB.flat_P(CFG, arm)


def apriori(arm, token, q_first_ln):
    key = (arm, token, round(q_first_ln, 12))
    if key not in _SCHED:
        sw, P = arm_def(arm)
        _SCHED[key] = PLB.ol_schedule(sw, get_tab(token), P, q_first_ln, P["dp_tie_eps"])[0]
    return _SCHED[key]


def bias_of(P, tab, j, xs, nuh):
    """planner_bias: slope of the residual of the realized ln q against table path + executed offset over the last update_window flown passes."""
    K = int(P["update_window"]); n = j - 1
    if n < K + 1:
        return 0.0
    idx = np.arange(n - K + 1, n + 1)
    p1 = tab["path"][0]
    r = np.array([xs[i - 1] - (xs[0] + PLB.tget(tab["path"], i) - p1 + nuh[i - 1]) for i in idx])
    xm = idx.mean(); den = ((idx - xm) ** 2).sum()
    return float(((idx - xm) * (r - r.mean())).sum() / den) if den > 0 else 0.0


# ------------------------------------------------------------------------------------------------ arrival prediction used by the guard

def predict_arrival(tab, kap, ra0, rp0, t0, x_plan, rp_plan, term=None):
    """Closing-apoapsis time (s) from the state after pass j-1 (ra0, rp0, periapsis time t0) for planned ln q x_plan[m] and periapsis rp_plan[m] of passes j+m (held at the last value beyond the plan).
    Returns (t_arr, N_steps, drops list) or (nan, N, drops) when not reached within 400 steps."""
    ra = ra0; rp_prev = rp0; t = t0; drops = []; nx = len(x_plan)
    for m in range(400):
        xm = x_plan[min(m, nx - 1)]; rpm = rp_plan[min(m, nx - 1)]
        t += 0.5 * Tk(ra, rp_prev) + 0.5 * Tk(ra, rpm)
        dr = kap * d0_at(tab, ra) * math.exp(xm) / tab["q_req"]
        if term is not None:                                                    # the terminal aim (OLRT3, Revision 14) takes over inside activate_steps x step: N equal drops d* = (R + aim) / N
            Reff = ra - RA_T + term["term_aim_below_m"]
            if 0 < Reff <= term["term_activate_steps"] * term["term_step_m"]:
                dr = Reff / max(1, math.ceil(Reff / term["term_step_m"]))
        drops.append(dr); ra = ra - dr
        if ra <= RA_T:
            return t + 0.5 * Tk(ra, rpm), m + 1, drops
        rp_prev = rpm
    return float("nan"), 400, drops


# ------------------------------------------------------------------------------------------------ the replay

def nominal_arrival_d(token, x0_ln=None, arm="OLR"):
    """Plan's nominal closing-apoapsis time (model world, calibrated): the open-loop plan from the first-pass level x0 on the ground table, Revision 19 B arrival with kappa."""
    x0_ln = math.log(QFIRST_CAL[token]) if x0_ln is None else x0_ln
    tab = get_tab(token); sw, P = arm_def(arm)
    inp = PLB.ol_inputs(sw, tab, P, x0_ln)
    pl = PLB.dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"])
    x = np.concatenate(([x0_ln], np.asarray(inp["level"]) + pl["nu"]))
    return O.arrival(tab, x, x0_ln, kappa=KAPPA[token])["days"]


def ic_shift(ic):
    """Shift u0 (ln q) of the first-pass level of an IC against the reference start: WP1 reports x0 against its own reference ln 0.212 (planner_rev20_ic.X_REF), the planner's calibrated q_first
    differs from it by 4e-4 ln (QFIRST_CAL), so the IC enters as a shift."""
    return 0.0 if ic is None else float(ic["x0"]) - math.log(Q_FIRST)


def replay(F, token, arm="OLRT3", ic=None, noise="pass", guard=None, tstar_d=None, h_des_m=H_DES_M, seed=0, trace=False, qs_cap=QS_LIMIT, term_below_m=None, eps_pool=None, eps_seed=0):
    """Closed-loop counterfactual replay of `arm` on the natural path of flight F (output of extract). guard=None: no time guard; guard=dict(GUARD_DEFAULT overrides): the guard with target tstar_d (days).
    ic=None: reference start (first pass at the flown level, x_ref = ln 0.212 in the plan); ic=dict(x0=..., nu_pre=..., dv_pre_mps=...) shifts the whole natural level by u0 = x0 - ln q_first_ref
    (the first-pass level becomes x0; dv_pre_mps is added to the total, not to the campaign dV). Returns the score dict."""
    tab = get_tab(token); kap = KAPPA[token]; sw, P = arm_def(arm); gd = None if guard is None else {**GUARD_DEFAULT, **guard}
    if term_below_m is not None:                                              # Revision 28B: the terminal-aim margin as a named parameter (ASSUMED 12.5 km in the config)
        P = dict(P, term_aim_below_m=float(term_below_m))
    n0 = F["n"]; pmax = n0 + P_MAX_EXTRA + 5
    F = extend(F, pmax, seed)
    x_ref = math.log(QFIRST_CAL[token]); u0 = ic_shift(ic); x0 = x_ref + u0
    dv_pre = 0.0 if ic is None else float(ic.get("dv_pre_mps", 0.0))
    sched = apriori(arm, token, x0)
    hi0 = math.log(P["q_cap_mean_w_cm2"]); W = tab["W"] if sw["slack"] else P["l0_width_ln"]; lo0f = hi0 - W
    Hloc = F["Hloc"]; ratio = P["burn_response_ratio"]
    path0 = float(tab["path"][0]); xn1 = F["xnat"][0]
    # state arrays (index p-1)
    xs = []; rps = []; ras = []; ts = []; nuh = []; qss = []; dvs = {}; c_hist = []; last_ra = float("nan")
    u = u0                                   # cumulative offset incl. the initial shift (true level = xnat + u)
    nu_plan = 0.0                            # the planner's own bookkeeping (executed burns after pass 1)
    last_g = (0.0, -10 ** 6); ra_prev = float(tab["ra"][0]); info = dict(engaged=0, infeasible=0, solves=0, decisions=0); lam_log = []
    reach = False; viol = False; ext_used = 0; tr = []
    if gd is not None:
        off = gd["t_offset_d"].get(token, 0.0) if isinstance(gd["t_offset_d"], dict) else float(gd["t_offset_d"])
        tstar_s = (tstar_d + off) * 86400.0                                       # aim point of the guard (round 2: shifted by the dev-measured bias); scoring always uses tstar_d
    for p in range(1, pmax + 1):
        i = p - 1
        # noise series of this pass: by pass index, or (noise="time") interpolated in time from the source flight at the burn-free time of this pass
        t_est = (ts[-1] + F["rhoT"][i - 1] * Tk(ras[-1], rps[-1])) if (p >= 2 and ts) else None
        inside = noise == "time" and t_est is not None and t_est <= F["t"][-1]
        Sd = {k: F[k][i] for k in TIME_SERIES}
        if inside:
            Sd = {k: float(np.interp(t_est, F["t"], F[k][:n0])) for k in TIME_SERIES}
        # natural level and periapsis
        if noise == "state" and p > 1:
            kb = kidx(tab, ra_prev)
            xnat_p = xn1 + tfun(tab, "path", kb) - path0 + F["rx"][i]
            rpnat_p = tfun(tab, "rp", kb) + F["rr"][i]
        elif inside:
            xnat_p = float(np.interp(t_est, F["t"], F["xnat"][:n0])); rpnat_p = float(np.interp(t_est, F["t"], F["rpnat"][:n0]))
        elif p > n0 or (noise == "time" and t_est is not None):
            # beyond the flown path (or beyond the source's time span): last smoothed level + the table's trend between the source end and now (by table time) + a bootstrapped zero-mean residual
            te = t_est if t_est is not None else (ts[-1] + F["rhoT"][i - 1] * Tk(ras[-1], rps[-1]))
            trend = lambda key, t: float(np.interp(t, tab["t"], tab[key]))
            xnat_p = F["xslow"][n0 - 1] + trend("path", te) - trend("path", F["t"][-1]) + F["hp"][i]
            rpnat_p = F["rpnat"][n0 - 1] + trend("rp", te) - trend("rp", F["t"][-1])
        else:
            xnat_p = F["xnat"][i]; rpnat_p = F["rpnat"][i]
        if p > n0:
            ext_used += 1
        dv = 0.0; w0 = float("nan")
        if p >= 2:
            ra_now = ras[-1]; rp_now = rps[-1]; a_now = 0.5 * (ra_now + rp_now); e_now = (ra_now - rp_now) / (ra_now + rp_now)
            w0 = Hloc / (ratio * PLB.lever_m_per_mps(MU, a_now, e_now) * F["lw"][i])
            z0 = xnat_p + u - Sd["e"]                                               # one-pass prediction before this pass (no further burn)
            info["decisions"] += 1
            # c_hist: observed drop / measured q of the previous pass (planner pushes from the second decision on)
            ra_osc = ra_now + F["dra"][i]
            if math.isfinite(last_ra):
                dprev = last_ra - ra_osc; qprev = math.exp(xs[-1])
                if dprev > 0 and qprev > 0:
                    c_hist.append(dprev / qprev)
            last_ra = ra_osc
            ta = PLB.terminal_aim(ra_now + F["dra"][i] - RA_T, c_hist, P) if sw["terminal"] else None
            if ta is not None:
                dv = -(ta["lnq"] - z0) * w0
                if abs(dv) < P["dp_min_burn_mps"]:
                    dv = 0.0
            else:
                outside = (z0 > hi0) or (z0 < lo0f)
                bias = bias_of(P, tab, p, xs, nuh) if (sw["mode"] == "olr" and sw["update"]) else 0.0
                inp = PLB.plan_inputs(sw, tab, P, p, z0, w0, bias, sw["lower_margin"])
                lvl = np.asarray(inp["level"], float); A_ = np.asarray(inp["lo"], float) - lvl; B_ = np.asarray(inp["hi"], float) - lvl; L_ = PLB.tie_weights(inp["w"], P["dp_tie_eps"])
                nu = dpf(A_, B_, L_, P["dp_burn_penalty_mps"])
                engaged = False
                if nu is not None and gd is not None and p >= gd["start_pass"]:
                    lvl_alt = None
                    if gd["anchor_k"] > 0 and p >= 6:
                        Kk = min(int(gd["anchor_k"]), p - 2); ii = np.arange(p - Kk, p)
                        cbar = float(np.mean([xs[q_ - 1] - (u0 + nuh[q_ - 1]) - PLB.tget(tab["path"], q_) for q_ in ii]))
                        lvl_alt = lvl.copy(); lvl_alt[1:] = u + cbar + np.array([PLB.tget(tab["slow"], p + m) for m in range(1, len(lvl))])
                    nu_ref = nu
                    tstar_g = tstar_s
                    if gd.get("aim_fn") is not None:                                       # Revision 26c node check v2 (planner_rev26_position): the aim time, clamped to the window
                        tstar_g = gd["aim_fn"](dict(tab=tab, kap=kap, gd=gd, info=info, j=p, lvl=lvl, P=P, nu_ref=nu, ra_now=ra_now, rp_now=rp_now, ts=list(ts), ras=list(ras), rps=list(rps), tstar_s=tstar_s, kap_meas=c_hist_ratio(xs, ras, tab, ra_prev0=float(tab["ra"][0]), k=gd["k_window"], stat=gd["kappa_stat"]), bias=bias, H=F["H"][i], lvl_alt=lvl_alt))
                        info["aim_shift_s"] = tstar_g - tstar_s
                    nu, engaged = guard_solve(tab, kap, gd, info, lam_log, p, lvl, A_, B_, L_, P, nu, ra_now, rp_now, ts[-1], tstar_g, c_hist_ratio(xs, ras, tab, ra_prev0=float(tab["ra"][0]), k=gd["k_window"], stat=gd["kappa_stat"]), bias, F["H"][i], lvl_alt)
                    if gd.get("pos_fn") is not None:                                       # Revision 26 position stage (planner_rev26_position): second objective on the accumulated apse rotation
                        nu, engaged = gd["pos_fn"](dict(tab=tab, kap=kap, gd=gd, info=info, lam_log=lam_log, j=p, lvl=lvl, A_=A_, B_=B_, L_=L_, P=P, nu_ref=nu_ref, nu_t=nu, eng_t=engaged, ra_now=ra_now, rp_now=rp_now, ts=list(ts), ras=list(ras), rps=list(rps),
                                                       tstar_s=tstar_g, kap_meas=c_hist_ratio(xs, ras, tab, ra_prev0=float(tab["ra"][0]), k=gd["k_window"], stat=gd["kappa_stat"]), bias=bias, H=F["H"][i], lvl_alt=lvl_alt))
                def finish(nu_, eng_):
                    d_ = -float(nu_[0]) * w0
                    if sw["dev_thr"] > 0 and sw["mode"] == "olr" and not outside and not eng_:
                        ap = sched.get(p, 0.0)
                        if abs(d_ - ap) <= sw["dev_thr"]:
                            d_ = ap
                    return 0.0 if abs(d_) < max(P["dp_min_burn_mps"], sw["min_burn_arm"]) else d_
                if nu is None:
                    dv = 0.0
                else:
                    dv = finish(nu, engaged)
                    if engaged and (gd["dev_thr_g"] > 0 or gd["min_burn_g"] > 0 or gd["rev_dead_g"] > 0):
                        dref = finish(nu_ref, False)
                        if gd["dev_thr_g"] > 0 and abs(dv - dref) <= gd["dev_thr_g"]:         # follow the planner's own (unguarded) burn unless the guard differs by more than the threshold
                            dv = dref
                        elif dv != dref and 0 < abs(dv) < gd["min_burn_g"]:                    # minimum guard burn: a smaller guard burn is replaced by the planner's own burn
                            dv = dref
                        elif (dv != dref and dv != 0 and gd["rev_dead_g"] > 0 and last_g[0] != 0 and np.sign(dv) != last_g[0] and p - last_g[1] <= gd["rev_win_g"] and abs(dv) < gd["rev_dead_g"]):
                            dv = dref                                                          # reversal deadband: a small guard burn against the previous guard burn within rev_win_g passes is not flown
                        if dv != dref and dv != 0:
                            last_g = (float(np.sign(dv)), p)
            if dv != 0.0:
                u += -dv / w0; nu_plan += -dv / w0; dvs[p] = dv
        # realized pass
        xp = xnat_p + u
        rpp = rpnat_p - F["H"][i] * u
        if p == 1:
            rpp = rpnat_p - F["H"][0] * u0
        epsv = Sd["eps"]
        if eps_pool is not None and ra_prev - RA_T <= 2.0e5:                    # Revision 28B: drop noise of the last passes drawn from a pool (a deterministic function of the pass index and the draw, common to all margins)
            epsv = float(eps_pool[(eps_seed * 7919 + p * 104729) % len(eps_pool)])
        kd = epsv * kap * d0_at(tab, ra_prev) * math.exp(xp) / tab["q_req"]
        ra_p = ra_prev - kd
        if p == 1:
            t_p = F["t"][0]
        else:
            t_p = ts[-1] + F["rhoT"][i - 1] * (0.5 * Tk(ras[-1], rps[-1]) + 0.5 * Tk(ras[-1], rpp))
        qs_p = math.exp(xp + Sd["s"])
        if qs_p > qs_cap:
            viol = True
        xs.append(xp); rps.append(rpp); ras.append(ra_p); ts.append(t_p); nuh.append(nu_plan); qss.append(qs_p)
        if trace:
            tr.append((p, xp, dv, ra_p, t_p, rpp, Sd["eps"]))
        ra_prev = ra_p
        if ra_p <= RA_T:
            reach = True; nfin = p
            break
        if t_p / 86400.0 > HORIZON_D:
            break
    out = dict(run_id=F["run_id"], arm=arm, reached=reach, hard_violation=bool(viol), qs_max=float(max(qss)), ext_used=int(ext_used), info=info)
    if reach:
        ra_f = ras[-1]; rp_f = rps[-1]
        t_arr = ts[-1] + F["rhoT"][nfin - 1] * 0.5 * Tk(ra_f, rp_f)
        dvc = float(sum(abs(v) for v in dvs.values()))
        Rpt = float(F["rp"][min(nfin, n0) - 1] - F["alt"][min(nfin, n0) - 1])
        fb = FL.first_burn_dv(ra_f, rp_f, Rpt + h_des_m)
        Tf = Tk(ra_f, rp_f)
        out.update(n=nfin, final_err_km=(ra_f - RA_T) / 1e3, t_arr_d=t_arr / 86400.0, T_final_s=Tf, w_max_d=W_MAX_PERIODS * Tf / 86400.0, dv_campaign=dvc, dv_pre=dv_pre, fb=fb, total=dvc + dv_pre + fb,
                   burns=len(dvs), net_dv=float(abs(sum(dvs.values()))), rp_f_km=rp_f / 1e3, alt_f_km=(rp_f - Rpt) / 1e3, ok_apo=bool(-TOL_KM <= (ra_f - RA_T) / 1e3 <= 0.0), dvs=dvs)
    else:
        out.update(n=len(ras), final_err_km=float("nan"), t_arr_d=float("nan"), ok_apo=False, dv_campaign=float(sum(abs(v) for v in dvs.values())), dv_pre=dv_pre, fb=float("nan"), total=float("nan"),
                   burns=len(dvs), dvs=dvs, w_max_d=float("nan"))
    if gd is not None:
        out["lam_log"] = lam_log
    if trace:
        out["trace"] = tr
    return out


def c_hist_ratio(xs, ras, tab, ra_prev0, k, stat="median"):
    """measured drop / model drop (state dependent, unit kappa) of the last k passes: the planner's estimator of kappa x noise (median)."""
    n = len(xs)
    if n < 4:
        return None
    vals = []
    for p in range(max(2, n - k + 1), n + 1):                        # passes 2..n (the first pass is excluded as in the terminal aim)
        rp_ = ras[p - 2]; dr = rp_ - ras[p - 1]; m = d0_at(tab, rp_) * math.exp(xs[p - 1]) / tab["q_req"]
        if dr > 0 and m > 0:
            vals.append(dr / m)
    return (float(np.mean(vals)) if stat == "mean" else float(np.median(vals))) if len(vals) >= 3 else None


def guard_solve(tab, kap, gd, info, lam_log, j, lvl, A_, B_, L_, P, nu0, ra_now, rp_now, t_now, tstar_s, kap_meas, bias, Hloc_rp, lvl_alt=None):
    """Revision 18 design (b): time-to-go guard around the DP of decision j. nu0 is the lambda = 0 solution of the planner's own DP (returned unchanged while the prediction is inside the band).
    The guard predicts with its own level forecast: the planner's (table slow increments, drift bias for `cap` steps) with the observed drift bias persisted for `bias_steps` steps (ASSUMED; = cap leaves
    the planner's forecast). Returns (nu, engaged)."""
    kh = kap if kap_meas is None else kap_meas
    n = len(lvl); cap = int(P["update_cap_passes"]); gain = float(P["update_gain"]); K2 = max(int(gd["bias_steps"]), cap)
    m_ = np.arange(n); corr = gain * bias * (np.minimum(m_, K2) - np.minimum(m_, cap))
    lvl0 = lvl if lvl_alt is None else lvl_alt
    lvl_g = lvl0 + corr
    if gd.get("fc_level_shift_ln"):                                            # v3 candidate: constant correction of the forecast level for steps >= 1
        lvl_g = lvl_g - gd["fc_level_shift_ln"] * (np.arange(n) >= 1)
    A2 = A_ + lvl - lvl_g; B2 = B_ + lvl - lvl_g                              # same tube [lo, hi] in the guard's level coordinates
    if gd.get("fc_reserve_scale") is not None:                                 # v3 candidate: scale the planner's reserve (half the fast reserve is kept free at each edge for steps >= 1)
        sh = (1.0 - gd["fc_reserve_scale"]) * float(P["fast_reserve_ln"]) / 2.0 * (np.arange(n) >= 1)
        A2 = A2 - sh; B2 = B2 + sh
    corr = lvl_g - lvl

    def pred(nu):
        x = lvl_g + nu
        k0 = min(j - 1, tab["n"])
        rp = rp_now + np.array([tfun(tab, "rp", min(j + m, tab["n"])) - tfun(tab, "rp", k0) for m in range(n)]) - Hloc_rp * nu
        return predict_arrival(tab, kh, ra_now, rp_now, t_now, x, rp, P if gd["term_model"] else None)

    W3 = W_MAX_PERIODS * Tk(RA_T, GUARD_W_RP_M)                                 # W_max basis: 3 periods of the nominal closing orbit
    nu0g = nu0 if not np.any(corr) else dpf(A2, B2, L_, P["dp_burn_penalty_mps"] if gd.get("pen_g") is None else gd["pen_g"])
    if nu0g is None:
        return nu0, False
    ta0, N0, drops0 = pred(nu0g)
    solves = 1
    if not math.isfinite(ta0):
        ta0 = t_now + 1.0e9
    ta_raw = ta0
    if gd["ema_beta"] > 0:                                                      # smooth the predicted arrival over decisions (the prediction is an absolute time: it does not drift with the clock)
        prev = info.get("ta_ema"); ta_ema = ta0 if prev is None else gd["ema_beta"] * prev + (1.0 - gd["ema_beta"]) * ta0
        info["ta_ema"] = ta_ema; ta0 = ta_ema
    dlt = ta_raw - ta0                                                          # raw minus smoothed: removed from every prediction of this decision
    s = gd["s_frac"] * W3 + gd["rel_slack"] * max(ta0 - t_now, 0.0); a_s = gd["aim_frac"] * s
    info.setdefault("pred_log", []).append((j, ta_raw / 86400.0, t_now / 86400.0))
    if gd.get("keep_fc"):
        info.setdefault("fc", {})[j] = (lvl_g + nu0g)[:60].copy(); info.setdefault("fc_lvl", {})[j] = lvl_g[:60].copy()
    if tstar_s - s <= ta0 <= tstar_s + s:
        lt = info.get("latch")
        if gd.get("latch") and lt is not None:                                  # v3 candidate: keep the last tightening while the prediction has not crossed the aim (hysteresis against re-engaging every few passes)
            if (lt[0] and ta0 <= tstar_s) or ((not lt[0]) and ta0 >= tstar_s):
                info["latch"] = None
            else:
                e_l = lt[1] * gd["eps_max"] * float(np.max(B2 - A2)); pg = P["dp_burn_penalty_mps"] if gd.get("pen_g") is None else gd["pen_g"]
                nu_l = dpf(np.minimum(A2 + e_l, B2), B2, L_, pg) if lt[0] else dpf(A2, np.maximum(B2 - e_l, A2), L_, pg)
                if nu_l is not None:
                    info["latch_steps"] = info.get("latch_steps", 0) + 1
                    return nu_l, True
        if gd.get("keep_fc"):
            info.setdefault("fcc", {})[j] = (lvl_g + nu0)[:60].copy()
        return nu0, False
    info["engaged"] += 1
    late = ta0 > tstar_s + s
    dr = np.array(list(drops0) + [0.0] * n)[:n]
    dN = drops0[N0 - 1] if N0 >= 1 else 1.0
    wT = np.where(np.arange(n) < N0, dr / max(dN, 1.0), 0.0)
    sign = 1.0 if late else -1.0
    mg = gd["edge_margin_ln"]
    Ag, Bg = A2, B2
    if mg > 0:
        Ag = A2 + mg; Bg = B2 - mg; bad = Ag > Bg; mid = 0.5 * (A2 + B2); Ag = np.where(bad, mid, Ag); Bg = np.where(bad, mid, Bg)

    pen_g = P["dp_burn_penalty_mps"] if gd.get("pen_g") is None else gd["pen_g"]
    tube = gd["actuator"] == "tube"
    width = float(np.max(B2 - A2))

    def solve(u):
        """u in [0, 1]: the guard's price. 'price': lambda = lam_min..lam_max on a log scale, linear term -lambda wT nu (Revision 18). 'tube': the tube is tightened by u x eps_max on the side that
        the arrival needs (late: lower edge up; early: upper edge down; the upper edge is never raised), the DP stays lazy."""
        if tube:
            e_ = u * gd["eps_max"] * width
            if late:
                A_t = np.minimum(Ag + e_, Bg); return dpf(A_t, Bg, L_, pen_g)
            B_t = np.maximum(Bg - e_, Ag); return dpf(Ag, B_t, L_, pen_g)
        lam = math.exp(math.log(gd["lam_min"]) + u * (math.log(gd["lam_max"]) - math.log(gd["lam_min"])))
        return dpf(Ag, Bg, L_, pen_g, lin=-sign * lam * wT)

    nu_hi = solve(1.0); solves += 1
    if nu_hi is None:
        info["infeasible"] += 1; info["solves"] += solves; return nu0, True
    th, _, _ = pred(nu_hi); th = (th - dlt) if math.isfinite(th) else t_now + 1.0e9
    if (late and th > tstar_s + s) or ((not late) and th < tstar_s - s):
        info["infeasible"] += 1; info["solves"] += solves; lam_log.append((j, "inf", 1.0, ta0 / 86400.0, th / 86400.0, s / 86400.0))
        return (nu_hi if gd["infeasible_mode"] == "pin" else nu0), True
    lo_u = 0.0; hi_u = 1.0; best = nu_hi; u_used = 1.0; tm_last = th
    while solves < gd["max_solves"]:
        mid = 0.5 * (lo_u + hi_u); nu_m = solve(mid); solves += 1
        if nu_m is None:
            break
        tm, _, _ = pred(nu_m); tm = (tm - dlt) if math.isfinite(tm) else t_now + 1.0e9; tm_last = tm
        if abs(tm - tstar_s) <= a_s:
            best = nu_m; u_used = mid; break
        too_late = tm > tstar_s
        if (late and too_late) or ((not late) and not too_late):
            lo_u = mid
        else:
            hi_u = mid; best = nu_m; u_used = mid
    info["solves"] += solves; lam_log.append((j, sign * u_used, solves, ta0 / 86400.0, tm_last / 86400.0 if solves > 2 else th / 86400.0, s / 86400.0))
    info["latch"] = (bool(late), float(u_used))
    if gd.get("keep_fc"):
        info.setdefault("fcc", {})[j] = (lvl_g + best)[:60].copy()
    return best, True


# ------------------------------------------------------------------------------------------------ flight index

TOKEN = {"ODY": "ODY", "G5863": "G1", "G0975": "G0975"}
OLR_FAMILY = ("OLR", "OLRZ", "OLRT", "OLRT1", "OLRT2", "OLRT3", "OLRT35", "OLRT4", "OLRT5")
_FJ = None


def split_id(run_id):
    g, rest = run_id.split("_PL_"); arm = rest.split("_pert_")[0]; tail = rest.split("_pert_ODY_s")[1]
    seed = int(tail.split("_")[0]); lvl = tail.split("_L")[1] if "_L" in tail else None
    return dict(geom=TOKEN[g], gid=g, arm=arm, seed=seed, level=(None if lvl is None else int("".join(ch for ch in lvl if ch.isdigit()))), level_tag=lvl)


def flights_index():
    """run_id -> corrected-dV row of the 366 flights (dv_correction/flights.json) plus every other flight directory found (rev17)."""
    global _FJ
    if _FJ is None:
        J = json.load(open(TD / "dv_correction" / "flights.json")); _FJ = {r["run_id"]: r for r in J if r.get("status") != "not_corrected"}
    return _FJ


def rev17_complete():
    """Rev 17 flight dirs whose summary row exists (complete). Summaries: summary_*rev17* anywhere under ROOT; a dir without a summary row is skipped."""
    rows = set()
    for f in glob.glob(str(ROOT / "summary_*r17*.csv")) + glob.glob(str(ROOT / "summary_rev17*.csv")) + glob.glob(str(ROOT / "rev17" / "summary*.csv")) + glob.glob(str(ROOT / "summary_*ws_r17*.csv")):
        try:
            d = pd.read_csv(f)
            for rid, st in zip(d.run_id, d.status):
                if st == "ok":
                    rows.add(rid)
        except Exception:
            pass
    return rows


_CACHE = {}


def get_flight(run_id):
    if run_id not in _CACHE:
        m = split_id(run_id); tab = get_tab(m["geom"])
        _CACHE[run_id] = extract(run_id, tab, KAPPA[m["geom"]])
    return _CACHE[run_id]


def dev_paths(geom):
    """DEV natural paths: devT OLRT flights, seeds 1-8 (tuning set; the OLRT arm is OLRT3 without the replan-deviation threshold, the closest flown relative that carries the terminal aim)."""
    pre = {"ODY": "ODY", "G1": "G5863"}[geom]
    return [f"{pre}_PL_OLRT_pert_ODY_s{s}" for s in range(1, 9)]


def dev_paths_ext(geom):
    """DEV paths for guard v2: the 8 devT OLRT paths plus the flown OLRT3 dev flights of the trx50 Revision 15 full-dev block (trx_pull4/planner/devR: G1 seeds 1-8, ODY seeds 4, 6, 8; dev seeds only)."""
    pre = {"ODY": "ODY", "G1": "G5863"}[geom]
    extra = [f"{pre}_PL_OLRT3_pert_ODY_s{s}" for s in range(1, 9) if (ROOT / "../trx_pull4/planner/devR" / f"{pre}_PL_OLRT3_pert_ODY_s{s}" / "passes.csv").exists()
             and not (ROOT / "devR" / f"{pre}_PL_OLRT3_pert_ODY_s{s}" / "passes.csv").exists()]
    return dev_paths(geom) + extra


def held_paths(geom):
    """HELD natural paths (TEST only): held OLR flights, seeds 13-28."""
    pre = {"ODY": "ODY", "G1": "G5863"}[geom]
    return [f"{pre}_PL_OLR_pert_ODY_s{s}" for s in range(13, 29) if run_dir(f"{pre}_PL_OLR_pert_ODY_s{s}") is not None]


def stat(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if len(v) == 0:
        return dict(n=0)
    return dict(n=int(len(v)), median=float(np.median(v)), q25=float(np.quantile(v, .25)), q75=float(np.quantile(v, .75)), min=float(v.min()), max=float(v.max()), mean=float(v.mean()))


# ------------------------------------------------------------------------------------------------ validation: own path (exact) and cross-arm (the test of the noise indexing)

EXCLUDE_ARMS = ("OLRT1", "OLRT4")          # remedies 1 and 4 are not implemented in the replay (not needed: only OLRT3 and its relatives are replayed)


def flown_summary(F):
    return dict(dv=F["dv_flown"], err_km=F["final_err_m"] / 1e3, t_arr_d=F["t_arr_d"], n=F["n"], hard=bool(F["qs_max"] > QS_LIMIT))


def own_path_validation(verbose=True):
    rows = []
    for rid in sorted(flights_index()):
        if "_PL_" not in rid:
            continue
        m = split_id(rid)
        if m["arm"] not in OLR_FAMILY or m["arm"] in EXCLUDE_ARMS or m["geom"] not in ("ODY", "G1"):
            continue
        if run_dir(rid) is None:
            continue
        F = get_flight(rid)
        if not F["reached"]:
            continue
        r = replay(F, m["geom"], m["arm"]); f = flown_summary(F)
        rows.append(dict(run_id=rid, arm=m["arm"], geom=m["geom"], seed=m["seed"], d_dv=r["dv_campaign"] - f["dv"], d_err_km=r["final_err_km"] - f["err_km"], d_t_d=r["t_arr_d"] - f["t_arr_d"], d_n=r["n"] - f["n"], hard_agree=bool(r["hard_violation"] == f["hard"])))
    df = pd.DataFrame(rows)
    out = dict(n=int(len(df)), max_abs_d_dv_mps=float(df.d_dv.abs().max()), max_abs_d_err_km=float(df.d_err_km.abs().max()), max_abs_d_t_d=float(df.d_t_d.abs().max()), max_abs_d_n=int(df.d_n.abs().max()),
               hard_agree_all=bool(df.hard_agree.all()), by_arm={a: int((df.arm == a).sum()) for a in sorted(df.arm.unique())})
    if verbose:
        print("own-path validation:", out)
    return out, df


def cross_arm_validation(noise_modes=("pass", "state"), geoms=("ODY", "G1"), targets=("OLRT3", "OLRT35", "OLRT", "OLRT5"), verbose=True):
    """Arm A (flown at seed s) replayed on the natural path of another flown arm B of the same seed and geometry (dev seeds 1-8 only: no held data used). Differences replay - flown for A."""
    idx = flights_index(); by = {}
    for rid in idx:
        if "_PL_" not in rid:
            continue
        m = split_id(rid)
        if m["geom"] in geoms and 1 <= m["seed"] <= 8 and m["level"] is None and run_dir(rid) is not None:
            by.setdefault((m["geom"], m["seed"]), {})[m["arm"]] = rid
    rows = []
    for (g, s), arms in sorted(by.items()):
        for a in targets:
            if a not in arms or a in EXCLUDE_ARMS:
                continue
            FA = get_flight(arms[a])
            if not FA["reached"]:
                continue
            fa = flown_summary(FA)
            for b, rb in sorted(arms.items()):
                if b == a:
                    continue
                try:
                    FB = get_flight(rb)
                except Exception:
                    continue
                if not FB["reached"] or not np.isfinite(FB["Hloc"]):
                    continue
                for nm in noise_modes:
                    try:
                        r = replay(FB, g, a, noise=nm, seed=s)
                    except Exception as ex:
                        continue
                    if not r["reached"]:
                        continue
                    rows.append(dict(geom=g, seed=s, target=a, source=b, noise=nm, d_dv=r["dv_campaign"] - fa["dv"], d_err_km=r["final_err_km"] - fa["err_km"], d_t_d=r["t_arr_d"] - fa["t_arr_d"], rel_t=r["t_arr_d"] / fa["t_arr_d"] - 1.0,
                                     flown_t=fa["t_arr_d"], flown_dv=fa["dv"], rep_dv=r["dv_campaign"], hard_rep=r["hard_violation"], hard_fl=fa["hard"], ext=r["ext_used"]))
    df = pd.DataFrame(rows)
    out = {}
    for nm in noise_modes:
        for g in geoms:
            d = df[(df.noise == nm) & (df.geom == g)]
            if len(d) == 0:
                continue
            out[f"{nm}@{g}"] = dict(n=int(len(d)), d_dv=stat(d.d_dv), abs_d_dv_median=float(d.d_dv.abs().median()), d_err_km=stat(d.d_err_km), abs_d_err_km_median=float(d.d_err_km.abs().median()), d_t_d=stat(d.d_t_d),
                                    abs_d_t_d_median=float(d.d_t_d.abs().median()), abs_rel_t_median=float(d.rel_t.abs().median()), rank_corr_dv=float(d[["rep_dv", "flown_dv"]].corr(method="spearman").iloc[0, 1]) if len(d) > 3 else float("nan"),
                                    hard_agree=float((d.hard_rep == d.hard_fl).mean()), hard_replay=int(d.hard_rep.sum()), hard_flown=int(d.hard_fl.sum()))
    if verbose:
        for k, v in out.items():
            print(k, "n", v["n"], "| |d dV| med %.3f | |d err| med %.2f km | |d t| med %.2f d (rel %.3f) | spearman dV %.2f | hard agree %.2f" % (v["abs_d_dv_median"], v["abs_d_err_km_median"], v["abs_d_t_d_median"], v["abs_rel_t_median"], v["rank_corr_dv"], v["hard_agree"]))
    return out, df


def noise_correlation():
    """Do two arms flown at the same seed meet the same noise? High-pass (15-pass moving mean removed) natural ln q of two arms: correlation by pass index against correlation after interpolating the second arm
    in TIME at the first arm's pass times; dev seeds 1-8 only, arms OLRT, OLRT3, OLRT35, OLRT5, OLRZ, L3 (pairs of different arms, same geometry and seed)."""
    idx = flights_index(); by = {}
    for rid in idx:
        if "_PL_" not in rid:
            continue
        m = split_id(rid)
        if m["geom"] in ("ODY", "G1") and 1 <= m["seed"] <= 8 and m["level"] is None and run_dir(rid) is not None and m["arm"] in ("OLRT", "OLRT3", "OLRT35", "OLRT5", "OLRZ", "L3"):
            by.setdefault((m["geom"], m["seed"]), {})[m["arm"]] = rid
    hp = lambda x, k=15: x - np.convolve(np.pad(x, (k // 2, k // 2), mode="edge"), np.ones(k) / k, mode="valid")[:len(x)]
    res = {"ODY": [], "G1": []}
    for (g, s_), arms in sorted(by.items()):
        ks = sorted(arms)
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                A = get_flight(arms[a]); B = get_flight(arms[b]); n = min(A["n"], B["n"]) - 1
                xa = A["xnat"][1:n + 1]; xb = B["xnat"][1:n + 1]; xbt = np.interp(A["t"][1:n + 1], B["t"], B["xnat"])
                res[g].append((float(np.corrcoef(hp(xa), hp(xb))[0, 1]), float(np.corrcoef(hp(xa), hp(xbt))[0, 1])))
    return {g: dict(pairs=len(v), corr_by_pass_median=float(np.median([x[0] for x in v])), corr_by_time_median=float(np.median([x[1] for x in v])),
                    corr_by_pass_q25_q75=[float(np.quantile([x[0] for x in v], .25)), float(np.quantile([x[0] for x in v], .75))], corr_by_time_q25_q75=[float(np.quantile([x[1] for x in v], .25)), float(np.quantile([x[1] for x in v], .75))]) for g, v in res.items()}


# ------------------------------------------------------------------------------------------------ scoring of a cell (a set of paths)

NOISE = "pass"            # set by `validate` (see results_rev20_replay.md); overridden by --noise


def tstar_for(token, ic=None):
    """t* (days) = the plan's nominal model-world calibrated arrival for the geometry and IC (WP1 ic_plan t_arr_d when given)."""
    return nominal_arrival_d(token, math.log(QFIRST_CAL[token]) + ic_shift(ic))


def run_cell(run_ids, token, arm, guard, tstar_d, ic=None, noise=None, h_des_m=H_DES_M):
    out = []
    for k, rid in enumerate(run_ids):
        F = get_flight(rid); m = split_id(rid)
        r = replay(F, token, arm, ic=ic, noise=noise or NOISE, guard=guard, tstar_d=tstar_d, h_des_m=h_des_m, seed=m["seed"] * 100 + (m["level"] or 0))
        r["tstar_d"] = tstar_d; r["seed"] = m["seed"]; r["level"] = m["level"]
        if r["reached"]:
            r["dt_d"] = r["t_arr_d"] - tstar_d; r["ok_win"] = bool(abs(r["dt_d"]) <= r["w_max_d"])
        else:
            r["dt_d"] = float("nan"); r["ok_win"] = False
        r["success"] = bool(r["reached"] and r["ok_apo"] and not r["hard_violation"] and r["ok_win"])
        out.append(r)
    return out


def summarize(res):
    n = len(res); ok = [r for r in res if r["reached"]]
    g = lambda k: [r[k] for r in ok]
    s = dict(n=n, reached=len(ok), hard=int(sum(r["hard_violation"] for r in res)), apo_ok=int(sum(r["ok_apo"] for r in res)), window_ok=int(sum(r["ok_win"] for r in res)), success=int(sum(r["success"] for r in res)),
             med_dv_campaign=float(np.median(g("dv_campaign"))) if ok else float("nan"), med_fb=float(np.median(g("fb"))) if ok else float("nan"), med_total=float(np.median(g("total"))) if ok else float("nan"),
             med_abs_dt_d=float(np.median(np.abs(g("dt_d")))) if ok else float("nan"), med_dt_d=float(np.median(g("dt_d"))) if ok else float("nan"), med_arrival_d=float(np.median(g("t_arr_d"))) if ok else float("nan"),
             tstar_d=float(np.mean([r["tstar_d"] for r in res])) if res else float("nan"), ext_max=int(max([r["ext_used"] for r in res] or [0])),
             infeasible_flights=int(sum(1 for r in res if r["info"].get("infeasible", 0) > 0)), engaged_decisions=int(sum(r["info"].get("engaged", 0) for r in res)), solves=int(sum(r["info"].get("solves", 0) for r in res)))
    ss = [r["total"] for r in res if r["success"]]
    s["med_total_success"] = float(np.median(ss)) if ss else float("nan")
    return s


def paired(base, test):
    """paired by path (same order): campaign dV differences test - base on paths both reached"""
    d = [t["dv_campaign"] - b["dv_campaign"] for b, t in zip(base, test) if b["reached"] and t["reached"]]
    dt = [t["total"] - b["total"] for b, t in zip(base, test) if b["reached"] and t["reached"]]
    return dict(n=len(d), med_extra_dv=float(np.median(d)) if d else float("nan"), q25=float(np.quantile(d, .25)) if d else float("nan"), q75=float(np.quantile(d, .75)) if d else float("nan"),
                max_extra_dv=float(np.max(d)) if d else float("nan"), med_extra_total=float(np.median(dt)) if dt else float("nan"))


# ------------------------------------------------------------------------------------------------ tuning on DEV paths only (seeds 1-8)

DEV_EXT = False          # guard v2 tuning sets this True: the dev set then includes the trx50 flown OLRT3 dev flights


def dev_cell(geom, arm="OLRT3", guard=None, ic=None, noise=None):
    ts = tstar_for(geom, ic)
    return run_cell((dev_paths_ext if DEV_EXT else dev_paths)(geom), geom, arm, guard, ts, ic, noise)


def tune_stage(configs, out_name, noise=None, geoms=("ODY", "G1")):
    """configs: list of guard-parameter dicts (overrides of GUARD_DEFAULT). Writes one JSON line per config to OUT/<out_name>.jsonl and returns the rows. DEV paths only."""
    OUT.mkdir(parents=True, exist_ok=True); rows = []
    base = {g: dev_cell(g, noise=noise) for g in geoms}
    sb = {g: summarize(base[g]) for g in geoms}
    for cfgd in configs:
        row = dict(guard=cfgd, per_geom={}); win = hard = apo = 0; extra = []
        t0 = time.time()
        for g in geoms:
            res = dev_cell(g, guard=cfgd, noise=noise); sm = summarize(res); pr = paired(base[g], res)
            row["per_geom"][g] = dict(summary=sm, paired=pr, baseline=sb[g])
            win += sm["window_ok"]; hard += sm["hard"]; apo += sm["apo_ok"]; extra += [t["dv_campaign"] - b["dv_campaign"] for b, t in zip(base[g], res) if b["reached"] and t["reached"]]
        row.update(window_ok=win, hard=hard, apo_ok=apo, apo_ok_base=sum(sb[g]["apo_ok"] for g in geoms), window_ok_base=sum(sb[g]["window_ok"] for g in geoms), n=sum(sb[g]["n"] for g in geoms),
                   med_extra_dv=float(np.median(extra)) if extra else float("nan"), wall_s=time.time() - t0)
        rows.append(row)
        with open(OUT / f"{out_name}.jsonl", "a") as fh:
            fh.write(json.dumps(row, default=float) + "\n")
        print("guard", cfgd, "| window %d/%d (base %d) hard %d apo %d (base %d) | median extra campaign dV %.2f | %.0fs" % (win, row["n"], row["window_ok_base"], hard, apo, row["apo_ok_base"], row["med_extra_dv"], row["wall_s"]), flush=True)
    return rows


# ------------------------------------------------------------------------------------------------ evaluation sets

GUARD_TUNED = dict(actuator="tube", pen_g=0.1, rel_slack=0.15, bias_steps=60, k_window=20, ema_beta=0.5)    # chosen on DEV paths only (seeds 1-8): tune_stage* jsonl in rev20_data/replay, rule in the results note

# Round 2 (POST HOC: defined after the first look at the held ODY arrival errors, which showed the same early bias that the dev paths show): the aim point is shifted by minus the median dev arrival error of the
# round-1 guard (evaluate_dev_ref.json, ODY -0.3043 d, G1 +0.0618 d). Dev-derived number, but the decision to add it followed a look at held results: held round-2 numbers are NOT independent; Rev 17 paths are.
GUARD_ROUND2 = dict(GUARD_TUNED, t_offset_d={"ODY": 0.3043, "G1": -0.0618})

# Guard v2 (Revision 20c): FROZEN after the dev-only v2 tuning (v2_stage*.jsonl; 30 configurations on 27 dev paths under the three indexings): round-1 parameters plus an ODY aim offset of +0.15 d
# (early-bias diagnosis, note results_rev20_guard_v2.md) and a guard deviation threshold of 0.05 m/s.
GUARD_V2 = dict(GUARD_TUNED, t_offset_d={"ODY": 0.15, "G1": 0.0}, dev_thr_g=0.05)


def toolbox_compare(geom, scenario, noise, guard):
    """Per held path: OLRT3 alone replayed with trace, then the cheapest end-of-campaign fix (planner_rev20_window.cheapest, W = 3 final periods, t* of the cell; None = window failure, late arrivals that
    phasing cannot absorb fail) against the guard's extra cost (guard total minus OLRT3-alone total, FB included)."""
    import planner_rev20_window as W
    ic = get_ic(geom, scenario); ts = tstar_for(geom, ic); rows = []
    for rid in held_paths(geom):
        F = get_flight(rid); m = split_id(rid); sd = m["seed"] * 100
        a = replay(F, geom, "OLRT3", ic=ic, noise=noise, seed=sd, trace=True); g = replay(F, geom, "OLRT3", ic=ic, noise=noise, guard=guard, tstar_d=ts, seed=sd)
        tr = a["trace"]; inc = float(np.median(pd.read_csv(run_dir(rid) / "passes.csv").i_rad.values))
        ra = np.array([t[3] for t in tr]); rp = np.array([t[5] for t in tr]); tc = a["t_arr_d"] * 86400.0; Wd = a["w_max_d"] * 86400.0
        st = W.prepare(ra[:-1], rp[:-1], np.full(len(ra) - 1, inc), ra[-1], rp[-1], a["alt_f_km"] * 1e3, inc)
        c = W.cheapest(tc, ts * 86400.0, Wd, st)
        rows.append(dict(run_id=rid, dt_alone_over_W=(a["t_arr_d"] - ts) / a["w_max_d"], dt_guard_over_W=(g["t_arr_d"] - ts) / g["w_max_d"], alone_win=abs(a["t_arr_d"] - ts) <= a["w_max_d"], guard_win=abs(g["t_arr_d"] - ts) <= g["w_max_d"],
                         tool_cost=None if c is None else c["cost"], tool_kind=None if c is None else c["kind"], guard_extra_total=g["total"] - a["total"], guard_extra_campaign=g["dv_campaign"] - a["dv_campaign"], hard_alone=a["hard_violation"], hard_guard=g["hard_violation"]))
    return rows


def get_ic(geom, scenario):
    """WP1 IC (planner_rev20_ic.select_ic) for a scenario 'S-fixed' | 'S-free', or None for the reference start / when the module is missing or its selftest fails."""
    if scenario in (None, "ref"):
        return None
    try:
        import planner_rev20_ic as IC
        return IC.select_ic(geom, scenario)
    except Exception as ex:                                  # pragma: no cover
        print("WP1 IC unavailable:", repr(ex)); return "unavailable"


def slim(r):
    keep = ("run_id", "reached", "hard_violation", "qs_max", "ext_used", "n", "final_err_km", "t_arr_d", "T_final_s", "w_max_d", "dv_campaign", "dv_pre", "fb", "total", "burns", "net_dv", "rp_f_km", "alt_f_km", "ok_apo",
            "tstar_d", "seed", "level", "dt_d", "ok_win", "success")
    d = {k: r.get(k) for k in keep}; d["info"] = {k: v for k, v in r["info"].items() if k != "pred_log"}
    return d


def evaluate_set(name, geoms, path_fn, guard, scenarios=("ref",), noise=None, tstar_shift_d=0.0, arm="OLRT3"):
    out = dict(name=name, guard=guard, noise=noise or NOISE, cells={})
    for g in geoms:
        ids = path_fn(g)
        if not ids:
            continue
        for sc in scenarios:
            ic = get_ic(g, sc)
            if ic == "unavailable" or (ic is not None and not ic.get("feasible", True)):
                out["cells"][f"{g}|{sc}"] = dict(status="IC unavailable or infeasible"); continue
            ts = tstar_for(g, ic) + tstar_shift_d
            base = run_cell(ids, g, arm, None, ts, ic, noise); gd = run_cell(ids, g, arm, guard, ts, ic, noise)
            out["cells"][f"{g}|{sc}"] = dict(geom=g, scenario=sc, tstar_d=ts, n=len(ids), baseline=summarize(base), guard=summarize(gd), paired=paired(base, gd), ic=None if ic is None else {k: ic.get(k) for k in ("x0", "nu_pre", "dv_pre_mps", "t_arr_d", "total_dv_mps", "campaign_dv_mps", "fb_mps")},
                                             rows_baseline=[slim(r) for r in base], rows_guard=[slim(r) for r in gd])
    return out


def acs(cell, level=None, n_paths=None):
    """AC2-AC4 verdicts for one cell (baseline = OLRT3 without the guard, guard = with)."""
    b, g = cell["baseline"], cell["guard"]; n = cell["n"]
    ac2 = (g["hard"] == 0) and (g["hard"] <= b["hard"])
    ac3 = g["apo_ok"] >= b["apo_ok"] - max(1, int(round(n / 16.0)))                   # not below OLRT3's by more than one path in 16
    need = None
    if level is None:
        need = math.ceil(14.0 / 16.0 * n - 1e-9)
    else:
        need = 3 if level in (0, 1, 2, 3, 4) and n == 4 and level != 0 else math.ceil(14.0 / 16.0 * n - 1e-9)
    ac4 = g["window_ok"] >= need
    return dict(AC2="PASS" if ac2 else "FAIL", AC3="PASS" if ac3 else "FAIL", AC4="PASS" if ac4 else "FAIL", need_window=need)


def dump(obj, name):
    OUT.mkdir(parents=True, exist_ok=True); p = OUT / name
    p.write_text(json.dumps(obj, indent=1, default=float)); return p


def sha16(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ Rev 17 paths (TEST only)

def rev17_paths(level, source_arm, geom="ODY"):
    """Rev 17 flown natural paths: arm `source_arm` (OL, OLR or OLRT3) at held seeds 13-16, level L<level>; only dirs that hold a complete flight (passes.csv ends at the target or the horizon, replan.csv has the same pass count)."""
    ids = []; done = rev17_complete()                       # summary rows (written when the block ends); empty until then, the passes.csv proxy below applies
    for s in range(13, 17):
        rid = f"ODY_PL_{source_arm}_pert_ODY_s{s}_L{level}"
        d = run_dir(rid)
        if d is None or (done and rid not in done):
            continue
        try:
            P = pd.read_csv(d / "passes.csv"); R_ = pd.read_csv(d / "replan.csv")
        except Exception:
            continue
        ra = P.ra_m.values.astype(float)
        if not (np.isfinite(ra).any() and (ra[np.isfinite(ra)] <= RA_T).any()) or len(R_) < len(P) - 1:
            continue
        ids.append(rid)
    return ids


def rev17_stopped(min_passes=190):
    """Rev 17 flights that stopped on safety_q (a flown hard violation) with at least min_passes passes (the rest of the campaign is then bootstrapped, ASSUMED): OLRT3 alone and with the guard replayed on their natural
    paths, hard-violation outcome only (the window is not meaningful on a truncated path). Needs the Rev 17 summary file."""
    out = []; ts = tstar_for("ODY"); gd = {**GUARD_DEFAULT, **GUARD_TUNED}
    for f in glob.glob(str(ROOT / "summary_rev17*.csv")):
        d = pd.read_csv(f)
        for rid, st, sr, npass in zip(d.run_id, d.status, d.stop_reason, d.passes):
            if st != "ok" or sr != "safety_q" or npass < min_passes:
                continue
            F = get_flight(rid); m = split_id(rid)
            a = replay(F, "ODY", "OLRT3", noise="pass", seed=m["seed"] * 100 + m["level"]); g = replay(F, "ODY", "OLRT3", noise="pass", guard=gd, tstar_d=ts, seed=m["seed"] * 100 + m["level"])
            out.append(dict(run_id=rid, source_arm=m["arm"], level=m["level"], seed=m["seed"], passes_flown=int(npass), qs_max_flown=F["qs_max"], ext_alone=a["ext_used"], hard_alone=a["hard_violation"], qs_max_alone=a["qs_max"],
                            hard_guard=g["hard_violation"], qs_max_guard=g["qs_max"], reached_alone=a["reached"], reached_guard=g["reached"]))
    return out


# ------------------------------------------------------------------------------------------------ selftest

def selftest():
    rng = np.random.default_rng(5)
    for n in (6, 40):                                   # fast DP = reference DPs
        for pen in (0.0, 0.02):
            lev = np.cumsum(rng.normal(0, .2, n)); w = rng.uniform(.5, 3, n); A = -.4 - lev; B = .4 - lev; L = PLB.tie_weights(w, 1e-6)
            assert np.abs(PLB.dp_exact(A, B, L, pen=pen)["nu"] - dpf(A, B, L, pen)).max() < 1e-12
            if pen == 0.0:
                lin = rng.normal(0, .3, n); assert np.abs(O.dp_lin(A, B, L, lin) - dpf(A, B, L, 0.0, lin)).max() < 1e-12
    for rid in ("ODY_PL_OLRT3_pert_ODY_s1", "G5863_PL_OLRT35_pert_ODY_s2"):       # own path: exact decision logic, state and clock
        if run_dir(rid) is None:
            continue
        m = split_id(rid); F = get_flight(rid); r = replay(F, m["geom"], m["arm"]); f = flown_summary(F)
        assert abs(r["dv_campaign"] - f["dv"]) < 1e-6 and abs(r["final_err_km"] - f["err_km"]) < 1e-4 and abs(r["t_arr_d"] - f["t_arr_d"]) < 1e-6 and r["n"] == f["n"], (rid, r["dv_campaign"], f["dv"])
    F = get_flight("ODY_PL_OLRT3_pert_ODY_s1")                      # guard with an unreachable-wide band never engages and leaves the plan unchanged
    r0 = replay(F, "ODY", "OLRT3"); r1 = replay(F, "ODY", "OLRT3", guard=dict(s_frac=1e6), tstar_d=70.0)
    assert r1["info"]["engaged"] == 0 and abs(r0["dv_campaign"] - r1["dv_campaign"]) < 1e-12 and abs(r0["t_arr_d"] - r1["t_arr_d"]) < 1e-12
    r2 = replay(F, "ODY", "OLRT3", guard=dict(), tstar_d=tstar_for("ODY"))        # guard: Q_s limit and tube untouched (hard violation flag computed), at most 12 solves per decision
    assert r2["info"]["solves"] <= 12 * r2["info"]["engaged"] + 1
    assert abs(tstar_for("ODY") - nominal_arrival_d("ODY")) < 1e-12
    print("planner_rev20_replay selftest ok: fast DP equals dp_exact/dp_lin, OLRT3/OLRT35 own-path replay exact, idle guard is the identity, solves <= 12 per engaged decision")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "toolbox", "rev17stopped", "noisecorr", "validate", "evaluate", "curve", "rev17", "gen1", "tune"]); ap.add_argument("--set", default="held"); ap.add_argument("--scenarios", default="ref")
    ap.add_argument("--noise", default=None); ap.add_argument("--tag", default=""); ap.add_argument("--shift", type=float, default=0.0); ap.add_argument("--round", type=int, default=1); a = ap.parse_args()
    global NOISE
    if a.noise:
        NOISE = a.noise
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "toolbox":
        gd = {**GUARD_DEFAULT, **GUARD_V2}
        out = {f"{g}|{sc}|{nm}": toolbox_compare(g, sc, nm, gd) for nm in ("pass", "time", "state") for g in ("ODY", "G1") for sc in a.scenarios.split(",")}
        dump(out, "toolbox_v2.json"); return
    if a.cmd == "rev17stopped":
        o = rev17_stopped(); dump(o, "rev17_stopped_replay.json")
        for r in o:
            print(r); 
        return
    if a.cmd == "noisecorr":
        o = noise_correlation(); dump(o, "noise_correlation.json"); print(o); return
    gd = {**GUARD_DEFAULT, **(GUARD_V2 if a.round == 3 else GUARD_ROUND2 if a.round == 2 else GUARD_TUNED)}
    if a.cmd == "validate":
        o, df = own_path_validation(); c, cdf = cross_arm_validation(noise_modes=("pass", "state", "time"))
        dump(dict(own_path=o, cross_arm=c), "validate.json"); cdf.to_csv(OUT / "validate_cross_arm_rows.csv", index=False); df.to_csv(OUT / "validate_own_path_rows.csv", index=False)
    elif a.cmd == "evaluate":
        fn = dev_paths if a.set == "dev" else held_paths
        res = evaluate_set(a.set, ("ODY", "G1"), fn, gd, a.scenarios.split(","), noise=a.noise, tstar_shift_d=a.shift)
        dump(res, f"evaluate_{a.set}_{a.scenarios.replace(',', '+')}{a.tag}.json")
        for k, c in res["cells"].items():
            if "baseline" in c:
                print(k, "t*=%.2f" % c["tstar_d"], "baseline", {x: c["baseline"][x] for x in ("hard", "apo_ok", "window_ok", "med_dv_campaign", "med_fb", "med_total")}, "guard", {x: c["guard"][x] for x in ("hard", "apo_ok", "window_ok", "med_dv_campaign", "med_fb", "med_total")}, c["paired"], acs(c))
    elif a.cmd == "curve":
        out = {}
        for g in ("ODY", "G1"):
            ids = held_paths(g); t0 = tstar_for(g); rows = []
            for d_ in (-8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0):
                res = run_cell(ids, g, "OLRT3", gd, t0 + d_, None, None)
                for r in res:                      # score against the SHIFTED target (the guard's own t*) and against the nominal t*
                    r["dt_nominal_d"] = r["t_arr_d"] - t0 if r["reached"] else float("nan")
                sm = summarize(res); rows.append(dict(shift_d=d_, tstar_d=t0 + d_, summary=sm, med_abs_dt_nominal_d=float(np.median([abs(r["dt_nominal_d"]) for r in res if r["reached"]])), med_dt_nominal_d=float(np.median([r["dt_nominal_d"] for r in res if r["reached"]])),
                                                      med_total=sm["med_total"], med_dv_campaign=sm["med_dv_campaign"]))
                print(g, d_, "window %d/%d" % (sm["window_ok"], sm["n"]), "dV %.2f" % sm["med_dv_campaign"], "dt_nom %.2f" % rows[-1]["med_dt_nominal_d"], flush=True)
            out[g] = rows
        dump(out, "curve_held.json")
    elif a.cmd == "rev17":
        out = {}
        for lvl in range(5):
            for src in ("OL", "OLR", "OLRT3"):
                ids = rev17_paths(lvl, src)
                if not ids:
                    continue
                ic = get_ic("ODY", a.scenarios.split(",")[0]); ts = tstar_for("ODY", ic); base = run_cell(ids, "ODY", "OLRT3", None, ts, ic); g_ = run_cell(ids, "ODY", "OLRT3", gd, ts, ic)
                out[f"L{lvl}|{src}"] = dict(level=lvl, source_arm=src, n=len(ids), ids=ids, baseline=summarize(base), guard=summarize(g_), paired=paired(base, g_), rows_baseline=[slim(r) for r in base], rows_guard=[slim(r) for r in g_])
                print(f"L{lvl} {src}", len(ids), "window base %d guard %d | hard %d/%d | apo %d/%d" % (out[f'L{lvl}|{src}']['baseline']['window_ok'], out[f'L{lvl}|{src}']['guard']['window_ok'], out[f'L{lvl}|{src}']['baseline']['hard'], out[f'L{lvl}|{src}']['guard']['hard'], out[f'L{lvl}|{src}']['baseline']['apo_ok'], out[f'L{lvl}|{src}']['guard']['apo_ok']), flush=True)
        dump(out, f"rev17_replay{a.tag}.json")
    elif a.cmd == "gen1":
        tab = get_tab("G0975")
        ids = [r for r in sorted(flights_index()) if r.startswith("G0975_PL_OLRT") and split_id(r)["arm"] in ("OLRT3", "OLRT35", "OLRT5")]
        ts = tstar_for("G0975"); base = run_cell(ids, "G0975", "OLRT3", None, ts, None); g_ = run_cell(ids, "G0975", "OLRT3", gd, ts, None)
        out = dict(n=len(ids), tstar_d=ts, baseline=summarize(base), guard=summarize(g_), paired=paired(base, g_), flown_arrival_median_d=float(np.median([get_flight(r)["t_arr_d"] for r in ids])))
        dump(out, "gen1_replay.json"); print(out)


if __name__ == "__main__":
    main()


# ------------------------------------------------------------------------------------------------ guard v2: tuning on DEV paths under all three noise indexings

def v2_tune(configs, out_name, indexings=("pass", "time", "state"), geoms=("ODY", "G1")):
    """Each config (guard-parameter overrides of GUARD_TUNED) is scored on the 16 DEV paths under each indexing. Row: per indexing window_ok, hard, apo_ok, median extra campaign dV; worst-of-three window = min over
    indexings; extra dV for the ranking = max over indexings. Held and Rev 17 paths are never touched here."""
    OUT.mkdir(parents=True, exist_ok=True); rows = []
    base = {(nm, g): dev_cell(g, noise=nm) for nm in indexings for g in geoms}
    for ov in configs:
        gdp = {**GUARD_TUNED, **ov}; t0 = time.time(); row = dict(guard=ov, per={}); 
        for nm in indexings:
            win = hard = apo = 0; extra = []; wins_g = {}
            for g in geoms:
                res = dev_cell(g, guard=gdp, noise=nm); sm = summarize(res); win += sm["window_ok"]; hard += sm["hard"]; apo += sm["apo_ok"]; wins_g[g] = sm["window_ok"]
                extra += [t["dv_campaign"] - b["dv_campaign"] for b, t in zip(base[(nm, g)], res) if b["reached"] and t["reached"]]
            row["per"][nm] = dict(window=win, window_by_geom=wins_g, hard=hard, apo=apo, med_extra=float(np.median(extra)))
        row["worst_window"] = min(v["window"] for v in row["per"].values()); row["hard_total"] = sum(v["hard"] for v in row["per"].values()); row["max_extra"] = max(v["med_extra"] for v in row["per"].values()); row["wall_s"] = time.time() - t0
        rows.append(row)
        with open(OUT / f"{out_name}.jsonl", "a") as fh:
            fh.write(json.dumps(row, default=float) + "\n")
        print("v2", ov, "| window pass/time/state %s | worst %d | hard %d | extra dV pass/time/state %s | %.0fs" % ("/".join(str(row["per"][n]["window"]) for n in indexings), row["worst_window"], row["hard_total"], "/".join("%.2f" % row["per"][n]["med_extra"] for n in indexings), row["wall_s"]), flush=True)
    return rows
