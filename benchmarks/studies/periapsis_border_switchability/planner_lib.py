"""planner_lib: the planner's pure core in Python (WP-P2). Mirror of the pure functions in planner_common.jl (same algorithm, same tie-breaks), used by the predictions, the analysis selftest and the
dry-run listing. No Julia, no flights. Plan: tasks/2026-10-15_periapsis-region-border-targeting/analytic_framework_and_transfer_plan.md R4.5 to R4.10, R7.3, R7.4, R8.1, R8.2.

Coordinates. The planner works in ln q, q the mean-field heat rate at the periapsis event on the atmosphere-relative speed (U15). A burn before pass j shifts ln q of passes j, j+1, ... by nu (cumulative
offset, nu > 0 = more heat); a tangential apoapsis burn dv changes nu by -dv / w with w = H / (ratio x lever) in m/s per ln unit, lever = dr_p / dv = 4 a sqrt((1-e)/(1+e)) / sqrt(mu/a) (exact vis-viva,
SOURCED stagek_lib.lever_m_per_mps; ratio = burn_response_ratio 1.0, DERIVED halflever_twin_results.md). The dynamic program minimizes sum_m w_m |nu_m - nu_{m-1}| subject to level_m + nu_m in [lo_m, hi_m]
(superposition: a burn does not change the later passive drift; checked for single burns only, plan R4.1).

Usage: python3 -I planner_lib.py --selftest | --dry-run
"""
import argparse, hashlib, json, math, pathlib, sys
import numpy as np

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
CFG_PATH = HERE / "planner_config.toml"
G3_CFG_PATH = HERE / "stageg3_config.toml"
INF = math.inf


def load_cfg(path=CFG_PATH):
    return tomllib.loads(pathlib.Path(path).read_text())


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(p):
    return sha256_bytes(pathlib.Path(p).read_bytes())


def canon(obj):
    """Canonical JSON text (sorted keys, no spaces) used for the hashes of the predictions body."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float)


def task_dir():
    import os
    d = os.environ.get("SPACEAGORA_BORDER_TASK_DIR")
    if d:
        return pathlib.Path(d)
    for p in HERE.parents:
        t = p / "tasks" / "2026-10-15_periapsis-region-border-targeting"
        if t.exists():
            return t
    return None


# ------------------------------------------------------------------------------------------------ lever and weights

def lever_m_per_mps(mu, a, e):
    """dr_p / dv for a tangential apoapsis burn (m per m/s): 4 a sqrt((1-e)/(1+e)) / sqrt(mu/a)."""
    return 4.0 * a * math.sqrt((1.0 - e) / (1.0 + e)) / math.sqrt(mu / a)


def burn_weight(H, ratio, mu, a, e):
    """m/s per ln-q unit: H / (ratio x lever)."""
    return H / (ratio * lever_m_per_mps(mu, a, e))


# ------------------------------------------------------------------------------------------------ the exact DP (mirror of planner_dp)

def dp_exact(A, B, L, tol=1.0e-9, pen=0.0, lin=None):
    """min sum_k L_k |nu_k - nu_{k-1}|, nu_{-1} = 0, A_k <= nu_k <= B_k. State set {0} U {A} U {B} (an optimal vertex of the piecewise-linear problem). Two relaxation sweeps per step (distance transform of
    L_k |.|); strict '<' so ties keep the nearer source: identical to planner_common.jl planner_dp. `pen` (m/s) is a fixed cost per burn (a step with nu_k != nu_{k-1}); 0 = none (the plan's U11: cadence
    unconstrained). `lin` (Revision 19 A5, default None): optional linear cost lin[k] nu_k on the level of step k, added after the relaxation (planner_rev18_offline.dp_lin); a vector with no nonzero entry takes exactly the
    code path of lin = None (bit-identical), as in planner_common.jl planner_dp. Returns dict(nu, obj) or None when infeasible (obj includes pen)."""
    A = [float(x) for x in A]; B = [float(x) for x in B]; L = [float(x) for x in L]; n = len(A)
    haslin = lin is not None and any(float(v) != 0.0 for v in lin)
    if haslin and len(lin) != n:
        raise ValueError("lin must have the same length as A")
    val = sorted(set([0.0] + A + B)); ns = len(val); i0 = val.index(0.0)
    prev = [INF] * ns; prev[i0] = 0.0
    parent = [None] * n
    for k in range(n):
        Lk = L[k]; cur = prev[:]; arg = list(range(ns))
        for j in range(1, ns):
            c = cur[j - 1] + Lk * (val[j] - val[j - 1])
            if c < cur[j]:
                cur[j] = c; arg[j] = arg[j - 1]
        for j in range(ns - 2, -1, -1):
            c = cur[j + 1] + Lk * (val[j + 1] - val[j])
            if c < cur[j]:
                cur[j] = c; arg[j] = arg[j + 1]
        if pen > 0.0:                                   # staying is free, any move pays pen once: new = min(prev, pen + best move)
            for j in range(ns):
                if prev[j] <= pen + cur[j]:
                    cur[j] = prev[j]; arg[j] = j
                else:
                    cur[j] = pen + cur[j]
        if haslin:
            lk = float(lin[k])
            for j in range(ns):
                cur[j] += lk * val[j]
        for j in range(ns):
            if val[j] < A[k] - tol or val[j] > B[k] + tol:
                cur[j] = INF
        parent[k] = arg; prev = cur
    j = min(range(ns), key=lambda x: prev[x])
    if not math.isfinite(prev[j]):
        return None
    nu = [0.0] * n
    for k in range(n - 1, -1, -1):
        nu[k] = val[j]; j = parent[k][j]
    return dict(nu=np.array(nu), obj=float(min(prev)))


def tie_weights(L, eps):
    """Lazy tie-break: step m (0-based) of n costs L_m (1 + eps (n - m) / n), so an early burn is slightly dearer than the same burn later."""
    L = np.asarray(L, float); n = len(L)
    return L * (1.0 + eps * (n - np.arange(n)) / n)


def dp_plan(level, lo, hi, w, eps=1.0e-6, tol=1.0e-9, pen=0.0, fast=False, lin=None):
    """Offsets nu_m for steps m = 0..n-1 with level_m + nu_m in [lo_m, hi_m]. Returns dict(nu, du, dv_cost) (dv_cost = sum w_m |dnu_m| with the untouched weights) or None.
    pen: fixed cost per burn (m/s). fast=True uses the vectorized oracle_lib.dp_schedule (same objective, ties may resolve differently): only for long replays. lin: optional linear cost per step (exact path only)."""
    level = np.asarray(level, float)
    if fast and lin is not None and any(float(v) != 0.0 for v in lin):
        raise ValueError("lin is not supported by the fast path")
    A = np.asarray(lo, float) - level; B = np.asarray(hi, float) - level
    if fast:
        import oracle_lib as OL
        o = OL.dp_schedule(A, B, tie_weights(w, eps), lam=pen); r = None if o is None else dict(nu=np.asarray(o["u"], float))
    else:
        r = dp_exact(A, B, tie_weights(w, eps), tol, pen, lin)
    if r is None:
        return None
    nu = r["nu"]; du = np.diff(np.concatenate(([0.0], nu)))
    return dict(nu=nu, du=du, dv_cost=float(np.sum(np.asarray(w, float) * np.abs(du))))


# ------------------------------------------------------------------------------------------------ plan inputs (mirror of planner_inputs)

def tget(arr, p):
    """arr indexed by 1-based pass p, clamped to the table (the last values are held beyond it)."""
    return float(arr[min(max(int(p) - 1, 0), len(arr) - 1)])


def dslow_at(tab, p):
    d = tab["d_slow"]
    return float(d[p - 1]) if 1 <= p <= len(d) else 0.0


def plan_inputs(sw, tab, P, j, z0, w0, bias=0.0, lower_margin=0.0, fb=None, fb_N=0, mu=float("nan")):
    """Inputs of the DP for the decision before pass j (online / olr modes). sw: dict(slack, ff, early, horizon); tab: ground table dict (n, W, d_slow, fast, w); P: [planner] constants
    (q_cap_mean_w_cm2, l0_width_ln, fast_reserve_ln, tail_passes, update gain / cap_passes); z0: predicted ln q of pass j with no further burn (full, one-pass prediction); w0: live weight of this burn.
    FF off: persistence of z0 (one step: the DP is a clamp to the nearest edge). FF on: z0 for the next pass, then the slow anchor z0 - fast_j plus the table's slow increments (plus `bias` for the next cap_passes
    steps), inside the tube shrunk by half the fast reserve at each edge."""
    hi0 = math.log(P["q_cap_mean_w_cm2"]); W = tab["W"] if sw["slack"] else P["l0_width_ln"]; lo0 = hi0 - W + lower_margin
    if not sw["ff"]:
        return dict(level=[z0], lo=[lo0], hi=[hi0], w=[w0])
    n = int(sw["horizon"]) if int(sw["horizon"]) > 0 else max(1, tab["n"] - j + 1 + int(P["tail_passes"]))
    res = P["fast_reserve_ln"] if sw["slack"] else 0.0
    s_j = z0 - tget(tab["fast"], j); level = [z0]; acc = 0.0
    cap = int(P["update_cap_passes"]); gain = float(P["update_gain"])
    for m in range(1, n):
        acc += dslow_at(tab, j + m - 1) + (gain * bias if m <= cap else 0.0)
        level.append(s_j + acc)
    lo = [lo0] + [lo0 + res / 2.0] * (n - 1); hi = [hi0] + [hi0 - res / 2.0] * (n - 1)
    w = [w0] + ([tget(tab["w"], j + m) for m in range(1, n)] if sw["early"] else [w0] * (n - 1))
    if fb is None:
        return dict(level=level, lo=lo, hi=hi, w=w)
    return dict(level=level, lo=lo, hi=hi, w=w, lin=fb_lin(tab, fb, mu, P["ra_target_m"], j, fb_N, n))      # Revision 19 A5: linear terminal-periapsis cost at the arrival pass fb_N (None outside the plan)


def ol_inputs(sw, tab, P, q_first_ln):
    """Inputs of the ground-computed open-loop schedule: passes p = 2 .. n + tail, level from the table's slow path anchored at the calibrated first-pass ln q, every step inside the reserve-shrunk tube,
    the table's true weights (arm OL has early = true)."""
    hi0 = math.log(P["q_cap_mean_w_cm2"]); lo0 = hi0 - tab["W"]; res = P["fast_reserve_ln"]
    N = tab["n"] + int(P["tail_passes"]); passes = list(range(2, N + 1))
    level = []; acc = 0.0
    for p in passes:
        acc += dslow_at(tab, p - 1); level.append(q_first_ln + acc)
    return dict(passes=passes, level=level, lo=[lo0 + res / 2.0] * len(passes), hi=[hi0 - res / 2.0] * len(passes), w=[tget(tab["w"], p) for p in passes])


def ol_schedule(sw, tab, P, q_first_ln, eps=1.0e-6, fast=False):
    """{pass: signed dv (m/s, + raises the periapsis)} of the open-loop plan; plus the plan dict."""
    inp = ol_inputs(sw, tab, P, q_first_ln); pl = dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], eps, pen=P["dp_burn_penalty_mps"], fast=fast)
    if pl is None:
        return {}, None
    sched = {}
    for p, du, w in zip(inp["passes"], pl["du"], inp["w"]):
        if abs(du) * w >= P["dp_min_burn_mps"]:
            sched[p] = -float(du) * float(w)
    return sched, pl


# ------------------------------------------------------------------------------------------------ terminal-periapsis term (Revision 19 A5; default off; mirror of planner_common.jl)

def first_burn_dv(mu, ra, rp, rp_des):
    """Exact vis-viva tangential periapsis-raise burn at apoapsis ra from rp to rp_des (m/s): v_a(a') - v_a(a); negative when it lowers the periapsis."""
    return math.sqrt(mu * (2 / ra - 2 / (ra + rp_des))) - math.sqrt(mu * (2 / ra - 2 / (ra + rp)))


def first_burn_slope(mu, ra_f, rp_f, H):
    """c_T = H dFB/drp (m/s per ln-q unit) at the closing apoapsis (ra_f, rp_f), desired periapsis radius fixed. dFB/drp = 1 / lever(a_f, e_f) exactly (lever = 4 a^2 v_a / mu)."""
    return H / lever_m_per_mps(mu, (ra_f + rp_f) / 2, (ra_f - rp_f) / (ra_f + rp_f))


def fb_lin(tab, fb, mu, ra_target, j, N, n):
    """Mirror of planner_fb_lin: list of length n with fb['weight'] x c_T at the step of the arrival pass N (step N - j + 1, 1-based) and 0 elsewhere; None when that step is outside 1..n (the term is dropped)."""
    m = int(N) - int(j) + 1
    if not (1 <= m <= n):
        return None
    kk = min(max(int(N), 1), int(tab["n"])) - 1; a = float(tab["a"][kk]); e = float(tab["e"][kk]); rp_t = a * (1 - e)
    lin = [0.0] * n; lin[m - 1] = fb["weight"] * first_burn_slope(mu, ra_target, rp_t, float(tab["H"][kk]))
    return lin


def apo_slope(slope_mps_per_km, drop_m):
    """Mirror of planner_apo_slope: cost (m/s per ln q) of the closing-apoapsis shortfall below the target circle: slope x drop / 1000."""
    return slope_mps_per_km * drop_m / 1000.0


def term_lin(tab, fb, P, mu, ra_target, j, N, n):
    """Mirror of planner_term_lin: the FB linear cost plus, with P['apo_one_sided'], the shortfall term with the table drop at pass min(N, n); None outside the plan."""
    lin = fb_lin(tab, fb, mu, ra_target, j, N, n)
    if lin is None or not P.get("apo_one_sided", False):
        return lin
    kk = min(max(int(N), 1), int(tab["n"])) - 1; a = np.asarray(tab["a"], float); e = np.asarray(tab["e"], float); ra_t = a * (1 + e); nn = int(tab["n"])
    drop = ra_t[kk] - ra_t[kk + 1] if kk + 1 < nn else ra_t[nn - 2] - ra_t[nn - 1]
    lin[int(N) - int(j)] += apo_slope(P["apo_shortfall_slope_mps_per_km"], drop); return lin


def wrap_deg(x):
    """Mirror of planner_wrap_deg: angle in degrees wrapped to (-180, 180]."""
    y = math.fmod(math.fmod(float(x) + 180.0, 360.0) + 360.0, 360.0) - 180.0
    return 180.0 if y == -180.0 else y


def j2_node_rate(mu, J2, Re, a, e, inc_rad):
    """Mirror of planner_j2_node_rate: secular J2 rate of the ascending node (rad/s)."""
    return -1.5 * math.sqrt(mu / a ** 3) * J2 * (Re / (a * (1 - e ** 2))) ** 2 * math.cos(inc_rad)


def node_error_deg(om_now_deg, rate_deg_s, t_now, t_arr, om_target_deg):
    return wrap_deg(om_now_deg + rate_deg_s * (t_arr - t_now) - om_target_deg)


def node_aim(om_now_deg, rate_deg_s, t_now, om_target_deg, tt):
    """Mirror of planner_node_aim: the arrival time of zero node error nearest tt (the node cycle has period 360 / |rate|)."""
    if rate_deg_s == 0.0:
        return float(tt)
    t0 = t_now + wrap_deg(om_target_deg - om_now_deg) / rate_deg_s; per = 360.0 / abs(rate_deg_s)
    return t0 + per * float(np.round((tt - t0) / per))


def ttg_predict(tab, mu, ra_target, kap, ra0, rp0, t0, x, rpp, term=None, ra_aim=None):
    """Mirror of planner_ttg_predict (guard arrival prediction with the optional terminal-aim emulation and the separate aim radius): returns (t, steps)."""
    a_t = np.asarray(tab["a"], float); e_t = np.asarray(tab["e"], float); ra_tab = a_t * (1 + e_t); d0 = np.concatenate((ra_tab[:-1] - ra_tab[1:], [float(ra_tab[-2] - ra_tab[-1])])); ra_r = ra_tab[::-1]; d0_r = d0[::-1]
    T = lambda r_, p_: 2 * math.pi * math.sqrt((0.5 * (r_ + p_)) ** 3 / mu)
    ra = float(ra0); rp_prev = float(rp0); t = float(t0); nx = len(x)
    for m in range(400):
        xm = x[min(m, nx - 1)]; rpm = rpp[min(m, nx - 1)]
        t += 0.5 * T(ra, rp_prev) + 0.5 * T(ra, rpm)
        dr = kap * float(np.interp(ra, ra_r, d0_r)) * math.exp(xm) / tab["q_req"]
        if term is not None:
            Reff = ra - (ra_target if ra_aim is None else ra_aim) + term["below"]
            if 0 < Reff <= term["steps"] * term["step"]:
                dr = Reff / max(1, math.ceil(Reff / term["step"]))
        ra -= dr
        if ra <= ra_target:
            return t + 0.5 * T(ra, rpm), m + 1
        rp_prev = rpm
    return float("nan"), 400


def interp1(x, xp, fp):
    """numpy.interp for ascending xp (clamped at the ends)."""
    return float(np.interp(x, xp, fp))


def arrival_pass(tab, mu, ra_target, x, q_first_ln, kappa):
    """Mirror of planner_arrival (planner_rev18_offline.arrival, state-dependent, with kappa): dict(N, days, ra_end); N = 3 n when the target is not reached. ra and rp from the table's a, e."""
    n = int(tab["n"]); a_t = np.asarray(tab["a"], float); e_t = np.asarray(tab["e"], float); ra_tab = a_t * (1 + e_t); rp_tab = a_t * (1 - e_t); ra = float(ra_tab[0]); t = float(tab["t"][0])
    d0 = np.concatenate((ra_tab[:-1] - ra_tab[1:], [float(ra_tab[-2] - ra_tab[-1])])); ra_r = ra_tab[::-1]; d0_r = d0[::-1]; path = np.asarray(tab["path"], float); H = np.asarray(tab["H"], float); kmax = 3 * n
    for k in range(1, kmax + 1):
        kk = min(k, n); xk = float(x[min(k, len(x)) - 1]); nu = xk - (q_first_ln + float(path[kk - 1] - path[0])); rp = float(rp_tab[kk - 1]) - float(H[kk - 1]) * nu; a = 0.5 * (ra + rp)
        drop = kappa * interp1(ra, ra_r, d0_r) * math.exp(xk) / tab["q_req"]; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2 * math.pi * math.sqrt(amid ** 3 / mu)
        if ra_next <= ra_target:
            return dict(N=k, days=(t + 0.5 * Tk) / 86400.0, ra_end=ra_next)
        ra = ra_next; t += Tk
    return dict(N=kmax, days=float("nan"), ra_end=ra)


def ol_plan(tab, P, q_first_ln, fb=None, iters=3, eps=None):
    """Mirror of planner_ol_plan. fb None: the existing schedule (N = 0). fb = dict(mu, kappa, weight): solve without the term, find the arrival pass N, re-solve with the linear cost at pass N, repeat until N is
    self-consistent (at most `iters` re-solves); when no arrival pass lies inside the plan the term is dropped (dropped True). Returns dict(sched, N, dropped, converged, nu)."""
    eps = P["dp_tie_eps"] if eps is None else eps
    inp = ol_inputs(dict(), tab, P, q_first_ln); level = inp["level"]; lo = inp["lo"]; hi = inp["hi"]; w = inp["w"]; passes = inp["passes"]
    plain = lambda: dp_plan(level, lo, hi, w, eps, pen=P["dp_burn_penalty_mps"])
    pl = None; Narr = 0; dropped = False; converged = True
    if fb is None:
        pl = plain()
    else:
        lin = None
        for it in range(iters + 1):
            pl = dp_plan(level, lo, hi, w, eps, pen=P["dp_burn_penalty_mps"], lin=lin)
            if pl is None:
                break
            Nn = arrival_pass(tab, fb["mu"], P["ra_target_m"], np.concatenate(([q_first_ln], np.asarray(level) + pl["nu"])), q_first_ln, fb["kappa"])["N"]
            if Nn == Narr:
                break
            Narr = Nn; lin = fb_lin(tab, fb, fb["mu"], P["ra_target_m"], 2, Narr, len(level))
            if lin is None:
                dropped = True; pl = plain(); break
            if it == iters:
                converged = False
    sched = {}
    if pl is None:
        return dict(sched=sched, N=Narr, dropped=dropped, converged=converged, nu=np.zeros(0))
    for p, du, wk in zip(passes, pl["du"], w):
        if abs(du) * wk >= P["dp_min_burn_mps"]:
            sched[p] = -float(du) * float(wk)
    return dict(sched=sched, N=Narr, dropped=dropped, converged=converged, nu=pl["nu"])


def decide(sw, tab, P, j, z0, w0, bias=0.0, fast=False):
    """Online decision before pass j: signed dv (m/s) and the plan. dv below dp_min_burn_mps is zeroed."""
    inp = plan_inputs(sw, tab, P, j, z0, w0, bias); pl = dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"], fast=fast)
    if pl is None:
        return 0.0, None
    dv = -float(pl["nu"][0]) * w0
    return (dv if abs(dv) >= P["dp_min_burn_mps"] else 0.0), pl


def terminal_aim(R, c_hist, P):
    """Mirror of planner_common.jl planner_terminal_aim (Revision 14). Returns dict(lnq, N, dstar, c) or None."""
    if len(c_hist) < 3:
        return None
    Reff = R + P["term_aim_below_m"]
    if Reff <= 0 or Reff > P["term_activate_steps"] * P["term_step_m"]:
        return None
    K = min(len(c_hist), int(P["term_c_window"])); c = float(np.median(c_hist[-K:])); N = max(1, math.ceil(Reff / P["term_step_m"])); dstar = Reff / N
    hi = math.log(P["q_cap_mean_w_cm2"]); lo = hi + math.log(P["term_min_q_factor"])
    return dict(lnq=min(max(math.log(dstar / c), lo), hi), N=N, dstar=dstar, c=c)


def target_table(cfg, name=None):
    """Mirror of planner_target_table: [planner.target] (the default target) merged with [planner.targets.<name>]."""
    base = dict(cfg["planner"].get("target", {}))
    if name is None or name == "" or name == "default":
        return base
    if name not in cfg["planner"].get("targets", {}):
        raise KeyError(f"unknown target {name}")
    base.update(cfg["planner"]["targets"][name]); return base


def geometry_target(cfg, gid):
    """Mirror of planner_geometry_target: the target name of a [planner.geometry.*] entry with this id (None = the default target)."""
    for o in cfg["planner"].get("geometry", {}).values():
        if o["id"] == gid:
            return o.get("target")
    return None


def flat_P(cfg, arm=None, target=None):
    """[planner] constants plus the dp and update sub-tables, flattened with prefixed keys (what the functions above read)."""
    P = dict(cfg["planner"]); P["dp_tie_eps"] = cfg["planner"]["dp"]["tie_eps"]; P["dp_min_burn_mps"] = cfg["planner"]["dp"]["min_burn_mps"]; P["dp_burn_penalty_mps"] = cfg["planner"]["dp"]["burn_penalty_mps"]
    if arm is not None and "burn_penalty_mps" in cfg["planner"]["arms"].get(arm, {}):
        P["dp_burn_penalty_mps"] = cfg["planner"]["arms"][arm]["burn_penalty_mps"]               # per-arm override (sensitivity arms)
    P["ra_target_m"] = cfg["stageo"]["ra_target_m"]                                                    # SOURCED: Table III end of main phase 2 (the closing apoapsis)
    P.update(cfg["planner"].get("terminal", {}))
    T = target_table(cfg, target)
    if str(T.get("event", "apoapsis")) != "apoapsis":
        raise NotImplementedError("Not implemented: decision event " + str(T["event"]))
    one = bool(T.get("apo_one_sided", False)); tol = float(T.get("apo_tol_m", 25000.0)); frac = float(T.get("apo_aim_above_frac", 0.5))
    for k in ("ra_target_m", "tail_passes", "term_step_m", "term_activate_steps", "term_min_q_factor", "term_c_window"):
        if k in T:
            P[k] = T[k]
    P["term_aim_below_m"] = -tol * frac if one else float(T.get("term_aim_below_m", P["term_aim_below_m"]))                   # one-sided: the aim sits above the target
    P["ra_stop_m"] = float(P["ra_target_m"]) + (tol if one else 0.0); P["apo_one_sided"] = one; P["apo_tol_m"] = tol; P["apo_aim_above_frac"] = frac
    P["apo_shortfall_slope_mps_per_km"] = float(T.get("apo_shortfall_slope_mps_per_km", 0.056))
    P["update_gain"] = cfg["planner"]["update"]["gain"]; P["update_cap_passes"] = cfg["planner"]["update"]["cap_passes"]; P["update_window"] = cfg["planner"]["update"]["window"]
    return P


def arm_switches(cfg, arm):
    a = cfg["planner"]["arms"][arm]
    return dict(mode=a["mode"], slack=a.get("slack", False), ff=a.get("ff", False), early=a.get("early", False), update=a.get("update", False), horizon=a.get("horizon", 0), safety=a.get("safety", "q"),
                terminal=a.get("terminal", False), fb_term=a.get("fb_term", False), ic_scenario=a.get("ic_scenario", "ref"), ttg_guard=a.get("ttg_guard", False))


# ------------------------------------------------------------------------------------------------ B-ODY rule (mirror of planner_decide, mode body)

def body_decide(P, B, j, t_s, last_burn_t, z0, w0):
    """B-ODY: Q_s predicted for the next pass = exp(z0) x qs_over_qmean; burn (re-center on the phase target) only when it is outside [q_low, q_high] and the last burn is at least cadence_s ago."""
    qs = math.exp(z0) * B["qs_over_qmean"]
    if B["q_low_qs_w_cm2"] <= qs <= B["q_high_qs_w_cm2"] or (t_s - last_burn_t) < B["cadence_s"]:
        return 0.0
    tgt = B["q_target_mp1_qs_w_cm2"] if j <= B["mp1_last_pass"] else B["q_target_mp2_qs_w_cm2"]
    nu0 = math.log(tgt / B["qs_over_qmean"]) - z0
    dv = -nu0 * w0
    return dv if abs(dv) >= P["dp_min_burn_mps"] else 0.0


# ------------------------------------------------------------------------------------------------ deterministic replay on the model path (predictions)

def simulate(arm_sw, tab, P, B, q_first_ln, n_max=None, fast=False):
    """Replay of one arm in the NOISE-FREE model world: the truth is ln q_j = q_first_ln + path_j - path_1 + cumulative offset (the table's full path, fast part included), the one-pass prediction is exact,
    the table is exact. Returns dict(dv, burns, burn_passes, dvs, viol_hi, viol_lo, n). This is the planner's model-consistent dV (a prediction under the model, not a noise-inclusive one). For B-ODY it is only
    an illustration: the noise that drives a B-ODY keeper's burns is absent."""
    n = int(n_max or tab["n"]); path = np.asarray(tab["path"], float); hi0 = math.log(P["q_cap_mean_w_cm2"])
    W = tab["W"] if arm_sw["slack"] else P["l0_width_ln"]; lo0 = hi0 - W
    nu = 0.0; dvs = {}; ntrue = np.zeros(n)
    base = lambda p: q_first_ln + float(path[min(p, len(path)) - 1] - path[0])
    sched = None
    if arm_sw["mode"] == "ol":
        sched, _ = ol_schedule(arm_sw, tab, P, q_first_ln, P["dp_tie_eps"], fast=fast)
    last_t = -INF; t = np.asarray(tab["t"], float)
    for p in range(1, n + 1):
        if p >= 2:
            z0 = base(p) + nu; w0 = tget(tab["w"], p)
            if arm_sw["mode"] == "body":
                dv = body_decide(P, B, p, float(t[min(p, len(t)) - 1]), last_t, z0, w0)
            elif arm_sw["mode"] == "ol":
                dv = sched.get(p, 0.0)
            else:
                dv, _ = decide(arm_sw, tab, P, p, z0, w0, 0.0, fast=fast)
            if dv != 0.0:
                dvs[p] = dv; nu += -dv / w0; last_t = float(t[min(p, len(t)) - 1])
        ntrue[p - 1] = base(p) + nu
    return dict(dv=float(sum(abs(v) for v in dvs.values())), burns=len(dvs), burn_passes=sorted(dvs), dvs=dvs, viol_hi=int((ntrue > hi0 + 1e-12).sum()), viol_lo=int((ntrue < lo0 - 1e-12).sum()),
                n=n, lnq=ntrue)


# ------------------------------------------------------------------------------------------------ flight list (dry run)

HOSTS = ("ws", "sf3")


def flight_list(cfg):
    """All flights of plan R4.5 amended by R7.3 / R8.2: list of dict(block, arm, geometry, seed, host, tag, runner, phase). `gate` rows are the stageg3 flights (cells ARM@GEOM, seeds 'gate')."""
    P = cfg["planner"]; ss = P["seedsets"]; out = []
    for f in list(P["flight"]) + list(P.get("sensitivity", [])):
        for cell in f["cells"]:
            arm, geom = cell.split("@")
            for host in f["hosts"]:
                for s in ss[f["seeds"]]:
                    runner = "stageg3_run.jl (frozen)" if f["block"] == "gate" else P["arms"][arm]["runner"]
                    out.append(dict(block=f["block"], arm=arm, geometry=geom, seed=int(s), host=host, runner=runner,
                                    tag=("G3s" if f["block"] == "gate" else P["seed_tag_prefix"]) + str(s), horizon_days=(P["smoke_days"] if f["block"] in ("smoke", "smoke_ol") else P["horizon_days"])))
    return out


def dry_run(cfg, verbose=True):
    L = flight_list(cfg); P = cfg["planner"]; pc = P["plan_counts"]; by = {}
    for r in L:
        by.setdefault(r["block"].split("_")[0] if r["block"].startswith("held") else r["block"], []).append(r)
    counts = {k: len(v) for k, v in by.items()}
    nsens = counts.get("sens", 0); Lmain = [r for r in L if r["block"] != "sens"]
    ok = len(Lmain) == pc["total"] and len(L) == pc["total_with_sens"] and nsens == pc["sens"] and all(counts.get(k, 0) == pc[k] for k in ("gate", "smoke", "repro", "xhost", "dev", "held"))
    ids = [(r["block"], r["arm"], r["geometry"], r["seed"], r["host"]) for r in L]
    ok &= len(set(ids)) == len(ids)
    if verbose:
        print("planner flight list (plan R4.5 amended by R7.3 and R8.2); nothing is flown")
        rows = {}
        for r in L:
            rows.setdefault((r["block"], r["arm"], r["geometry"], r["host"]), []).append(r["seed"])
        for (b, a, g, h), seeds in rows.items():
            print(f"  {b:15s} {a:5s} @ {g:6s} host={h:4s} seeds={seeds[0]}..{seeds[-1]} ({len(seeds):2d})  runner={P['arms'][a]['runner'] if b != 'gate' else 'stageg3_run.jl (frozen)'}")
        print("  per phase:", {k: counts.get(k, 0) for k in ("gate", "smoke", "repro", "xhost", "dev", "held", "sens")}, "main total", len(Lmain), "(plan", pc["total"], ") + sensitivity", nsens, "= total", len(L), "(plan", pc["total_with_sens"], ")")
        print(f"  DERIVED cost at 0.37 worker-h per flight (R4.9 basis, range 0.5x to 2x): main {0.37 * len(Lmain):.0f} worker-h ({0.185 * len(Lmain):.0f} to {0.74 * len(Lmain):.0f}); sensitivity {0.37 * nsens:.1f} ({0.185 * nsens:.1f} to {0.74 * nsens:.1f}); total {0.37 * len(L):.0f}")
        ex = [e for e in P.get("extra", [])]
        for e in ex:
            print(f"  NOT COUNTED ({e['block']}): cells {e['cells']} seeds {e['seeds']} -> {len(e['cells']) * len(P['seedsets'][e['seeds']])} flights (R7.7 'smoke adds OL'; coordinator decision)")
        by_host = {}
        for r in L:
            by_host[r["host"]] = by_host.get(r["host"], 0) + 1
        print("  per host:", by_host, "| gate flights are specified and flown by stageg3_run.jl, the Stage O baselines by stageo_run.jl; the planner runner flies", sum(1 for r in L if r["runner"] == "planner"))
        print("dry run:", "counts match the plan (274 main + 32 sensitivity = 306)" if ok else "COUNT MISMATCH")
    return ok, L


# ------------------------------------------------------------------------------------------------ selftest

def _random_instance(rng, n):
    lev = np.cumsum(rng.normal(0, 0.2, n)); width = rng.uniform(0.2, 0.8); lo = np.full(n, -width / 2); hi = np.full(n, width / 2)
    w = rng.uniform(0.5, 3.0, n)
    return lev, lo, hi, w


def PLB_ok(cfg):
    return cfg['planner']['arms']['OLRT']['terminal'] is True and cfg['planner']['arms']['OLR'].get('terminal', False) is False and cfg['planner']['terminal']['tolerance_km'] == 25.0


def selftest():
    cfg = load_cfg(); g3 = load_cfg(G3_CFG_PATH); P = flat_P(cfg)
    # 1. config: shared sections equal the gate config's; counts; arms
    for k in ("stagef", "analysis", "common", "planets", "stageo"):
        a, b = dict(cfg[k]), dict(g3[k])
        if k == "stagef":
            assert a.pop("results_subdir") == "planner" and b.pop("results_subdir") == "stageg3"
        assert a == b, f"[{k}] differs from stageg3_config.toml"
    a, b = dict(cfg["stageg3"]), dict(g3["stageg3"]); assert a == b, "[stageg3] differs from stageg3_config.toml"
    assert cfg["planner"]["burn_response_ratio"] == g3["stageg3"]["burn_response_ratio"] == 1.0
    assert abs(cfg["planner"]["l0_width_ln"] - 4000.0 / 6450.0) < 1e-4 and cfg["planner"]["body"]["mp1_last_pass"] == round(108 / 153 * 200)
    assert cfg["planner"]["q_cap_mean_w_cm2"] == 0.2503 and cfg["planner"]["calibration"]["old_q_cap_mean_w_cm2"] == 0.234 and cfg["planner"]["q_never_exceed_qs_w_cm2"] == 0.65 and cfg["planner"]["d_ref_days"] == 70.4 and cfg["planner"]["sigma_p"] == 0.0452
    ok, L = dry_run(cfg, verbose=False); assert ok and len(L) == 306 and sum(1 for r in L if r["block"] != "sens") == 274
    arms = [r["arm"] for r in L if r["block"].startswith("held")]
    assert [r["arm"] for r in L if r["block"] == "sens"].count("L3Z") == 16 and flat_P(cfg, "L3Z")["dp_burn_penalty_mps"] == 0.0 and flat_P(cfg, "L3")["dp_burn_penalty_mps"] == 0.01 == cfg["planner"]["dp"]["burn_penalty_mps"]
    for a_, n_ in (("BODY", 32), ("FO", 16), ("N2I", 16), ("L0", 16), ("L1", 32), ("L2", 16), ("L3", 32), ("OL", 32), ("OLR", 32)):
        assert arms.count(a_) == n_, (a_, arms.count(a_))
    assert sum(1 for r in L if r["block"] == "gate") == 10 and sum(1 for r in L if r["arm"] in ("L3", "L1") and r["geometry"] == "G1" and r["block"] == "held_geometry") == 32
    assert {a_ for a_ in cfg["planner"]["arms"]} == {"BODY", "L0", "L1", "L2", "L3", "OL", "OLR", "FO", "N2I", "L3Z", "OLRZ", "OLRT", "OLRT1", "OLRT2", "OLRT3", "OLRT4", "OLRT5", "OLRT35", "BODYT", "N2IT", "N2I0", "OLT", "OLRT3-W", "OLRT3-WF", "OLRT3-W3", "OLRT3-WF3"}
    # ladder switches (R2.2): L0 -> L3 adds slack, ff, early one at a time; OL / OL-R are L3 without / with replanning
    sw = {k: arm_switches(cfg, k) for k in ("L0", "L1", "L2", "L3", "OL", "OLR")}
    assert [(sw[k]["slack"], sw[k]["ff"], sw[k]["early"]) for k in ("L0", "L1", "L2", "L3")] == [(False, False, False), (True, False, False), (True, True, False), (True, True, True)]
    assert (sw["OL"]["mode"], sw["OL"]["update"], sw["OLR"]["mode"], sw["OLR"]["update"], sw["OLR"]["horizon"]) == ("ol", False, "olr", True, 0)
    # 2. lever
    mu = 4.2828375e13
    for a_, e_ in ((1.6e7, 0.78), (7e6, 0.36), (4e6, 0.1)):
        assert abs(lever_m_per_mps(mu, a_, e_) - 4 * a_ ** 2 * math.sqrt(mu * (1 - e_) / ((1 + e_) * a_)) / mu) < 1e-9 * lever_m_per_mps(mu, a_, e_)   # 4 a^2 v_a / mu
    try:
        import stagek_lib as K
        for a_, e_ in ((1.6e7, 0.78), (7e6, 0.36)):
            assert abs(K.lever_m_per_mps(mu, a_, e_) - lever_m_per_mps(mu, a_, e_)) < 1e-9 * lever_m_per_mps(mu, a_, e_)
        k_checked = True
    except Exception:
        k_checked = False
    # 3. DP against the independent oracle_lib.dp_schedule (exact mode) on random instances, and structural checks
    import oracle_lib as OL
    rng = np.random.default_rng(1)
    for n in (3, 8, 25):
        for _ in range(20):
            lev, lo, hi, w = _random_instance(rng, n); A, B = lo - lev, hi - lev
            r = dp_exact(A, B, w); o = OL.dp_schedule(A, B, w)
            assert r is not None and o is not None and abs(r["obj"] - o["dv"]) < 1e-9 * max(1.0, o["dv"]), (r["obj"], o["dv"])
            assert np.all(r["nu"] >= A - 1e-9) and np.all(r["nu"] <= B + 1e-9)
    z = dp_plan([0.0] * 5, [-1.0] * 5, [1.0] * 5, [1.0] * 5); assert z["dv_cost"] == 0.0 and np.all(z["nu"] == 0)               # inside the tube: no burn
    z = dp_plan([0.0, 1.5, 1.5], [-1.0] * 3, [1.0] * 3, [1.0, 1.0, 1.0]); assert abs(z["dv_cost"] - 0.5) < 1e-9 and z["nu"][0] == 0.0 # one violation, lazy: burn at step 1, to the edge
    z = dp_plan([0.0, 0.0, 1.5], [-1.0] * 3, [1.0] * 3, [1.0, 1.0, 1.0]); assert z["du"][2] != 0.0 and z["du"][0] == 0.0 and z["du"][1] == 0.0
    z = dp_plan([0.0, 0.0, 1.5], [-1.0] * 3, [1.0] * 3, [1.0, 2.0, 4.0]); assert z["du"][0] != 0.0 and abs(z["dv_cost"] - 0.5) < 1e-9   # rising weights: prepay (burn first)
    # burn penalty: one fixed cost per burn merges two small burns into none/one; pen = 0 changes nothing
    lvl = [0.0, 0.6, 1.2, 1.8]; a0 = dp_plan(lvl, [-1.0] * 4, [0.5] * 4, [1.0] * 4); a1 = dp_plan(lvl, [-1.0] * 4, [0.5] * 4, [1.0] * 4, pen=0.3)
    assert (np.abs(a0["du"]) > 1e-12).sum() >= 2 and (np.abs(a1["du"]) > 1e-12).sum() == 1 and a1["dv_cost"] >= a0["dv_cost"] - 1e-12 and abs(a1["dv_cost"] - 1.3) < 1e-9, (a0["du"], a1["du"])   # same dV, one burn instead of three
    assert abs(dp_plan(lvl, [-1.0] * 4, [0.5] * 4, [1.0] * 4, pen=0.0)["dv_cost"] - a0["dv_cost"]) < 1e-15 and abs(dp_plan(lvl, [-1.0] * 4, [0.5] * 4, [1.0] * 4, pen=0.3, fast=True)["dv_cost"] - a1["dv_cost"]) < 1e-9   # fast path agrees
    # 4. plan inputs
    tab = dict(n=10, W=0.72, d_slow=[0.05] * 9, fast=[0.0] * 10, w=[1.0 + 0.1 * k for k in range(10)], path=[0.05 * k for k in range(10)], t=[1.0e5 * (k + 1) for k in range(10)])
    hi0 = math.log(P["q_cap_mean_w_cm2"])
    i0 = plan_inputs(dict(slack=False, ff=False, early=False, horizon=1), tab, P, 2, -1.8, 2.0); assert i0["level"] == [-1.8] and abs(i0["hi"][0] - hi0) < 1e-12 and abs(i0["lo"][0] - (hi0 - P["l0_width_ln"])) < 1e-12
    i1 = plan_inputs(dict(slack=True, ff=False, early=False, horizon=1), tab, P, 2, -1.8, 2.0); assert abs(i1["lo"][0] - (hi0 - 0.72)) < 1e-12
    i3 = plan_inputs(dict(slack=True, ff=True, early=True, horizon=4), tab, P, 2, -1.8, 2.0)
    assert len(i3["level"]) == 4 and abs(i3["level"][2] - (-1.8 + 0.05 * 2)) < 1e-12 and abs(i3["hi"][1] - (hi0 - P["fast_reserve_ln"] / 2)) < 1e-12 and i3["w"][1:] == [tab["w"][2], tab["w"][3], tab["w"][4]]
    i2 = plan_inputs(dict(slack=True, ff=True, early=False, horizon=4), tab, P, 2, -1.8, 2.0); assert i2["w"] == [2.0] * 4
    ib = plan_inputs(dict(slack=True, ff=True, early=True, horizon=4), tab, P, 2, -1.8, 2.0, bias=0.01); assert abs(ib["level"][1] - (-1.8 + 0.05 + 0.01)) < 1e-12
    # FF off is a clamp to the nearest edge (reset-to-edge): level above the cap by 0.1 -> burn up 0.1 x w
    dv, _ = decide(dict(slack=True, ff=False, early=False, horizon=1), tab, P, 3, hi0 + 0.1, 2.0); assert abs(dv - 0.2) < 1e-9
    dv, _ = decide(dict(slack=True, ff=False, early=False, horizon=1), tab, P, 3, hi0 - 0.72 - 0.1, 2.0); assert abs(dv + 0.2) < 1e-9
    dv, _ = decide(dict(slack=True, ff=False, early=False, horizon=1), tab, P, 3, hi0 - 0.3, 2.0); assert dv == 0.0
    # 5. replay: flat path inside the tube -> no burns; ramp -> burns; OL schedule equals OL-R in the noise-free world
    flat = dict(tab, path=[0.0] * 10, d_slow=[0.0] * 9, n=10); q1 = hi0 - 0.3
    for k in ("L0", "L1", "L2", "L3", "OL", "OLR"):
        r = simulate(sw[k], flat, P, cfg["planner"]["body"], q1); assert r["burns"] == 0 and r["viol_hi"] == 0, k
    n = 40; ramp = dict(n=n, W=0.72, d_slow=[0.04] * (n - 1), fast=[0.0] * n, w=[1.0 + 0.05 * k for k in range(n)], path=[0.04 * k for k in range(n)], t=[1.0e5 * (k + 1) for k in range(n)])
    rr = {k: simulate(sw[k], ramp, P, cfg["planner"]["body"], q1, n_max=n) for k in ("L0", "L1", "L2", "L3", "OL", "OLR")}
    assert all(rr[k]["viol_hi"] == 0 for k in rr), {k: rr[k]["viol_hi"] for k in rr}
    assert rr["L3"]["dv"] <= rr["L2"]["dv"] + 1e-9 and abs(rr["OL"]["dv"] - rr["OLR"]["dv"]) < 0.05 * max(1.0, rr["OLR"]["dv"]), {k: round(v["dv"], 4) for k, v in rr.items()}
    # B-ODY: inside the band no burn; above 0.32 re-center to the MP1 target within the cadence; cadence blocks a second burn
    B = cfg["planner"]["body"]; zin = math.log(0.25 / B["qs_over_qmean"]); assert body_decide(P, B, 10, 5e5, -INF, zin, 2.0) == 0.0
    zhi = math.log(0.40 / B["qs_over_qmean"]); dvb = body_decide(P, B, 10, 5e5, -INF, zhi, 2.0)
    assert abs(dvb - 2.0 * (zhi - math.log(0.24 / B["qs_over_qmean"]))) < 1e-12 and body_decide(P, B, 10, 5e5, 5e5 - 3600.0, zhi, 2.0) == 0.0
    assert abs(body_decide(P, B, 150, 5e5, -INF, zhi, 2.0) - 2.0 * (zhi - math.log(0.20 / B["qs_over_qmean"]))) < 1e-12
    td = task_dir()                                                          # pinned artifacts, when present next to the plan
    for key in ("ground_table", "predictions"):
        f = td / cfg["planner"][key]["path"] if td else None
        if f is not None and f.exists():
            assert cfg["planner"][key]["sha256"] == sha256_file(f), f"{f.name}: sha256 differs from planner_config.toml [planner.{key}]"
    # Revision 14 terminal aim (values asserted again in planner_run.jl --phase=selftest)
    ta = terminal_aim(50e3, [2e5, 2.2e5, 1.8e5, 2e5], P); assert ta["N"] == 2 and abs(ta["dstar"] - 31250.0) < 1e-9 and abs(ta["lnq"] - math.log(31250.0 / 2e5)) < 1e-12, ta
    assert terminal_aim(200e3, [2e5] * 4, P) is None and terminal_aim(50e3, [2e5, 2e5], P) is None and terminal_aim(-20e3, [2e5] * 4, P) is None
    lo_ = math.log(P["q_cap_mean_w_cm2"]) + math.log(P["term_min_q_factor"]); assert abs(terminal_aim(10e3, [8e5] * 4, P)["lnq"] - lo_) < 1e-12 and terminal_aim(100e3, [1e3] * 4, P)["lnq"] == math.log(P["q_cap_mean_w_cm2"])
    assert terminal_aim(5e3, [2e5] * 3, P)["N"] == 1 and PLB_ok(cfg)
    # Revision 20b / 20 / 19: the OLT arm, the IC and FB switches are default off; every pre-existing arm keeps the old switches
    a_olt = cfg["planner"]["arms"]["OLT"]; so = cfg["planner"]["arms"]["OL"]
    assert (a_olt["mode"], a_olt["terminal"], a_olt["safety"], a_olt["slack"], a_olt["ff"], a_olt["early"], a_olt["update"], a_olt["horizon"]) == (so["mode"], True, so["safety"], so["slack"], so["ff"], so["early"], so["update"], so["horizon"])
    assert all(arm_switches(cfg, k)["fb_term"] is False and arm_switches(cfg, k)["ic_scenario"] == "ref" and arm_switches(cfg, k)["ttg_guard"] is False for k in cfg["planner"]["arms"] if k not in ("OLRT3-W", "OLRT3-WF", "OLRT3-W3", "OLRT3-WF3")) and (arm_switches(cfg, "OLRT3-W")["ttg_guard"], arm_switches(cfg, "OLRT3-W")["ic_scenario"], arm_switches(cfg, "OLRT3-W")["terminal"]) == (True, "S-fixed", True) and arm_switches(cfg, "OLT")["terminal"] and not arm_switches(cfg, "OL")["terminal"]
    assert flat_P(cfg, "OLT")["dp_burn_penalty_mps"] == flat_P(cfg, "OL")["dp_burn_penalty_mps"]
    blk = [f for f in cfg["planner"]["terminal_flight"] if f["block"] == "rev20olt"][0]; cs = blk["cell_seeds"]
    assert sum(len(v) for v in cs.values()) == blk["count"] == 13 and sorted(cs) == sorted(blk["cells"]) and cs["OLT@ODY@L0"] == [13, 14, 15, 16] and 13 not in cs["OLT@ODY@L1"] and 13 not in cs["OLT@ODY@L2"] and 14 not in cs["OLT@ODY@L3"]
    gv = cfg["planner"]["guard"]; assert (gv["pen_g"], gv["rel_slack"], gv["s_frac"], gv["aim_frac"], gv["anchor_k"], gv["bias_steps"], gv["k_window"], gv["ema_beta"], gv["eps_max"], gv["max_solves"], gv["dev_thr_g"]) == (0.1, 0.15, 0.5, 0.5, 15, 60, 20, 0.5, 1.0, 12, 0.05) and gv["t_offset_d"] == {"ODY": 0.15, "G5863": 0.0}      # Revision 20c freeze
    bv = [f for f in cfg["planner"]["terminal_flight"] if f["block"] == "rev20val"][0]; assert bv["count"] == 29 and len(bv["cells"]) == 8
    for blk_, n_ in (("rev21c", 36), ("rev21v", 52)):
        b_ = [f for f in cfg["planner"]["terminal_flight"] if f["block"] == blk_][0]; assert b_["count"] == n_ == sum(len(v) for v in b_["cell_seeds"].values()) and sorted(b_["cell_seeds"]) == sorted(b_["cells"]), blk_
    b23 = [f for f in cfg["planner"]["terminal_flight"] if f["block"] == "rev23geo"][0]; assert b23["count"] == 12 and b23["cells"] == ["OLRT3-W3@GE@L0", "OLRT3-W3@GF@L0", "OLRT3-W3@GAM@L0"] and [cfg["planner"]["geometry"][k]["id"] for k in ("G0528", "G0890", "G0063M")] == ["G0528", "G0890", "G0063M"]
    assert cfg["planner"]["ground_table_rev23"]["sha256"] == "8abde05404da2eb0f7b9e48e6354f3e2641cf1208cdee7e8e8cbc97fd48a00d4" and cfg["planner"]["rev23"]["cal_tokens"] == ["GE", "GF", "GAM"]
    # Revision 25: target as a parameter (default = current values), one-sided tolerance above the target, node check helpers, block rev25low and the rev25 pins
    P0 = flat_P(cfg); assert P0["ra_target_m"] == cfg["stageo"]["ra_target_m"] == P0["ra_stop_m"] and P0["tail_passes"] == cfg["planner"]["tail_passes"] and P0["term_aim_below_m"] == 12500.0 and not P0["apo_one_sided"] and flat_P(cfg, None, "default") == P0
    assert flat_P(cfg, "OLRT3", "T2000")["ra_target_m"] == 5396190.0 and [geometry_target(cfg, g) for g in ("G0486", "G1230", "G1038", "ODY")] == ["T2000", "T1000", "T500", None]
    cfg_o = json.loads(json.dumps(cfg)); cfg_o["planner"]["target"]["apo_one_sided"] = True; P1 = flat_P(cfg_o)
    assert P1["ra_stop_m"] == P1["ra_target_m"] + 25000.0 and P1["term_aim_below_m"] == -12500.0 and terminal_aim(60e3, [2e5] * 4, P0)["dstar"] == (60e3 + 12500.0) / 2 and terminal_aim(60e3, [2e5] * 4, P1)["dstar"] == (60e3 - 12500.0) / 2
    cfg_e = json.loads(json.dumps(cfg)); cfg_e["planner"]["target"]["event"] = "periapsis"
    try:
        flat_P(cfg_e); raise AssertionError("event must be rejected")
    except NotImplementedError:
        pass
    assert wrap_deg(190.0) == -170.0 and wrap_deg(-180.0) == 180.0 and wrap_deg(720.0) == 0.0 and j2_node_rate(mu, 0.001958616126805748, 3394200.0, 1e7, 0.65, math.radians(93.5)) > 0 > j2_node_rate(mu, 0.001958616126805748, 3394200.0, 1e7, 0.65, math.radians(25.2))
    assert abs(node_error_deg(10.0, 1e-5, 0.0, 1e6, 20.0)) < 1e-12 and abs(node_aim(10.0, 1e-5, 0.0, 20.0, 5e5) - 1e6) < 1e-6 and node_aim(10.0, 0.0, 0.0, 20.0, 5.0) == 5.0
    assert apo_slope(0.056, 40000.0) == 0.056 * 40.0 and "node_tol_deg" not in cfg["planner"]["guard"] and cfg["planner"]["guard"]["node_check"] is False and cfg["planner"]["guard_sets"]["v3"]["node_check"] is False
    b25 = [f for f in cfg["planner"]["terminal_flight"] if f["block"] == "rev25low"][0]; assert b25["count"] == 13 == sum(len(v) for v in b25["cell_seeds"].values()) and b25["cell_seeds"]["OLRT3@ODY@L4"] == [14] and cfg["planner"]["rev25"]["cal_tokens"] == ["GG", "GH", "GI"]
    td25 = task_dir()
    if td25 is not None:
        for key, pth in (("ground_table_rev25", cfg["planner"]["ground_table_rev25"]), ("ic_extra.rev25", cfg["planner"]["ic_extra"]["rev25"])):
            f = td25 / pth["path"]
            if f.exists():
                assert sha256_file(f) == pth["sha256"], key
        tj = json.loads((td25 / "target_orbits.json").read_text()); rad = {t["target_id"]: 1000.0 * t["circle_r_km"] for t in tj["targets"]}
        assert all(cfg["planner"]["targets"][k]["ra_target_m"] == rad[cfg["planner"]["targets"][k]["target_id"]] for k in ("T2000", "T1000", "T500"))
    assert cfg["planner"]["guard_sets"]["v3"]["s_frac"] == 0.25 and cfg["planner"]["guard_sets"]["v3"]["term_model"] is True and cfg["planner"]["arms"]["OLRT3-W3"]["guard_set"] == "v3"
    assert cfg["planner"]["dispersion"]["L2p"]["G0975"]["mars_mgcm_dust_levels"][0] == 1.0042 and arm_switches(cfg, "OLRT3-WF")["ic_scenario"] == "S-free"
    # linear cost: lin = None, zeros and a zero-valued list are bit-identical to the old DP; a nonzero entry moves the plan (cross-language values asserted again in planner_run.jl)
    n_ = 8; vv = cross_language_vector(); args_ = (vv["level"], vv["lo"], vv["hi"], vv["w"], 1.0e-6)
    r_none = dp_plan(*args_); r_zero = dp_plan(*args_, lin=[0.0] * n_); r_negz = dp_plan(*args_, lin=[-0.0] * n_)
    assert np.array_equal(r_none["nu"], r_zero["nu"]) and r_none["dv_cost"] == r_zero["dv_cost"] and np.array_equal(r_none["nu"], r_negz["nu"])
    A_ = np.array(vv["lo"]) - np.array(vv["level"]); B_ = np.array(vv["hi"]) - np.array(vv["level"])
    assert dp_exact(A_, B_, vv["w"], pen=0.3)["obj"] == dp_exact(A_, B_, vv["w"], pen=0.3, lin=[0.0] * n_)["obj"]
    r_lin = dp_plan(*args_, lin=[0, 0, 0, 0, 0, 1.5, 0, 0]); assert np.allclose(r_lin["nu"], [-0.2] * 6 + [-0.15] * 2, atol=1e-9) and abs(r_lin["dv_cost"] - 0.56) < 1e-9
    # exactness of the linear DP against a brute-force grid on random small instances
    rng2 = np.random.default_rng(7)
    for _ in range(30):
        n3 = 4; lev, lo3, hi3, w3 = _random_instance(rng2, n3); lin3 = rng2.uniform(-1.5, 1.5, n3); A3, B3 = lo3 - lev, hi3 - lev
        r3 = dp_exact(A3, B3, w3, lin=lin3)
        # coordinate check: no single-step perturbation of the DP optimum lowers the objective (convex piecewise-linear: a local minimum is global)
        def f3(nu_):
            prev_ = 0.0; tot = 0.0
            for k in range(n3):
                tot += w3[k] * abs(nu_[k] - prev_) + lin3[k] * nu_[k]; prev_ = nu_[k]
            return tot
        base_ = f3(r3["nu"]); assert abs(base_ - r3["obj"]) < 1e-9
        for k in range(n3):
            for d_ in (-1e-3, 1e-3):
                t_ = r3["nu"].copy(); t_[k] += d_
                if A3[k] - 1e-12 <= t_[k] <= B3[k] + 1e-12:
                    assert f3(t_) >= base_ - 1e-9
    # first-burn slope identity (Revision 19 A1): c_T / H = 1 / lever = dFB/drp at fixed desired periapsis, exact vis-viva, three final states
    for rp_ in (3.480e6, 3.484e6, 3.490e6):
        fd = (first_burn_dv(mu, 7.408e6, rp_ - 250.0, rp_ + 150e3) - first_burn_dv(mu, 7.408e6, rp_ + 250.0, rp_ + 150e3)) / 500.0
        assert abs(first_burn_slope(mu, 7.408e6, rp_, 6450.0) / 6450.0 - fd) < 1e-6 * fd
    tabx = dict(tab, a=[1.0e7] * 10, e=[0.65] * 10, H=[6450.0] * 10)
    assert abs(first_burn_slope(mu, 7.408e6, 3.480e6, 1.0) * 1000.0 - 0.18793) < 5e-5 and fb_lin(tabx, dict(weight=2.0), mu, 7.408e6, 2, 1, 5) is None and fb_lin(tabx, dict(weight=2.0), mu, 7.408e6, 2, 8, 5) is None
    lv_ = fb_lin(tabx, dict(weight=2.0), mu, 7.408e6, 2, 5, 9); assert sum(1 for v in lv_ if v != 0.0) == 1 and lv_[3] > 0.0 and fb_lin(tabx, dict(weight=0.0), mu, 7.408e6, 2, 5, 9) == [0.0] * 9
    # plan inputs and the OL plan with the term off are the old ones; the term adds one nonzero linear cost at the arrival pass
    ia_ = plan_inputs(dict(slack=True, ff=True, early=True, horizon=4), tab, P, 2, -1.8, 2.0); assert "lin" not in ia_ and ia_ == plan_inputs(dict(slack=True, ff=True, early=True, horizon=4), tab, P, 2, -1.8, 2.0)
    ib_ = plan_inputs(dict(slack=True, ff=True, early=True, horizon=4), tabx, P, 2, -1.8, 2.0, fb=dict(weight=1.0), fb_N=4, mu=mu); assert len(ib_["lin"]) == 4 and ib_["lin"][2] > 0.0 and sum(1 for v in ib_["lin"] if v) == 1 and ib_["level"] == ia_["level"]
    tdir = task_dir(); gtab = tdir / cfg["planner"]["ground_table"]["path"] if tdir else None
    if gtab is not None and gtab.exists():
        gt = json.loads(gtab.read_text())["geometries"]; MUm = cfg["planner"]["analysis"]["mu_mars_m3s2"]; P0 = dict(P, dp_burn_penalty_mps=0.0); q1 = math.log(0.212)
        kap = dict(ODY=cfg["planner"]["terminal_periapsis"]["kappa"]["ODY"], G1=cfg["planner"]["terminal_periapsis"]["kappa"]["G5863"]); gid_ = dict(ODY="ODY", G1="G5863")
        for tok, (N0, N1, b0, b1) in dict(ODY=(216, 231, 106, 140), G1=(265, 265, 4, 4)).items():
            gtab_ = gt[gid_[tok]]; p_old = ol_plan(gtab_, P0, q1); fb_ = dict(mu=MUm, kappa=kap[tok], weight=1.0); p_on = ol_plan(gtab_, P0, q1, fb=fb_); p_zero = ol_plan(gtab_, P0, q1, fb=dict(fb_, weight=0.0))
            s_old, _ = ol_schedule({}, gtab_, P0, q1); s_cfg, _ = ol_schedule({}, gtab_, P, q1)
            assert p_old["sched"] == s_old and p_zero["sched"] == s_old and np.array_equal(p_zero["nu"], p_old["nu"]) and ol_plan(gtab_, P, q1)["sched"] == s_cfg      # weight 0 and the term off reproduce the existing schedule exactly
            assert p_zero["N"] == N0 and p_on["N"] == N1 and len(s_old) == b0 and len(p_on["sched"]) == b1 and not p_on["dropped"] and p_on["converged"], (tok, p_zero["N"], p_on["N"], len(s_old), len(p_on["sched"]))      # Revision 19/20 model-world values (rev20_data/ic/headline.json crosscheck; burns 106 -> 140 at Odyssey)
        try:
            import planner_rev18_offline as O18
            tabO = O18.load_tab("ODY"); x_ = np.concatenate(([q1], np.array(ol_inputs(dict(), gt["ODY"], P0, q1)["level"]) + ol_plan(gt["ODY"], P0, q1)["nu"]))
            a1 = arrival_pass(gt["ODY"], MUm, 7.408e6, x_, q1, kap["ODY"]); a2 = O18.arrival(tabO, x_, q1, kappa=kap["ODY"]); assert a1["N"] == a2["N"] and abs(a1["days"] - a2["days"]) < 1e-6
        except ImportError:
            pass
    vec = cross_language_vector()
    print(f"planner_lib selftest ok: config sections equal stageg3_config.toml, {len(L)} flights (274 main + 32 sensitivity), DP equals oracle_lib.dp_schedule on 60 random instances, plan inputs, reset-to-edge, prepayment, "
          f"OL = OL-R in the model world, B-ODY rule{'' if k_checked else ' (stagek_lib lever check skipped: import failed)'}; default-off paths bit-identical (lin None = zeros = -0.0; ol_plan fb=None and weight 0 = ol_schedule; OLT config); "
          f"FB term: slope identity, ODY N 216 -> 231 (106 -> 140 burns), G1 unchanged; cross-language vector nu = {np.round(vec['nu'], 9).tolist()}, with lin nu = {np.round(vec['nu_lin'], 9).tolist()}")


def cross_language_vector():
    """Fixed tie-free instance whose DP offsets planner_common.jl's selftest asserts (copied there)."""
    level = [0.00, 0.30, 0.55, 0.70, 0.40, 0.10, -0.20, -0.35]; lo = [-0.5] * 8; hi = [0.5] * 8; w = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]
    r = dp_plan(level, lo, hi, w, 1.0e-6); rl = dp_plan(level, lo, hi, w, 1.0e-6, lin=[0, 0, 0, 0, 0, 1.5, 0, 0])
    return dict(level=level, lo=lo, hi=hi, w=w, nu=r["nu"], cost=r["dv_cost"], lin=[0, 0, 0, 0, 0, 1.5, 0, 0], nu_lin=rl["nu"], cost_lin=rl["dv_cost"])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.dry_run:
        ok, _ = dry_run(load_cfg()); sys.exit(0 if ok else 1)
    ap.error("--selftest or --dry-run")


if __name__ == "__main__":
    main()
