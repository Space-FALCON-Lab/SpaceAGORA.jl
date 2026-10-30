"""Revision 18 offline computations (analysis only, no Julia, no flights): (1) the apoapsis-burn timing lever dT = 3 T sqrt((1-e)/(1+e)) dv / V_c, (2) the in-model arrival time of a planned ln-q path and the dV-versus-arrival-time curve obtained
by sweeping a price lambda (m/s per orbit of arrival time) inside the existing exact DP (a linear term on the level), per geometry from planner_ground_tables.json, with D_min / D_max at the tube edges.
Model (DERIVED, see results_rev18_design.md section 3): drop per pass = d0_k q_k / q_req (the table is the q_req schedule; solve_q_req uses time ~ 1/q), period from the osculating a of each orbit, burn-induced periapsis shift -H nu.
Usage: python3 -I planner_rev18_offline.py --out F.json"""
import argparse, json, math, pathlib, sys
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_analysis as A
MU = A.MU; RA_T = A.RA_TARGET; TD = pathlib.Path(PLB.task_dir()) if hasattr(PLB, "task_dir") else None
CFG = PLB.load_cfg(); P = PLB.flat_P(CFG)
Q_FIRST = 0.212      # MEASURED: calibration target / q_first of the cal flights (planner_config cal q target 0.212)


def load_tab(token):
    T = json.load(open(pathlib.Path(PLB.task_dir()) / "planner_ground_tables.json"))["geometries"]
    g = T[{"ODY": "ODY", "G1": "G5863"}[token]]; tab = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in g.items()}; tab["token"] = token; return tab


def lever_T_s_per_mps(a, e):
    """DERIVED: dT/dv (s per m/s) of a tangential apoapsis burn: 3 T sqrt((1-e)/(1+e)) / V_c, V_c = sqrt(mu/a), T = 2 pi sqrt(a^3/mu)."""
    T = 2 * math.pi * math.sqrt(a ** 3 / MU); return 3 * T * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)


def arrival_index(tab, x, q_first_ln=math.log(Q_FIRST)):
    """In-model closing-apoapsis time (days) and pass count of a planned ln-q path x_k (k = 1.. ; held at the last value beyond its length). Returns dict(days, N, ra_end)."""
    n = tab["n"]; q_req = tab["q_req"]; ra = float(tab["ra"][0]); t = float(tab["t"][0]); path = tab["path"]; kmax = int(n * 3)
    d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])]))      # drop of pass k at q_req
    for k in range(1, kmax + 1):
        kk = min(k, n); xk = float(x[min(k, len(x)) - 1]); nu = xk - (q_first_ln + float(path[kk - 1] - path[0]))
        rp = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nu; a = 0.5 * (ra + rp)
        drop = float(d0[kk - 1]) * math.exp(xk) / q_req; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2 * math.pi * math.sqrt(amid ** 3 / MU)
        if ra_next <= RA_T:
            return dict(days=(t + 0.5 * Tk) / 86400.0, N=k, ra_end=ra_next)
        ra = ra_next; t += Tk
    return dict(days=float("nan"), N=kmax, ra_end=ra)


KAPPA = {"ODY": 1.0, "G1": 1.0}      # drop-per-q calibration factor per geometry; set from the flights by planner_rev18_calibrate.py (MEASURED), default 1 = the uncalibrated table


def arrival(tab, x, q_first_ln=math.log(Q_FIRST), kappa=None):
    """State-dependent version (Revision 18b): the drop per unit q is read at the CURRENT apoapsis radius (d0 interpolated on the table's ra), not at the pass index, so a path that runs slower or faster than the table schedule is
    treated correctly; kappa scales the drop (calibration, MEASURED per geometry). Returns dict(days, N, ra_end)."""
    kap = KAPPA.get(tab.get("token"), 1.0) if kappa is None else kappa
    n = tab["n"]; q_req = tab["q_req"]; ra = float(tab["ra"][0]); t = float(tab["t"][0]); path = tab["path"]; kmax = int(n * 3)
    d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); ra_t = tab["ra"][::-1]; d0_r = d0[::-1]
    for k in range(1, kmax + 1):
        kk = min(k, n); xk = float(x[min(k, len(x)) - 1]); nu = xk - (q_first_ln + float(path[kk - 1] - path[0]))
        rp = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nu; a = 0.5 * (ra + rp)
        drop = kap * float(np.interp(ra, ra_t, d0_r)) * math.exp(xk) / q_req; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2 * math.pi * math.sqrt(amid ** 3 / MU)
        if ra_next <= RA_T:
            return dict(days=(t + 0.5 * Tk) / 86400.0, N=k, ra_end=ra_next)
        ra = ra_next; t += Tk
    return dict(days=float("nan"), N=kmax, ra_end=ra)



def dp_lin(A_, B_, L, lin, tol=1.0e-9):
    """Exact DP of planner_lib.dp_exact with an added linear cost lin_k * nu_k (vertex set {0} U {A} U {B} stays optimal for a piecewise-linear objective). Returns nu or None."""
    A_ = np.asarray(A_, float); B_ = np.asarray(B_, float); L = np.asarray(L, float); lin = np.asarray(lin, float); n = len(A_)
    val = np.array(sorted(set([0.0] + list(A_) + list(B_)))); ns = len(val); i0 = int(np.argmin(np.abs(val))); INF = float("inf")
    prev = np.full(ns, INF); prev[i0] = 0.0; parent = np.zeros((n, ns), dtype=np.int32)
    for k in range(n):
        cur = prev.copy(); arg = np.arange(ns); Lk = L[k]
        for j in range(1, ns):
            c = cur[j - 1] + Lk * (val[j] - val[j - 1])
            if c < cur[j]: cur[j] = c; arg[j] = arg[j - 1]
        for j in range(ns - 2, -1, -1):
            c = cur[j + 1] + Lk * (val[j + 1] - val[j])
            if c < cur[j]: cur[j] = c; arg[j] = arg[j + 1]
        cur = cur + lin[k] * val
        cur[(val < A_[k] - tol) | (val > B_[k] + tol)] = INF; parent[k] = arg; prev = cur
    j = int(np.argmin(prev))
    if not np.isfinite(prev[j]): return None
    nu = np.zeros(n)
    for k in range(n - 1, -1, -1): nu[k] = val[j]; j = parent[k][j]
    return nu


def plan_for_price(tab, lam_orb, eps=1.0e-6, pinned=None):
    """OL-style schedule on the ground table with a price lam_orb (m/s per orbit of arrival time, + rewards a shorter arrival) on the level. pinned = 'hi' or 'lo' pins the level to a tube edge. Returns dict(nu, dv, x, burns)."""
    inp = PLB.ol_inputs(dict(), tab, P, math.log(Q_FIRST)); lvl = np.array(inp["level"]); lo = np.array(inp["lo"]); hi = np.array(inp["hi"]); w = np.array(inp["w"])
    if pinned == "hi": lo = hi.copy()
    if pinned == "lo": hi = lo.copy()
    # reference path drop weights: orbits of arrival per ln-q unit at pass k = drop_k / drop_N on the nominal (lambda = 0) path
    ref = plan_for_price(tab, 0.0) if (lam_orb != 0.0 and pinned is None) else None
    if lam_orb != 0.0 and pinned is None:
        xr = ref["x_full"]; d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); n = tab["n"]
        dk = np.array([d0[min(k, n) - 1] * math.exp(xr[min(k, len(xr)) - 1]) / tab["q_req"] for k in range(1, len(xr) + 1)]); wT = dk[1:1 + len(lvl)] / dk[min(len(dk) - 1, ref["N"] - 1)]
        lin = -lam_orb * wT
    else:
        lin = np.zeros(len(lvl))
    A_ = lo - lvl; B_ = hi - lvl; tw = PLB.tie_weights(w, eps)
    nu = dp_lin(A_, B_, tw, lin)
    if nu is None: return None
    du = np.diff(np.concatenate(([0.0], nu))); dv = float(np.sum(w * np.abs(du))); x = lvl + nu; xfull = np.concatenate(([math.log(Q_FIRST)], x))
    r = arrival(tab, xfull)
    return dict(nu=nu, dv=dv, x_full=xfull, days=r["days"], N=r["N"], burns=int((np.abs(du) * w >= P["dp_min_burn_mps"]).sum()))


def curve(token, lams):
    tab = load_tab(token); rows = []
    for lam in lams:
        r = plan_for_price(tab, lam)
        if r: rows.append(dict(lam_orb=lam, dv=r["dv"], days=r["days"], N=r["N"], burns=r["burns"]))
    return tab, rows


def selftest():
    # derivation: dT = (3/2) T da / a, da = drp/2, drp = 4 a sqrt((1-e)/(1+e)) dv / Vc  ->  dT = 3 T sqrt((1-e)/(1+e)) dv / Vc
    for a, e in ((1.6e7, 0.78), (1.0e7, 0.65), (5.4e6, 0.36)):
        T = 2 * math.pi * math.sqrt(a ** 3 / MU); Vc = math.sqrt(MU / a); dv = 1.0; drp = PLB.lever_m_per_mps(MU, a, e) * dv; dT = 1.5 * T * (drp / 2) / a
        assert abs(dT - lever_T_s_per_mps(a, e) * dv) < 1e-9 * dT, (a, e)
        T1 = 2 * math.pi * math.sqrt(((a + 1e-3 * 0 + drp / 2) ** 3) / MU); assert abs((T1 - T) - dT) < 0.01 * dT, ((T1 - T), dT)    # exact period change of da = drp/2 agrees to 1% (first order in da/a)
    tab = load_tab("ODY"); r = arrival(tab, np.full(300, math.log(tab["q_req"])), math.log(tab["q_req"]) - float(tab["path"][0]) * 0 )
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out"); a = ap.parse_args(); print("selftest arrival at constant q_req (ODY):", selftest())
