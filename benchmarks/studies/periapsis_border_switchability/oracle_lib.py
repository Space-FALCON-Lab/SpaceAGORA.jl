"""Oracle-schedule library: lever, dynamic program over passes, and reference policies.

SUPERPOSITION ASSUMPTION (unverified; see the results note): a tangential burn at apoapsis shifts the next
periapsis radius by lever * dv and does not change the later passive drift, so the periapsis history under
a schedule is  x_k + u_k  with x_k the passive history and u_k the cumulative offset.

Units: lengths in m, speeds in m/s, lever L in m/s per m of periapsis shift (so dV = L * |du|).
No third-party package other than numpy is used.
"""
import math
import numpy as np

MU = 3.986004418e14   # Earth GM, m^3/s^2 (same value as analyze_static.py)
INF = float("inf")


def apoapsis_speed(a, e):
    ra = a * (1.0 + e)
    return np.sqrt(MU * (2.0 / ra - 1.0 / a))


def lever_kepler(a, e):
    """m/s per m. Tangential impulse dv at apoapsis: da = 2 a^2 v_a dv / mu, r_a fixed, so
    dr_p = 2 da = 4 a^2 v_a dv / mu  (km per m/s = 1e-3 * 4 a^2 v_a / mu); the inverse is returned."""
    a = np.asarray(a, float); e = np.asarray(e, float)
    return MU / (4.0 * a * a * apoapsis_speed(a, e))


def shift_per_dv_km(a, e):
    """Periapsis shift in km per m/s (the direct form of the same relation)."""
    return 4.0 * a * a * apoapsis_speed(a, e) / MU / 1000.0


class MeasuredLever:
    """Keplerian lever times a correction ratio R(a) = measured / Keplerian, interpolated in log a from the
    finite-difference points (a, e, measured m/s per km), held constant outside the measured a range."""

    def __init__(self, points):
        pts = sorted(points)             # (a_m, e, measured_mps_per_km)
        self.a = np.array([p[0] for p in pts])
        self.ratio = np.array([p[2] / (1000.0 * float(lever_kepler(p[0], p[1]))) for p in pts])
        self.points = pts

    def __call__(self, a, e):
        r = np.interp(np.log(np.asarray(a, float)), np.log(self.a), self.ratio)
        return lever_kepler(a, e) * r


# ------------------------------------------------------------------------------------------------ DP

def _run_argmin(c, idx):
    cm = np.minimum.accumulate(c)
    ia = np.maximum.accumulate(np.where(c <= cm, idx, -1))
    return cm, ia


def dp_schedule(A, B, L, res=None, lam=0.0, dvmin=0.0, tie=0.0, want_path=True):
    """Minimize  sum_k [ L_k |u_k - u_{k-1}| + (lam + tie) * 1(u_k != u_{k-1}) ]  over offsets u_k, u_0 = 0,
    subject to A_k <= u_k <= B_k (offset interval that keeps the periapsis inside the band at pass k), one burn
    opportunity per pass. Returns dict(dv, burns, u, du, obj) or None when infeasible; `dv` and `burns` are
    recomputed from the path (the tie penalty is not in dv).

    res=None: EXACT mode. The state set is the candidate offsets {0} U {A_j} U {B_j}. With cost piecewise linear
    in u, a path with a fixed burn pattern is an LP whose optimal vertex has every constant segment at one of
    the interval bounds A_j or B_j, so this set contains an optimum for any L_k and any burn penalty lam >= 0.
    res=float: uniform grid of spacing res (m) anchored at u_0 = 0 (restricts the path, so dV >= exact).
    `dvmin` (forbid burns with L_k |du| < dvmin, m/s) is supported on the grid only."""
    A = np.asarray(A, float); B = np.asarray(B, float); L = np.asarray(L, float)
    n = len(A)
    if res is None:
        if dvmin > 0: raise ValueError("dvmin needs a grid (res)")
        val = np.unique(np.concatenate(([0.0], A, B)))
    else:
        g0 = min(0, int(math.floor(A.min() / res)) - 1); g1 = max(0, int(math.ceil(B.max() / res)) + 1)
        val = (np.arange(g0, g1 + 1)) * res
    ns = len(val)
    idx = np.arange(ns)
    i0 = int(np.argmin(np.abs(val)))
    prev = np.full(ns, INF); prev[i0] = 0.0
    parent = np.empty((n, ns), dtype=np.int32) if want_path else None
    pen = lam + tie
    tol = 1e-7
    for k in range(n):
        Lk = L[k]
        w = max(1, int(math.ceil(dvmin / (Lk * res) - 1e-12))) if dvmin > 0 else 1
        conv = np.full(ns, INF); arg = np.zeros(ns, dtype=np.int64)
        if w < ns:
            c1 = prev - Lk * val
            cm, ia = _run_argmin(c1, idx)
            fwd = np.full(ns, INF); af = np.zeros(ns, dtype=np.int64)
            fwd[w:] = Lk * val[w:] + cm[:-w]; af[w:] = ia[:-w]
            c2 = (prev + Lk * val)[::-1]
            cm2r, ia2r = _run_argmin(c2, idx)
            cm2 = cm2r[::-1]; ab = (ns - 1 - ia2r)[::-1]
            bwd = np.full(ns, INF); ab2 = np.zeros(ns, dtype=np.int64)
            bwd[:-w] = -Lk * val[:-w] + cm2[w:]; ab2[:-w] = ab[w:]
            useb = bwd < fwd
            conv = np.where(useb, bwd, fwd); arg = np.where(useb, ab2, af)
        cand = conv + pen
        burn = cand < prev
        new = np.where(burn, cand, prev)
        if want_path:
            parent[k] = np.where(burn, arg, idx)
        new[(val < A[k] - tol) | (val > B[k] + tol)] = INF
        prev = new
    j = int(np.argmin(prev))
    if not np.isfinite(prev[j]):
        return None
    out = dict(obj=float(prev[j]))
    if want_path:
        gs = np.empty(n, dtype=np.int64)
        for k in range(n - 1, -1, -1):
            gs[k] = j; j = int(parent[k, j])
        u = val[gs]
        du = np.diff(np.concatenate(([0.0], u)))
        out.update(u=u, dv=float(np.sum(L * np.abs(du))), burns=int((np.abs(du) > 1e-9).sum()), du=du)
    return out


def dp_budget(A, B, L, cmax):
    """Exact trade-off curve: best[c] = minimum dV using at most c burns (c = 0..cmax; inf if infeasible), over the
    candidate-offset state set of the exact DP. State = (burns used, offset); O(cmax * n_states) per pass."""
    A = np.asarray(A, float); B = np.asarray(B, float); L = np.asarray(L, float)
    val = np.unique(np.concatenate(([0.0], A, B))); ns = len(val)
    prev = np.full((cmax + 1, ns), INF); prev[0, int(np.argmin(np.abs(val)))] = 0.0
    tol = 1e-7
    for k in range(len(A)):
        Lk = L[k]
        cm = np.minimum.accumulate(prev - Lk * val, axis=1)
        cm2 = np.minimum.accumulate((prev + Lk * val)[:, ::-1], axis=1)[:, ::-1]
        conv = np.full((cmax + 1, ns), INF)
        conv[:, 1:] = Lk * val[1:] + cm[:, :-1]
        conv[:, :-1] = np.minimum(conv[:, :-1], -Lk * val[:-1] + cm2[:, 1:])
        new = prev.copy()
        new[1:] = np.minimum(prev[1:], conv[:-1])
        new[:, (val < A[k] - tol) | (val > B[k] + tol)] = INF
        prev = new
    return np.minimum.accumulate(prev.min(axis=1))


def path_cost(u, L):
    du = np.diff(np.concatenate(([0.0], np.asarray(u, float))))
    return float(np.sum(np.asarray(L) * np.abs(du))), int((np.abs(du) > 1e-9).sum()), du


# ------------------------------------------------------------------------------------------ policies

def _enforce_dvmin(pos, tgt, lo, hi, L, dvmin):
    """If the burn from pos to tgt is cheaper than dvmin, lengthen it in the same direction (clamped to the band)."""
    d = tgt - pos
    if dvmin > 0 and L * abs(d) < dvmin:
        d = math.copysign(dvmin / L, d)
        tgt = min(max(pos + d, lo), hi)
    return tgt


def lazy(x, lo, hi, L, s=0.0, dvmin=0.0):
    """On a forecast exit (the periapsis at this pass would leave [lo, hi]) burn to the violated edge moved
    `s` meters inside, capped at the band center. s = 0: step just enough. s >= half width: near_passive."""
    u = 0.0; cost = 0.0; n = 0; steps = np.zeros(len(x))
    for k in range(len(x)):
        pos = x[k] + u
        if pos > hi[k]:
            tgt = max(hi[k] - s, 0.5 * (lo[k] + hi[k]))
        elif pos < lo[k]:
            tgt = min(lo[k] + s, 0.5 * (lo[k] + hi[k]))
        else:
            continue
        tgt = _enforce_dvmin(pos, tgt, lo[k], hi[k], L[k], dvmin)
        du = tgt - pos; u += du; steps[k] = du; cost += L[k] * abs(du); n += 1
    return cost, n, steps


def near_passive(x, lo, hi, L, dvmin=0.0):
    return lazy(x, lo, hi, L, s=1e18, dvmin=dvmin)


def reset_opposite(x, lo, hi, L, dvmin=0.0):
    u = 0.0; cost = 0.0; n = 0; steps = np.zeros(len(x))
    for k in range(len(x)):
        pos = x[k] + u
        if pos > hi[k]: tgt = lo[k]
        elif pos < lo[k]: tgt = hi[k]
        else: continue
        tgt = _enforce_dvmin(pos, tgt, lo[k], hi[k], L[k], dvmin)
        du = tgt - pos; u += du; steps[k] = du; cost += L[k] * abs(du); n += 1
    return cost, n, steps


def hold(x, lo, hi, L, mode, max_span=None):
    """Greedy farthest-hold with the perfect forecast: on a forced burn, find the longest run of passes that one
    offset can satisfy (intersection of the offset intervals) and pick a point in it:
    'nearest' (least dV for this burn), 'center' (most margin), 'far' (the end of the interval farthest from
    the current offset, i.e. aim at the upstream edge). With unlimited forecast all three have the same burn passes
    (the minimum burn count). `max_span` limits the forecast to that many passes including the current one
    (max_span=1 is lazy-to-edge)."""
    A = lo - x; B = hi - x; nn = len(x)
    u = 0.0; cost = 0.0; n = 0; k = 0; steps = np.zeros(nn)
    while k < nn:
        if A[k] <= u <= B[k]:
            k += 1; continue
        mA, mB, m = A[k], B[k], k
        while m + 1 < nn and (max_span is None or m + 1 - k < max_span) and max(mA, A[m + 1]) <= min(mB, B[m + 1]):
            m += 1; mA = max(mA, A[m]); mB = min(mB, B[m])
        if mode == "nearest": tgt = min(max(u, mA), mB)
        elif mode == "center": tgt = 0.5 * (mA + mB)
        elif mode == "far": tgt = mA if abs(u - mA) > abs(u - mB) else mB
        else: raise ValueError(mode)
        du = tgt - u; u = tgt; steps[k] = du; cost += L[k] * abs(du); n += 1; k += 1
    return cost, n, steps
