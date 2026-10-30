"""Checks for part2_derivations.md (tasks/2026-10-15_periapsis-region-border-targeting): burn-timing optimality and the lever, the slack-band DP structure (clamp optimality and the
prepayment bound), the harmonic time law, the reflected-walk (limit-cycle) dV, the pass-maximum statistic, the time-risk frontier, and the post hoc ladder replays.
Python only, single process. No Julia, no flights. Every input file is read-only; the only files written are analytic_data/part2_results.json and analytic_data/cache/part2_dev_tables.json
in the task directory (the dev-config ground table, rebuilt in memory from the frozen ranking code with the pre-calibration values q_cap 0.234 and d_ref 64.2 d, about 35 s once).

Usage: python3 -I part2_derivations_check.py --selftest            synthetic checks only (under 10 s)
       python3 -I part2_derivations_check.py --data [--results-root R]   flight-table checks (Stage O dev passes, planner dev replan/passes), writes the JSON
"""
import argparse, glob, json, math, pathlib, sys, time
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB

MU = 4.2828370371e13          # SOURCED: Mars mu (half-lever note, planner tables)
RHO_S = 0.5826                # -zeta(1/2)/sqrt(2 pi): expected overshoot / sigma of a Gaussian random walk (Siegmund); verified numerically below
import os
RESULTS_ROOT = os.environ.get("PLANNER_DATA_ROOT", ".")


# ------------------------------------------------------------------------------------------------ item 1: lever and burn timing

def lever(a, e):
    return PLB.lever_m_per_mps(MU, a, e)


def rp_after(a, e, f, dv, alpha, beta=0.0):
    """Exact periapsis radius after an impulse dv at true anomaly f; alpha = in-plane angle from the velocity (0 = along v), beta = out-of-plane angle."""
    p = a * (1 - e * e); r = p / (1 + e * math.cos(f)); k = math.sqrt(MU / p)
    vr, vt = k * e * math.sin(f), k * (1 + e * math.cos(f)); V = math.hypot(vr, vt)
    tr, tt = vr / V, vt / V; nr, nt = -tt, tr                                 # in-plane unit normal (rotated +90 deg)
    ip = dv * math.cos(beta); vr2 = vr + ip * (math.cos(alpha) * tr + math.sin(alpha) * nr); vt2 = vt + ip * (math.cos(alpha) * tt + math.sin(alpha) * nt)
    w = dv * math.sin(beta); V2 = vr2 ** 2 + vt2 ** 2 + w * w
    a2 = 1.0 / (2.0 / r - V2 / MU); h2 = r * math.hypot(vt2, w); e2 = math.sqrt(max(0.0, 1.0 - h2 * h2 / (MU * a2)))
    return a2 * (1 - e2)


def lever_first_order(a, e, f):
    """(A, B): d r_p / d v_T and d r_p / d v_N (in-plane normal) to first order (Gauss); |optimal lever| = hypot(A, B)."""
    p = a * (1 - e * e); r = p / (1 + e * math.cos(f)); V = math.sqrt(MU * (2 / r - 1 / a))
    return 2 * a * (1 - e) * (1 - math.cos(f)) / ((1 + e) * V), r * math.sin(f) / V


def P_poly(c, e):
    return 16 * (1 + 2 * e * c + e * e) * (1 + e * c) ** 2 - (1 - e) ** 2 * (4 * (1 - c) ** 2 * (1 + e * c) ** 2 + (1 - c * c) * (1 + e) ** 4)


def K_poly(c, e):
    return -(2 * e * (1 - e) ** 2 * c * c - (e ** 4 + 4 * e ** 3 + 2 * e * e + 12 * e - 3) * c + (e ** 4 + 2 * e ** 3 - 2 * e * e - 6 * e - 11))


def dlnlever_dk(e, kQ):
    """Per-pass relative change of the lever at fixed r_p under the Laplace pass law d(e^(3/2))/dk = -3 kQ."""
    return -(1 + 3 * e) * kQ / (math.sqrt(e) * (1 - e * e))


# ------------------------------------------------------------------------------------------------ item 2: clamp (play operator) and DP

def clamp_cost(A, B, w, u0=0.0):
    u = u0; c = 0.0; nb = 0
    for Ak, Bk, wk in zip(A, B, w):
        v = min(max(u, Ak), Bk)
        if v != u:
            c += wk * abs(v - u); nb += 1
        u = v
    return c, nb


# ------------------------------------------------------------------------------------------------ items 4 and 5: reflected walk

def tv_rate(d, s, W, rho=RHO_S):
    """Expected |burn| per pass of the two-sided clamp on a walk with drift d, step sd s, band W (BM formula d coth(d We / s^2), We = W + 2 rho s)."""
    We = W + 2 * rho * s
    return s * s / We if abs(d) < 1e-12 else abs(d) / math.tanh(abs(d) * We / (s * s))


def reaim_rate(s, W, g, rho=RHO_S):
    """Re-aim keeper (burn to g inside the edge it crossed), driftless: expected |burn| per pass and burns per pass."""
    return s * s / (W - g + rho * s), s * s / ((g + rho * s) * (W - g + rho * s))


def sim_clamp(W, s, d, n, rng):
    e = rng.normal(d, s, n); y = -W / 2; c = 0.0
    for v in e:
        y += v
        if y > 0:
            c += y; y = 0.0
        elif y < -W:
            c += -W - y; y = -W
    return c / n


def sim_reaim(W, s, g, n, rng):
    e = rng.normal(0, s, n); y = -g; c = 0.0; nb = 0
    for v in e:
        y += v
        if y > 0:
            c += y + g; y = -g; nb += 1
        elif y < -W:
            c += (-W + g) - y; y = -W + g; nb += 1
    return c / n, nb / n


# ------------------------------------------------------------------------------------------------ item 3: pass maximum of an OU log-perturbation under a Gaussian drag window

def ou_paths(n, tau, dt, T, rng):
    t = np.arange(-T, T + dt / 2, dt); K = len(t); a = math.exp(-dt / tau)
    Z = np.empty((n, K)); Z[:, 0] = rng.standard_normal(n)
    for k in range(1, K):
        Z[:, k] = a * Z[:, k - 1] + math.sqrt(1 - a * a) * rng.standard_normal(n)
    return t, Z


def pass_max(sig_t, X, t, G):
    """M = max_t [X(t) - (t - t*)^2 / (2 sig_t^2) + G], t* = sig_t sqrt(2G) so that the profile is 0 at the event t = 0 and peaks at G."""
    ts = sig_t * math.sqrt(2 * G) if np.isscalar(sig_t) else sig_t * np.sqrt(2 * G)
    st = np.asarray(sig_t, float).reshape(-1, 1) if not np.isscalar(sig_t) else sig_t
    tsc = np.asarray(ts, float).reshape(-1, 1) if not np.isscalar(ts) else ts
    prof = -(t[None, :] - tsc) ** 2 / (2 * st ** 2) + G
    return (X + prof).max(1)


# ------------------------------------------------------------------------------------------------ selftest

def selftest():
    t0 = time.time(); rng = np.random.default_rng(20261027); out = {}
    # S1a the lever and its closed forms
    a0, e0 = 16017e3, 0.7814; rp = a0 * (1 - e0); L0 = lever(a0, e0)
    ra1 = 7408e3; a1 = (ra1 + rp) / 2; e1 = (ra1 - rp) / (ra1 + rp); L1 = lever(a1, e1)
    Vcp = math.sqrt(MU / rp)
    for a, e in ((a0, e0), (a1, e1)):
        ra = a * (1 + e)
        assert abs(lever(a, e) - 4 * rp / (Vcp * (1 - e) * math.sqrt(1 + e))) < 1e-9 * lever(a, e)
        assert abs(lever(a, e) - math.sqrt(2) * (ra + rp) ** 1.5 / (Vcp * math.sqrt(ra))) < 1e-9 * lever(a, e)
    assert 13.6e3 < L0 < 13.8e3 and 5.3e3 < L1 < 5.4e3
    out["lever_start_m_per_mps"], out["lever_end_m_per_mps"], out["lever_ratio"] = L0, L1, L1 / L0
    # S1b exact finite impulses: the best |dr_p| per dv over f, alpha, beta is at apoapsis along +-v; first order equals the exact derivative
    best = {}
    for e in (0.05, 0.36, 0.6, 0.7814, 0.9):
        a = rp / (1 - e); dv = 0.01; m = (0.0, None)
        for f in np.linspace(0, math.pi, 91):
            for al in np.linspace(-math.pi, math.pi, 73):
                g = abs(rp_after(a, e, f, dv, al) - a * (1 - e)) / dv
                if g > m[0]:
                    m = (g, (round(math.degrees(f), 1), round(math.degrees(al), 1)))
        Lx = lever(a, e); assert m[1][0] == 180.0 and m[1][1] in (0.0, 180.0, -180.0) and abs(m[0] / Lx - 1) < 2e-4, (e, m, Lx)
        for f in (0.7, 1.6, 2.5):
            A, B = lever_first_order(a, e, f); h = 1e-4
            dA = (rp_after(a, e, f, h, 0.0) - rp_after(a, e, f, -h, 0.0)) / (2 * h); dB = (rp_after(a, e, f, h, math.pi / 2) - rp_after(a, e, f, -h, math.pi / 2)) / (2 * h)
            assert abs(dA - A) < 1e-5 * Lx and abs(abs(dB) - abs(B)) < 1e-5 * Lx
        assert abs(rp_after(a, e, math.pi, dv, 0.0, math.radians(30)) - a * (1 - e) - math.cos(math.radians(30)) * Lx * dv) < 2e-3 * Lx * dv   # out-of-plane part is wasted (first order)
        best[str(e)] = m
    out["burn_timing_grid_best"] = best
    # S1c the inequality: P(c,e) = (1+c)(1+2ec+e^2) K(c,e) and K >= 0 on c in [-1,1], e in [0,1)
    for _ in range(200):
        c, e = rng.uniform(-1, 1), rng.uniform(0, 0.99)
        assert abs(P_poly(c, e) - (1 + c) * (1 + 2 * e * c + e * e) * K_poly(c, e)) < 1e-9 * (1 + abs(P_poly(c, e)))
    C, E = np.meshgrid(np.linspace(-1, 1, 401), np.linspace(0, 0.999, 400)); assert K_poly(C, E).min() >= 0
    # S1d per-pass lever shrink: closed form vs finite difference of the pass law (kQ 1.01e-3 is the review's approximate Mars value at Q 0.224)
    kQ = 1.01e-3
    for e in (0.78, 0.6, 0.36):
        e_next = (e ** 1.5 - 3 * kQ) ** (2 / 3); fd = math.log(lever(rp / (1 - e_next), e_next) / lever(rp / (1 - e), e))
        assert abs(fd / dlnlever_dk(e, kQ) - 1) < 0.02, (e, fd, dlnlever_dk(e, kQ))
    out["dlnlever_per_pass"] = {str(e): dlnlever_dk(e, kQ) for e in (0.78, 0.36)}
    # S2 clamp optimality and the prepayment bound (random instances, exact DP of planner_lib)
    worst = 0.0
    for trial in range(60):
        n = int(rng.integers(8, 40)); lev = np.cumsum(rng.normal(rng.normal(0, 0.05), 0.15, n)); W = rng.uniform(0.2, 0.8)
        A = -lev - W / 2; B = -lev + W / 2
        wc = np.full(n, rng.uniform(0.3, 1.2)); assert abs(clamp_cost(A, B, wc)[0] - PLB.dp_exact(A, B, wc)["obj"]) < 1e-9
        wd = np.sort(rng.uniform(0.3, 1.2, n))[::-1]; assert abs(clamp_cost(A, B, wd)[0] - PLB.dp_exact(A, B, wd)["obj"]) < 1e-9
        wi = np.sort(rng.uniform(0.3, 1.2, n)); gap = clamp_cost(A, B, wi)[0] - PLB.dp_exact(A, B, wi)["obj"]
        assert -1e-9 <= gap <= W * (wi[-1] - wi[0]) + 1e-9; worst = max(worst, gap / (W * (wi[-1] - wi[0])))
    n = 400; d = 0.01; W = 0.5; A = d * np.arange(1, n + 1) - 0.25; B = A + W; wi = np.linspace(0.5, 1.2, n)
    tight = (clamp_cost(A, B, wi)[0] - PLB.dp_exact(A, B, wi)["obj"]) / (W * (wi[-1] - wi[0]))
    assert tight > 0.9, tight
    out["prepay_bound_random_worst_fraction"], out["prepay_bound_ramp_fraction"] = worst, tight
    # S3 harmonic time law: integrate the energy decay with a wandering q, compare with (1/c) int T/(tau q) dx and the period-weighted harmonic mean on the reference grid
    x0, x1, cdec = 1 / 16e6, 1 / 5.5e6, 2.0e-12; xs = np.linspace(x0, x1, 20001)
    Tf = lambda x: 2 * math.pi * math.sqrt((1 / x) ** 3 / MU); tauf = lambda x: 60.0 * (1 + 0.5 * (x - x0) / (x1 - x0))
    qf = lambda x: 0.2 * math.exp(0.35 * math.sin(40 * (x - x0) / (x1 - x0)))
    x = x0; t = 0.0
    while x < x1:
        x += cdec * qf(x) * tauf(x); t += Tf(x)
    I = np.trapz([Tf(v) / (tauf(v) * qf(v)) for v in xs], xs) / cdec
    assert abs(t / I - 1) < 0.01, (t, I)
    for Wb in (0.72, 1.034):
        zz = np.linspace(-Wb / 2, Wb / 2, 200001); assert abs(math.log(np.mean(np.exp(-zz))) - math.log(math.sinh(Wb / 2) / (Wb / 2))) < 1e-6
    out["jensen_penalty_uniform"] = {str(Wb): math.log(math.sinh(Wb / 2) / (Wb / 2)) for Wb in (0.7186, 1.0335)}
    # S4 reflected walk: clamp cost vs d coth(d We/s^2); one-sided Lindley E[burn] = d; re-aim keeper
    rw = []
    for Wb, s, dd in ((0.72, 0.066, 0.0), (0.46, 0.066, 0.0), (0.72, 0.066, 0.003)):
        sim = sim_clamp(Wb, s, dd, 120000, rng); th = tv_rate(dd, s, Wb); rw.append((Wb, s, dd, sim, th)); assert abs(sim / th - 1) < 0.06, (Wb, s, dd, sim, th)
    e = rng.normal(0.004, 0.069, 120000); y = 0.0; cb = 0.0
    for v in e:
        y += v
        if y > 0:
            cb += y; y = 0.0
    assert abs(cb / len(e) / 0.004 - 1) < 0.05
    for Wb, s, g in ((0.22, 0.069, 0.10), (0.72, 0.066, 0.36)):
        sc, sf = sim_reaim(Wb, s, g, 120000, rng); tc, tf = reaim_rate(s, Wb, g); assert abs(sc / tc - 1) < 0.08 and abs(sf / tf - 1) < 0.1, (sc, tc, sf, tf)
    out["reflected_walk_checks"] = rw
    # S5 pass maximum: exact lower bound, the two limits, and the consistency of the SOURCED moments
    t, Z = ou_paths(3000, 1e9, 1.0, 120.0, rng); M = pass_max(55.0, 0.277 * Z, t, 0.0); X0 = 0.277 * Z[:, len(t) // 2]
    assert np.allclose(M, X0)                                               # tau -> inf, G = 0: the maximum is the event value
    t, Z = ou_paths(3000, 22.2, 1.0, 200.0, rng); X = 0.277 * Z; M = pass_max(55.0, X, t, 0.0555); X0 = X[:, len(t) // 2]
    assert np.all(M >= X0 - 1e-12) and np.corrcoef(M - X0, X0)[0, 1] < -0.5
    sdM = math.sqrt(0.210 ** 2 + 0.277 ** 2 + 2 * (-0.73) * 0.210 * 0.277); assert abs(sdM - 0.188) < 0.004
    out["sourced_moment_consistency_sd_M"] = sdM
    # S6 frontier slope with an exponential tail: d ln q / d ln E = beta / (1 - beta)
    beta, q99, Nref, qref = 0.045, 0.707, 200.0, 0.163
    def qE(Ev):
        lo, hi = 0.05, 0.54
        for _ in range(80):
            q = math.sqrt(lo * hi); x = math.log(0.54 / q); N = Nref * qref / q
            if 0.01 * math.exp(-(x - q99) / beta) * N > Ev: hi = q
            else: lo = q
        return q
    sl = (math.log(qE(1.0)) - math.log(qE(0.5))) / math.log(2.0); assert abs(sl - beta / (1 - beta)) < 1e-4
    # S7 waiting versus prepaying an uncertain need: cost formulas vs Monte Carlo
    from scipy.stats import norm
    mu_D, sd_D, wk, wj = 0.3, 0.2, 0.5, 0.6; D = rng.normal(mu_D, sd_D, 400000)
    lazy = wj * np.mean(np.abs(D)); early = wk * mu_D + wj * np.mean(np.abs(D - mu_D))
    lazy_th = wj * (sd_D * math.sqrt(2 / math.pi) * math.exp(-mu_D ** 2 / (2 * sd_D ** 2)) + mu_D * (1 - 2 * norm.cdf(-mu_D / sd_D)))
    early_th = wk * mu_D + wj * sd_D * math.sqrt(2 / math.pi)
    assert abs(lazy / lazy_th - 1) < 0.01 and abs(early / early_th - 1) < 0.01
    # S8 least-squares slope variance over N passes: 12 s^2 / (N (N^2 - 1))
    N = 15; xx = np.arange(N); sl = [np.polyfit(xx, rng.normal(0, 1, N), 1)[0] for _ in range(20000)]
    assert abs(np.var(sl) / (12 / (N * (N * N - 1))) - 1) < 0.05
    out["selftest_s"] = time.time() - t0
    print("part2 selftest ok:", json.dumps({k: v for k, v in out.items() if k != "burn_timing_grid_best"}, default=float))
    print("  burn-timing grid optimum (deg f, alpha) per e:", best)
    assert out["selftest_s"] < 10.0
    return out


# ------------------------------------------------------------------------------------------------ data mode

def _task():
    t = PLB.task_dir(); assert t is not None, "task directory not found"; return t


def dev_tables(task):
    """The pre-calibration ground tables the dev flights flew (q_cap 0.234, d_ref 64.2 d), rebuilt in memory with the frozen code; cached."""
    p = task / "analytic_data" / "cache" / "part2_dev_tables.json"
    if p.exists():
        return json.loads(p.read_text())
    import copy
    import planner_predictions as PP, stageg3_analysis as G3A
    cfg = PLB.load_cfg(); cfg3 = G3A.load3(); cfgd = copy.deepcopy(cfg); cfgd["planner"]["q_cap_mean_w_cm2"] = 0.234; cfgd["planner"]["d_ref_days"] = 64.2
    ctx = PP.make_ctx_cal(cfg3, cfgd); T = PP.build_tables(cfgd, cfg3, ctx, PP.geometries(cfgd, cfg3))
    p.parent.mkdir(parents=True, exist_ok=True); p.write_text(json.dumps(T)); return T


def stageo_passes(task):
    import pandas as pd
    fr = []
    for f in sorted(glob.glob(str(task / "stageo_data" / "dev" / "runs" / "*" / "passes.csv"))):
        P = pd.read_csv(f, usecols=["pass", "vrel_mps", "rho", "e", "h_local_m", "rp_m", "a_m", "r_pert", "rho_peak_kg_m3"]); P["run"] = pathlib.Path(f).parent.name; fr.append(P)
    P = pd.concat(fr); return P[np.isfinite(P.rho_peak_kg_m3) & (P.rho_peak_kg_m3 > 0)]


def flight(root, geom, arm, s):
    import pandas as pd
    d = pathlib.Path(root) / "planner" / "dev" / f"{geom}_PL_{arm}_pert_ODY_s{s}"
    Rp = pd.read_csv(d / "replan.csv"); P = pd.read_csv(d / "passes.csv")
    lnq = np.log(P.q_w_cm2.values); n = len(lnq); m = dict(zip(Rp["pass"], Rp["nu_cum"])); nu = np.zeros(n); cur = 0.0
    for p in range(1, n + 1):
        cur = m.get(p, cur); nu[p - 1] = cur
    x = lnq - nu; flown_all = float(np.abs(Rp.executed_dv_mps).sum()); Rp = Rp[Rp["pass"] <= n]
    e = {int(p): float(lnq[p - 1] - (z - dv / w)) for p, z, dv, w in zip(Rp["pass"], Rp.z0_pred_lnq, Rp.executed_dv_mps, Rp.w0)}
    w = {int(p): float(v) for p, v in zip(Rp["pass"], Rp.w0)}
    T = 2 * np.pi * np.sqrt(P.a_m.values ** 3 / MU)
    return dict(x=x, e=e, w=w, n=n, lnq=lnq, T=T, days=float(P.time_s.values[-1] / 86400), flown=float(np.abs(Rp.executed_dv_mps).sum()), flown_all=flown_all,
                a=Rp.a_m.values, ecc=Rp.e.values, w0=Rp.w0.values, z0=Rp.z0_pred_lnq.values)


def Pcfg(cfg, qcap, pen=None, res=None):
    P = PLB.flat_P(cfg); P["q_cap_mean_w_cm2"] = qcap
    if pen is not None: P["dp_burn_penalty_mps"] = pen
    if res is not None: P["fast_reserve_ln"] = res
    return P


def replay(F, tab, P, sw, olr=False, fast_zero=False, oracle=None):
    """Re-run the planner's decision rule (planner_lib.decide, the Python mirror of the flight code) on the flight's realized natural path x (burn-removed ln q_mean) and its recorded
    one-pass prediction errors. With the flown config it reproduces the flown burns exactly (checked). Counterfactual arms assume the natural path does not depend on the policy (superposition,
    the planner's own assumption). oracle='res' / 'full': the feed-forward level is the realized future natural path (perfect forecast), tube shrunk by half the reserve or not."""
    x = F["x"]; u = 0.0; dvs = []; z = np.empty(F["n"]); z[0] = x[0]
    t = dict(tab, fast=[0.0] * len(tab["fast"])) if fast_zero else tab
    for p in range(2, F["n"] + 1):
        if p in F["w"]:
            z0 = x[p - 1] + u - F["e"][p]; w0 = F["w"][p]
            if oracle:
                hi0 = math.log(P["q_cap_mean_w_cm2"]); lo0 = hi0 - tab["W"]; r = 0.0 if oracle == "full" else P["fast_reserve_ln"] / 2; m_end = min(F["n"], p + int(sw["horizon"]) - 1)
                lev = [z0] + [x[q - 1] + u for q in range(p + 1, m_end + 1)]; k = len(lev)
                pl = PLB.dp_plan(lev, [lo0] + [lo0 + r] * (k - 1), [hi0] + [hi0 - r] * (k - 1), [w0] + [PLB.tget(tab["w"], q) for q in range(p + 1, m_end + 1)], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"])
                dv = 0.0 if pl is None else -float(pl["nu"][0]) * w0; dv = dv if abs(dv) >= P["dp_min_burn_mps"] else 0.0
            else:
                bias = 0.0
                if olr and p - 1 >= int(P["update_window"]) + 1:
                    K = int(P["update_window"]); idx = np.arange(p - K, p)
                    r = np.array([x[i - 1] - (x[0] + PLB.tget(tab["path"], i) - tab["path"][0]) for i in idx]); bias = float(np.polyfit(idx, r, 1)[0])
                dv, _ = PLB.decide(sw, t, P, p, z0, w0, bias, fast=olr)
            if dv != 0.0:
                dvs.append(dv); u += -dv / w0
        z[p - 1] = x[p - 1] + u
    dvs = np.array(dvs); sg = np.sign(dvs)
    return dict(dv=float(np.abs(dvs).sum()), burns=int(len(dvs)), rev=int((np.diff(sg) != 0).sum()) if len(dvs) > 1 else 0, up=float(dvs[dvs > 0].sum()), down=float(-dvs[dvs < 0].sum()), z=z)


def med(v):
    return float(np.median(v))


def data(root):
    import pandas as pd
    from scipy.stats import norm
    t0 = time.time(); task = _task(); cfg = PLB.load_cfg(); R = {}
    # ---------------- D1/D2/D3 the pass-maximum statistic
    P = stageo_passes(task); vp = np.sqrt(MU * (2 / P.rp_m - 1 / P.a_m))
    es_vp = np.log(P.rho_peak_kg_m3 / P.rho * (vp / P.vrel_mps) ** 3).values       # R1.1 definition (Q_s on the inertial V_p, q_mean on V_rel)
    es = np.log(P.rho_peak_kg_m3 / P.rho).values                                    # U15 definition (both on V_rel)
    lr = np.log(P.r_pert).values; dd = es - lr
    sig_t = (P.rp_m * np.sqrt(P.h_local_m / (MU * P.e))).values
    from scipy.stats import skew, kurtosis
    R["D1_measured"] = dict(n=int(len(P)), runs=int(P.run.nunique()), eps_s_vp=dict(mean=float(es_vp.mean()), median=med(es_vp), sd=float(es_vp.std()), q99=float(np.quantile(es_vp, .99)), max=float(es_vp.max())),
                            eps_s_vrel=dict(mean=float(es.mean()), median=med(es), sd=float(es.std()), q99=float(np.quantile(es, .99)), q999=float(np.quantile(es, .999)), max=float(es.max())),
                            ln_rpert=dict(mean=float(lr.mean()), sd=float(lr.std()), q99=float(np.quantile(lr, .99)), skew=float(skew(lr)), excess_kurtosis=float(kurtosis(lr))),
                            ln_Qs_over_qeff=dict(mean=float(dd.mean()), median=med(dd), ratio_median=float(math.exp(med(dd))), sd=float(dd.std()), q05=float(np.quantile(dd, .05)), min=float(dd.min()), corr_with_ln_rpert=float(np.corrcoef(dd, lr)[0, 1])),
                            sigma_t_s=dict(median=med(sig_t), min=float(sig_t.min()), max=float(sig_t.max())))
    rng = np.random.default_rng(20261027); t, Z = ou_paths(len(P), 22.2, 1.0, 200.0, rng); i0 = len(t) // 2; srt = np.sort(lr)
    models = {}
    for lab, X in (("gaussian_m-0.005_s0.277", -0.005 + 0.277 * Z), ("gaussian_twin_m-0.019_s0.298", -0.019 + 0.298 * Z), ("measured_marginal_copula", np.quantile(srt, norm.cdf(Z)))):
        for G in (0.0, 0.0554):
            M = pass_max(sig_t, X, t, G); D = M - X[:, i0]
            bins = {}
            for lo_, hi_ in ((0.3, 0.5), (0.5, 0.7), (0.7, 0.8)):
                sel = (P.e.values >= lo_) & (P.e.values < hi_); bins[f"e{lo_}-{hi_}"] = dict(model=float(M[sel].mean()), measured=float(es[sel].mean()), n=int(sel.sum()))
            models[f"{lab}_G{G}"] = dict(mean=float(M.mean()), median=med(M), sd=float(M.std()), q99=float(np.quantile(M, .99)), max=float(M.max()), ratio_median=float(math.exp(med(D))),
                                         sd_Qs_over_qeff=float(D.std()), corr=float(np.corrcoef(D, X[:, i0])[0, 1]), by_e=bins)
    tau_s = {}
    for tau in (11.1, 22.2, 44.4, 60.0):                                          # sensitivity: Stage F's "decorrelates over about 60 s" against the twin's tau_c 22.2 s
        tt, Zt = ou_paths(len(P), tau, 1.0, 200.0, np.random.default_rng(7)); Xt = np.quantile(srt, norm.cdf(Zt)); M = pass_max(sig_t, Xt, tt, 0.0); D = M - Xt[:, len(tt) // 2]
        tau_s[str(tau)] = dict(mean=float(M.mean()), sd=float(M.std()), ratio_median=float(math.exp(med(D))), sd_Qs_over_qeff=float(D.std()), corr=float(np.corrcoef(D, Xt[:, len(tt) // 2])[0, 1]))
    models["copula_G0_tau_sensitivity"] = tau_s
    R["D2_pass_max_models"] = models
    # ---------------- D4 the time-risk frontier (R1.2 reproduction)
    q99 = float(np.quantile(es_vp, .99)); beta = 0.045; body = es_vp[es_vp <= q99]; fb = len(body) / len(es_vp); yy = np.linspace(0, 12 * beta, 400); wy = np.exp(-yy / beta) / beta
    def tail(x, sp):
        if sp == 0:
            return float(np.mean(es_vp > x)) if x <= q99 else 0.01 * math.exp(-(x - q99) / beta)
        return float(np.mean(norm.sf((x - body) / sp)) * fb + np.trapz(norm.sf((x - q99 - yy) / sp) * wy, yy) * (1 - fb))
    def qplan(E, sp, N0=200.0, qref=0.224 / 1.372):
        lo, hi = 0.05, 0.54
        for _ in range(60):
            q = math.sqrt(lo * hi); N = N0 * qref / q
            if tail(math.log(0.54 / q), sp) * N > E: hi = q
            else: lo = q
        return q
    fr = {}
    for sp in (0.0, 0.0452, 0.10, 0.15):
        fr[str(sp)] = {str(E): dict(q_plan=qplan(E, sp), days=64.2 * (0.224 / 1.372) / qplan(E, sp), passes=200 * (0.224 / 1.372) / qplan(E, sp)) for E in (0.05, 0.5, 1, 2, 3)}
    R["D4_frontier"] = fr
    # dead-band frontier at the dev cap (0.234): occupancy of a driftless (or drifting) reflected walk in [hi - W, hi]; time factor E[e^-y], loss E[-y], risk per pass, dV per campaign
    hi_d = math.log(0.234); sig, sumw, npass = 0.0628, 160.9, 195                    # sig = sqrt(var rate at lag 20, dev L1 median 0.00394); sumw, passes: dev L1 medians (filled from D6 below if available)
    db = {}
    for W in (0.2, 0.4, 0.6202, 0.7186, 1.0335):
        for dlt in (0.0, 0.0019):
            th = 2 * dlt / sig ** 2; ys = np.linspace(-W, 0, 201)
            p = np.full_like(ys, 1.0 / W) if th == 0 else th * np.exp(th * ys) / (1 - math.exp(-th * W)); p /= np.trapz(p, ys)
            risk = np.trapz(p * np.array([tail(math.log(0.54) - hi_d - y, 0.0452) for y in ys]), ys)
            db[f"W{W}_d{dlt}"] = dict(time_factor=float(np.trapz(p * np.exp(-ys), ys)), mean_loss=float(np.trapz(-ys * p, ys)), risk_per_pass=float(risk), E=float(risk * npass),
                                      dv_campaign=float(sumw * tv_rate(dlt, sig, W)), exit_time_from_mid=float((W / 2) ** 2 / sig ** 2))
    R["D4_deadband"] = dict(sigma=sig, sumw=sumw, passes=npass, rows=db, retarget_every_pass_dv=float(sumw * math.sqrt(2 / math.pi) * 0.114))
    # ---------------- planner dev flights (read-only): lever, time law, natural paths, replays
    TD = dev_tables(task); TC = json.loads((task / "planner_ground_tables.json").read_text())["geometries"]
    Pd = Pcfg(cfg, 0.234); Pc = Pcfg(cfg, 0.2503); SW = {a: PLB.arm_switches(cfg, a) for a in ("L0", "L1", "L2", "L3", "OLR")}
    # D9 lever along the flights
    F1 = flight(root, "ODY", "L1", 1)
    R["D9_lever"] = dict(first=dict(a=float(F1["a"][0]), e=float(F1["ecc"][0]), lever=lever(F1["a"][0], F1["ecc"][0]), H_implied=float(lever(F1["a"][0], F1["ecc"][0]) * F1["w0"][0])),
                         last=dict(a=float(F1["a"][-1]), e=float(F1["ecc"][-1]), lever=lever(F1["a"][-1], F1["ecc"][-1]), H_implied=float(lever(F1["a"][-1], F1["ecc"][-1]) * F1["w0"][-1])),
                         w0_first=float(F1["w0"][0]), w0_last=float(F1["w0"][-1]))
    # D5 time law and K
    tl = []
    for d in sorted(glob.glob(str(pathlib.Path(root) / "planner" / "dev" / "ODY_PL_*"))):
        Pp = pd.read_csv(pathlib.Path(d) / "passes.csv", usecols=["time_s", "q_w_cm2", "a_m"]); T = 2 * np.pi * np.sqrt(Pp.a_m.values ** 3 / MU); q = Pp.q_w_cm2.values
        qH = np.sum(T) / np.sum(T / q); tl.append(dict(arm=pathlib.Path(d).name.split("_PL_")[1].split("_")[0], days=Pp.time_s.values[-1] / 86400, qH=qH, qA=q.mean(), qG=math.exp(np.mean(np.log(q)))))
    tl = pd.DataFrame(tl); tl["days_qH"] = tl.days * tl.qH; tl["days_qA"] = tl.days * tl.qA; tl["pen"] = np.log(tl.qG / tl.qH)
    Kd = float(tl.days_qH.median())
    R["D5_time_law"] = dict(K_days_W_cm2=Kd, by_arm=tl.groupby("arm")[["days", "qH", "qA", "days_qH", "days_qA", "pen"]].median().to_dict(orient="index"),
                            spread_days_qH=float(tl.groupby("arm").days_qH.median().max() / tl.groupby("arm").days_qH.median().min() - 1),
                            spread_days_qA=float(tl.groupby("arm").days_qA.median().max() / tl.groupby("arm").days_qA.median().min() - 1))
    # D6/D7/D8 per flight
    rows = []; fe_rows = []
    TDo, TCo = TD["ODY"], TC["ODY"]; fast = np.asarray(TDo["fast"]); path = np.asarray(TDo["path"]); slow = np.asarray(TDo["slow"])
    for arm in ("L1", "L3"):
        for s in range(1, 9):
            F = flight(root, "ODY", arm, s); x = F["x"]; n = F["n"]; hi = math.log(0.234); W = TDo["W"]
            wv = np.r_[F["w"][2], [F["w"][p] for p in range(2, n + 1)]][:n]
            A, B = hi - W - x[1:], hi - x[1:]
            same = replay(F, TDo, Pd, SW[arm])
            ex = PLB.dp_exact(A, B, wv[1:]); import oracle_lib as OLB; h = W / 40; gr = OLB.dp_schedule(A, B, wv[1:], res=h)
            bnd = h * float(np.sum(wv[1:][np.abs(np.diff(np.r_[0.0, ex["nu"]])) > 1e-12]))
            r = dict(arm=arm, seed=s, n=n, flown=F["flown"], flown_all=F["flown_all"], replay_same=same["dv"], days=F["days"],
                     floor_full=ex["obj"], floor_grid_W40=float(gr["obj"]), grid_bound=bnd, floor_reserve=PLB.dp_exact(A + 0.1285, B - 0.1285, wv[1:])["obj"], clamp_truth=clamp_cost(A, B, wv[1:])[0],
                     ep_sd=float(np.std(list(F["e"].values()))), ep_mean=float(np.mean(list(F["e"].values()))))
            vg = {L: float(np.var(x[L:] - x[:-L]) / L) for L in (1, 10, 20, 30)}; dr = (x[-1] - x[1]) / (n - 2); sumw = float(np.sum(wv[1:]))
            r.update(var_rate_1=vg[1], var_rate_10=vg[10], var_rate_20=vg[20], var_rate_30=vg[30], drift=dr, sumw=sumw,
                     pred_clamp_bm=sumw * tv_rate(dr, math.sqrt(vg[20]), W, rho=0.0), pred_clamp_bm_siegmund=sumw * tv_rate(dr, math.sqrt(vg[20]), W),
                     pred_clamp_w10=sumw * tv_rate(dr, math.sqrt(vg[10]), W), pred_clamp_w30=sumw * tv_rate(dr, math.sqrt(vg[30]), W),
                     pred_floor_reserve_bm=sumw * tv_rate(dr, math.sqrt(vg[20]), W - 0.257, rho=0.0))
            for a2 in ("L1", "L2", "L3"):
                rr = replay(F, TDo, Pd, SW[a2]); r[f"dev_{a2}"] = rr["dv"]; r[f"dev_{a2}_burns"] = rr["burns"]; r[f"dev_{a2}_rev"] = rr["rev"]; r[f"dev_{a2}_down"] = rr["down"]
            r["dev_L3_nofast"] = replay(F, TDo, Pd, SW["L3"], fast_zero=True)["dv"]; r["dev_L3_pen0"] = replay(F, TDo, Pcfg(cfg, 0.234, pen=0.0), SW["L3"])["dv"]
            r["dev_L3_res0"] = replay(F, TDo, Pcfg(cfg, 0.234, res=0.0), SW["L3"])["dv"]; r["dev_L3_oracle_res"] = replay(F, TDo, Pd, SW["L3"], oracle="res")["dv"]; r["dev_L3_oracle_full"] = replay(F, TDo, Pd, SW["L3"], oracle="full")["dv"]
            r["dev_L3_e0"] = replay(dict(F, e={k: 0.0 for k in F["e"]}), TDo, Pd, SW["L3"])["dv"]; r["dev_L1_e0"] = replay(dict(F, e={k: 0.0 for k in F["e"]}), TDo, Pd, SW["L1"])["dv"]
            r["dev_OLR"] = replay(F, TDo, Pd, SW["OLR"], olr=True)["dv"]
            for a2 in ("L0", "L1", "L2", "L3"):
                rr = replay(F, TCo, Pc, SW[a2]); r[f"cal_{a2}"] = rr["dv"]; r[f"cal_{a2}_burns"] = rr["burns"]; r[f"cal_{a2}_rev"] = rr["rev"]
                qH = np.sum(F["T"]) / np.sum(F["T"] / np.exp(rr["z"])); r[f"cal_{a2}_days"] = Kd / qH
            r["cal_L3Z"] = replay(F, TCo, Pcfg(cfg, 0.2503, pen=0.0), SW["L3"])["dv"]
            rr = replay(F, TCo, Pc, SW["OLR"], olr=True); r["cal_OLR"] = rr["dv"]; r["cal_OLR_burns"] = rr["burns"]; r["cal_OLR_rev"] = rr["rev"]
            r["cal_OLRZ"] = replay(F, TCo, Pcfg(cfg, 0.2503, pen=0.0), SW["OLR"], olr=True)["dv"]
            rr = replay(F, TDo, Pd, SW["L1"]); qH = np.sum(F["T"]) / np.sum(F["T"] / np.exp(rr["z"])); r["dev_L1_days_pred"] = Kd / qH
            rows.append(r)
            # D8 forecast error of the table and the anchor jitter
            N = min(n, len(path)); zp = dict(zip(range(2, n + 2), F["z0"]))
            nu = F["lnq"] - x; anc = np.array([zp[p] - nu[p - 2] - fast[p - 1] for p in range(2, N + 1) if p in zp]); zz = np.array([zp[p] - nu[p - 2] for p in range(2, N + 1) if p in zp])
            fe = {tau: float(np.std([(x[k + tau] - x[k]) - (path[k + tau] - path[k]) for k in range(1, N - tau)])) for tau in (1, 5, 10, 30)}
            fe30 = [x[k + 30] - (anc[k - 1] + slow[k + 30] - slow[k]) for k in range(1, min(N - 30, len(anc) + 1))]
            fe_rows.append(dict(arm=arm, seed=s, corr_dx_dfast=float(np.corrcoef(np.diff(x[1:N]), np.diff(fast[1:N]))[0, 1]), anchor_step_sd=float(np.std(np.diff(anc) - np.diff(slow[1:len(anc) + 1]))),
                                z_step_sd=float(np.std(np.diff(zz))), **{f"fe_path_{k}": v for k, v in fe.items()}, fe_anchor30_sd=float(np.std(fe30)), fe_anchor30_mean=float(np.mean(fe30))))
    df = pd.DataFrame(rows); fdf = pd.DataFrame(fe_rows)
    assert np.allclose(df.replay_same, df.flown, atol=1e-9), "replay does not reproduce the flown burns"
    R["D6_D7_per_flight"] = df.to_dict(orient="records"); R["D6_D7_medians"] = df.groupby("arm").median(numeric_only=True).to_dict(orient="index")
    R["D8_forecast"] = dict(per_flight=fdf.to_dict(orient="records"), medians=fdf.groupby("arm").median(numeric_only=True).to_dict(orient="index"))
    R["D6_corr_pred_vs_flown_L1"] = float(np.corrcoef(df[df.arm == "L1"].pred_clamp_bm_siegmund, df[df.arm == "L1"].flown)[0, 1])
    # G1 dev L3 flights
    g = []
    for s in range(1, 9):
        F = flight(root, "G5863", "L3", s); r = dict(seed=s, flown=F["flown"], replay_same=replay(F, TD["G5863"], Pd, SW["L3"])["dv"])
        for a2 in ("L1", "L2", "L3"):
            r[f"dev_{a2}"] = replay(F, TD["G5863"], Pd, SW[a2])["dv"]; r[f"cal_{a2}"] = replay(F, TC["G5863"], Pc, SW[a2])["dv"]
        r["cal_L3Z"] = replay(F, TC["G5863"], Pcfg(cfg, 0.2503, pen=0.0), SW["L3"])["dv"]; r["cal_OLR"] = replay(F, TC["G5863"], Pc, SW["OLR"], olr=True)["dv"]
        r["cal_OLRZ"] = replay(F, TC["G5863"], Pcfg(cfg, 0.2503, pen=0.0), SW["OLR"], olr=True)["dv"]
        g.append(r)
    g = pd.DataFrame(g); assert np.allclose(g.replay_same, g.flown, atol=1e-9)
    R["D7_G1"] = dict(per_flight=g.to_dict(orient="records"), medians=g.median(numeric_only=True).to_dict())
    # D10 N2-I / F-O lever weight sums (Stage O dev) for the re-aim heuristic
    sw_ = {}
    for arm, pat in (("N2I", "N2I_*"), ("FO", "FO_*")):
        sums = []
        for f in sorted(glob.glob(str(task / "stageo_data" / "dev" / "runs" / pat / "passes.csv"))):
            Q = pd.read_csv(f, usecols=["a_m", "e", "h_local_m"]); sums.append(float(np.sum(Q.h_local_m / [lever(a, e) for a, e in zip(Q.a_m, Q.e)])))
        sw_[arm] = dict(median_sumw=med(sums), n=len(sums))
    s_ = 0.069; W_, g_ = 0.827 - 0.61, 0.10; rate, freq = reaim_rate(s_, W_, g_)
    R["D10_reaim_heuristic"] = dict(sumw=sw_, assumed_W=W_, assumed_g=g_, sigma=s_, rate_ln_per_pass=rate, burns_per_pass=freq, pred_dv=rate * sw_["N2I"]["median_sumw"],
                                    pred_burns_208=freq * 208, flown_dv=4.21, flown_burns=26)
    R["runtime_s"] = time.time() - t0
    out = task / "analytic_data" / "part2_results.json"; out.write_text(json.dumps(R, indent=1, default=float))
    print(f"wrote {out} ({R['runtime_s']:.0f} s)")
    pd.set_option("display.width", 250)
    print(df.groupby("arm").median(numeric_only=True).T.round(3).to_string()); print(fdf.groupby("arm").median(numeric_only=True).T.round(3).to_string()); print(g.median().round(3).to_string())
    return R


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--data", action="store_true"); ap.add_argument("--results-root", default=RESULTS_ROOT)
    a = ap.parse_args()
    if not (a.selftest or a.data):
        ap.error("--selftest | --data")
    if a.selftest:
        selftest()
    if a.data:
        data(a.results_root)


if __name__ == "__main__":
    main()
