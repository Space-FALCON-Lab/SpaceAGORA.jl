"""Stage M library: apoapsis-targeting controller formulas as pure functions (Python only, no Julia).

Plan: tasks/2026-10-15_periapsis-region-border-targeting/apoapsis_controller_plan.md (local note), sections 2.1, 2.4, 3.4, A.2 item 3, A.7.
SI units throughout, except where a name says otherwise: heat rates enter drop_per_pass in W/m2 (Q = 0.5 rho V^3) and the end-game rule in
W/cm2 (Q_main = 0.15), joined by drop_unit_m, the apoapsis drop (m) per pass for Q = 1 W/cm2 (= drop_per_pass(.., q_w_m2=1e4, ..)).
The Julia mirror (stagem_common.jl) must reproduce these functions to 1e-12.
"""
import math
from statistics import NormalDist

import numpy as np

Z95 = NormalDist().inv_cdf(0.95)   # two-piece-normal 5/95% quantile factor, 1.6449 (the planning arithmetic used the 4-digit value)


# ---- decay map (empirical form of Pi_E) -----------------------------------------------------------

def drop_per_pass(a, e, rp, q_w_m2, K, mu):
    """Apoapsis drop [m] per pass at fixed r_p: K a^2 Q r_p / (sqrt(mu e) mu) (= -2 da, since r_a + r_p = 2a)."""
    return K * a * a * q_w_m2 * rp / (math.sqrt(mu * e) * mu)


def k_theory(H, beta, f_rot=1.0):
    """Theory gain of the decay map: K = 4 f_rot sqrt(2 pi H) / beta (H [m], beta = m / (C_D A) [kg/m2]); from -2 da = 4 a^2 f_rot Q tau / (mu beta)."""
    return 4 * f_rot * math.sqrt(2 * math.pi * H) / beta


def k_update(k_prev, obs_drop_m, a, e, rp, q_eff_w_m2, mu, span=10):
    """EWMA of the observed gain: k_obs = obs_drop / (a^2 q_eff r_p / (sqrt(mu e) mu)), alpha = 2 / (span + 1).

    obs_drop_m = -(change of r_a + r_p over the pass) >= 0 (a, e, rp: elements at the start of the pass; q_eff: measured peak heat rate).
    Unusable observations (non-finite, no drag, non-positive drop) leave K unchanged.
    """
    unit = a * a * q_eff_w_m2 * rp / (math.sqrt(mu * e) * mu)
    if not (math.isfinite(obs_drop_m) and math.isfinite(unit)) or unit <= 0 or obs_drop_m <= 0:
        return k_prev
    alpha = 2.0 / (span + 1)
    return k_prev + alpha * (obs_drop_m / unit - k_prev)


# ---- end game -------------------------------------------------------------------------------------

def q_end(delta_goal_m, e_hi, e_lo, drop_unit_m, q_min_abs, q_main):
    """Q_end [W/cm2] = clip(delta_goal / (u90 drop_unit), q_min_abs, q_main), u90 = max(exp(e_hi) - 1, 1 - exp(e_lo)).

    e_hi > 0 and e_lo < 0 are the planet's 95% / 5% one-pass ln-density error quantiles (Stage I refit).
    """
    u90 = max(math.exp(e_hi) - 1.0, 1.0 - math.exp(e_lo))
    return min(max(delta_goal_m / (u90 * drop_unit_m), q_min_abs), q_main)


def endgame_decision(D_m, drop_unit_m, q_main, q_end, q_min):
    """Plan 2.1 decision rule at apoapsis k; D = r_a,k - r_a,t, Q_need = D / drop_unit [W/cm2]. Returns (action, q_plan, phase).

    STOP if D < drop(q_min) / 2 (q_plan 0, phase "last"); else fly:
      LAST        Q_need <= q_end         q_plan = clip(Q_need, q_min, q_main)
      PENULTIMATE Q_need <= q_main + q_end q_plan = clip(Q_need - q_end, q_end, q_main)   (leaves at most q_end for the last pass)
      MAIN        otherwise                q_plan = q_main
    Deviation from the plan text, which reads "LAST if Q_need <= Q_main": that order can fly the last pass at up to Q_main (e.g. Q_need 0.20
    -> MAIN leaves 0.05 -> LAST at 0.05, not <= q_end) and cannot give the expected last-pass Q of q_end. Gating LAST on q_end is the
    reading that matches the stated goal (the last pass removes at most d_last = delta_goal / u90). Logged in the project log.
    """
    if D_m < drop_unit_m * q_min / 2.0:
        return "stop", 0.0, "last"
    q_need = D_m / drop_unit_m
    if q_need <= q_end:
        return "fly", min(max(q_need, q_min), q_main), "last"
    if q_need <= q_main + q_end:
        return "fly", min(max(q_need - q_end, q_end), q_main), "penultimate"
    return "fly", q_main, "main"


def stop_nearest(D_m, d_obs_last_m):
    """Baseline stop S0 (arms without Q control): stop when D < d_obs,last / 2 (persistence of the last observed drop)."""
    return D_m < d_obs_last_m / 2.0


# ---- final raise and trim (impulsive, analytic) ---------------------------------------------------

def _v_apo(ra, rp, mu):
    """Speed at apoapsis of the ellipse (ra, rp), vis-viva."""
    return math.sqrt(2 * mu * rp / (ra * (ra + rp)))


def _v_peri(ra, rp, mu):
    return math.sqrt(2 * mu * ra / (rp * (ra + rp)))


def raise_dv(ra, rp, rp_final, mu):
    """Apoapsis burn [m/s] raising r_p to rp_final at fixed r_a: v_a(rp_final) - v_a(rp)."""
    return _v_apo(ra, rp_final, mu) - _v_apo(ra, rp, mu)


def rp_final_from_raise(ra, rp, dv, mu):
    """Inverse of raise_dv in rp_final (bisection on [rp, ra]; raise_dv increases with rp_final)."""
    lo, hi = rp, ra
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if raise_dv(ra, rp, mid, mu) < dv: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def trim_dv(ra_end, ra_t, rp_final, mu):
    """Periapsis burn [m/s, magnitude] at rp_final moving the apoapsis from ra_end to ra_t."""
    return abs(_v_peri(ra_t, rp_final, mu) - _v_peri(ra_end, rp_final, mu))


def trim_lever(ra, rp, mu):
    """Apoapsis change [m] per m/s of periapsis burn: 4 a sqrt((1 + e) / (1 - e)) / sqrt(mu / a)."""
    a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
    return 4 * a * math.sqrt((1 + e) / (1 - e)) / math.sqrt(mu / a)


# ---- campaign map ---------------------------------------------------------------------------------

def campaign_map(ra0, ra_t, rp, q_w_m2, K, mu, max_passes=100000):
    """Pass-by-pass map at constant Q and r_p: returns (passes, days, drop_at_target [m]); days sum the orbital periods."""
    ra, n, t = ra0, 0, 0.0
    while ra > ra_t and n < max_passes:
        a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
        ra -= drop_per_pass(a, e, rp, q_w_m2, K, mu); t += 2 * math.pi * math.sqrt(a ** 3 / mu); n += 1
    a = 0.5 * (ra_t + rp); e = (ra_t - rp) / (ra_t + rp)
    return n, t / 86400.0, drop_per_pass(a, e, rp, q_w_m2, K, mu)


# ---- final-apoapsis precision model (3.4, A.7) ----------------------------------------------------

def _two_piece(rng, q05, q95, n):
    """Two-piece normal with the given 5% (<= 0) and 95% (>= 0) quantiles."""
    g = rng.standard_normal(n)
    return g * np.where(g < 0, -q05 / Z95, q95 / Z95)


def precision_mc(d_last_m, e_q05, e_q95, t_q05_m, t_q95_m, n, seed, s0=False):
    """Monte Carlo of the pre-trim final-apoapsis error: err = d - d exp(eps) + t (positive = apoapsis still too high).

    eps: two-piece normal on the ln-density error (e_q05 < 0 < e_q95); t: two-piece normal tesseral floor [m] (t_q05_m <= 0 <= t_q95_m).
    s0=True models the baseline stop without Q control: the distance left before the last decision is U(0, d); no pass is flown if it is < d / 2.
    Returns dict(within10, within5, mean_abs [m], p90_abs [m]).
    """
    rng = np.random.default_rng(seed)
    eps = _two_piece(rng, e_q05, e_q95, n); t = _two_piece(rng, t_q05_m, t_q95_m, n)
    err = d_last_m - d_last_m * np.exp(eps) + t
    if s0:
        D = rng.random(n) * d_last_m
        err = np.where(D < d_last_m / 2, D, D - d_last_m * np.exp(eps) + t)
    ab = np.abs(err)
    return dict(within10=float(np.mean(ab <= 10e3)), within5=float(np.mean(ab <= 5e3)), mean_abs=float(ab.mean()), p90_abs=float(np.quantile(ab, 0.9)))


# ---- J2 apsidal rotation and periapsis altitude (A.2 item 3, closed form) ------------------------

def apsidal_rotation_per_pass(a, e, inc, J2, R, mu):
    """Secular J2 apsidal rotation per orbit [rad]: 3 pi J2 (R / p)^2 (2 - 5/2 sin^2 i), p = a (1 - e^2) (angle of AOP change per pass)."""
    p = a * (1 - e * e)
    return 3 * math.pi * J2 * (R / p) ** 2 * (2 - 2.5 * math.sin(inc) ** 2)


def periapsis_geodetic_altitude(rp, inc, omega, R_eq, f):
    """Altitude of the periapsis point over the oblate ellipsoid [m]: rp - R_eq (1 - f sin^2(lat)), sin(lat) = sin i sin(omega).

    First order in f and geocentric latitude: the exact conversion differs by meters (Addendum B uses the exact one).
    """
    s = math.sin(inc) * math.sin(omega)
    return rp - R_eq * (1 - f * s * s)
