"""Apoapsis-targeting Pi parameters: closed forms and reference quadratures.

Derivation: tasks/2026-10-15_periapsis-region-border-targeting/apoapsis_pi_derivation.md (local note).
SI units throughout. beta = m / (C_D A) [kg/m2], Q = 0.5 rho V^3 [W/m2], H scale height [m].
"""
import numpy as np
from scipy.integrate import quad


# ---- closed forms -----------------------------------------------------------------------------

def tau_eff(rp, e, H, mu):
    """Laplace pass duration: int (Q/Q_p) dt = r_p sqrt(2 pi H / (mu e)) [s] (leading order in H/(r_p e))."""
    return rp * np.sqrt(2 * np.pi * H / (mu * e))


def delta_a_pass(a, e, rp, Q, H, beta, mu, f_rot=1.0):
    """Change of semi-major axis over one drag pass [m] (negative). Energy identity d(eps)/dt = -f_rot Q / beta."""
    return -2 * a**2 * f_rot * Q * tau_eff(rp, e, H, mu) / (mu * beta)


def pi_E(a, e, rp, Q, H, beta, mu, f_rot=1.0):
    """Energy depletion per pass |d eps|/|eps| = |da|/a."""
    return -delta_a_pass(a, e, rp, Q, H, beta, mu, f_rot) / a


def pi_q(rp, Q, beta, mu):
    """Heat-rate drag group Q r_p / (beta V_c^3), V_c = sqrt(mu / r_p)."""
    return Q * rp / (beta * (mu / rp) ** 1.5)


def k_Q(rp, Q, H, beta, mu, f_rot=1.0):
    """Per-pass decrement of e^(3/2) is 3 k_Q at constant Q and r_p; k_Q = Pi_q sqrt(2 pi H / r_p)."""
    return f_rot * pi_q(rp, Q, beta, mu) * np.sqrt(2 * np.pi * H / rp)


def passes_to_go(e, e_t, kq):
    """N = (e^1.5 - e_t^1.5) / (3 k_Q)."""
    return (e**1.5 - e_t**1.5) / (3 * kq)


def g_time(e):
    """Time primitive: t = T_cp [g(e0) - g(et)] / k_Q, g(e) = sqrt(e/(1-e)) - asin(sqrt(e))."""
    return np.sqrt(e / (1 - e)) - np.arcsin(np.sqrt(e))


def time_to_go(rp, e, e_t, kq, mu):
    """Campaign time [s] at constant Q and r_p; T_cp = 2 pi sqrt(r_p^3/mu) the circular period at r_p."""
    return 2 * np.pi * np.sqrt(rp**3 / mu) * (g_time(e) - g_time(e_t)) / kq


def lam_rp(e, piE):
    """Drag-driven periapsis drop per pass in scale heights (positive = drop), leading order:
    |dr_p| / H = Pi_E (1 - e) / (2 e (1 + e))."""
    return piE * (1 - e) / (2 * e * (1 + e))


def n_life(lam, Q, Q_hard):
    """Passes until the peak heat rate reaches Q_hard with no maneuver (drag-only periapsis decay, dQ/dk = Lambda Q)."""
    return (1 - Q / Q_hard) / lam


def lever_rp(a, e):
    """d r_p / d Pi_dv for a tangential apoapsis burn, Pi_dv = dv / sqrt(mu / a): 4 a sqrt((1 - e)/(1 + e))."""
    return 4 * a * np.sqrt((1 - e) / (1 + e))


def corridor_width_H(q_lo, q_hi):
    """Width of a heat-rate corridor in scale heights at fixed V: ln(Q_hi / Q_lo)."""
    return np.log(q_hi / q_lo)


def f_rot(r, a, vrel, lat, omega, mu):
    """Drag power factor (v . v_rel) / V_rel^2 at the periapsis point: power per unit mass = f_rot Q / beta."""
    v2 = mu * (2 / r - 1 / a); w = omega * r * np.cos(lat)
    return (v2 + vrel**2 - w**2) / (2 * vrel**2)


def elements_from_x(rp, x):
    """a, e from x = r_a / r_p."""
    return rp * (1 + x) / 2, (x - 1) / (x + 1)


# ---- reference quadratures (Kepler ellipse, exponential density, non-rotating atmosphere) ------

def delta_a_quad(a, e, H, rho_p, beta, mu):
    """Exact Kepler-path integral of da/dt = -(a^2/(mu beta)) rho V^3 over one orbit (first order in drag)."""
    p = a * (1 - e**2); rp = a * (1 - e); h = np.sqrt(mu * p)
    def f(nu):
        r = p / (1 + e * np.cos(nu)); v3 = (mu * (2 / r - 1 / a)) ** 1.5
        return rho_p * np.exp(-(r - rp) / H) * v3 * r**2 / h
    w = 60 * np.sqrt(H / (rp * e)) if e > 0 else np.pi
    w = min(w, np.pi)
    I = 2 * quad(f, 0, w, limit=400, epsabs=0, epsrel=1e-12)[0]
    return -(a**2 / (mu * beta)) * I


def gauss_pass_quad(a, e, H, rho_p, beta, mu):
    """First-order Gauss integrals (tangential drag) over one orbit: returns (da, de, drp)."""
    p = a * (1 - e**2); rp = a * (1 - e); h = np.sqrt(mu * p)
    def acc(nu):
        r = p / (1 + e * np.cos(nu)); v2 = mu * (2 / r - 1 / a); v = np.sqrt(v2)
        return rho_p * np.exp(-(r - rp) / H) * v2 / (2 * beta), r, v
    def fa(nu):
        aT, r, v = acc(nu); return -(2 * a**2 * v / mu) * aT * r**2 / h
    def fe(nu):
        aT, r, v = acc(nu); return -(2 / v) * (e + np.cos(nu)) * aT * r**2 / h
    w = min(60 * np.sqrt(H / (rp * e)), np.pi)
    da = 2 * quad(fa, 0, w, limit=400, epsabs=0, epsrel=1e-12)[0]
    de = 2 * quad(fe, 0, w, limit=400, epsabs=0, epsrel=1e-12)[0]
    return da, de, (1 - e) * da - a * de
