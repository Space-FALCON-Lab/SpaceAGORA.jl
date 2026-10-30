"""Closed-form per-pass terms for the periapsis altitude (WP-A2/A3 of analytic_framework_and_transfer_plan.md; Python only, SI units).

Terms (plan 3.1): T2 oblate-surface geometry, T4b/T4c tesseral energy and resonance calendar, T5 drag density factor, T6 drag at small c = a e / H
(modified Bessel functions), T8 third-body tide, T9 SRP bound, Pi groups (3.2). Existing libraries are imported, not copied: apoapsis_pi_lib (drag Laplace),
stagek_lib (zonal rates and maps), stagem_lib (apsidal rotation, geodetic altitude), tesseral_lib (field, Gauss rates, Kaula series), kaula_lib.
Labels: SOURCED (file named), DERIVED (derived here, checked numerically in analytic_selftest.py), ASSUMED (named choice), RECALLED (from memory, verified or listed).
"""
import math, pathlib, sys
import numpy as np
from scipy.special import ive

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import apoapsis_pi_lib as AP
import stagem_lib as SM
import stagek_lib as SK
import tesseral_lib as TL

# SOURCED src/environment/ephemerides/planets.jl (Mars, Earth, Venus structs) and pck00011.tpc radii; H = exponential-atmosphere scale height of the struct
PLANETS = {
    "mars": dict(Re=3396190.0, Rp=3376200.0, w_rot=7.08823596e-5, H=6308.278108),
    "earth": dict(Re=6378136.6, Rp=6356751.9, w_rot=7.2921158553e-5, H=8500.0),
    "venus": dict(Re=6051800.0, Rp=6051800.0, w_rot=-2.99e-7, H=15900.0),
}
for _p in PLANETS.values(): _p["f"] = (_p["Re"] - _p["Rp"]) / _p["Re"]
SRP_P1AU = 1361.0 / 299792458.0     # N/m2, radiation pressure at 1 AU for TSI 1361 W/m2 (ASSUMED nominal TSI, not a repo value)
AU = 1.495978707e11                  # m (IAU 2012 exact)


# ---- T2: surface geometry (DERIVED) ----------------------------------------------------------------------------------------------------------
def surface_radius(inc, om, planet, first_order=False):
    """Geocentric radius of the reference ellipsoid under the periapsis point, sin(phi) = sin i sin omega:
    R_s = a b / sqrt(b^2 cos^2 phi + a^2 sin^2 phi) (exact; the Stage O engine altitude equals r - R_s to 0.02-1.5 m), or the first-order form R_eq (1 - f sin^2 phi)
    (differs by up to 34 m at Mars).  Periapsis height over the surface: h_p = r_p - R_s."""
    P = PLANETS[planet]; s = (math.sin(inc) * math.sin(om)) ** 2
    if first_order: return P["Re"] * (1 - P["f"] * s)
    a, b = P["Re"], P["Rp"]; return a * b / math.sqrt(b * b + (a * a - b * b) * s)


def dh_geom(inc, om, d_om, planet, d_inc=0.0, first_order=False):
    """Per-pass change of the geometry term of the periapsis height, h_geom = -R_s: -(R_s(i + di, omega + d_om) - R_s(i, omega)).  Small-step form
    R_eq f sin^2 i sin(2 omega) d_om (first order in f and d_om); the exact difference is used here (no singularity at phi = +-90 deg)."""
    return -(surface_radius(inc + d_inc, om + d_om, planet, first_order) - surface_radius(inc, om, planet, first_order))


def zonal_rates_per_pass(P, a, e, inc, om):
    """T1/T3: first-order orbit-averaged zonal (J2..J20) rates times the Kepler period: returns (d_e, d_omega, d_Omega, d_rp) per pass (rp = -a de, da = 0)."""
    T = 2 * math.pi * math.sqrt(a ** 3 / P["mu"]); ed, wd, Od = SK.zonal_secular(P, a, e, inc, om)
    return ed * T, wd * T, Od * T, -a * ed * T


# ---- T4b/T4c: tesseral energy and crossing calendar --------------------------------------------------------------------------------------------
def theta_p_from_pass(raan, om, inc, lon_body):
    """Planet rotation angle at the periapsis event: inertial longitude of the periapsis point minus its body-fixed longitude (as stagek_tess.py)."""
    Px = math.cos(raan) * math.cos(om) - math.sin(raan) * math.sin(om) * math.cos(inc)
    Py = math.sin(raan) * math.cos(om) + math.cos(raan) * math.sin(om) * math.cos(inc)
    return (math.atan2(Py, Px) - lon_body) % (2 * math.pi)


def tess_arc(cs, Ct, St, a, e, inc, raan, om, theta_p, w_rot, nf=128, with_boundary=False, f_lo=0.0, f_hi=2 * math.pi):
    """First-order change of (a, e, i, Omega, omega, rp) over the arc f in [f_lo, f_hi] (true anomaly; default periapsis to periapsis, f_lo = -pi, f_hi = pi is apoapsis to
    apoapsis) of the Kepler ellipse (a, e, inc, raan, om), from the body-fixed field (Ct, St) with the planet turning at w_rot (theta = theta_p + w_rot t, t = 0 at periapsis).
    Also returns r_a = a(1+e) kick 'ra', dE (specific Kepler energy change, mu da / (2 a^2)) and dhz (change of the z component of the specific angular momentum).
    Identity (DERIVED, checked in analytic_selftest): dE - w_rot dhz = U(end) - U(start) (potential of the perturbing field at the two arc ends, body frame);
    with_boundary adds that dU.  At apoapsis ends U is small, so dE ~ w_rot dhz there (Jacobi)."""
    mu = cs["mu"]; x, w = np.polynomial.legendre.leggauss(nf); half = 0.5 * (f_hi - f_lo); f = f_lo + half * (x + 1.0); w = half * w
    p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f)); N = math.sqrt(mu / a ** 3)
    tau_of = lambda ff: ((lambda EE: (EE - e * np.sin(EE)) / N)(2 * np.arctan2(np.sqrt(1 - e) * np.sin(ff / 2), np.sqrt(1 + e) * np.cos(ff / 2))))
    tau = tau_of(f); Q = TL.pqw_matrix(inc, raan, om)
    r = Q @ np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); v = Q @ (math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f]))
    th = theta_p + w_rot * tau; rb = TL.rot_z(-th, r); ai = TL.rot_z(th, TL.accel_body(rb[0], rb[1], rb[2], cs, Ct, St))
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    aR, aS, aW = (ai * rhat).sum(0), (ai * shat).sum(0), (ai * hhat).sum(0)
    rates = TL.gauss_rates(a, e, inc, om, f, rr, h, aR, aS, aW, mu); k = dict(zip(TL.KICK_NAMES, rates @ (w * rr ** 2 / h)))
    k["ra"] = (1 + e) * k["a"] + a * k["e"]; k["dE"] = mu * k["a"] / (2 * a * a)
    dh = 0.5 * h * (k["a"] / a - 2 * e * k["e"] / (1 - e * e)); k["dhz"] = dh * math.cos(inc) - h * math.sin(inc) * k["i"]
    if with_boundary:
        def U_at(ff):
            rf = p / (1 + e * math.cos(ff)); pos = Q @ np.array([[rf * math.cos(ff)], [rf * math.sin(ff)], [0.0]]); thf = theta_p + w_rot * float(tau_of(np.array(ff)))
            return float(TL.potential_body(*TL.rot_z(-thf, pos), cs, Ct, St)[0])
        k["dU"] = U_at(f_hi) - U_at(f_lo)
    return k


def a_osc_minus_mean(a, e, inc, om, J2, R):
    """First-order J2 short-period difference a_osc - a_mean at periapsis (f = 0):
    (J2 R^2 / a) { [(a/r)^3 - (1-e^2)^(-3/2)] (1 - (3/2) sin^2 i) + (3/2) (a/r)^3 sin^2 i cos(2 omega) }, r = a (1-e)   (RECALLED Kozai 1959 form;
    verified numerically: it takes the median error of the periapsis-to-periapsis commensurability nu on the Stage O tables from 1.4% to 0.08%)."""
    r = a * (1 - e); s2 = math.sin(inc) ** 2
    return (J2 * R * R / a) * (((a / r) ** 3 - (1 - e * e) ** -1.5) * (1 - 1.5 * s2) + 1.5 * (a / r) ** 3 * s2 * math.cos(2 * om))


def nu_commens(a, mu, w_rot):
    """Commensurability nu = w_rot T / (2 pi) with the Kepler period of the given semi-major axis (a MEAN value; see nu_from_osc for osculating periapsis values)."""
    return w_rot * math.sqrt(a ** 3 / mu)


def nu_from_osc(a0, e0, om0, a1, e1, om1, inc, mu, w_rot, J2, R):
    """nu of the interval between two periapsis events from the osculating elements at its two ends: mean a at each end (Kozai short-period removal), then the midpoint
    (drag changes a by ~1.5% per pass, so each event value is mid-drag and the interval mean is between the two)."""
    am0 = a0 - a_osc_minus_mean(a0, e0, inc, om0, J2, R); am1 = a1 - a_osc_minus_mean(a1, e1, inc, om1, J2, R)
    return nu_commens(0.5 * (am0 + am1), mu, w_rot)


def nu_of_e(rp, e, mu, w_rot):
    """nu as a function of e at fixed r_p: Pi_omega (1-e)^(-3/2), Pi_omega = w_rot sqrt(r_p^3/mu)."""
    return w_rot * math.sqrt(rp ** 3 / mu) * (1 - e) ** -1.5


def crossing_calendar(rp, e_series, mu, w_rot, mmax=4):
    """Resonance-crossing calendar from a decay history e_k (one value per pass): returns list of (m, j, pass index k_c (fractional), nu'(per pass)) where m nu(e_k) = j integer."""
    nu = np.array([nu_of_e(rp, e, mu, w_rot) for e in e_series]); out = []
    for m in range(1, mmax + 1):
        x = m * nu
        for j in range(int(math.ceil(min(x[0], x[-1]))), int(math.floor(max(x[0], x[-1]))) + 1):
            s = np.where(np.diff(np.sign(x - j)) != 0)[0]
            for i in s:
                frac = (j - x[i]) / (x[i + 1] - x[i]); out.append((m, j, i + frac, (nu[i + 1] - nu[i])))
    return sorted(out, key=lambda t: t[2])


# ---- T5: drag with a lognormal density error (DERIVED) ----------------------------------------------------------------------------------------
def drag_density_factor(sigma):
    """E[exp(y)] for y ~ N(mu_y, sigma^2) relative to exp(mu_y): exp(sigma^2 / 2); the mean decay is this factor times the median-density decay."""
    return math.exp(0.5 * sigma ** 2)


def drag_drop_ra_m(a, e, rp, Q_w_m2, H, beta, mu, f_rot=1.0):
    """T5: apoapsis (r_a + r_p = 2a) change per pass [m] from the Laplace form: 2 delta_a (negative)."""
    return 2 * AP.delta_a_pass(a, e, rp, Q_w_m2, H, beta, mu, f_rot)


def drag_drop_rp_m(a, e, rp, Q_w_m2, H, beta, mu, f_rot=1.0):
    """T5 periapsis change per pass [m] (negative): Lambda H, Lambda = Pi_E (1-e) / (2 e (1+e))."""
    pe = AP.pi_E(a, e, rp, Q_w_m2, H, beta, mu, f_rot); return -AP.lam_rp(e, pe) * H


# ---- T6: drag at small c, modified Bessel forms (DERIVED, checked against AP.gauss_pass_quad) -------------------------------------------------------
def _series_coeffs(p, q, N):
    """Coefficients b_k of y^k in (1+y)^p (1-y)^q, k = 0..N."""
    def binom(al, n):
        c = np.ones(n + 1)
        for k in range(1, n + 1): c[k] = c[k - 1] * (al - k + 1) / k
        return c
    A, B = binom(p, N), binom(q, N) * (-1.0) ** np.arange(N + 1)
    return np.convolve(A, B)[: N + 1]


def _bessel_integral(poly_in_x, c):
    """(1/(2 pi)) e^{-c} int_0^{2 pi} P(cos E) e^{c cos E} dE = sum_j chebcoef_j ive(j, c) (poly2cheb: P(x) = sum chebcoef_j T_j(x), T_j(cos E) = cos jE)."""
    ch = np.polynomial.chebyshev.poly2cheb(poly_in_x)
    return float(sum(cj * ive(j, c) for j, cj in enumerate(ch)))


def bessel_drag_pass(a, e, H, rho_p, beta, mu, N=8, f_rot=1.0):
    """Per-pass first-order drag changes (da, de, drp) for tangential drag on an exponential atmosphere, truncated at order N in e (error O(e^(N+1))):
        da = -2 pi f_rot (a^2 rho_p / beta) e^{-c} sum_j g_j I_j(c),   g(x) = (1+e x)^{3/2} (1-e x)^{-1/2}
        de = -2 pi f_rot (a rho_p (1-e^2) / beta) e^{-c} sum_j h_j I_j(c),   h(x) = x (1+e x)^{1/2} (1-e x)^{-1/2},  x = cos E,  c = a e / H
    where g_j, h_j are the Chebyshev coefficients of the e-series polynomials.  e^{-c} I_j(c) is evaluated as scipy ive (stable for large c).
    Reduces to the Laplace form for c >> 1 and to -2 pi a^2 rho/beta (da) for e -> 0.  Exact Kepler ellipse (E-variable), non-rotating atmosphere at f_rot = 1."""
    c = a * e / H
    pg = _series_coeffs(1.5, -0.5, N) * e ** np.arange(N + 1)                       # polynomial in x: sum_k b_k e^k x^k
    ph = np.concatenate([[0.0], (_series_coeffs(0.5, -0.5, N - 1) * e ** np.arange(N))])   # x * series
    da = -2 * math.pi * f_rot * (a * a * rho_p / beta) * _bessel_integral(pg, c)
    de = -2 * math.pi * f_rot * (a * rho_p * (1 - e * e) / beta) * _bessel_integral(ph, c)
    return da, de, (1 - e) * da - a * de


# ---- T8: third-body tide (DERIVED; closed form singly averaged, direct quadrature as the check) ---------------------------------------------------
def tide_rp_averaged(a, e, mu, mu_b, rb_vec, P_hat, Q_hat):
    """Singly averaged quadrupole: dr_p per orbit = 15 pi a e sqrt(1-e^2) Pi_3B (u.P)(u.Q), Pi_3B = (mu_b/mu)(a/r_b)^3, u = unit vector to the perturber
    (frozen over the orbit), P, Q the perihelion and semilatus-rectum unit vectors.  Derivation (DERIVED): R = (K a^2/2)[(3/2)(1+4e^2) b_P^2 + (3/2)(1-e^2) b_Q^2 - 1 - 3e^2/2],
    dR/domega = (15/2) K a^2 e^2 b_P b_Q, de/dt = -(sqrt(1-e^2)/(n a^2 e)) dR/domega, da = 0, dr_p = -a de."""
    rb = np.linalg.norm(rb_vec); u = rb_vec / rb; Pi = (mu_b / mu) * (a / rb) ** 3
    return 15 * math.pi * a * e * math.sqrt(1 - e * e) * Pi * (u @ P_hat) * (u @ Q_hat)


def tide_arc_direct(mu, mu_b, a, e, inc, raan, om, rb_of_tau, nf=64):
    """Direct first-order Gauss quadrature of the exact third-body tidal acceleration along the Kepler arc between two periapsis passages (f from 0 to 2 pi):
    returns dict (a, e, i, raan, omega, rp).  rb_of_tau(tau) -> perturber position vector (3,) or (3,K) in the same inertial frame, tau = time since periapsis."""
    x, w = np.polynomial.legendre.leggauss(nf); f = math.pi * (x + 1.0); w = math.pi * w
    p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f)); N = math.sqrt(mu / a ** 3)
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f / 2), np.sqrt(1 + e) * np.cos(f / 2)); tau = (E - e * np.sin(E)) / N
    Q = TL.pqw_matrix(inc, raan, om)
    r = Q @ np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); v = Q @ (math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f]))
    rb = rb_of_tau(tau); d = rb - r
    ai = mu_b * (d / np.linalg.norm(d, axis=0) ** 3 - rb / np.linalg.norm(rb, axis=0) ** 3)
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    aR, aS, aW = (ai * rhat).sum(0), (ai * shat).sum(0), (ai * hhat).sum(0)
    rates = TL.gauss_rates(a, e, inc, om, f, rr, h, aR, aS, aW, mu)
    return dict(zip(TL.KICK_NAMES, rates @ (w * rr ** 2 / h)))


# ---- T9: SRP bound (DERIVED) ----------------------------------------------------------------------------------------------------------------
def srp_accel(area_m2, mass_kg, cr, r_sun_m):
    """Radiation acceleration [m/s2] (ASSUMED TSI 1361 W/m2 via SRP_P1AU): cr (A/m) P_1AU (AU/r)^2."""
    return cr * area_m2 / mass_kg * SRP_P1AU * (AU / r_sun_m) ** 2


def srp_drp_bound(a, e, mu, a_srp):
    """|dr_p| per orbit <= 3 pi a_srp a^3 sqrt(1-e^2) / mu for a fixed force direction without eclipse (dr_p = 3 pi eta F_Q a^3/mu; F_Q = F cos(angle to Q) <= F).
    The first-order orbit-averaged change of a is zero for a fixed-direction force (conservative potential -F.r closes over the Kepler orbit)."""
    return 3 * math.pi * a_srp * a ** 3 * math.sqrt(1 - e * e) / mu


# ---- Pi groups (plan 3.2) ----------------------------------------------------------------------------------------------------------------------
def pi_f(planet, H=None):
    """Pi_f = f R_eq / H, flattening in scale heights (amplitude of the T2 geometry term times sin^2 i)."""
    P = PLANETS[planet]; return P["f"] * P["Re"] / (P["H"] if H is None else H)


def pi_omega(rp, mu, w_rot):
    """Pi_omega = w_rot sqrt(r_p^3/mu)."""
    return w_rot * math.sqrt(rp ** 3 / mu)


def pi_srp(a_srp, a, mu):
    """Pi_SRP = a_srp / (mu / a^2)."""
    return a_srp / (mu / a ** 2)


def pi_3b(mu_b, mu, a, r_b):
    """Pi_3B = (mu_b / mu)(a / r_b)^3 (spec section 2)."""
    return (mu_b / mu) * (a / r_b) ** 3
