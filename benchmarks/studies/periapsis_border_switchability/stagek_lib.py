"""Stage K library: natural periapsis drift of a burn-free orbit about an oblate planet (Python only, no Julia).

Facts used (DERIVED in tasks/.../stagek_plan.md section 2, checked by stagek_selftest.py):
  * r_p of the osculating orbit at the pass (f = 0, rdot = 0) is a (1 - e) exactly (the pass tables carry radius_m = a (1 - e) to 1e-8 m).
  * Secular J2 gives da = de = di = 0; the J2 contribution to the pass-to-pass change of r_p is second order in J2 (the first-order
    short-period theory is ill-conditioned for r_p at e ~ 0.8: terms of 1e5 m cancel to 1e3 m and the neglected second-order terms are of
    the size of the answer), so the zonal prediction is a numerical one-orbit map: integrate the Keplerian + zonal-harmonic (J2..Jn)
    equations of motion from the periapsis state at pass k to the next periapsis and read a (1 - e).
  * First-order orbit-averaged zonal rates (kaula_lib.zonal_rates) give the apsidal rate and, for degrees >= 3 (and 4 ..), the secular de/dt;
    they are used for the apsidal rotation per pass and as an analytic cross-check of the numerical map.
"""
import math, pathlib, sys
import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kaula_lib as KL


def planet(name, nmax=20):
    df, R, mu, meta = KL.read_gravity(name)
    J = KL.zonal_J(df)
    return dict(name=name, R=R, mu=mu, J={n: J[n] for n in range(2, nmax + 1) if n in J})


def elements(mu, r, v):
    """(a, e, inc, Omega, omega) from a Cartesian state (equatorial frame of the planet pole)."""
    rn = np.linalg.norm(r); h = np.cross(r, v); hn = np.linalg.norm(h); a = -mu / (2 * (v @ v / 2 - mu / rn))
    ev = np.cross(v, h) / mu - r / rn; e = np.linalg.norm(ev); inc = math.acos(h[2] / hn)
    node = np.array([-h[1], h[0], 0.0]); node /= np.linalg.norm(node)
    Om = math.atan2(node[1], node[0]); w = math.atan2(np.cross(node, ev) @ (h / hn), node @ ev)
    return a, e, inc, Om % (2 * math.pi), w % (2 * math.pi)


def periapsis_state(mu, a, e, inc, Om, w):
    """Cartesian state at periapsis from osculating elements."""
    rp = a * (1 - e); vp = math.sqrt(mu * (1 + e) / rp)
    cO, sO, cw, sw, ci, si = math.cos(Om), math.sin(Om), math.cos(w), math.sin(w), math.cos(inc), math.sin(inc)
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si]); Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    return rp * P, vp * Q


def make_rhs(P, nmax):
    mu, R = P["mu"], P["R"]; Js = [(n, P["J"][n]) for n in range(2, nmax + 1) if n in P["J"]]

    def rhs(t, u):
        x, y, z, vx, vy, vz = u; r2 = x * x + y * y + z * z; r = math.sqrt(r2); s = z / r
        # Legendre P_n(s) and P_n'(s) by recurrence
        p0, p1 = 1.0, s; Pn = {1: s};
        for k in range(2, nmax + 1):
            p0, p1 = p1, ((2 * k - 1) * s * p1 - (k - 1) * p0) / k; Pn[k] = p1
        # Phi = -mu/r + sum mu J_n R^n P_n / r^(n+1);  dP_n/ds = n (s P_n - P_{n-1})/(s^2 - 1)  (n >= 1)
        dFr = 0.0; dFs = 0.0   # dFr = sum (n+1) c_n P_n,  dFs = sum c_n dP_n/ds   with c_n = mu J_n (R/r)^n / r^2 ... scaled below
        for n, Jn in Js:
            c = mu * Jn * (R / r) ** n / r2
            Pm = Pn[n - 1] if n > 1 else 1.0
            dP = n * (s * Pn[n] - Pm) / (s * s - 1.0) if abs(s) < 1 - 1e-12 else 0.5 * n * (n + 1) * (1.0 if s > 0 else (-1.0) ** (n + 1))
            dFr += (n + 1) * c * Pn[n]; dFs += c * dP
        # force = -grad Phi; Phi_r = mu/r^2 - sum (n+1) c_n P_n / r ... (c already divided by r^2, so Phi_r = mu/r^2 - dFr / r * r ... see below)
        # Phi = -mu/r + sum mu J_n R^n P_n r^-(n+1):  dPhi/dr = mu/r^2 - sum (n+1) mu J_n R^n P_n r^-(n+2) = mu/r^2 - dFr
        # grad of P_n(z/r) term: dPhi/dx_i = (dPhi/dr) x_i/r + sum mu J_n R^n r^-(n+1) P_n'(s) (delta_iz - s x_i/r)/r
        Phir = mu / r2 - dFr
        gx = Phir * x / r + dFs * (-s * x / r) / 1.0 * 1.0
        gy = Phir * y / r + dFs * (-s * y / r)
        gz = Phir * z / r + dFs * (1.0 - s * s)
        return [vx, vy, vz, -gx, -gy, -gz]
    return rhs


def zonal_map(P, a, e, inc, Om, w, nmax=20, rtol=1e-11, atol=1e-8):
    """One-orbit map: from the periapsis state with osculating elements (a, e, inc, Om, w) integrate the zonal (J2..J_nmax) problem to the next
    periapsis; return the osculating elements there (a, e, inc, Om, w) and the interval (s)."""
    mu = P["mu"]; r0, v0 = periapsis_state(mu, a, e, inc, Om, w); T = 2 * math.pi * math.sqrt(a ** 3 / mu)
    rhs = make_rhs(P, nmax)
    ev = lambda t, u: u[0] * u[3] + u[1] * u[4] + u[2] * u[5]; ev.direction = 1
    # start slightly after periapsis so that the first event is the NEXT one
    t0 = 1e-3 * T
    s0 = solve_ivp(rhs, [0, t0], np.concatenate([r0, v0]), method="DOP853", rtol=rtol, atol=atol)
    sol = solve_ivp(rhs, [t0, 1.5 * T], s0.y[:, -1], method="DOP853", rtol=rtol, atol=atol, events=ev)
    te = sol.t_events[0]
    if len(te) == 0: raise RuntimeError("no periapsis found")
    u = sol.y_events[0][0]
    return (*elements(mu, np.array(u[:3]), np.array(u[3:])), te[0])


def zonal_secular(P, a, e, inc, w, nlo=2, nhi=20, nq=1024):
    """First-order orbit-averaged (de/dt, dw/dt, dOmega/dt) summed over degrees nlo..nhi (per second), mean motion sqrt(mu/a^3)."""
    N = math.sqrt(P["mu"] / a ** 3); p = a * (1 - e * e); ed = wd = Od = 0.0
    for n, Jn in P["J"].items():
        if n < nlo or n > nhi: continue
        r = KL.zonal_rates(n, e, inc, w, Jn * (P["R"] / p) ** n, N, nq=nq)
        ed += r[0]; wd += r[1]; Od += r[2]
    return ed, wd, Od


def lever_m_per_mps(mu, a, e):
    """d r_p / d(dv) for a tangential apoapsis burn: 4 a sqrt((1-e)/(1+e)) / sqrt(mu/a)  (m per m/s)."""
    return 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(mu / a)


def pi_drift(a, e, dr):
    """Dimensionless drift cost Pi_drift = (dr/a) / (4 sqrt((1-e)/(1+e))) = dv_cancel / sqrt(mu/a) (signed: > 0 means r_p rising)."""
    return (dr / a) / (4 * np.sqrt((1 - e) / (1 + e)))


def burn_apsidal_gain(e, f):
    """|d omega| per Pi_dv (dv / sqrt(mu/a)) of the best-direction in-plane impulse at true anomaly f (Gauss):
    d omega = (sqrt(p/mu)/e) sqrt(cos^2 f + (1 + r/p)^2 sin^2 f) dv,  r/p = 1/(1 + e cos f)."""
    return math.sqrt(1 - e * e) / e * np.sqrt(np.cos(f) ** 2 + (1 + 1 / (1 + e * np.cos(f))) ** 2 * np.sin(f) ** 2)
