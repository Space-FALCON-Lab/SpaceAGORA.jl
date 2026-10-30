"""Conservative zonal Cowell propagator (point mass + zonals J2..Jnmax about the body z axis) for checking the analytic zonal geometry.

No drag, no third bodies, no tesserals. Returns the osculating periapsis radius at every periapsis passage (radial-velocity zero crossing, dense output).
Units inside: length R, time sqrt(R^3/mu).
"""
import numpy as np
from scipy.integrate import solve_ivp


def _leg(nmax, s):
    P = [1.0, s]
    for k in range(2, nmax + 1): P.append(((2 * k - 1) * s * P[-1] - (k - 1) * P[-2]) / k)
    dP = [0.0] + [n * (s * P[n] - P[n - 1]) / (s * s - 1.0) if abs(s) < 1 - 1e-14 else 0.5 * n * (n + 1) * np.sign(s) ** (n + 1) for n in range(1, nmax + 1)]
    return P, dP


def make_rhs(Js, mu=1.0, R=1.0):
    """Js: dict n -> J_n. Returns rhs(t, y) for y = (x,y,z,vx,vy,vz)."""
    nmax = max(Js) if Js else 2; items = [(n, J) for n, J in sorted(Js.items())]
    def rhs(t, y):
        x, yy, z = y[0], y[1], y[2]; r2 = x * x + yy * yy + z * z; r = np.sqrt(r2); s = z / r
        ax = -mu * x / (r2 * r); ay = -mu * yy / (r2 * r); az = -mu * z / (r2 * r)
        P, dP = _leg(nmax, s)
        for n, J in items:
            f = -mu * J * R ** n / r ** (n + 2)
            rad = -(n + 1) * P[n]; zc = dP[n]; rr = -dP[n] * s       # grad = r^-(n+2)[ rad*rhat + dP*(zhat - s*rhat) ]
            ax += f * ((rad + rr) * x / r); ay += f * ((rad + rr) * yy / r); az += f * ((rad + rr) * z / r + zc)
        return np.array([y[3], y[4], y[5], ax, ay, az])
    return rhs


def apoapsis_state(a, e, inc, raan, omega, mu=1.0):
    """Cartesian state at true anomaly 180 deg."""
    nu = np.pi; p = a * (1 - e * e); r = p / (1 + e * np.cos(nu))
    rp_ = np.array([r * np.cos(nu), r * np.sin(nu), 0.0]); vp_ = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    def R3(t): c, s = np.cos(t), np.sin(t); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    def R1(t): c, s = np.cos(t), np.sin(t); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    Q = R3(raan) @ R1(inc) @ R3(omega)
    return np.concatenate([Q @ rp_, Q @ vp_])


def periapsis_series(Js, a, e, inc, raan, omega, t_end_s, mu, R, rtol=1e-11, atol_len=1e-9):
    """Return (t_s array, rp_m array) at periapsis passages for the zonal set Js. a in m, angles in rad."""
    T0 = np.sqrt(R ** 3 / mu); y0 = apoapsis_state(a / R, e, inc, raan, omega)
    rhs = make_rhs(Js)
    ev = lambda t, y: y[0] * y[3] + y[1] * y[4] + y[2] * y[5]; ev.direction = 1.0
    sol = solve_ivp(rhs, (0.0, t_end_s / T0), y0, method="DOP853", rtol=rtol, atol=atol_len, events=ev)
    ts, rps = [], []
    for t, y in zip(sol.t_events[0], sol.y_events[0]):
        ts.append(t * T0); rps.append(np.linalg.norm(y[:3]) * R)
    return np.array(ts), np.array(rps)
