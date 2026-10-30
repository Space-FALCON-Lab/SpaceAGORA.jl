"""Tesseral/sectorial gravity helpers for the per-pass derivation (tesseral_*.py). Python only, SI units.

Contents
  coefficients(planet, L)       unnormalized C, S (L+1 x L+1) from the repo gravity file (kaula_lib reader); C[0,0] = 1
  vw(xb, yb, zb, L, R)          Cunningham V_nm, W_nm (Montenbruck & Gill 2000, sec. 3.2; recalled, verified numerically in tesseral_selftest.py)
  accel_body(...)               body-fixed acceleration of the degree >= 2 part (vectorized over K positions, per-position coefficient sets allowed)
  potential_body(...)           matching potential (for the gradient check and the Kaula reconstruction)
  kepler_arc(...)               Kepler arc from apoapsis to apoapsis (true-anomaly Gauss-Legendre nodes), times, inertial states
  gauss_kick(...)               first-order per-pass change of (a, e, i, Omega, omega, r_p) by Gauss quadrature of a body-fixed field along the arc,
                                with the planet rotating (theta = theta_p + w_rot (t - t_p)); the "direct" per-pass kick
  incl_fn(n, m, inc)            complex inclination function D_nmk(i), k = -n..n:  P_nm(sin phi) e^{i m (alpha - Omega)} = sum_k D_nmk e^{i k u}
  hansen(n, k, e, NM)           X_j^{-(n+1),k}(e) for j = -NM/2 .. NM/2-1 by FFT:  (a/r)^{n+1} e^{i k f} = sum_j X_j e^{i j M}
  kaula_F(n, m, p, inc)         Kaula's closed-form inclination function F_nmp(i) (recalled, Kaula 1966; checked against incl_fn)
  series_kick(...)              per-pass kick from the Lagrange planetary equations applied to the (n, m) series, summed over k and |j| <= J,
                                with the commensurability filter W_j = 2 sin(pi s)/s, s = j - m nu
"""
import math
import numpy as np
import kaula_lib as kl

# ---------------------------------------------------------------- coefficients ----------------------------------------------------------


def norm_factor(n, m):
    """Fully normalized -> unnormalized factor N_nm = sqrt((2 - delta_m0)(2n+1)(n-m)!/(n+m)!)."""
    return math.sqrt((1 if m == 0 else 2) * (2 * n + 1) * math.factorial(n - m) / math.factorial(n + m))


def coefficients(planet, L, gdir=None):
    df, R, mu, meta = kl.read_gravity(planet, gdir)
    Cb = np.zeros((L + 1, L + 1)); Sb = np.zeros((L + 1, L + 1))
    for r in df.itertuples():
        n, m = int(r.degree), int(r.order)
        if n <= L and m <= n: Cb[n, m], Sb[n, m] = r.C, r.S
    C = np.zeros_like(Cb); S = np.zeros_like(Sb)
    for n in range(L + 1):
        for m in range(n + 1):
            f = norm_factor(n, m); C[n, m] = Cb[n, m] * f; S[n, m] = Sb[n, m] * f
    C[0, 0] = 1.0; C[1, :] = 0.0; S[1, :] = 0.0
    return dict(C=C, S=S, Cbar=Cb, Sbar=Sb, R=R, mu=mu, L=L, meta=meta)


def select(C, S, keep):
    """Copy of (C, S) keeping only entries where keep(n, m) is True (degree >= 2)."""
    C2 = np.zeros_like(C); S2 = np.zeros_like(S); L = C.shape[0] - 1
    for n in range(2, L + 1):
        for m in range(n + 1):
            if keep(n, m): C2[n, m], S2[n, m] = C[n, m], S[n, m]
    return C2, S2

# ---------------------------------------------------------------- Cunningham V, W ---------------------------------------------------------


def vw(xb, yb, zb, L, R):
    """V, W of shape (L+3, L+3, K):  V_nm + i W_nm = (R/r)^{n+1} P_nm(sin phi) e^{i m lambda} (P_nm without Condon-Shortley phase)."""
    xb, yb, zb = np.atleast_1d(xb), np.atleast_1d(yb), np.atleast_1d(zb)
    K = xb.shape[0]; N = L + 3
    r2 = xb * xb + yb * yb + zb * zb; rho = R * R / r2; x0, y0, z0 = R * xb / r2, R * yb / r2, R * zb / r2
    V = np.zeros((N, N, K)); W = np.zeros((N, N, K))
    V[0, 0] = R / np.sqrt(r2)
    ms = np.arange(N)
    for n in range(1, N):
        # sectoral
        V[n, n] = (2 * n - 1) * (x0 * V[n - 1, n - 1] - y0 * W[n - 1, n - 1])
        W[n, n] = (2 * n - 1) * (x0 * W[n - 1, n - 1] + y0 * V[n - 1, n - 1])
        # n, m < n: V_nm = ((2n-1) z0 V_{n-1,m} - (n+m-1) rho V_{n-2,m}) / (n-m)
        m = ms[:n][:, None]
        a1 = ((2 * n - 1) / (n - m)) * z0
        V[n, :n] = a1 * V[n - 1, :n]; W[n, :n] = a1 * W[n - 1, :n]
        if n >= 2:
            a2 = ((n + m - 1) / (n - m)) * rho
            V[n, :n] -= a2 * V[n - 2, :n]; W[n, :n] -= a2 * W[n - 2, :n]
    return V, W


def accel_body(xb, yb, zb, cs, C=None, S=None):
    """Acceleration (body frame) of degrees >= 2 of the field (C, S default: cs's full set). Returns (3, K)."""
    C = cs["C"] if C is None else C; S = cs["S"] if S is None else S
    L = C.shape[0] - 1; R, mu = cs["R"], cs["mu"]
    V, W = vw(xb, yb, zb, L, R)
    Cp = C.copy(); Sp = S.copy(); Cp[:2] = 0.0; Sp[:2] = 0.0
    if Cp.ndim == 2: Cp, Sp = Cp[..., None], Sp[..., None]
    nnn = np.arange(L + 1)[:, None, None]; mm = np.arange(L + 1)[None, :, None]
    half = np.where(mm == 0, 1.0, 0.5); fac = (nnn - mm + 2) * (nnn - mm + 1) * (mm > 0)
    nn = np.arange(L + 1)
    Vn1 = V[nn + 1]; Wn1 = W[nn + 1]                    # (L+1, L+3, K): index n+1
    Vp = Vn1[:, 1:L + 2]; Wp = Wn1[:, 1:L + 2]          # m+1
    Vm = np.concatenate([np.zeros_like(Vn1[:, :1]), Vn1[:, :L]], axis=1); Wm = np.concatenate([np.zeros_like(Wn1[:, :1]), Wn1[:, :L]], axis=1)  # m-1
    V0 = Vn1[:, :L + 1]; W0 = Wn1[:, :L + 1]
    ax = (half * (-Cp * Vp - Sp * Wp) + 0.5 * fac * (Cp * Vm + Sp * Wm)).sum(axis=(0, 1))
    ay = (half * (-Cp * Wp + Sp * Vp) + 0.5 * fac * (-Cp * Wm + Sp * Vm)).sum(axis=(0, 1))
    az = ((nnn - mm + 1) * (-Cp * V0 - Sp * W0)).sum(axis=(0, 1))
    return (mu / R ** 2) * np.array([ax, ay, az])


def potential_body(xb, yb, zb, cs, C=None, S=None):
    """Potential of degrees >= 2 (U = (mu/R) sum C V + S W), shape (K,)."""
    C = cs["C"] if C is None else C; S = cs["S"] if S is None else S
    L = C.shape[0] - 1; V, W = vw(xb, yb, zb, L, cs["R"])
    Cp = C.copy(); Sp = S.copy(); Cp[:2] = 0.0; Sp[:2] = 0.0
    if Cp.ndim == 2: Cp, Sp = Cp[..., None], Sp[..., None]
    return (cs["mu"] / cs["R"]) * (Cp * V[:L + 1, :L + 1] + Sp * W[:L + 1, :L + 1]).sum(axis=(0, 1))

# ---------------------------------------------------------------- frames, Kepler -----------------------------------------------------------


def rot_z(th, v):
    """Rotate vectors v (3, K) by angle th (K,) about z (inertial = rot_z(theta) body)."""
    c, s = np.cos(th), np.sin(th)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]])


def pqw_matrix(inc, raan, om):
    c, s = np.cos, np.sin
    R3 = lambda t: np.array([[c(t), -s(t), 0], [s(t), c(t), 0], [0, 0, 1]])
    R1 = lambda t: np.array([[1, 0, 0], [0, c(t), -s(t)], [0, s(t), c(t)]])
    return R3(raan) @ R1(inc) @ R3(om)


def elements(r, v, mu):
    """Osculating (a, e, i, raan, omega, f) from inertial r, v (3,) arrays."""
    h = np.cross(r, v); rn = np.linalg.norm(r); vn2 = v @ v
    a = 1.0 / (2.0 / rn - vn2 / mu)
    ev = np.cross(v, h) / mu - r / rn; e = np.linalg.norm(ev)
    hn = np.linalg.norm(h); inc = math.acos(h[2] / hn)
    nvec = np.array([-h[1], h[0], 0.0]); nn = np.linalg.norm(nvec)
    raan = math.atan2(nvec[1], nvec[0]) % (2 * math.pi)
    om = math.atan2((np.cross(nvec, ev) @ h) / hn, nvec @ ev) % (2 * math.pi)
    f = math.atan2((np.cross(ev, r) @ h) / hn, ev @ r) % (2 * math.pi)
    return a, e, inc, raan, om, f


def kepler_arc(a, e, inc, raan, om, mu, nf=256):
    """Gauss-Legendre nodes in true anomaly over [-pi, pi]. Returns f, w (weights in f), t - t_p, r (3,nf), v (3,nf), r scalar, h."""
    x, w = np.polynomial.legendre.leggauss(nf); f = math.pi * x; w = math.pi * w
    p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f))
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f / 2)); M = E - e * np.sin(E); N = math.sqrt(mu / a ** 3)
    Q = pqw_matrix(inc, raan, om)
    rp = np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); vp = math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f])
    return f, w, M / N, Q @ rp, Q @ vp, rr, h


def gauss_rates(a, e, inc, om, f, rr, h, aR, aS, aW, mu):
    """Gauss variational equations (RSW accelerations), rates per unit time of a, e, i, Omega, omega, r_p."""
    p = a * (1 - e * e); u = om + f
    da = 2 * a * a / h * (e * np.sin(f) * aR + p / rr * aS)
    de = (p * np.sin(f) * aR + ((p + rr) * np.cos(f) + rr * e) * aS) / h
    di = rr * np.cos(u) * aW / h
    dO = rr * np.sin(u) * aW / (h * np.sin(inc))
    dw = (-p * np.cos(f) * aR + (p + rr) * np.sin(f) * aS) / (h * e) - rr * np.sin(u) * np.cos(inc) * aW / (h * np.sin(inc))
    return np.array([da, de, di, dO, dw, (1 - e) * da - a * de])


KICK_NAMES = ("a", "e", "i", "raan", "omega", "rp")


def gauss_kick(a, e, inc, raan, om, theta_p, w_rot, cs, C=None, S=None, nf=256, fmask=None):
    """First-order change of (a, e, i, raan, omega, rp) over the apoapsis-to-apoapsis arc of the osculating ellipse, with the planet's
    rotation angle theta_p at periapsis. C, S: the field to apply (degree >= 2). fmask(f) -> bool restricts the arc (localization test)."""
    mu = cs["mu"]; f, w, tau, r, v, rr, h = kepler_arc(a, e, inc, raan, om, mu, nf)
    th = theta_p + w_rot * tau
    rb = rot_z(-th, r); ab = accel_body(rb[0], rb[1], rb[2], cs, C, S); ai = rot_z(th, ab)
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    aR, aS, aW = (ai * rhat).sum(0), (ai * shat).sum(0), (ai * hhat).sum(0)
    rates = gauss_rates(a, e, inc, om, f, rr, h, aR, aS, aW, mu)
    wt = w * rr ** 2 / h                                  # dt = r^2/h df
    if fmask is not None: wt = wt * fmask(f)
    return dict(zip(KICK_NAMES, rates @ wt))

def term_kicks(a, e, inc, raan, om, theta_p, w_rot, mu, R, L, nf=128, half=False, mmin=1, fmask=None):
    """Per-term first-order periapsis-radius kick over one apoapsis-to-apoapsis arc for a UNIT fully normalized coefficient:
    Y[n, m] complex with  kick_rp(n, m) = Re[(Cbar_nm - i Sbar_nm) * Y[n, m]] = Cbar Re Y + Sbar Im Y.  half=True: inbound half only
    (apoapsis to periapsis; the physical periapsis radius at this pass = osculating r_p at the previous apoapsis + this half kick).
    One spherical-harmonic recursion per node for all (n, m); m < mmin are zeroed."""
    f, w, tau, r, v, rr, h = kepler_arc(a, e, inc, raan, om, mu, nf)
    if half:
        x, w = np.polynomial.legendre.leggauss(nf); f = 0.5 * math.pi * (x - 1.0); w = 0.5 * math.pi * w
        p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f))
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f / 2)); tau = (E - e * np.sin(E)) / math.sqrt(mu / a ** 3)
        Q = pqw_matrix(inc, raan, om)
        r = Q @ np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); v = Q @ (math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f]))
    th = theta_p + w_rot * tau; rb = rot_z(-th, r)
    V, W = vw(rb[0], rb[1], rb[2], L, R)
    nn = np.arange(L + 1)[:, None, None]; mm = np.arange(L + 1)[None, :, None]
    halfw = np.where(mm == 0, 1.0, 0.5); fac = (nn - mm + 2) * (nn - mm + 1) * (mm > 0)
    idx = np.arange(L + 1); Vn1 = V[idx + 1]; Wn1 = W[idx + 1]
    Vp, Wp = Vn1[:, 1:L + 2], Wn1[:, 1:L + 2]
    Vm = np.concatenate([np.zeros_like(Vn1[:, :1]), Vn1[:, :L]], axis=1); Wm = np.concatenate([np.zeros_like(Wn1[:, :1]), Wn1[:, :L]], axis=1)
    V0, W0 = Vn1[:, :L + 1], Wn1[:, :L + 1]
    Nf = np.array([[norm_factor(n, m) if m <= n else 0.0 for m in range(L + 1)] for n in range(L + 1)])[..., None]
    Nf[:2] = 0.0; Nf[:, :mmin] = 0.0
    g = mu / R ** 2 * Nf
    aC = g * np.array([-halfw * Vp + 0.5 * fac * Vm, -halfw * Wp - 0.5 * fac * Wm, -(nn - mm + 1) * V0])      # (3, n, m, nodes) body frame
    aS = g * np.array([-halfw * Wp + 0.5 * fac * Wm, halfw * Vp + 0.5 * fac * Vm, -(nn - mm + 1) * W0])
    c, s = np.cos(th), np.sin(th)
    def to_rsw(ab):
        ax = c * ab[0] - s * ab[1]; ay = s * ab[0] + c * ab[1]; az = ab[2]
        rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
        return [ax * rhat[0] + ay * rhat[1] + az * rhat[2], ax * shat[0] + ay * shat[1] + az * shat[2], ax * hhat[0] + ay * hhat[1] + az * hhat[2]]
    p = a * (1 - e * e); wt = w * rr ** 2 / h * (1.0 if fmask is None else fmask(f))
    def drp(ab):
        aR, aS_, aW = to_rsw(ab)
        da = 2 * a * a / h * (e * np.sin(f) * aR + p / rr * aS_); de = (p * np.sin(f) * aR + ((p + rr) * np.cos(f) + rr * e) * aS_) / h
        return (((1 - e) * da - a * de) * wt).sum(-1)
    return drp(aC) + 1j * drp(aS)


# ---------------------------------------------------------------- Kaula functions ----------------------------------------------------------


def incl_fn(n, m, inc, nu_pts=None):
    """D_nmk(i) for k = -n..n (array index k+n):  P_nm(sin phi) e^{i m (alpha - Omega)} = sum_k D_nmk(i) e^{i k u}  (trig polynomial of degree n)."""
    nu_pts = nu_pts or 4 * (n + 1)
    u = 2 * np.pi * np.arange(nu_pts) / nu_pts
    x, y, z = np.cos(u), np.sin(u) * np.cos(inc), np.sin(u) * np.sin(inc)
    V, W = vw(x, y, z, n, 1.0)                            # r = R = 1: V + iW = P_nm e^{i m alpha}
    Y = V[n, m] + 1j * W[n, m]
    c = np.fft.fft(Y) / nu_pts
    ks = np.arange(-n, n + 1)
    return c[ks % nu_pts]


def hansen(n, k, e, NM=8192):
    """X_j^{-(n+1),k}(e), j = -NM/2 .. NM/2 - 1 (array in np.fft order: index j % NM)."""
    M = 2 * np.pi * np.arange(NM) / NM
    E = M.copy()
    for _ in range(60): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    f = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    a_r = 1.0 / (1 - e * np.cos(E))
    return np.fft.fft(a_r ** (n + 1) * np.exp(1j * k * f)) / NM


def kaula_F(n, m, p, inc):
    """Kaula's inclination function (as recalled from Kaula 1966, eq. 3.61; checked numerically against incl_fn in tesseral_selftest.py)."""
    kk = (n - m) // 2; s_i, c_i = math.sin(inc), math.cos(inc); tot = 0.0
    for t in range(0, min(p, kk) + 1):
        pre = math.factorial(2 * n - 2 * t) / (math.factorial(t) * math.factorial(n - t) * math.factorial(n - m - 2 * t) * 2 ** (2 * n - 2 * t))
        sm = 0.0
        for s in range(0, m + 1):
            inner = 0.0
            for c in range(0, p - t + 1):
                if n - m - 2 * t + s < c or m - s < p - t - c or p - t - c < 0: continue
                inner += math.comb(n - m - 2 * t + s, c) * math.comb(m - s, p - t - c) * (-1) ** (c - kk)
            sm += math.comb(m, s) * c_i ** s * inner
        tot += pre * s_i ** (n - m - 2 * t) * sm
    return tot


class SeriesTerm:
    """Cached Kaula-series pieces of one (n, m) term at (e, i): D_k, dD_k/di, X_j^k, dX_j^k/de (finite differences)."""
    def __init__(self, n, m, e, inc, NM=8192, de=1e-6, di=1e-6):
        self.n, self.m, self.e, self.inc, self.NM = n, m, e, inc, NM
        self.ks = np.arange(-n, n + 1, 1)
        self.D = incl_fn(n, m, inc); self.dD = (incl_fn(n, m, inc + di) - incl_fn(n, m, inc - di)) / (2 * di)
        self.X = np.array([hansen(n, k, e, NM) for k in self.ks])
        self.dX = np.array([(hansen(n, k, e + de, NM) - hansen(n, k, e - de, NM)) / (2 * de) for k in self.ks])
        self.j = np.fft.fftfreq(NM, 1.0 / NM)


def series_kick(term, a, raan, om, theta_p, nu, cs, J=None, only_j=None):
    """Per-pass change of (a, e, i, raan, omega, rp) from the (n, m) term by the Lagrange planetary equations, integrated over
    M in [-pi, pi] (apoapsis to apoapsis) with theta = theta_p + nu M, elements frozen. J: keep |j| <= J; only_j: keep listed j only."""
    n, m, e, inc = term.n, term.m, term.e, term.inc; mu, R = cs["mu"], cs["R"]
    N = math.sqrt(mu / a ** 3); eta = math.sqrt(1 - e * e)
    amp = (mu / a) * (R / a) ** n * (cs["C"][n, m] - 1j * cs["S"][n, m])
    cols = np.arange(term.NM) if only_j is None else (np.asarray(only_j, int) % term.NM)
    j = term.j[cols]; s = j - m * nu
    Wj = np.where(np.abs(s) < 1e-12, 2 * np.pi, 2 * np.sin(np.pi * s) / np.where(np.abs(s) < 1e-12, 1.0, s))
    if J is not None: Wj = Wj * (np.abs(j) <= J)
    Wj = Wj / N                                          # time integral factor: int e^{i s M} dt over the orbit
    ph = np.exp(1j * (term.ks * om + m * (raan - theta_p)))   # e^{i psi_p} without the j part (M = 0 at periapsis)
    # sums S_x = sum_k sum_j coef_x(k, j) A_kj e^{i psi_p} W_j
    X, dX = term.X[:, cols], term.dX[:, cols]
    A = amp * (term.D * ph)[:, None] * X                 # (k, j)
    Ai = amp * (term.dD * ph)[:, None] * X
    Ae = amp * (term.D * ph)[:, None] * dX
    SM = (1j * j[None, :] * A * Wj).sum(); Sw = (1j * term.ks[:, None] * A * Wj).sum(); SO = (1j * m * A * Wj).sum()
    Si = (Ai * Wj).sum(); Se = (Ae * Wj).sum()
    da = 2 / (N * a) * SM
    de = (eta ** 2 / (N * a * a * e)) * SM - eta / (N * a * a * e) * Sw
    di = (math.cos(inc) * Sw - SO) / (N * a * a * eta * math.sin(inc))
    dO = Si / (N * a * a * eta * math.sin(inc))
    dw = eta / (N * a * a * e) * Se - math.cos(inc) / (N * a * a * eta * math.sin(inc)) * Si
    drp = (1 - e) * da - a * de
    return dict(zip(KICK_NAMES, [float(np.real(x)) for x in (da, de, di, dO, dw, drp)]))
