"""Kaula-rule helpers shared by the kaula_*.py scripts (gravity-file reader, degree amplitudes, K fit, zonal secular geometry).

Conventions. Files hold fully normalized coefficients (columns degree, order, C, S). The zonal J_n = -sqrt(2n+1) Cbar_n0.
Degree amplitude (rms per coefficient) A_n = sqrt( sum_m (C_nm^2 + S_nm^2) / (2n+1) ); Kaula's rule: A_n ~ K / n^2.
"""
import os, pathlib
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
FILES = {"Earth": "EarthGGM05C.csv", "Mars": "Mars50c.csv", "Venus": "MGNP180U.csv"}


def gravity_dir(cli=None):
    if cli: return pathlib.Path(cli)
    env = os.environ.get("SPACEAGORA_GRAVITY_DIR")
    if env: return pathlib.Path(env)
    # worktrees do not carry data/: fall back to the main checkout through the git common dir
    for q in [REPO, *REPO.parents]:
        if (q / "data" / "Gravity_harmonics_data").is_dir(): return q / "data" / "Gravity_harmonics_data"
    raise FileNotFoundError("data/Gravity_harmonics_data not found; set SPACEAGORA_GRAVITY_DIR")


def read_gravity(name, gdir=None):
    path = gravity_dir(gdir) / FILES.get(name, name)
    meta = {}
    for line in open(path):
        if not line.startswith("#"): break
        if ":" in line:
            k, v = line[1:].split(":", 1); meta[k.strip()] = v.strip()
    df = pd.read_csv(path, comment="#")
    return df, float(meta["reference_radius_m"]), float(meta["gm_m3s2"]), meta


def degree_amplitude(df, nmax=None):
    nmax = int(df.degree.max()) if nmax is None else nmax
    n = np.arange(2, nmax + 1); A = np.full(n.shape, np.nan)
    for j, d in enumerate(n):
        g = df[df.degree == d]
        if len(g): A[j] = np.sqrt((g.C ** 2 + g.S ** 2).sum() / (2 * d + 1))
    return n, A


def zonal_J(df):
    """dict n -> J_n (unnormalized) for the zonal rows of the file."""
    z = df[df.order == 0]
    return {int(r.degree): -np.sqrt(2 * r.degree + 1) * r.C for r in z.itertuples() if r.degree >= 2}


def fit_K(n, A, nlo, nhi):
    """Geometric-mean fit of K_n = A_n n^2 over nlo <= n <= nhi. Returns K, scatter in dex (std of log10 K_n), n_used, log-log free slope."""
    m = (n >= nlo) & (n <= nhi) & np.isfinite(A) & (A > 0)
    lk = np.log10(A[m] * n[m] ** 2)
    slope = np.polyfit(np.log10(n[m]), np.log10(A[m]), 1)[0]
    return 10 ** lk.mean(), lk.std(ddof=1), int(m.sum()), slope


# ---- orbit-averaged zonal geometry (single average over the mean anomaly) ----------------------------------
# <R_n> = -(mu/a) J_n sqrt(1-e^2) Phi_n(e,i,w),  Phi_n = (R/p)^n (1/2pi) Int (1+e cos nu)^(n-1) P_n(sin i sin(w+nu)) d nu,
# R_disturbing = U - mu/r = -(mu/r) sum J_n (R/r)^n P_n(sin phi).  Rates with Pi_n = J_n (R/p)^n and N = sqrt(mu/a^3):
#   de/dt = N Pi_n Gamma_n,  Gamma_n = (1-e^2)/e * (1/2pi) Int (1+e cos nu)^(n-1) P_n'(x) sin i cos(w+nu) d nu,  x = sin i sin(w+nu)
from numpy.polynomial import legendre as _L

_NQ = 4096


def _legendre_and_derivative(n, x):
    """P_n(x), P_n'(x) by the three-term recurrence (stable for n up to several hundred)."""
    p0, p1 = np.ones_like(x), x.copy()
    if n == 0: return p0, np.zeros_like(x)
    for k in range(2, n + 1): p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
    dp = n * (x * p1 - p0) / (x * x - 1.0) if n >= 1 else 0 * x
    # fix |x| = 1 (not hit for sin i < 1 off-pole nodes except i = 90 deg); use the closed form n(n+1)/2 * x^(n+1)
    bad = np.abs(np.abs(x) - 1.0) < 1e-12
    if bad.any(): dp = np.where(bad, 0.5 * n * (n + 1) * np.sign(x) ** (n + 1), dp)
    return p1, dp


def zonal_averages(n, e, inc, w, nq=_NQ):
    """Return dict with Phi, Phi_e, Phi_i, Phi_w (scaled by (R/p)^n, see header). inc, w in radians."""
    nu = (np.arange(nq) + 0.5) * (2 * np.pi / nq); c = np.cos(nu); u = w + nu
    S = np.exp(n * np.log(1 + 0 * c) + (n - 1) * np.log1p(e * c))   # (1+e cos)^(n-1), the (R/p)^n factor is carried outside
    x = np.sin(inc) * np.sin(u); P, dP = _legendre_and_derivative(n, x)
    m = lambda f: f.mean()
    return dict(Phi=m(S * P), Phi_e=m((n - 1) * np.exp((n - 2) * np.log1p(e * c)) * c * P) if n >= 2 else 0.0,
                Phi_i=m(S * dP * np.cos(inc) * np.sin(u)), Phi_w=m(S * dP * np.sin(inc) * np.cos(u)))


def zonal_rates(n, e, inc, w, Pi_n, N, nq=_NQ):
    """Secular rates (e, omega, Omega) in rad/s or 1/s caused by degree n, given Pi_n = J_n (R/p)^n and mean motion N."""
    A = zonal_averages(n, e, inc, w, nq); q = 1 - e * e
    edot = N * Pi_n * q / e * A["Phi_w"]
    wdot = -N * Pi_n * (q / e * ((2 * n - 1) * e / q * A["Phi"] + A["Phi_e"]) - np.cos(inc) / np.sin(inc) * A["Phi_i"])
    Odot = -N * Pi_n * A["Phi_i"] / np.sin(inc)
    return edot, wdot, Odot


def gamma_max(n, e, inc, nw=24, nq=_NQ):
    """max over omega of |Gamma_n| (so |de/dt|max = N |Pi_n| Gamma_max); grid over omega in [0, pi) then refine is unnecessary (harmonics of w)."""
    ws = np.linspace(0, 2 * np.pi, nw, endpoint=False)
    q = 1 - e * e
    return max(abs(q / e * zonal_averages(n, e, inc, w, nq)["Phi_w"]) for w in ws)
