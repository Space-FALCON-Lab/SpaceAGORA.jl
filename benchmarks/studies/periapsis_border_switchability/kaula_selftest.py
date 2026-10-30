"""Self-check of the zonal geometry (python3 -I kaula_selftest.py): textbook J2 rates, hand-derived J3 and sympy J4 closed forms, the Kozai third-body
coefficient 15/8, and the analytic J3 periapsis drift against the Cowell propagator (Mars-like orbit, 20 days)."""
import pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import kaula_lib as kl, kaula_cowell as kc
e, inc, w = 0.78, np.radians(37.0), np.radians(40.0); s2 = np.sin(inc) ** 2
ed, wd, Od = kl.zonal_rates(2, e, inc, w, 1.0, 1.0)
assert abs(ed) < 1e-12 and abs(wd - 0.75 * (4 - 5 * s2)) < 1e-12 and abs(Od + 1.5 * np.cos(inc)) < 1e-12          # J2: edot = 0, omega_dot, Omega_dot (textbook)
assert abs(kl.zonal_rates(3, e, inc, w, 1.0, 1.0)[0] - 3 / 8 * (1 - e * e) * np.sin(inc) * (5 * s2 - 4) * np.cos(w)) < 1e-12   # J3 closed form
assert abs(kl.zonal_rates(4, e, inc, w, 1.0, 1.0)[0] - 15 / 32 * e * (1 - e * e) * s2 * (7 * s2 - 6) * np.sin(2 * w)) < 1e-12  # J4 closed form (sympy)
# Kozai quadrupole third-body coefficient, perturber on a circular orbit in the reference plane
def Rhat(e, inc, w, nq=4096):
    nu = (np.arange(nq) + .5) * 2 * np.pi / nq; r = (1 - e * e) / (1 + e * np.cos(nu)); f = r ** 2 * (0.75 * (1 - (np.sin(inc) * np.sin(w + nu)) ** 2) - 0.5); d = r ** 2
    return np.sum(f * d) / np.sum(d)
h = 1e-6; num = -(np.sqrt(1 - e * e) / e) * (Rhat(e, inc, w + h) - Rhat(e, inc, w - h)) / (2 * h)
assert abs(num - 15 / 8 * e * np.sqrt(1 - e * e) * s2 * np.sin(2 * w)) < 1e-8
df, R, mu, _ = kl.read_gravity("Mars"); J = kl.zonal_J(df); a = 16.03e6; T = 20 * 86400.0; i0 = np.radians(25.2)
t2, r2 = kc.periapsis_series({2: J[2]}, a, e, i0, 0.0, 0.0, T, mu, R); t3, r3 = kc.periapsis_series({2: J[2], 3: J[3]}, a, e, i0, 0.0, 0.0, T, mu, R)
n = min(len(r2), len(r3)); num = (r3[:n] - r2[:n])[-1]
p = a * (1 - e * e); N = np.sqrt(mu / a ** 3); wdot = kl.zonal_rates(2, e, i0, 0.0, J[2] * (R / p) ** 2, N)[1]
tt = np.linspace(0, t3[n - 1], 3000); rate = np.array([-a * kl.zonal_rates(3, e, i0, wdot * t, J[3] * (R / p) ** 3, N, nq=1024)[0] for t in tt])
an = np.sum(0.5 * (rate[1:] + rate[:-1]) * np.diff(tt))
assert abs(num - an) < 0.02 * abs(num), (num, an)
print("kaula_selftest OK; J3 drift Cowell %.1f m vs analytic %.1f m" % (num, an))
