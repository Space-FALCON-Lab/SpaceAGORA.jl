"""Self-tests of tesseral_lib.py (python3 -I tesseral_selftest.py). Prints PASS/FAIL per test; exits nonzero on a failure.

1 acceleration = -grad(potential) by central differences (Mars50c 20x20, random body points at 1.0-1.5 R)
2 zonal part = kaula_cowell.make_rhs zonal acceleration (J_n = -C_n0), independent implementation
3 Kaula's closed-form F_nmp(i) (recalled) against the defining Fourier coefficients D_{nm,n-2p}(i) (x 1 for n-m even, x (-i) for n-m odd)
4 Hansen + inclination-function reconstruction of the (n, m) potential at orbit points against the direct potential
5 Lagrange-series per-pass kick (all j) against the Gauss quadrature kick along the same Kepler arc, (n, m) = (2,2), (3,1), (5,3), (12,7)
"""
import math, sys, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import tesseral_lib as tl
import kaula_cowell as kc

fails = []
def check(name, val, tol):
    ok = val <= tol; print(f"{'PASS' if ok else 'FAIL'}  {name}: {val:.3e} (tol {tol:.0e})");  (None if ok else fails.append(name))

cs = tl.coefficients("Mars", 20); R, mu = cs["R"], cs["mu"]
rng = np.random.default_rng(1)
P = rng.normal(size=(3, 6)); P = P / np.linalg.norm(P, axis=0) * R * rng.uniform(1.0, 1.5, 6)
acc = tl.accel_body(P[0], P[1], P[2], cs)
h = 1.0; g = np.zeros_like(acc)
for d in range(3):
    dp = np.zeros((3, 1)); dp[d] = h
    up = tl.potential_body(*(P + dp), cs); um = tl.potential_body(*(P - dp), cs); g[d] = (up - um) / (2 * h)
check("1 accel vs grad U (relative, max)", np.max(np.abs(acc - g)) / np.max(np.abs(acc)), 1e-6)

Cz, Sz = tl.select(cs["C"], cs["S"], lambda n, m: m == 0)
az = tl.accel_body(P[0], P[1], P[2], cs, Cz, Sz)
Js = {n: -cs["C"][n, 0] for n in range(2, 21)}
rhs = kc.make_rhs(Js, mu=mu, R=R)
ref = np.array([rhs(0, np.concatenate([P[:, k], [0, 0, 0]]))[3:] for k in range(P.shape[1])]).T + mu * P / np.linalg.norm(P, axis=0) ** 3
check("2 zonal accel vs kaula_cowell (relative, max)", np.max(np.abs(az - ref)) / np.max(np.abs(az)), 1e-11)

worst = 0.0
for (n, m) in [(2, 1), (2, 2), (3, 1), (3, 3), (4, 2), (5, 3), (7, 6), (10, 4), (12, 7)]:
    for inc in (0.44, 1.1, 1.63, 2.03):
        D = tl.incl_fn(n, m, inc)
        for p in range(n + 1):
            F = tl.kaula_F(n, m, p, inc); Dk = D[(n - 2 * p) + n] * (1 if (n - m) % 2 == 0 else 1j)
            worst = max(worst, abs(Dk - F) / max(1e-300, np.max(np.abs(D))))
check("3 Kaula F_nmp vs Fourier D (max abs / max|D|)", worst, 1e-9)

worst = 0.0
a, e, inc, raan, om = 1.6e7, 0.78, 0.44, 0.7, 2.1
for (n, m) in [(2, 2), (3, 1), (6, 5)]:
    T = tl.SeriesTerm(n, m, e, inc, NM=8192)
    Cn, Sn = tl.select(cs["C"], cs["S"], lambda nn, mm: (nn, mm) == (n, m))
    for M in (0.0, 0.3, -1.0, 3.0):
        E = M
        for _ in range(60): E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        f = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2)); r = a * (1 - e * math.cos(E))
        Q = tl.pqw_matrix(inc, raan, om); x = Q @ np.array([r * math.cos(f), r * math.sin(f), 0.0]); th = 0.9
        xb = tl.rot_z(np.array([-th]), x[:, None])
        Ud = tl.potential_body(xb[0], xb[1], xb[2], cs, Cn, Sn)[0]
        amp = (mu / a) * (R / a) ** n * (cs["C"][n, m] - 1j * cs["S"][n, m])
        ph = np.exp(1j * (T.ks * om + m * (raan - th)))[:, None] * np.exp(1j * T.j * M)[None, :]
        Us = np.real(amp * ((T.D[:, None] * T.X) * ph).sum())
        worst = max(worst, abs(Us - Ud) / abs(amp) / (a / r) ** (n + 1))
check("4 series reconstruction of U_nm (relative to the term amplitude)", worst, 1e-6)

worst = 0.0
w_rot = 7.08823596e-5
for (n, m) in [(2, 2), (3, 1), (5, 3), (12, 7)]:
    T = tl.SeriesTerm(n, m, e, inc, NM=8192)
    Cn, Sn = tl.select(cs["C"], cs["S"], lambda nn, mm: (nn, mm) == (n, m))
    N = math.sqrt(mu / a ** 3); nu = w_rot / N
    for thp in (0.0, 1.3):
        gk = tl.gauss_kick(a, e, inc, raan, om, thp, w_rot, cs, Cn, Sn, nf=512)
        sk = tl.series_kick(T, a, raan, om, thp, nu, cs)
        for key, scale in (("rp", None), ("a", None), ("i", None), ("raan", None), ("omega", None)):
            ref = max(abs(gk[key]), 1e-30)
            worst = max(worst, abs(sk[key] - gk[key]) / max(abs(gk["rp"]) / a if key in ("i", "raan", "omega") else abs(gk["rp"]), ref))
        print(f"   ({n},{m}) theta_p={thp}: rp kick quadrature {gk['rp']:.4f} m, series {sk['rp']:.4f} m;  i {gk['i']:.3e} / {sk['i']:.3e} rad")
check("5 series kick vs Gauss quadrature kick (relative)", worst, 1e-5)

Ct, St = tl.select(cs["C"], cs["S"], lambda n, m: m >= 1)
worst = 0.0
for thp in (0.0, 2.2):
    Y = tl.term_kicks(a, e, inc, raan, om, thp, w_rot, mu, R, 20, nf=256)
    k1 = float(np.real(((cs["Cbar"] - 1j * cs["Sbar"]) * Y).sum())); k2 = tl.gauss_kick(a, e, inc, raan, om, thp, w_rot, cs, Ct, St, nf=256)["rp"]
    worst = max(worst, abs(k1 - k2) / abs(k2)); print(f"   term_kicks sum {k1:.4f} m, gauss_kick {k2:.4f} m")
check("6 per-term unit kicks summed with the file coefficients vs gauss_kick (relative)", worst, 1e-10)

print("FAILED:" if fails else "ALL PASS", fails)
sys.exit(1 if fails else 0)
