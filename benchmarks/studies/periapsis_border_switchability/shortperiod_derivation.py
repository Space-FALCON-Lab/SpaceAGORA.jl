"""Symbolic derivation of the first-order short-period and slow periapsis-radius terms of a third-body tide (sympy), plus a numerical check
against direct integration of the equations of motion (scipy DOP853, frozen perturber, quadrupole force).

Usage: python3 -I shortperiod_derivation.py [--results-root DIR] [--skip-harmonics]
Writes <root>/estimator_v2/shortperiod_derivation.txt (the closed forms printed by sympy and the numerical check).

Setup (all in the orbit frame: x along the apse line, y in the plane toward the semilatus rectum, K = mu_b / r_b^3, Pi = K a^3 / mu = Pi_3B):
  quadrupole force  f = K (3 b (b.r) - r),  potential R2 = (K/2) [3 (b.r)^2 - r^2];   X = cos E - e, Y = eta sin E (position / a);   eta = sqrt(1 - e^2).
Semi-major axis. With the perturber frozen the energy is conserved, v^2/2 - mu/r - R2 = const, hence exactly
  -mu/(2a_osc) = C + R2(r)   =>   da = (2 a^2 / mu) (R2 - <R2>_M) = a Pi (q - <q>),   q = 3 (bx X + by Y)^2 - (X^2 + Y^2).
Angular momentum. d h_z / dt = 3 K a^2 G,  G = bx by (X^2 - Y^2) + (by^2 - bx^2) X Y  (a trigonometric polynomial in E), so
  d h_z / h = 3 Pi I(E) / eta,   I(E) = int (G - <G>)(1 - e cos E) dE + const,  <I>_M = 0.
Eccentricity from h^2 = mu a (1 - e^2):  d e = -(eta^2 / 2e) (2 dh/h - da/a);   periapsis radius d r_p = (1 - e) da - a de.
Slow (orbit-averaged) rate by Lagrange's equation:  d r_p / dt = (eta / (n a e)) dR_bar/d omega with d bx / d omega = by, d by / d omega = -bx.
"""
import argparse, pathlib, sys, time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--skip-harmonics", action="store_true")
args = ap.parse_args()
out = pathlib.Path(results_root(args.results_root)) / "estimator_v2"; out.mkdir(parents=True, exist_ok=True)
log = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)

E = sp.symbols("E", real=True); e = sp.symbols("e", positive=True); eta = sp.symbols("eta", positive=True)
bx, by = sp.symbols("bx by", real=True)
def poly_eta(expr):
    expr = sp.expand(expr); res = 0
    for t in sp.Add.make_args(expr):
        c, pw = t.as_coeff_exponent(eta)
        res += c * (1 - e**2)**(pw // 2) * eta**(pw % 2)
    return sp.expand(res)
X = sp.cos(E) - e; Y = eta * sp.sin(E)
def avgM(f):
    g = sp.expand(poly_eta(f) * (1 - e * sp.cos(E)))
    return sp.simplify(sp.integrate(g, (E, 0, 2 * sp.pi)) / (2 * sp.pi))

q = poly_eta(3 * (bx * X + by * Y)**2 - (X**2 + Y**2)); qa = avgM(q)
G = poly_eta(bx * by * (X**2 - Y**2) + (by**2 - bx**2) * X * Y); Ga = avgM(G)
say("<q>_M =", sp.factor(qa)); say("<G>_M =", sp.factor(Ga))
Iint = sp.expand(sp.integrate(sp.expand(poly_eta((G - Ga) * (1 - e * sp.cos(E)))), E)); Iint = sp.expand(Iint - avgM(Iint))
def at(Ev, expr): return sp.simplify(poly_eta(sp.expand(expr.subs(E, Ev))))
for nm, Ev in (("periapsis", 0), ("apoapsis", sp.pi)):
    dq = at(Ev, q - qa); Iv = at(Ev, Iint)
    da = dq                                                   # da / (a Pi)
    de = sp.simplify(-(eta**2 / (2 * e)) * (6 * Iv / eta - dq))   # de / Pi
    drp = sp.simplify((1 - e) * da - de)                      # d r_p / (a Pi)
    say(f"{nm}: da/(a Pi) =", sp.collect(sp.expand(da.subs(eta, sp.sqrt(1 - e**2))), [bx, by]))
    say(f"{nm}: de/Pi     =", sp.collect(sp.expand(de.subs(eta, sp.sqrt(1 - e**2))), [bx, by]))
    say(f"{nm}: dr_p/(a Pi) =", sp.collect(sp.expand(drp.subs(eta, sp.sqrt(1 - e**2))), [bx, by]))

# Slow rate: averaged quadrupole and octupole disturbing functions
bdr = bx * X + by * Y; r2 = X**2 + Y**2
R2a = sp.expand(avgM(sp.Rational(1, 2) * (3 * bdr**2 - r2)))
R3a = sp.expand(avgM(sp.Rational(1, 2) * (5 * bdr**3 - 3 * r2 * bdr)))
dw = lambda R: sp.diff(R, bx) * by - sp.diff(R, by) * bx
D2 = sp.factor(sp.simplify(eta / e * dw(R2a))); D3 = sp.factor(sp.simplify(eta / e * dw(R3a)))
say("<R2>/(mu_b a^2/r_b^3) =", sp.collect(R2a, [bx, by])); say("<R3>/(mu_b a^3/r_b^4) =", sp.collect(R3a, [bx, by]))
say("d r_p/dt / (n a Pi): quadrupole D2 =", D2, "   octupole (times a/r_b) D3 =", D3)

if not args.skip_harmonics:
    lam, om, inc = sp.symbols("lambda omega i", real=True)
    bxl = sp.cos(lam) * sp.cos(om) + sp.sin(lam) * sp.sin(om) * sp.cos(inc)
    byl = -sp.cos(lam) * sp.sin(om) + sp.sin(lam) * sp.cos(om) * sp.cos(inc)
    for nm, D in (("quadrupole", D2), ("octupole", D3)):
        Dt = sp.expand(sp.expand_trig(D.subs({bx: bxl, by: byl})))
        for k in range(0, 4):
            ca = sp.simplify(sp.integrate(sp.expand(Dt * sp.cos(k * lam)), (lam, 0, 2 * sp.pi)) / (sp.pi if k > 0 else 2 * sp.pi))
            sa = sp.simplify(sp.integrate(sp.expand(Dt * sp.sin(k * lam)), (lam, 0, 2 * sp.pi)) / sp.pi) if k > 0 else 0
            if ca != 0 or sa != 0:
                say(f"{nm}: harmonic k={k}: cos coeff = {sp.factor(ca)};  sin coeff = {sp.factor(sa)}")

# ---- numerical check against direct integration (quadrupole force, frozen perturber) ----
mu = 3.986004415e14; a = 3.68e7; ee = 0.8233; et = np.sqrt(1 - ee * ee); n = np.sqrt(mu / a**3); Pi = 1e-4; K = Pi * mu / a**3
rng = np.random.default_rng(1); b = rng.normal(size=3); b /= np.linalg.norm(b); bxn, byn, _ = b
def rhs(t, s):
    r = s[:3]; rn = np.linalg.norm(r)
    return np.concatenate([s[3:], -mu * r / rn**3 + K * (3 * b * (b @ r) - r)])
E0 = 2.0
r0 = np.array([a * (np.cos(E0) - ee), a * et * np.sin(E0), 0.0]); v0 = (n * a / (1 - ee * np.cos(E0))) * np.array([-np.sin(E0), et * np.cos(E0), 0.0])
T = 2 * np.pi / n; ts = np.linspace(0, T, 41)
sol = solve_ivp(rhs, [0, T], np.concatenate([r0, v0]), t_eval=ts, rtol=1e-13, atol=1e-6, method="DOP853")
def elems(s):
    r, v = s[:3], s[3:]; rn = np.linalg.norm(r); aa = 1 / (2 / rn - v @ v / mu); h = np.cross(r, v); ev = np.cross(v, h) / mu - r / rn
    return aa, np.linalg.norm(ev)
el = np.array([elems(sol.y[:, k]) for k in range(len(ts))])
Mt = E0 - ee * np.sin(E0) + n * ts; Et = Mt.copy()
for _ in range(60): Et -= (Et - ee * np.sin(Et) - Mt) / (1 - ee * np.cos(Et))
Eg = np.linspace(0, 2 * np.pi, 400001)
def parts(Eg):
    Xg = np.cos(Eg) - ee; Yg = et * np.sin(Eg)
    qg = 3 * (bxn * Xg + byn * Yg)**2 - (Xg**2 + Yg**2) - (12 * bxn**2 * ee**2 + 3 * bxn**2 - 3 * byn**2 * ee**2 + 3 * byn**2 - 3 * ee**2 - 2) / 2
    Gg = bxn * byn * (Xg**2 - Yg**2) + (byn**2 - bxn**2) * Xg * Yg - 2.5 * bxn * byn * ee**2
    f = Gg * (1 - ee * np.cos(Eg)); I = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(Eg))])
    return qg, I
qg, Ig = parts(Eg); Ig = Ig - getattr(np, "trapezoid", getattr(np, "trapz", None))(Ig * (1 - ee * np.cos(Eg)), Eg) / (2 * np.pi)
Em = np.mod(Et, 2 * np.pi); dq = np.interp(Em, Eg, qg); dI = np.interp(Em, Eg, Ig)
dA = a * Pi * dq; dE = -(et**2 / (2 * ee)) * (6 * Pi * dI / et - Pi * dq) + (-(15 / 2) * Pi * n * ee * et * bxn * byn) * ts
pa = dA - dA[0]; pe = dE - dE[0]; prp = (1 - ee) * pa - a * pe
oa = el[:, 0] - el[0, 0]; oe = el[:, 1] - el[0, 1]; orp = el[:, 0] * (1 - el[:, 1]) - el[0, 0] * (1 - el[0, 1])
say(f"direct integration vs closed form over one orbit (Pi = {Pi}, e = {ee}): max |a error| = {np.abs(pa - oa).max():.2f} m of a range {np.ptp(oa):.0f} m; "
    f"max |e error| = {np.abs(pe - oe).max():.2e} of range {np.ptp(oe):.2e}; max |r_p error| = {np.abs(prp - orp).max():.2f} m of range {np.ptp(orp):.0f} m "
    f"(second order in Pi: relative {np.abs(prp - orp).max() / np.ptp(orp):.1e})")
(out / "shortperiod_derivation.txt").write_text("\n".join(log) + "\n")
