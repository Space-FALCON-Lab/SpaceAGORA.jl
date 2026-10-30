"""Stage K selftest (synthetic orbit, no flight data): checks stagek_lib against an independent integration and against derived identities.

Usage: python3 -I stagek_selftest.py       (exit code 1 on failure; prints each check)
  a  the one-orbit zonal map (hand-written Legendre-recurrence force) reproduces a sympy-gradient integration of the same potential (degrees 2..6)
     to 0.05 m in r_p after 1 and after 5 chained orbits;
  b  J2 only, Mars-like orbit: the pass-to-pass change of r_p is small (|dr| < 5 m; informational magnitude printed), and J3 added changes it
     by tens of m (the J2 first-order short-period theory is NOT used: see stagek_lib docstring);
  c  apsidal rate: first-order zonal average (kaula_lib) vs the numerical map within 5 % (J2..J20);
  d  analytic J3-and-higher secular term vs the numerical map (informational ratio, bound 0.5..1.5);
  e  identities: lever and Pi_drift round trip; burn apsidal gain at e = 0.82 against the spec's (2/e): ratio at f = 0 and its maximum.
"""
import pathlib, sys, math
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagek_lib as K

P = K.planet("Mars"); mu, R = P["mu"], P["R"]
fails = []
def check(name, ok, msg):
    print(("PASS " if ok else "FAIL ") + name + ": " + msg)
    if not ok: fails.append(name)

# --- independent sympy force
x, y, z = sp.symbols("x y z", real=True); rr = sp.sqrt(x * x + y * y + z * z)
def sym_accel(nmax):
    V = -mu / rr * (1 - sum(P["J"][n] * (R / rr) ** n * sp.legendre(n, z / rr) for n in range(2, nmax + 1)))
    return sp.lambdify((x, y, z), [-sp.diff(V, s) for s in (x, y, z)], "numpy")

def sym_passes(nmax, a, e, inc, Om, w, npass):
    g = sym_accel(nmax); r0, v0 = K.periapsis_state(mu, a, e, inc, Om, w); T = 2 * math.pi * math.sqrt(a ** 3 / mu)
    f = lambda t, u: np.concatenate([u[3:], np.array(g(*u[:3]))])
    ev = lambda t, u: u[:3] @ u[3:]; ev.direction = 1
    s0 = solve_ivp(f, [0, 1e-3 * T], np.concatenate([r0, v0]), method="DOP853", rtol=1e-12, atol=1e-9)
    sol = solve_ivp(f, [1e-3 * T, (npass + 0.5) * T], s0.y[:, -1], method="DOP853", rtol=1e-12, atol=1e-9, events=ev)
    return [K.elements(mu, u[:3], u[3:]) for u in sol.y_events[0]]

a0, e0, i0, O0, w0 = 1.6e7, 0.78, math.radians(40), math.radians(70), math.radians(60)
ref = sym_passes(6, a0, e0, i0, O0, w0, 5)
st = (a0, e0, i0, O0, w0); rp_chain = []
for k in range(5):
    out = K.zonal_map(P, *st, nmax=6); st = out[:5]; rp_chain.append(st[0] * (1 - st[1]))
err = [abs(rp_chain[k] - ref[k][0] * (1 - ref[k][1])) for k in range(5)]
check("a_map_vs_sympy", err[0] < 0.05 and err[4] < 0.05, f"|dr_p| after 1 and 5 orbits {err[0]:.4f} m, {err[4]:.4f} m")

# --- J2 only and J2 + J3 per-pass change from the map
def dr_map(nmax, a=a0, e=e0, i=i0, O=O0, w=w0):
    o = K.zonal_map(P, a, e, i, O, w, nmax=nmax); return o[0] * (1 - o[1]) - a * (1 - e), o
d2, o2 = dr_map(2); d3, o3 = dr_map(3)
check("b_J2_small", abs(d2) < 5.0 and abs(d3 - d2) > 10.0, f"J2-only dr_p = {d2:+.3f} m/pass; J3 adds {d3 - d2:+.2f} m/pass")

# --- apsidal rate
wd_map = ((K.zonal_map(P, a0, e0, i0, O0, w0, nmax=20)[4] - w0 + math.pi) % (2 * math.pi)) - math.pi
T = 2 * math.pi * math.sqrt(a0 ** 3 / mu); _, wd, _ = K.zonal_secular(P, a0, e0, i0, w0)
check("c_apsidal_rate", abs(wd * T / wd_map - 1) < 0.05, f"map {math.degrees(wd_map):.4f} deg/pass, first-order average {math.degrees(wd * T):.4f} deg/pass")

# --- analytic secular (n >= 3) vs numerical map (J3 only difference); evaluated at mid-pass omega
ed, _, _ = K.zonal_secular(P, a0, e0, i0, w0 + 0.5 * wd * T, nlo=3, nhi=3)
ana = -a0 * ed * T; ratio = (d3 - d2) / ana
check("d_J3_analytic", 0.5 < ratio < 1.5, f"numerical/analytic = {ratio:.3f} (numerical {d3 - d2:+.2f} m, analytic {ana:+.2f} m per pass)")

# --- identities
lev = K.lever_m_per_mps(mu, a0, e0)
check("e_lever_pi", abs(K.pi_drift(a0, e0, lev) * math.sqrt(mu / a0) - 1) < 1e-12, f"lever {lev:.1f} m per m/s")
e_ = 0.82; fs = np.linspace(0, math.pi, 2001); g_ = K.burn_apsidal_gain(e_, fs)
check("e_burn_gain", abs((2 / e_) / g_[0] - 3.4906) < 0.01, f"spec 2/e = {2 / e_:.3f}; gain at f=0 {g_[0]:.3f} (ratio {(2 / e_) / g_[0]:.2f}); max {g_.max():.3f} at f = {math.degrees(fs[g_.argmax()]):.1f} deg (ratio {(2 / e_) / g_.max():.2f})")
print("SELFTEST", "FAILED: " + ",".join(fails) if fails else "OK")
sys.exit(1 if fails else 0)
