"""Selftest of analytic_lib.py (WP-A2). Reproduces the DERIVED (plan) numbers of analytic_framework_and_transfer_plan.md that depend only on constants
(tolerance = the plan's quoted precision, stated per line; mismatches are printed, never hidden), and checks each new closed form against an independent numeric.
Usage: python3 -I analytic_selftest.py   (exit 1 on any hard failure)
"""
import math, pathlib, sys
import numpy as np
from scipy.integrate import solve_ivp
from statistics import NormalDist

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analytic_lib as L, apoapsis_pi_lib as AP, stagem_lib as SM, tesseral_lib as TL

MU = 42828.375e9; R = 3396.19e3; H = 6308.278108; BETA = 20.2402; RP = 3504.4675e3     # SOURCED: R1.2 mu, planets.jl R/H, dev_pairs beta, Stage F truth pass-1 radius
RA0, RAT = 28533e3, 7408e3                                                             # SOURCED Odyssey main phase (plan 2.1)
bad = []


def check(name, got, want, tol, note=""):
    ok = abs(got - want) <= tol * abs(want) if want != 0 else abs(got) <= tol
    print(f"{'ok  ' if ok else 'MISMATCH'} {name}: got {got:.6g} want {want:.6g} (rel tol {tol:g}) {note}")
    if not ok: bad.append(name)


def plan_numbers():
    check("K(f_rot=1.05)", SM.k_theory(H, 20.24, 1.05), 41.31, 1e-3, "plan 2 'Drop law verified'")
    check("K(f_rot=1.00)", SM.k_theory(H, 20.24, 1.0), 39.35, 1e-3)
    ra, rp = 3900e3, 3501e3; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
    check("drop per W/cm2 at 3900/3501 km [km]", SM.drop_per_pass(a, e, rp, 1e4, SM.k_theory(H, 20.24, 1.05), MU) / 1e3, 304.3, 2e-3)
    e0 = (RA0 - RP) / (RA0 + RP); et = (RAT - RP) / (RAT + RP); kq = AP.k_Q(RP, 0.224e4, H, 20.24, MU, 1.0)
    check("passes at constant Q 0.224", AP.passes_to_go(e0, et, kq), 164.5, 2e-3, "plan 5.1 (r_p = Stage F pass-1 radius; the plan leaves r_p unstated)")
    check("days at constant Q 0.224", AP.time_to_go(RP, e0, et, kq, MU) / 86400, 52.9, 2e-3)
    kq54 = AP.k_Q(RP, 0.54e4, H, 20.24, MU, 1.0)
    check("passes at 0.54", AP.passes_to_go(e0, et, kq54), 68.2, 3e-3); check("days at 0.54", AP.time_to_go(RP, e0, et, kq54, MU) / 86400, 21.9, 3e-3)
    for lab, ra_, want in (("start", RA0, 1984), ("7408 km", RAT, 309), ("3900 km", 3900e3, 31)):
        a_ = 0.5 * (ra_ + 3501e3); e_ = (ra_ - 3501e3) / (ra_ + 3501e3); check(f"c = a e / H at {lab}", a_ * e_ / H, want, 5e-3 if want > 100 else 3e-2, "(plan quotes 31 for 31.6: two-digit rounding)" if want < 100 else "")
    lev0 = SM.__dict__.get("lever_m_per_mps", None)
    import stagek_lib as SK
    a0 = 0.5 * (RA0 + 3501e3); a1 = 0.5 * (RAT + 3501e3); e0b = (RA0 - 3501e3) / (RA0 + 3501e3); e1b = (RAT - 3501e3) / (RAT + 3501e3)
    check("lever start [km per m/s]", SK.lever_m_per_mps(MU, a0, e0b) / 1e3, 13.7, 5e-3); check("lever 7408 km [km per m/s]", SK.lever_m_per_mps(MU, a1, e1b) / 1e3, 5.4, 1e-2)
    check("heat load total [MJ/m2], f_rot 1", 20.24 * MU / 2 * (1 / a1 - 1 / a0) / 1e6, 52.4, 3e-3, "plan 5.2(e): (mu beta / 2 f_rot)(1/a_t - 1/a_0)")
    tau = lambda e, rp: AP.tau_eff(rp, e, H, MU)
    check("per-pass load at start [kJ/m2] (Q 0.224)", 0.224e4 * tau(e0b, 3501e3) / 1e3, 270, 1e-2); check("per-pass load at 7408 km [kJ/m2]", 0.224e4 * tau(e1b, 3501e3) / 1e3, 399, 1e-2)
    for lab, a_, e_, want in (("start", a0, e0b, 0.694), ("7408 km", a1, e1b, 0.138)):
        check(f"nu at {lab}", L.nu_commens(a_, MU, L.PLANETS["mars"]["w_rot"]), want, 3e-2, "plan 'nu 0.694 to 0.138'")
    # hybrid marginal value d t / d dv (plan 5.4): T(e) sqrt(e (1+e)) / (k_Q V_c), per m/s, at the start
    Vc = math.sqrt(MU / RP); T0 = 2 * math.pi * math.sqrt((RP / (1 - e0)) ** 3 / MU)
    check("hybrid d t/d dv at start [d per m/s]", T0 * math.sqrt(e0 * (1 + e0)) / (kq * Vc) / 86400, 0.249, 5e-3, "plan 5.4")
    # SRP
    a_srp = L.srp_accel(2.2 * 1.7 + 2 * 1.945 * 1.7, 461.0, 1.3, 1.524 * L.AU)
    check("a_SRP [m/s2]", a_srp, 5.7e-8, 1e-2, "plan T9 (A/m from the drag area, C_R 1.3 ASSUMED, TSI 1361 ASSUMED)")
    check("SRP |dr_p| per orbit at start [m]", L.srp_drp_bound(a0, e0b, MU, a_srp), 32.0, 3e-2); check("SRP |dr_p| per orbit at 7408 km [m]", L.srp_drp_bound(a1, e1b, MU, a_srp), 2.0, 0.15)


def e_table():
    """Plan 5.2 E-table: uniform eps, Gaussian quantile on a lognormal sigma 0.32, realized mean Q = Q_plan exp(sigma^2/2); passes and days from the constant-Q law."""
    nd = NormalDist(); sig = 0.32; qa = 0.54; e0 = (RA0 - RP) / (RA0 + RP); et = (RAT - RP) / (RAT + RP)
    def run(E):
        q = 0.25
        for _ in range(200):
            kq = AP.k_Q(RP, q * math.exp(0.5 * sig ** 2) * 1e4, H, 20.24, MU, 1.0); N = AP.passes_to_go(e0, et, kq)
            z = nd.inv_cdf(1 - min(E / N, 0.5)); q = 0.5 * q + 0.5 * qa * math.exp(-z * sig)
        kq = AP.k_Q(RP, q * math.exp(0.5 * sig ** 2) * 1e4, H, 20.24, MU, 1.0); return q, AP.passes_to_go(e0, et, kq), AP.time_to_go(RP, e0, et, kq, MU) / 86400
    for E, (qw, dw, nw) in {0.05: (0.177, 63.6, 198), 0.5: (0.226, 50.0, 155), 1.0: (0.246, 45.9, 143), 1.8: (0.267, 42.4, 132)}.items():
        q, N, d = run(E); check(f"E={E} Q_plan", q, qw, 3e-2, "plan 5.2 table"); check(f"E={E} days", d, dw, 3e-2); check(f"E={E} passes", N, nw, 3e-2)


def bessel_vs_quadrature():
    worst = 0.0
    for a, e, rho in [(3.7e6, 0.05, 1e-8), (3.7e6, 0.02, 1e-8), (4.2e6, 0.1, 1e-8), (5e6, 0.3, 2e-9), (8e6, 0.5, 1e-9), (1.6e7, 0.78, 1e-10)]:
        da, de, drp = L.bessel_drag_pass(a, e, H, rho, BETA, MU, N=30); q = AP.gauss_pass_quad(a, e, H, rho, BETA, MU)
        err = max(abs(da / q[0] - 1), abs(de / q[1] - 1)); worst = max(worst, err)
    check("Bessel (N=30) vs Gauss quadrature, worst of da/de relative error (e <= 0.78)", worst, 0.0, 1e-3, "(truncation O(e^31); 1.3e-4 at e = 0.78, < 1e-9 for e <= 0.3)")
    # circular limit: da -> -2 pi a^2 rho / beta
    da, de, _ = L.bessel_drag_pass(7e6, 1e-9, H, 1e-9, BETA, MU, N=6); check("circular limit da", da, -2 * math.pi * 7e6 ** 2 * 1e-9 / BETA, 1e-5)
    # Laplace limit at c >> 1
    a, e = 1.0e7, 0.4; rho = 1e-9; da, de, _ = L.bessel_drag_pass(a, e, H, rho, BETA, MU, N=30); rp = a * (1 - e)
    Q = 0.5 * rho * (MU * (1 + e) / rp) ** 1.5; check("Laplace limit da (c = %.0f)" % (a * e / H), da, AP.delta_a_pass(a, e, rp, Q, H, BETA, MU), 2e-3)


def tide_vs_numeric():
    """Direct tide quadrature and the averaged closed form against a DOP853 integration of Kepler + tidal acceleration (frozen perturber)."""
    mu = 3.98600436233e14; mub = 4.9028e15; rb = np.array([2.0e8, 2.5e8, 0.8e8]) * 10.0;   # far perturber (a/r_b = 0.007): quadrupole regime; mu_b scaled to keep the tide measurable
    a, e = 2.4e7, 0.8; inc, raan, om = 0.8, 0.4, 1.1
    Q = TL.pqw_matrix(inc, raan, om); r0 = Q @ np.array([a * (1 - e), 0, 0]); v0 = Q @ np.array([0, math.sqrt(mu * (1 + e) / (a * (1 - e))), 0])
    def acc(t, s):
        r = s[:3]; d = rb - r; return np.concatenate([s[3:], -mu * r / np.linalg.norm(r) ** 3 + mub * (d / np.linalg.norm(d) ** 3 - rb / np.linalg.norm(rb) ** 3)])
    T = 2 * math.pi * math.sqrt(a ** 3 / mu); ev = lambda t, s: s[:3] @ s[3:]; ev.direction = 1
    sol = solve_ivp(acc, [0, 1.2 * T], np.concatenate([r0, v0]), method="DOP853", rtol=1e-12, atol=1e-6, events=ev, first_step=1.0)
    te = [t for t in sol.t_events[0] if t > 0.5 * T]; s1 = sol.sol(te[0]) if sol.sol is not None else sol.y_events[0][-1]
    s1 = sol.y_events[0][np.argmax(np.array(sol.t_events[0]) > 0.5 * T)]; rp1 = TL.elements(s1[:3], s1[3:], mu)
    rp_num = rp1[0] * (1 - rp1[1]) - a * (1 - e)
    k = L.tide_arc_direct(mu, mub, a, e, inc, raan, om, lambda tau: np.repeat(rb[:, None], np.size(tau), axis=1), 96)
    avg = L.tide_rp_averaged(a, e, mu, mub, rb, Q[:, 0], Q[:, 1])
    check("T8 direct quadrature vs DOP853, one orbit r_p change", k["rp"], rp_num, 2e-2, "first order in the tide")
    check("T8 averaged closed form vs DOP853", avg, rp_num, 2e-2, "(quadrupole; octupole error O(a/r_b))")
    check("T8 averaged vs direct", avg, k["rp"], 2e-2)
    rb2 = rb / 10.0; mub2 = mub / 1000.0; k2 = L.tide_arc_direct(mu, mub2, a, e, inc, raan, om, lambda tau: np.repeat(rb2[:, None], np.size(tau), axis=1), 96)
    avg2 = L.tide_rp_averaged(a, e, mu, mub2, rb2, Q[:, 0], Q[:, 1]); print(f"info octupole size at a/r_b = {a / np.linalg.norm(rb2):.3f}, e = {e}: averaged/direct - 1 = {avg2 / k2['rp'] - 1:+.3f}")


def jacobi_identity():
    """dE_K - w_rot dh_z = U(end) - U(start) for the first-order Gauss kick of a body-fixed field over one periapsis-to-periapsis arc (derived in analytic_lib.tess_arc)."""
    cs = TL.coefficients("Mars", 20); Ct, St = TL.select(cs["C"], cs["S"], lambda n, m: m >= 1); W = L.PLANETS["mars"]["w_rot"]
    worst = 0.0
    for a, e, inc, raan, om, th in [(1.6e7, 0.78, 1.63, 0.5, 1.9, 0.3), (1.0e7, 0.65, 0.44, 6.2, 0.05, 2.0), (8e6, 0.5, 2.0, 1.0, 4.0, 5.0)]:
        for lo, hi in ((0.0, 2 * math.pi), (-math.pi, math.pi)):
          k = L.tess_arc(cs, Ct, St, a, e, inc, raan, om, th, W, 128, with_boundary=True, f_lo=lo, f_hi=hi)
        # first-order identity; note: the Gauss kick is first order, the boundary term is exact U(r_p, theta1) - U(r_p, theta0): they agree to O(kick^2)
          lhs = k["dE"] - W * k["dhz"]; worst = max(worst, abs(lhs - k["dU"]) / max(abs(k["dU"]), abs(k["dE"]), 1e-9))
    check("Jacobi: (dE - w dhz) vs dU, worst relative", worst, 0.0, 5e-2, "(first order; Gauss kick vs exact boundary potential)")


def t2_check():
    """dh_geom against the small-step first-order rate R_eq f sin^2 i sin(2 omega) d_om, and the surface radius at the poles/equator."""
    P = L.PLANETS["mars"]; inc, om, dw = 1.632, 1.9, 1e-4
    check("T2 dh_geom vs R_eq f sin^2 i sin(2 om) d_om", L.dh_geom(inc, om, dw, "mars", first_order=True), P["Re"] * P["f"] * math.sin(inc) ** 2 * math.sin(2 * (om + 0.5 * dw)) * dw, 1e-6)
    check("T2 surface radius at equator (i = 0)", L.surface_radius(0.0, 1.0, "mars"), P["Re"], 1e-12); check("T2 surface radius at the pole (i = 90, omega = 90)", L.surface_radius(math.pi / 2, math.pi / 2, "mars"), P["Rp"], 1e-12)
    d = abs(L.surface_radius(1.0, 1.2, "mars") - L.surface_radius(1.0, 1.2, "mars", True)); print(f"info exact vs first-order surface radius at i = 1.0, omega = 1.2: {d:.1f} m (second order in f)")
    if not d < 50.0: bad.append("T2 second-order size")
    J2, R = 1.958616126805748e-3, 3394200.0; print(f"info Kozai a_osc - a_mean at a = 16000 km, e = 0.78, i = 93.5 deg, omega = 1.9: {L.a_osc_minus_mean(1.6e7, 0.78, 1.632, 1.9, J2, R) / 1e3:.1f} km (checked through nu on the Stage O tables)")


def planet_pi_groups():
    for pl in ("mars", "earth", "venus"):
        print(f"info Pi_f {pl}: {L.pi_f(pl):.4g}")
    check("Pi_f Mars", L.pi_f("mars"), 3.17, 5e-3, "plan 3.2: 19.99 km / 6.308 km")
    # spec Pi_drag vs AP.pi_q:  Pi_q = Pi_drag sqrt(1+e)  (Pi_drag = 0.5 rho_p beta^-1 r_p (1+e), Q = 0.5 rho V_p^3, V_p = V_c sqrt(1+e))
    rho, rp, e = 1e-8, 3.5e6, 0.7; Q = 0.5 * rho * (MU * (1 + e) / rp) ** 1.5
    e0, e1 = 0.7, 0.6; P0 = 1.0 / (rp * (1 + e0)) ** 2; P1 = 1.0 / (rp * (1 + e1)) ** 2
    check("spec 5.2 Pi_J2,k+1 = Pi_J2,k ((1+e_k)/(1+e_k+1))^2 at fixed r_p (Pi_J2 ~ p^-2, p = r_p (1+e))", P1 / P0, ((1 + e0) / (1 + e1)) ** 2, 1e-12)
    check("spec 5.2 Pi_3B,k+1 = Pi_3B,k ((1-e_k)/(1-e_k+1))^3 at fixed r_p (a = r_p/(1-e))", L.pi_3b(1.0, 1.0, rp / (1 - e1), 1.0) / L.pi_3b(1.0, 1.0, rp / (1 - e0), 1.0), ((1 - e0) / (1 - e1)) ** 3, 1e-12)
    check("Pi_q = Pi_drag sqrt(1+e) (spec 2.4 vs apoapsis-Pi)", AP.pi_q(rp, Q, BETA, MU), 0.5 * rho / BETA * rp * (1 + e) * math.sqrt(1 + e), 1e-12)


if __name__ == "__main__":
    plan_numbers(); e_table(); bessel_vs_quadrature(); tide_vs_numeric(); jacobi_identity(); t2_check(); planet_pi_groups()
    print("\nSELFTEST", "OK" if not bad else f"MISMATCHES: {bad}"); sys.exit(1 if bad else 0)
