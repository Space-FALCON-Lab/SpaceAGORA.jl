"""Checks for the four derivations in tasks/.../analytic_derivations.md (Python only, no flights, no Julia).

  D1 drag-induced apse rotation      d_omega = d_e_drag * f_star / (e (1 + e)),  f_star = (1 + e) (R_s' + H kappa) / (r_p e)
  D2 pass-integrated vs event drag   F = sqrt(H_loc / H_0) * exp(r'' f_star^2 / (2 H)) * (1 + 3 H'/8) * F_pert
  D3 octupole of the averaged tide   dr_p,oct = -(5 pi / 8) a sqrt(1-e^2) (mu_b/mu) (a/r_b)^4 b_Q [(15 + 90 e^2) b_P^2 + 15 (1-e^2) b_Q^2 - 12 - 9 e^2]
  D4 all-e Bessel drag and dH/dh     exact Fourier coefficients of g, h as series in beta = e / (1 + sqrt(1-e^2)); H' enters as the polynomial factor 1 + (H' c^2 / 2)(1-x)^2

Usage: python3 -I analytic_derivations_check.py --selftest       (synthetic checks only, about 10 s)
       python3 -I analytic_derivations_check.py --data           (selftest + checks against the flight files on disk; first run builds a cache, a few minutes)
Labels as in the note: SOURCED (file named), DERIVED, ASSUMED.
"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.special import ive

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analytic_lib as L, apoapsis_pi_lib as AP, tesseral_lib as TL

MARS = L.PLANETS["mars"]; MU_MARS = 42828370371000.0   # SOURCED as analytic_pairs.MU_M (Mars50c.csv GM)


# ---- D1: drag-induced apse rotation -------------------------------------------------------------------------------------------------------
def surface_slope(inc, om, planet="mars"):
    """R_s' = dR_s/du at u = omega (first order in f): -R_eq f sin^2 i sin 2 omega [m/rad]."""
    P = L.PLANETS[planet]; return -P["Re"] * P["f"] * math.sin(inc) ** 2 * math.sin(2 * om)


def f_star(e, rp, inc, om, H=0.0, kappa=0.0, planet="mars"):
    """True anomaly of the drag-weight centroid: minimum of (h(f) - h_p)/H - kappa f with h = r(f) - R_s(omega + f), r'' = r_p e / (1 + e)."""
    return (surface_slope(inc, om, planet) + H * kappa) * (1 + e) / (rp * e)


def apse_rotation_drag(de_drag, e, rp, inc, om, H=0.0, kappa=0.0, planet="mars"):
    """Secular apse rotation of one drag pass [rad]: de_drag f_star / (e (1 + e)) = de_drag (R_s' + H kappa) / (r_p e^2).  de_drag < 0."""
    return de_drag * f_star(e, rp, inc, om, H, kappa, planet) / (e * (1 + e))


def kappa_trailing(om, rho_mean, alt, H_loc, W=40, nmin=8):
    """Along-track log-density gradient at fixed altitude, causal: slope of ln rho_mean + alt/H_loc on omega over the passes k-W..k (event points move
    along track by d omega each pass).  Returns NaN where fewer than nmin passes are available."""
    om = np.unwrap(np.asarray(om)); y = np.log(np.asarray(rho_mean, float)) + np.asarray(alt) / np.asarray(H_loc); out = np.full(len(om), np.nan)
    for k in range(len(om)):
        lo = max(0, k - W)
        if k + 1 - lo >= nmin: out[k] = np.polyfit(om[lo:k + 1], y[lo:k + 1], 1)[0]
    return out


# ---- D2: pass-integrated vs event drag -------------------------------------------------------------------------------------------------
def profile_factor(e, rp, H_loc, H0, fstar, Hprime=0.0):
    """Mean-field drag over the Laplace prediction made with the event density and H0: sqrt(H_loc/H0) exp(r'' f*^2 / (2 H_loc)) (1 + 3 H'/8)."""
    rpp = rp * e / (1 + e); return math.sqrt(H_loc / H0) * math.exp(rpp * fstar ** 2 / (2 * H_loc)) * (1 + 3 * Hprime / 8)


def pert_factor_mean(m, s):
    """Expected pass-integrated perturbation factor for a stationary lognormal factor exp(X), X ~ N(m, s^2), averaged over many correlation times."""
    return math.exp(m + 0.5 * s * s)


def pert_factor_event(sig_t, s, tau_c, n=241):
    """Median of (pass-integrated factor)/(event factor) for an Ornstein-Uhlenbeck log factor (sd s, correlation exp(-|t|/tau_c)) under a Gaussian drag
    window of sd sig_t centered on the event: int w(t) exp(s^2 (1 - rho(t)^2) / 2) dt / int w."""
    t = np.linspace(-6, 6, n); w = np.exp(-0.5 * t * t)
    return float(np.trapz(w * np.exp(0.5 * s * s * (1 - np.exp(-2 * np.abs(t * sig_t) / tau_c))), t) / np.trapz(w, t))


# ---- D3: octupole -----------------------------------------------------------------------------------------------------------------------
def tide_rp_octupole(a, e, mu, mu_b, rb_vec, P_hat, Q_hat):
    rb = np.linalg.norm(rb_vec); u = rb_vec / rb; bP, bQ = u @ P_hat, u @ Q_hat
    return -(5 * math.pi / 8) * a * math.sqrt(1 - e * e) * (mu_b / mu) * (a / rb) ** 4 * bQ * ((15 + 90 * e * e) * bP ** 2 + 15 * (1 - e * e) * bQ ** 2 - 12 - 9 * e * e)


# ---- D4: all-e Bessel drag, dH/dh ---------------------------------------------------------------------------------------------------------
def _binom_series(p, sgn, beta, K):
    c = np.ones(K + 1)
    for k in range(1, K + 1): c[k] = c[k - 1] * (p - k + 1) / k
    return c * (sgn * beta) ** np.arange(K + 1)


def fourier_cheb(p, q, e, K=None):
    """Chebyshev (= cosine) coefficients c_j of (1 + e cos E)^p (1 - e cos E)^q, exact up to the truncation of a geometric series in beta^2:
    1 +- e cos E = (1 +- beta z)(1 +- beta/z) / (1 + beta^2), z = exp(iE), beta = e / (1 + sqrt(1 - e^2)), so the function is (1+beta^2)^-(p+q) |F(z)|^2 with
    F(z) = (1 + beta z)^p (1 - beta z)^q = sum f_k z^k, and c_0 = (1+beta^2)^-(p+q) sum f_k^2, c_j = 2 (1+beta^2)^-(p+q) sum_k f_k f_(k+j)."""
    beta = e / (1 + math.sqrt(1 - e * e))
    if K is None: K = int(min(4000, max(20, math.ceil(math.log(1e-17) / math.log(max(beta, 1e-12))))))
    f = np.convolve(_binom_series(p, 1, beta, K), _binom_series(q, -1, beta, K))[: K + 1]
    full = np.correlate(f, f, mode="full")[K:]; c = 2 * full; c[0] = full[0]
    return c * (1 + beta * beta) ** (-(p + q))


def bessel_drag_pass_alle(a, e, H, rho_p, beta_b, mu, Hprime=0.0, f_rot=1.0):
    """Per-pass first-order (da, de, drp), tangential drag, any e < 1 and any c = a e / H (no e-series): da = -2 pi f (a^2 rho_p / beta) sum g_j e^-c I_j(c),
    de = -2 pi f (a rho_p (1-e^2)/beta) sum h_j e^-c I_j(c), g = (1+ex)^1.5 (1-ex)^-0.5, h = x (1+ex)^0.5 (1-ex)^-0.5.  Non-isothermal layer H(h) = H + H' (h - h_p):
    first order in H' the density gains the factor 1 + (H' c^2 / 2)(1 - x)^2, multiplied into g and h as a Chebyshev polynomial (1.5 T0 - 2 T1 + 0.5 T2)."""
    C = np.polynomial.chebyshev; c = a * e / H
    g = fourier_cheb(1.5, -0.5, e); hh = C.chebmul([0, 1], fourier_cheb(0.5, -0.5, e))
    if Hprime:
        corr = np.array([1.0, 0, 0]) + 0.5 * Hprime * c * c * np.array([1.5, -2.0, 0.5]); g = C.chebmul(g, corr); hh = C.chebmul(hh, corr)
    S = lambda co: float(np.sum(co * ive(np.arange(len(co)), c)))
    da = -2 * math.pi * f_rot * (a * a * rho_p / beta_b) * S(g); de = -2 * math.pi * f_rot * (a * rho_p * (1 - e * e) / beta_b) * S(hh)
    return da, de, (1 - e) * da - a * de


def gauss_quad_layer(a, e, H, rho_p, beta_b, mu, Hprime=0.0):
    """Reference: first-order Gauss integrals with the exact non-isothermal profile rho = rho_p (1 + H' dh/H)^(-1/H') (H(h) = H + H' dh)."""
    p = a * (1 - e * e); rp = a * (1 - e); h = math.sqrt(mu * p)
    def acc(nu):
        r = p / (1 + e * math.cos(nu)); v2 = mu * (2 / r - 1 / a); dh = r - rp
        rho = rho_p * (math.exp(-dh / H) if Hprime == 0 else (1 + Hprime * dh / H) ** (-1 / Hprime)); return rho * v2 / (2 * beta_b), r, math.sqrt(v2)
    fa = lambda nu: (lambda A, r, v: -(2 * a * a * v / mu) * A * r * r / h)(*acc(nu)); fe = lambda nu: (lambda A, r, v: -(2 / v) * (e + math.cos(nu)) * A * r * r / h)(*acc(nu))
    w = min(80 * math.sqrt(H / (rp * e)), math.pi) if e > 0 else math.pi
    da = 2 * quad(fa, 0, w, limit=800, epsabs=0, epsrel=1e-12)[0]; de = 2 * quad(fe, 0, w, limit=800, epsabs=0, epsrel=1e-12)[0]
    return da, de


# ---- selftest ----------------------------------------------------------------------------------------------------------------------------
FAILS = []
def check(name, got, want, rtol, note=""):
    ok = abs(got - want) <= rtol * abs(want); print(f"{'OK ' if ok else 'BAD'} {name}: {got:.6g} vs {want:.6g} (rel {abs(got / want - 1):.2e}, tol {rtol:g}) {note}")
    if not ok: FAILS.append(name)


def selftest():
    out = {}
    # D1 synthetic: Gauss quadrature of a drag pass whose density is exp(-(h(f) - h_p)/H + kappa f), oblate surface, non-rotating atmosphere
    mu, H, beta = MU_MARS, 6500.0, 20.0; rows = []
    for e, om, kap in ((0.78, 1.9, 3.0), (0.6, 1.6, 2.0), (0.4, 1.1, 0.0), (0.4, 1.1, -2.0), (0.35, 1.0, 3.0)):
        rp = 3.50e6; a = rp / (1 - e); inc = 1.632; p = a * (1 - e * e); hm = math.sqrt(mu * p)
        def dens(f):
            r = p / (1 + e * math.cos(f)); return math.exp(-((r - L.surface_radius(inc, om + f, "mars")) - (rp - L.surface_radius(inc, om, "mars"))) / H + kap * f)
        def rates(f):
            r = p / (1 + e * math.cos(f)); V = math.sqrt(mu * (2 / r - 1 / a)); D = 1e-9 * dens(f) * V * V / (2 * beta); dt = r * r / hm
            return -2 * D * math.sin(f) / (e * V) * dt, -2 * D * (e + math.cos(f)) / V * dt
        w = 0.6; dom = quad(lambda f: rates(f)[0], -w, w, limit=800, epsrel=1e-11, points=[0])[0]; de = quad(lambda f: rates(f)[1], -w, w, limit=800, epsrel=1e-11, points=[0])[0]
        pred = apse_rotation_drag(de, e, rp, inc, om, H, kap); rows.append((e, om, kap, dom, pred))
        scale = abs(de) * (abs(surface_slope(inc, om)) + H * abs(kap)) / (rp * e * e)    # size of the larger of the two terms (they can cancel)
        ok = abs(pred - dom) <= 0.05 * scale; print(f"{'OK ' if ok else 'BAD'} D1 synthetic apse rotation e={e} omega={om} kappa={kap}: {pred:.4e} vs {dom:.4e}, error/term scale {abs(pred - dom) / scale:.2e} (tol 0.05)")
        if not ok: FAILS.append("D1 synthetic")
    out["D1_synthetic"] = rows
    # D2: Laplace with H' -> factor 1 + 3H'/8 (large c), and the OU event factor limits
    for Hp in (0.03, 0.1):
        a, e = 9.4e6, 0.63; num = gauss_quad_layer(a, e, 6500.0, 1e-9, 20.0, mu, Hp)[0] / gauss_quad_layer(a, e, 6500.0, 1e-9, 20.0, mu, 0.0)[0]
        check(f"D2/D4 Laplace dH/dh factor 1 + 3H'/8, H'={Hp}", 1 + 3 * Hp / 8, num, 3e-3 if Hp < 0.05 else 1.5e-2, "(first order in H')")
    check("D2 OU event factor, tau_c -> 0 limit is exp(s^2/2)", pert_factor_event(50.0, 0.3, 1e-6), math.exp(0.045), 1e-3)
    check("D2 OU event factor, tau_c -> inf limit is 1", pert_factor_event(50.0, 0.3, 1e9), 1.0, 1e-6)
    # D3: octupole vs direct quadrature (frozen perturber), error must fall as (a/r_b)^2
    mu_e = 3.98600436233e14; rb0 = np.array([2.0e8, 2.5e8, 0.8e8]); errs = []
    for sc in (3.0, 1.0):
        rb = rb0 * sc; mub = 4.9028e15 * sc ** 3 / 1000; a, e = 2.4e7, 0.8; inc, raan, om = 0.8, 0.4, 1.1; Q = TL.pqw_matrix(inc, raan, om)
        k = L.tide_arc_direct(mu_e, mub, a, e, inc, raan, om, lambda tau: np.repeat(rb[:, None], np.size(tau), axis=1), 128)["rp"]
        q = L.tide_rp_averaged(a, e, mu_e, mub, rb, Q[:, 0], Q[:, 1]); o = tide_rp_octupole(a, e, mu_e, mub, rb, Q[:, 0], Q[:, 1])
        errs.append((a / np.linalg.norm(rb), q / k - 1, (q + o) / k - 1)); print(f"info D3 a/r_b={errs[-1][0]:.3f}: quadrupole {errs[-1][1]:+.4f}, quadrupole+octupole {errs[-1][2]:+.5f}")
    check("D3 a/r_b=0.073 e=0.8: (quad+oct)/direct - 1 within 2.5%", 1 + errs[1][2], 1.0, 0.025)
    check("D3 residual scales as (a/r_b)^2", errs[1][2] / errs[0][2], (errs[1][0] / errs[0][0]) ** 2, 0.25)
    out["D3_synthetic"] = errs
    # D4: all-e Bessel vs quadrature, including e beyond the e-series range; and the e-series at e=0.78 for reference
    worst = 0.0
    for e in (0.01, 0.3, 0.78, 0.9, 0.95, 0.99):
        for c in (2, 25, 300, 3000):
            a = c * 6308.278108 / e
            if not (3.0e6 < a < 2e9): continue
            da, de, _ = bessel_drag_pass_alle(a, e, 6308.278108, 1e-9, 20.24, MU_MARS); qa, qe = gauss_quad_layer(a, e, 6308.278108, 1e-9, 20.24, MU_MARS)
            worst = max(worst, abs(da / qa - 1), abs(de / qe - 1))
    check("D4 all-e Bessel vs quadrature, worst of da, de (e 0.01..0.99, c 2..3000)", 1 + worst, 1.0, 1e-6)
    for Hp in (0.03, 0.1):
        worstH = 0.0
        for e, c in ((0.78, 25), (0.78, 1000), (0.36, 300), (0.05, 10), (0.95, 5000)):
            a = c * 6500.0 / e; d1 = bessel_drag_pass_alle(a, e, 6500.0, 1e-9, 20.0, MU_MARS, Hp); d0 = gauss_quad_layer(a, e, 6500.0, 1e-9, 20.0, MU_MARS, Hp)
            worstH = max(worstH, abs(d1[0] / d0[0] - 1), abs(d1[1] / d0[1] - 1))
        check(f"D4 Bessel with dH/dh vs exact layer quadrature, H'={Hp}", 1 + worstH, 1.0, 0.25 * Hp * Hp + 2e-4, "(error O(H'^2))")
    out["D4_worst_alle"] = worst; return out


# ---- data checks ----------------------------------------------------------------------------------------------------------------------------
def _elems(r, v, mu):
    rn = np.linalg.norm(r, axis=1); h = np.cross(r, v); hn = np.linalg.norm(h, axis=1); ev = np.cross(v, h) / mu - r / rn[:, None]; e = np.linalg.norm(ev, axis=1)
    p = hn ** 2 / mu; inc = np.arccos(h[:, 2] / hn); nv = np.stack([-h[:, 1], h[:, 0], 0 * hn], 1); nn = np.linalg.norm(nv, axis=1)
    om = np.arccos(np.clip((nv * ev).sum(1) / (nn * e), -1, 1)); om = np.where(ev[:, 2] < 0, 2 * np.pi - om, om)
    f = np.arccos(np.clip((ev * r).sum(1) / (e * rn), -1, 1)); f = np.where((r * v).sum(1) < 0, -f, f)
    return e, inc, om, f, p, hn, rn, h


def twin_passes(tasks, runs):
    """Per logged pass of the half-lever twin flights (SOURCED halflever_twin_*/<run>/steps.csv, passes.csv): Gauss-integrated drag d omega and d e (perturbed and
    mean-field), fitted kappa (ln rho_mean on altitude, altitude^2 and f within 3 H of the minimum), drag-work d a, perturbation factors."""
    rows = []
    for run in runs:
        s = pd.read_csv(tasks / run / "steps.csv"); ps = pd.read_csv(tasks / run / "passes.csv"); br = np.r_[0, np.where(s.t_s.diff() > 100)[0], len(s)]
        for j in range(min(len(br) - 1, len(ps) - 1)):
            g = s.iloc[br[j]:br[j + 1]]; r = g[["rx", "ry", "rz"]].values; v = g[["vx", "vy", "vz"]].values; ad = g[["adx", "ady", "adz"]].values; rp_ = g.r_pert.values
            e, inc, om, f, p, hn, rn, h = _elems(r, v, MU_MARS); rh = r / rn[:, None]; hh = h / hn[:, None]; sh = np.cross(hh, rh); t = g.t_s.values
            R, S_, Wc = (ad * rh).sum(1), (ad * sh).sum(1), (ad * hh).sum(1); u = om + f
            dom_in = (1 / (hn * e)) * (-p * np.cos(f) * R + (p + rn) * np.sin(f) * S_); dom_W = -(rn * np.sin(u) * np.cos(inc) / (hn * np.sin(inc))) * Wc
            de_r = (1 / hn) * (p * np.sin(f) * R + ((p + rn) * np.cos(f) + rn * e) * S_); I = lambda y: float(np.trapz(y, t))
            cum = np.r_[0, np.cumsum(0.5 * (dom_in[1:] + dom_in[:-1] + dom_W[1:] + dom_W[:-1]) * np.diff(t))]; k0 = int(np.argmax(f >= 0)); w0 = -f[k0 - 1] / (f[k0] - f[k0 - 1])
            first_half = cum[k0 - 1] + w0 * (cum[k0] - cum[k0 - 1])
            alt = g.alt_m.values; hp = alt.min(); m = (alt - hp) < 3 * 6500; lr = np.log(g.rho_mean.values)
            cf = np.linalg.lstsq(np.c_[np.ones(m.sum()), alt[m] - hp, (alt[m] - hp) ** 2, f[m]], lr[m], rcond=None)[0]; Hf = -1 / cf[1]
            k = int(np.argmin(np.abs(f))); x = ps.iloc[j]
            work = I((ad * v).sum(1)); work_m = I((ad * v).sum(1) / rp_)
            rows.append(dict(run=run, j=j, e=e[k], inc=inc[k], om=om[k], rp=rn[k], kappa=cf[3], H_fit=Hf, Hprime=2 * cf[2] * Hf * Hf,
                             dom=I(dom_in + dom_W), dom_W=I(dom_W), dom_mean=I((dom_in + dom_W) / rp_), de=I(de_r), de_mean=I(de_r / rp_), dom_first=first_half, dom_second=I(dom_in + dom_W) - first_half,
                             da_obs=2 * x.a_m ** 2 * work / MU_MARS, da_mean=2 * x.a_m ** 2 * work_m / MU_MARS, t_ev=x.time_s, r0=x.r_pert, a=x.a_m, rp_ev=x.radius_m, e_ev=x.e,
                             q_mean=x.q_w_cm2, q_eff=x.q_eff_w_cm2, H_loc=x.h_local_m, mass=x.mass_kg, vrel=x.vrel_mps, lat=x.lat_rad))
    return pd.DataFrame(rows)


def data_checks(tasks):
    import analytic_pairs as PP, stagek_lib as SK, analytic_ephem as E
    from scipy.interpolate import CubicSpline
    R = {}; H0 = MARS["H"]; late = ["halflever_twin_late/s1_none", "halflever_twin_late/s2_none"]; early = [f"halflever_twin_data/s{i}_none" for i in range(1, 5)]
    tw = twin_passes(tasks, late + early); P = SK.planet("Mars")
    # D1a: kinematic identity and the f* decomposition on the logged passes (mean-field drag, so the perturbation does not enter)
    tw["rpp"] = tw.rp * tw.e / (1 + tw.e); tw["Rs1"] = [surface_slope(i, o) for i, o in zip(tw.inc, tw.om)]
    tw["pred_geo"] = tw.de_mean * tw.Rs1 / (tw.rp * tw.e ** 2); tw["pred_full"] = tw.de_mean * (tw.Rs1 + tw.H_fit * tw.kappa) / (tw.rp * tw.e ** 2)
    sl = lambda x, y: float((x * y).sum() / (x * x).sum())
    R["D1_twin_mean_field"] = dict(n=int(len(tw)), sum_dom_gauss=float(tw.dom_mean.sum()), ratio_gauss_over_geometry_only=float(tw.dom_mean.sum() / tw.pred_geo.sum()),
                                   ratio_gauss_over_geometry_plus_kappa=float(tw.dom_mean.sum() / tw.pred_full.sum()), slope_full=sl(tw.pred_full, tw.dom_mean), corr_full=float(np.corrcoef(tw.pred_full, tw.dom_mean)[0, 1]),
                                   geometry_share_of_f_star=float((tw.Rs1 / (tw.Rs1 + tw.H_fit * tw.kappa)).median()), W_term_share=float(tw.dom_W.sum() / tw.dom.sum()))
    # D1b: event-to-event omega residual of the late twins (no burns) after numeric zonal map + tesseral + Sun, against the half-pass drag integrals
    cs = TL.coefficients("Mars", 20); Ct, St = TL.select(cs["C"], cs["S"], lambda n, m: m >= 1); et0 = E.utc_to_et(2001, 11, 6, 19, 0, 32); clos = {}
    for run in late:
        ps = pd.read_csv(tasks / run / "passes.csv"); omu = np.unwrap(ps.omega_rad.values); g = tw[tw.run == run].reset_index(drop=True); res_s = []; res_n = []
        for k in range(len(g) - 1):
            r = ps.iloc[k]; th = L.theta_p_from_pass(r.raan_rad, r.omega_rad, r.i_rad, r.lon_rad)
            tk = L.tess_arc(cs, Ct, St, r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, th, MARS["w_rot"], 128)["omega"]
            zs = L.zonal_rates_per_pass(P, r.a_m, r.e, r.i_rad, r.omega_rad)[1]; zn = SK.zonal_map(P, r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, nmax=20, rtol=1e-10, atol=1e-6)[4] - r.omega_rad
            Tp = ps.time_s[k + 1] - r.time_s; ts = np.linspace(0, Tp, 5); pos = np.array([E.third_body_in_frame("mars", 10, et0 + r.time_s + x) for x in ts]).T; cs_ = CubicSpline(ts, pos, axis=1)
            sun = L.tide_arc_direct(MU_MARS, E.GM[10], r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, lambda tau: cs_(np.minimum(tau, Tp)), 32)["omega"]
            obs = omu[k + 1] - omu[k]; res_s.append(obs - zs - tk); res_n.append(obs - zn - tk - sun)
        drag = (g.dom_second.values[:-1] + g.dom_first.values[1:])
        clos[run] = dict(sum_resid_secular_zonal=float(np.sum(res_s)), sum_resid_numeric_zonal_sun=float(np.sum(res_n)), sum_drag_gauss=float(drag.sum()),
                         drag_over_resid=float(drag.sum() / np.sum(res_n)), corr_per_pair=float(np.corrcoef(drag, res_n)[0, 1]), mean_numeric_minus_secular_zonal=float(np.mean(np.array(res_s) - np.array(res_n))))
    R["D1_twin_event_closure"] = clos
    # D1c: Stage O T2 cumulative (SOURCED stageo_data/dev/runs; pair cache of analytic_pairs)
    aug = PP.CACHE / "derivations_stageo_aug.csv"
    if aug.exists(): n = pd.read_csv(aug)
    else:
        d = PP.stageo_pairs(); n = d[~d.burn].copy(); zn = []; sun = []
        for r in n.itertuples():
            zn.append((SK.zonal_map(P, r.a, r.e, r.inc, r.raan, r.om, nmax=20, rtol=1e-10, atol=1e-6)[4] - r.om + math.pi) % (2 * math.pi) - math.pi)
            ts = np.linspace(0, r.dt, 5); pos = np.array([E.third_body_in_frame("mars", 10, et0 + r.t + x) for x in ts]).T; cs_ = CubicSpline(ts, pos, axis=1)
            sun.append(L.tide_arc_direct(MU_MARS, E.GM[10], r.a, r.e, r.inc, r.raan, r.om, lambda tau, cs_=cs_, T=r.dt: cs_(np.minimum(tau, T)), 32)["omega"])
        n["znum"] = zn; n["sun_dom"] = sun
        da = np.array([L.drag_drop_ra_m(r.a, r.e, r.rp, r.q_eff * 1e4, H0, r.beta, PP.MU_M, r.frot) / 2 for r in n.itertuples()])
        drp = np.array([L.drag_drop_rp_m(r.a, r.e, r.rp, r.q_eff * 1e4, H0, r.beta, PP.MU_M, r.frot) for r in n.itertuples()]); n["de_drag"] = ((1 - n.e) * da - drp) / n.a
        n.to_csv(aug, index=False)
    n["held"] = n.seed > 6; n["Rs1"] = [surface_slope(i, o) for i, o in zip(n.inc, n.om)]; n["rpp"] = n.rp * n.e / (1 + n.e)
    n["dom_geo"] = n.de_drag * n.Rs1 / (n.rp * n.e ** 2); runs = {d.name: pd.read_csv(d / "passes.csv") for d, _, _ in PP.stageo_runs()}
    def cum(dom):
        g = np.array([L.dh_geom(r.inc, r.om, x, "mars", r.t_di) for r, x in zip(n.itertuples(), dom)])
        c = (pd.Series(g, index=n.index).groupby(n.run).sum() - n.groupby("run").o_dg.sum()); h = n.groupby("run").held.first()
        return dict(all_median=float(c.abs().median()), all_max=float(c.abs().max()), dev_max=float(c.abs()[~h].max()), held_median=float(c.abs()[h].median()), held_max=float(c.abs()[h].max()), signed_median=float(c.median()))
    base = n.z_dom + n.t_dom; m1 = n.znum + n.t_dom + n.sun_dom; T2 = dict(M0_secular_zonal_tesseral=cum(base), M1_numeric_zonal_tesseral_sun=cum(m1), M2_plus_geometry=cum(m1 + n.dom_geo))
    for W in (10, 20, 30, 40):
        ks = pd.concat([pd.DataFrame(dict(run=r, k=ps["pass"].values, ks=kappa_trailing(ps.omega_rad.values, ps.rho.values, ps.alt_m.values, ps.h_local_m.values, W))) for r, ps in runs.items()])
        kap = n[["run", "k"]].merge(ks, on=["run", "k"], how="left").ks.fillna(0.0).values; n[f"kap{W}"] = kap
        T2[f"M3_plus_geometry_plus_kappa_W{W}"] = cum(m1 + n.dom_geo + n.de_drag * n.h_local * kap / (n.rp * n.e ** 2))
    Wbest = min((10, 20, 30, 40), key=lambda W: T2[f"M3_plus_geometry_plus_kappa_W{W}"]["dev_max"]); T2["W_selected_on_dev_max"] = Wbest; R["D1_stageO_T2_cumulative"] = T2
    # D2: decomposition on the twin passes and the Stage O prediction by a-bin
    tw["beta"] = tw.mass / (PP.CD * PP.AREA); tw["frot"] = [AP.f_rot(r.rp_ev, r.a, r.vrel, r.lat, MARS["w_rot"], MU_MARS) for r in tw.itertuples()]
    tw["lap_mean"] = [L.drag_drop_ra_m(r.a, r.e_ev, r.rp_ev, r.q_mean * 1e4, H0, r.beta, MU_MARS, r.frot) / 2 for r in tw.itertuples()]
    tw["B2"] = tw.da_mean / tw.lap_mean; tw["B3"] = tw.da_obs / tw.da_mean; tw["A_eff"] = tw.da_obs / (tw.lap_mean * tw.q_eff / tw.q_mean)
    tw["fstar"] = (tw.Rs1 + tw.H_fit * tw.kappa) / tw.rpp; tw["B2_pred"] = [profile_factor(r.e, r.rp, r.H_fit, H0, r.fstar, r.Hprime) for r in tw.itertuples()]
    X = [];
    for run in late + early:
        s = pd.read_csv(tasks / run / "steps.csv"); ps = pd.read_csv(tasks / run / "passes.csv"); br = np.r_[0, np.where(s.t_s.diff() > 100)[0], len(s)]
        for j in range(min(len(br) - 1, len(ps))):
            g = s.iloc[br[j]:br[j + 1]]; tt = np.arange(-150, 151, 1.0); X.append(np.interp(tt, g.t_s.values - ps.time_s.iloc[j], np.log(g.r_pert.values)))
    X = np.array(X); mX, vX = X.mean(), X.var(); acf = {L_: float(np.mean((X[:, L_:] - mX) * (X[:, :-L_] - mX)) / vX) for L_ in (5, 10, 20, 30, 60)}
    tau_c = float(np.median([-L_ / math.log(r) for L_, r in acf.items() if L_ <= 30])); s_near = float(X[:, 100:201].std())
    sig_t = tw.rp_ev * np.sqrt(tw.H_loc / (MU_MARS * tw.e_ev)); Fpe = np.array([pert_factor_event(x, s_near, tau_c) for x in sig_t])
    R["D2_twin"] = dict(A_eff_median=float(tw.A_eff.median()), B2_profile_median=float(tw.B2.median()), B3_pert_median=float(tw.B3.median()), B3_pert_mean=float(tw.B3.mean()),
                        B3_over_event_median=float((tw.B3 / tw.r0).median()), B2_pred_median=float(tw.B2_pred.median()), B2_over_pred_median=float((tw.B2 / tw.B2_pred).median()),
                        B2_corr=float(np.corrcoef(tw.B2, tw.B2_pred)[0, 1]), shift_factor_median=float(np.exp(tw.rpp * tw.fstar ** 2 / (2 * tw.H_fit)).median()),
                        Hprime_median=float(tw.Hprime.median()), H_fit_median=float(tw.H_fit.median()), ln_r_mean=float(mX), ln_r_sd_within_100s=s_near, acf=acf, tau_c_s=tau_c,
                        F_pert_mean_pred=pert_factor_mean(float(X[:, 100:201].mean()), s_near), F_pert_event_pred_median=float(np.median(Fpe)))
    Hp = float(tw.Hprime.median()); W = Wbest; k = n[f"kap{W}"].values
    pm = np.array([L.drag_drop_ra_m(r.a, r.e, r.rp, r.q_mean * 1e4, H0, r.beta, PP.MU_M, r.frot) for r in n.itertuples()])
    n["A_mean"] = 2 * n.o_da / pm; n["A_eff"] = n.A_mean * n.q_mean / n.q_eff; fs = (n.Rs1 + n.h_local * k) / n.rpp
    n["B2p"] = [profile_factor(r.e, r.rp, r.h_local, H0, f, Hp) for r, f in zip(n.itertuples(), fs)]
    n["Fpe"] = [pert_factor_event(x, s_near, tau_c) for x in n.rp * np.sqrt(n.h_local / (PP.MU_M * n.e))]
    bins = pd.cut(n.a, [3e6, 1e7, 1.3e7, 1.6e7, 3e7]); Fm = R["D2_twin"]["F_pert_mean_pred"]
    R["D2_stageO_by_a_bin"] = {str(b): dict(n=int(len(g)), obs_over_laplace_qmean=float(g.A_mean.median()), obs_over_laplace_qeff=float(g.A_eff.median()),
                                            pred_qmean=float((g.B2p * Fm).median()), pred_qeff=float((g.B2p * g.Fpe).median()), obs_over_pred_qmean=float((g.A_mean / (g.B2p * Fm)).median()),
                                            obs_over_pred_qeff=float((g.A_eff / (g.B2p * g.Fpe)).median())) for b, g in n.groupby(bins, observed=True)}
    # D3: octupole on the Stage K Earth and Venus pairs (cached tide tables; perturber at the arc midpoint, as the averaged form)
    EVf = lambda y, p: float(1 - ((y - p) ** 2).sum() / (y ** 2).sum()); oc = {}
    for pl in ("earth", "venus"):
        d = PP.tide_pairs(pl); Pp = SK.planet(pl.capitalize()); mu = Pp["mu"]; e0 = {ep: E.utc_to_et(*PP.EPOCHS[pl][ep]) for ep in PP.EPOCHS[pl]}; o = []
        for r in d.itertuples():
            T = 2 * math.pi * math.sqrt(r.a ** 3 / mu); Q = TL.pqw_matrix(r.i, r.raan, r.w)
            o.append(sum(tide_rp_octupole(r.a, r.e, mu, E.GM[b], E.third_body_in_frame(pl, b, e0[int(r.epoch)] + r.t0 + 0.5 * T), Q[:, 0], Q[:, 1]) for b in PP.BODIES[pl]))
        d["oct"] = o; m = d.split == "heldout"
        oc[pl] = dict(EV_held_direct=EVf(d.dr[m], (d.dr_M + d.dr_T8)[m]), EV_held_quadrupole=EVf(d.dr[m], (d.dr_M + d.dr_T8avg)[m]), EV_held_quad_oct=EVf(d.dr[m], (d.dr_M + d.dr_T8avg + d.oct)[m]),
                      rms_quad_minus_direct_m=float(np.sqrt(((d.dr_T8avg - d.dr_T8) ** 2).mean())), rms_quad_oct_minus_direct_m=float(np.sqrt(((d.dr_T8avg + d.oct - d.dr_T8) ** 2).mean())))
    R["D3_stageK"] = oc
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--data", action="store_true"); ap.add_argument("--out", default=None)
    A = ap.parse_args(); res = {"selftest": selftest()}
    if A.data:
        import analytic_pairs as PP
        res["data"] = data_checks(PP.TASKS); out = pathlib.Path(A.out) if A.out else PP.TASKS / "analytic_data" / "derivations_results.json"
        out.write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res["data"], indent=1, default=float)); print("wrote", out)
    print("FAILURES:", FAILS if FAILS else "none"); sys.exit(1 if FAILS else 0)
