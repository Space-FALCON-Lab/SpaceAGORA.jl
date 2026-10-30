"""Mechanism and formula checks against the vectorized Cowell runs (tesseral_cowell.py --set mech / drag).

Usage: python3 -I tesseral_analyze.py --set mech|drag [--results-root DIR]
mech writes <root>/tesseral/{mech_summary.csv, mech_terms.csv, mech_truncation.csv, mech_localization.csv, mech_passes.csv}:
  effect      Cowell variant minus the zonal-only run, osculating elements at each apoapsis (k = 0..40)
  P0          orbit-averaged, non-resonant: zero for every m != 0
  P_slow      the slow (near-resonant) Kaula terms only: for each (n, m) the single j nearest m*nu (series_kick, only_j)
  P_map       per-pass map: Gauss-quadrature kick of the variant's tesseral field along the Kepler arc of the zonal reference's apoapsis
              elements, with the rotation angle at the reference's periapsis; cumulative sum (first order)
  P_ser(J)    the same map from the Lagrange series, keeping the 2J+1 values of j nearest m*nu (truncation study)
drag writes <root>/tesseral/drag_passes.csv (map along the zonal reference with drag, nu per pass, resonance crossings).
"""
import argparse, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import tesseral_lib as tl
from tesseral_cowell import W_ROT

KEYS = ["a_m", "e", "i_rad", "raan_rad", "omega_rad", "rp_m"]; KK = dict(zip(KEYS, tl.KICK_NAMES))


def rms(x): x = np.asarray(x, float); return float(np.sqrt(np.mean(x ** 2))) if len(x) else float("nan")


def nu_eff(ref):
    """Commensurability nu = (w_rot T_anom - dOmega)/(2 pi) per pass from the zonal reference (T_anom = periapsis-to-periapsis)."""
    tp = ref.t_peri_s.values; Om = np.unwrap(ref.raan_rad.values)
    T = np.diff(tp); dO = np.diff(Om)[:len(T)]
    nu = (W_ROT * T - dO) / (2 * math.pi); ok = np.isfinite(nu)
    return np.where(ok, nu, np.interp(np.arange(len(nu)), np.nonzero(ok)[0], nu[ok]))


def map_prediction(ref, cs, C, S, nf=192, fmask=None):
    """Cumulative first-order kicks along the reference (rows k = 0..K-1 used; returns array (K+1, 6) with row 0 = 0)."""
    out = [np.zeros(6)]
    for r in ref.itertuples():
        if not np.isfinite(r.t_peri_s): break
        k = tl.gauss_kick(r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, r.theta_peri, W_ROT, cs, C, S, nf=nf, fmask=fmask)
        out.append(out[-1] + np.array([k[n] for n in tl.KICK_NAMES]))
    return np.array(out)


def mech(root):
    d = pd.read_csv(root / "tesseral" / "cowell_mech.csv"); cs = tl.coefficients("Mars", 20); L = 20
    rows, terms, trunc, loc, passes = [], [], [], [], []
    for o, g in d.groupby("orbit"):
        ref = g[g.variant == "zonal"].sort_values("pass").reset_index(drop=True); K = len(ref) - 1
        nu = nu_eff(ref); nu0 = float(np.mean(nu)); inc = ref.inc_deg.iloc[0]; om0 = ref.omega_deg.iloc[0]
        e0, i0 = ref.e.iloc[0], ref.i_rad.iloc[0]
        eff = {v: (gv.sort_values("pass")[KEYS].values - ref[KEYS].values) for v, gv in g.groupby("variant") if v != "zonal"}
        # map predictions per variant (tesseral part only)
        Cz, Sz = tl.select(cs["C"], cs["S"], lambda n, m: m == 0)
        def tess(v):
            if v == "full": keep = lambda n, m: m > 0
            elif v.startswith("m"): mm = int(v[1:]); keep = lambda n, m: m == mm
            else: nn = int(v[1:]); keep = lambda n, m: n == nn and m > 0
            return tl.select(cs["C"], cs["S"], keep), keep
        refk = ref.iloc[:K]
        for v in eff:
            (C, S), keep = tess(v); P = map_prediction(refk, cs, C, S)
            E = eff[v][:K + 1]
            rp_eff, rp_map = E[1:, 5], P[1:, 5]
            terms.append(dict(orbit=o, inc_deg=inc, omega_deg=om0, variant=v, nu=nu0, rp_eff_rms_m=rms(rp_eff), rp_eff_max_m=float(np.max(np.abs(rp_eff))),
                              rp_map_err_rms_m=rms(rp_map - rp_eff), frac_map=rms(rp_map - rp_eff) / max(rms(rp_eff), 1e-12),
                              **{f"{KK[key]}_frac_map": rms(P[1:, j] - E[1:, j]) / max(rms(E[1:, j]), 1e-300) for j, key in enumerate(KEYS)},
                              **{f"{KK[key]}_eff_rms": rms(E[1:, j]) for j, key in enumerate(KEYS)}))
            if v == "full":
                for k in range(K + 1):
                    passes.append(dict(orbit=o, inc_deg=inc, k=k, rp_eff_m=E[k, 5], rp_map_m=P[k, 5], i_eff=E[k, 2], i_map=P[k, 2], a_eff=E[k, 0], a_map=P[k, 0],
                                       theta_peri=ref.theta_peri.iloc[k] if k < K else np.nan))
        # slow-term-only and truncated-series predictions for the full tesseral field (frozen e, i; per-pass a, raan, omega, theta)
        E = eff["full"][:K + 1]
        series = {(n, m): tl.SeriesTerm(n, m, e0, i0, NM=4096) for n in range(2, L + 1) for m in range(1, n + 1)}
        Js = [0, 1, 2, 5, 10, 20, 40, None]; preds = {J: np.zeros((K + 1, 6)) for J in Js}
        for k, r in enumerate(refk.itertuples()):
            nuk = nu[min(k, len(nu) - 1)]
            for J in Js:
                tot = np.zeros(6)
                for (n, m), T in series.items():
                    if J is None: sk = tl.series_kick(T, r.a_m, r.raan_rad, r.omega_rad, r.theta_peri, nuk, cs)
                    else:
                        jc = np.round(m * nuk); js = np.arange(jc - J - 1, jc + J + 2); js = js[np.argsort(np.abs(js - m * nuk))][:2 * J + 1]
                        sk = tl.series_kick(T, r.a_m, r.raan_rad, r.omega_rad, r.theta_peri, nuk, cs, only_j=js)
                    tot += np.array([sk[x] for x in tl.KICK_NAMES])
                preds[J][k + 1] = preds[J][k] + tot
        for J in Js:
            trunc.append(dict(orbit=o, inc_deg=inc, omega_deg=om0, J=("all" if J is None else J), rp_err_rms_m=rms(preds[J][1:, 5] - E[1:, 5]),
                              rp_frac=rms(preds[J][1:, 5] - E[1:, 5]) / rms(E[1:, 5]), i_frac=rms(preds[J][1:, 2] - E[1:, 2]) / rms(E[1:, 2]),
                              a_frac=rms(preds[J][1:, 0] - E[1:, 0]) / rms(E[1:, 0])))
        # localization: fraction of the predicted per-pass rp kick collected within |f| < fc
        (C, S), _ = tess("full"); full = map_prediction(refk, cs, C, S)[1:, 5]
        for fc in (30, 45, 60, 90, 120):
            part = map_prediction(refk, cs, C, S, fmask=lambda f, fc=fc: (np.abs(f) < math.radians(fc)).astype(float))[1:, 5]
            loc.append(dict(orbit=o, inc_deg=inc, fc_deg=fc, frac_err=rms(part - full) / rms(full)))
        kicks = np.diff(np.r_[0.0, full])
        rows.append(dict(orbit=o, inc_deg=inc, raan_deg=ref.raan_deg.iloc[0], omega_deg=om0, nu=nu0, nu_min=float(nu.min()), nu_max=float(nu.max()),
                         rp_eff_rms_m=rms(E[1:, 5]), rp_eff_max_m=float(np.max(np.abs(E[1:, 5]))), kick_rms_m=rms(kicks),
                         P0_frac=1.0, slow_frac=trunc[-len(Js)]["rp_frac"], map_frac=[t for t in terms if t["orbit"] == o and t["variant"] == "full"][0]["frac_map"]))
        print(f"orbit {o} inc {inc} nu {nu0:.4f}: rp effect rms {rows[-1]['rp_eff_rms_m']:.0f} m, kick rms {rows[-1]['kick_rms_m']:.0f} m, map err frac {rows[-1]['map_frac']:.3f}, slow-only err frac {rows[-1]['slow_frac']:.3f}", flush=True)
    out = root / "tesseral"
    for name, x in (("mech_summary", rows), ("mech_terms", terms), ("mech_truncation", trunc), ("mech_localization", loc), ("mech_passes", passes)):
        pd.DataFrame(x).to_csv(out / f"{name}.csv", index=False)


def drag(root):
    d = pd.read_csv(root / "tesseral" / "cowell_drag.csv"); cs = tl.coefficients("Mars", 20); rows = []
    (C, S) = tl.select(cs["C"], cs["S"], lambda n, m: m > 0)
    for o, g in d.groupby("orbit"):
        ref = g[g.variant == "zonal"].sort_values("pass").reset_index(drop=True); full = g[g.variant == "full"].sort_values("pass").reset_index(drop=True)
        K = min(len(ref), len(full)) - 1; ref = ref.iloc[:K + 1]; full = full.iloc[:K + 1]
        nu = nu_eff(ref); P = map_prediction(ref.iloc[:K], cs, C, S); Pf = map_prediction(full.iloc[:K], cs, C, S)
        for k in range(K + 1):
            rows.append(dict(orbit=o, inc_deg=ref.inc_deg.iloc[0], k=k, t_d=ref.t_apo_s.iloc[k] / 86400, a_ref_m=ref.a_m.iloc[k], nu=nu[k] if k < len(nu) else np.nan,
                             rp_eff_m=full.rp_m.iloc[k] - ref.rp_m.iloc[k], rp_map_ref_m=P[k, 5] if k < len(P) else np.nan, rp_map_full_m=Pf[k, 5] if k < len(Pf) else np.nan))
    pd.DataFrame(rows).to_csv(root / "tesseral" / "drag_passes.csv", index=False); print("wrote drag_passes.csv", len(rows))


def loc(root):
    """Where along the arc the per-pass r_p kick is collected: unit-coefficient kicks Y_nm restricted to |f| < fc against the full arc,
    pass-1 elements (A0, E0), rms over m for each degree; and the share of the arc time spent there."""
    from tesseral_cowell import A0, E0
    cs = tl.coefficients("Mars", 20); rows = []
    for inc in (25.2, 63.4, 93.522):
        for om in (0.0, 90.0):
            Yf = tl.term_kicks(A0, E0, math.radians(inc), 0.0, math.radians(om), 0.0, W_ROT, cs["mu"], cs["R"], 20, nf=512)
            for fc in (10, 20, 30, 45, 60, 90, 120, 150):
                Yr = tl.term_kicks(A0, E0, math.radians(inc), 0.0, math.radians(om), 0.0, W_ROT, cs["mu"], cs["R"], 20, nf=512,
                                   fmask=lambda f, fc=fc: (np.abs(f) < math.radians(fc)).astype(float))
                for n in (2, 3, 4, 6, 10, 15, 20):
                    rows.append(dict(inc_deg=inc, omega_deg=om, fc_deg=fc, n=n, frac_missing=float(np.sqrt((np.abs(Yr[n, 1:n + 1] - Yf[n, 1:n + 1]) ** 2).sum() / (np.abs(Yf[n, 1:n + 1]) ** 2).sum()))))
    df = pd.DataFrame(rows); df.to_csv(root / "tesseral" / "localization_by_degree.csv", index=False)
    print(df.pivot_table(index="fc_deg", columns="n", values="frac_missing", aggfunc="median").round(3).to_string())
    p = A0 * (1 - E0 ** 2)
    for fc in (30, 60, 90, 120):
        f = math.radians(fc); E = 2 * math.atan(math.sqrt((1 - E0) / (1 + E0)) * math.tan(f / 2)); print(f"|f| < {fc} deg: r <= {p / (1 + E0 * math.cos(f)) / cs['R']:.2f} R, time share {(E - E0 * math.sin(E)) / math.pi:.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--set", required=True, choices=["mech", "drag", "loc"])
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); {"mech": mech, "drag": drag, "loc": loc}[a.set](root)
