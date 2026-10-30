"""Feasibility test of a per-pass (apoapsis-to-apoapsis) map that carries the tesseral kick, against the drag Cowell set.

Usage: python3 -I tesseral_propagator.py [--results-root DIR] [--nf 4097]
The map is a minimal stand-in for an estimator step (it is not estimator v3): from the osculating elements at apoapsis k it integrates
the Gauss equations of drag (same atmosphere and beta as tesseral_cowell.py) and of the chosen gravity field (zonal 20x0 or full 20x20)
along the Kepler arc with the planet rotating, to first order (elements frozen on the arc; cumulative trapezoid in true anomaly), and
times the next apoapsis from the mean-anomaly equation. Runs, per orbit of cowell_drag.csv:
  zonal map vs Cowell zonal   (the map's own error without tesserals)
  full map  vs Cowell full    (the map with the tesseral kick: what an estimator could reach)
  zonal map vs Cowell full    (the omitted-tesseral error, the v3 situation)
Writes <root>/tesseral/{propagator_passes.csv, propagator_bins.csv}.
"""
import argparse, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import tesseral_lib as tl
from tesseral_cowell import W_ROT, RHO_REF, H_REF, SCALE_H, RE, FLAT, BETA


def step(el, t_a, cs, C, S, nf):
    a, e, inc, O, w = el; mu = cs["mu"]; N0 = math.sqrt(mu / a ** 3); T0 = 2 * math.pi / N0
    f = np.linspace(-math.pi, math.pi, nf); p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f))
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f / 2)); tau = (E - e * np.sin(E)) / N0
    Q = tl.pqw_matrix(inc, O, w)
    r = Q @ np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); v = Q @ (math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f]))
    th = W_ROT * (t_a + T0 / 2 + tau)
    rb = tl.rot_z(-th, r); acc = tl.rot_z(th, tl.accel_body(rb[0], rb[1], rb[2], cs, C, S))
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    dtdf = rr ** 2 / h; df = f[1] - f[0]
    cum = lambda g: np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * df)])
    # radial offset of the real path from the Kepler arc near periapsis = the accumulated gravity change of r_p (the physical periapsis
    # radius is the osculating r_p at the previous apoapsis plus the inbound kick); without it the density is off by exp(dr/H)
    gR, gS, gW = (acc * rhat).sum(0), (acc * shat).sum(0), (acc * hhat).sum(0)
    dr = cum(tl.gauss_rates(a, e, inc, w, f, rr, h, gR, gS, gW, mu)[5] * dtdf)
    vrel = v - np.array([-W_ROT * r[1], W_ROT * r[0], 0 * f]); vr = np.linalg.norm(vrel, axis=0)
    hgt = rr + dr - RE * (1 - FLAT * (r[2] / rr) ** 2); acc = acc - 0.5 * RHO_REF * np.exp(-(hgt - H_REF) / SCALE_H) * vr * vrel / BETA
    aR, aS, aW = (acc * rhat).sum(0), (acc * shat).sum(0), (acc * hhat).sum(0)
    rates = tl.gauss_rates(a, e, inc, w, f, rr, h, aR, aS, aW, mu)[:5]
    b = a * math.sqrt(1 - e * e)
    dM = b / (a * h * e) * ((p * np.cos(f) - 2 * e * rr) * aR - (p + rr) * np.sin(f) * aS)   # Gauss mean-anomaly equation minus N (Vallado form; recalled, checked by the zonal run)
    da_cum = cum(rates[0] * dtdf); dx = np.array([cum(q * dtdf)[-1] for q in rates]); dMp = cum(dM * dtdf)[-1]
    int_da = cum(da_cum * dtdf)[-1]                       # integral of the accumulated Delta a over time
    # a and e from the accumulated energy and angular momentum (both first-order exact for a frozen arc); the linearized
    # (1-e) da - a de of r_p cancels to a residual smaller than its second-order error once drag takes 1-2% of a per pass
    dE = cum((v * acc).sum(0) * dtdf)[-1]; dH = cum(rr * aS * dtdf)[-1]
    En = -mu / (2 * a) + dE; Hn = h + dH; a1 = -mu / (2 * En); e1 = math.sqrt(max(0.0, 1 + 2 * En * Hn ** 2 / mu ** 2))
    dx[0], dx[1] = a1 - a, e1 - e; N1 = math.sqrt(mu / a1 ** 3)
    dt_ = T0 + (1.5 * N0 / a * int_da - dMp) / N1
    return np.array(el) + dx, t_a + dt_


def run(root, nf):
    d = pd.read_csv(root / "tesseral" / "cowell_drag.csv"); cs = tl.coefficients("Mars", 20)
    fields = {"zonal": tl.select(cs["C"], cs["S"], lambda n, m: m == 0), "full": tl.select(cs["C"], cs["S"], lambda n, m: True)}
    rows = []
    for o, g in d.groupby("orbit"):
        truth = {v: g[g.variant == v].sort_values("pass").reset_index(drop=True) for v in ("zonal", "full")}
        r0 = truth["zonal"].iloc[0]; el0 = (r0.a_m, r0.e, r0.i_rad, r0.raan_rad, r0.omega_rad)
        for fld in ("zonal", "full"):
            el, t = np.array(el0), 0.0; C, S = fields[fld]; K = max(len(truth["zonal"]), len(truth["full"]))
            for k in range(K):
                row = dict(orbit=o, inc_deg=r0.inc_deg, field=fld, k=k, t_d=t / 86400, rp_m=el[0] * (1 - el[1]), a_m=el[0])
                for tv in ("zonal", "full"):
                    T = truth[tv]
                    if k < len(T): row[f"err_vs_{tv}_m"] = row["rp_m"] - T.rp_m.iloc[k]; row[f"dt_vs_{tv}_s"] = t - T.t_apo_s.iloc[k]; row[f"truth_{tv}_rp_m"] = T.rp_m.iloc[k]
                rows.append(row)
                if el[0] * (1 - el[1]) - RE < 80e3: break
                el, t = step(el, t, cs, C, S, nf)
        print(f"orbit {o}: done", flush=True)
    P = pd.DataFrame(rows); P.to_csv(root / "tesseral" / "propagator_passes.csv", index=False)
    edges = [0, 5, 10, 20, 30, 40, 50, 60]; out = []
    for (fld, tv), lab in {("zonal", "zonal"): "zonal map vs zonal truth", ("full", "full"): "full map vs full truth",
                           ("zonal", "full"): "zonal map vs full truth"}.items():
        for lo, hi in zip(edges[:-1], edges[1:]):
            b = P[(P.field == fld) & (P.t_d >= lo) & (P.t_d < hi)].dropna(subset=[f"err_vs_{tv}_m"])
            if len(b): out.append(dict(case=lab, bin=f"{lo}-{hi}", n=len(b), rms_km=math.sqrt(np.mean(b[f"err_vs_{tv}_m"] ** 2)) / 1e3,
                                       timing_rms_s=math.sqrt(np.mean(b[f"dt_vs_{tv}_s"] ** 2))))
    B = pd.DataFrame(out); B.to_csv(root / "tesseral" / "propagator_bins.csv", index=False)
    print(B.pivot_table(index="case", columns="bin", values="rms_km").round(3).to_string())
    print(B.pivot_table(index="case", columns="bin", values="timing_rms_s").round(0).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--nf", type=int, default=4097)
    a = ap.parse_args(); run(pathlib.Path(results_root(a.results_root)), a.nf)
