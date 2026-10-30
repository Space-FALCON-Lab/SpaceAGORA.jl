"""Stage K supplement (UNFROZEN, added after the pre-registered criteria were evaluated; changes no verdict): does the non-zonal (tesseral and sectorial)
part of the 20x20 field explain the Mars pass-to-pass periapsis-radius change that the zonal map and drag do not?

Method (first order, DERIVED, same quadrature as tesseral_lib.gauss_kick but over the periapsis-to-periapsis arc, true anomaly 0 to 2 pi):
  dr_T = integral over the arc of d r_p / dt = (1 - e) da/dt - a de/dt  (Gauss equations, RSW components of the order m >= 1 acceleration of the
  body-fixed field rotated into the inertial frame, planet turning at w_rot), evaluated on the Kepler ellipse of the state at pass k.
Rotation phase at pass k: theta_k = lon_inertial(periapsis) - lon_rad (pass table, body-fixed longitude), checked for consistency with w_rot dt between passes.
w_rot = 7.08823596e-5 rad/s (SOURCED: planets.jl Mars struct, as tesseral_cowell.py). Mars only (Earth and Venus drift is tide-dominated).

Usage: python3 -I stagek_tess.py --f20 <stagef20 dir> --out <dir with pairs_full_mars.csv>
"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import tesseral_lib as tl

W_ROT = 7.08823596e-5


def tess_kick(cs, Ct, St, a, e, inc, raan, om, theta_p, nf=128):
    mu = cs["mu"]; x, w = np.polynomial.legendre.leggauss(nf); f = math.pi * (x + 1.0); w = math.pi * w
    p = a * (1 - e * e); h = math.sqrt(mu * p); rr = p / (1 + e * np.cos(f)); N = math.sqrt(mu / a ** 3)
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f / 2), np.sqrt(1 + e) * np.cos(f / 2)); tau = (E - e * np.sin(E)) / N
    Q = tl.pqw_matrix(inc, raan, om)
    r = Q @ np.array([rr * np.cos(f), rr * np.sin(f), 0 * f]); v = Q @ (math.sqrt(mu / p) * np.array([-np.sin(f), e + np.cos(f), 0 * f]))
    th = theta_p + W_ROT * tau; rb = tl.rot_z(-th, r); ai = tl.rot_z(th, tl.accel_body(rb[0], rb[1], rb[2], cs, Ct, St))
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    aR, aS, aW = (ai * rhat).sum(0), (ai * shat).sum(0), (ai * hhat).sum(0)
    rates = tl.gauss_rates(a, e, inc, om, f, rr, h, aR, aS, aW, mu)
    return float(rates[5] @ (w * rr ** 2 / h))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--f20", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    out = pathlib.Path(a.out); df = pd.read_csv(out / "pairs_full_mars.csv"); root = pathlib.Path(a.f20) / "mars"
    cs = tl.coefficients("Mars", 20); Ct, St = tl.select(cs["C"], cs["S"], lambda n, m: m >= 1)
    lon = {}; lat_err = []
    for tag in df.tag.unique():
        d = pd.read_csv(root / "runs" / f"open_mean_{tag}" / "passes.csv").set_index("pass"); lon[tag] = d
    th = []; dth_err = []
    for r in df.itertuples():
        d = lon[r.tag]; P = np.array([math.cos(r.raan) * math.cos(r.w) - math.sin(r.raan) * math.sin(r.w) * math.cos(r.i),
                                       math.sin(r.raan) * math.cos(r.w) + math.cos(r.raan) * math.sin(r.w) * math.cos(r.i)])
        lam = math.atan2(P[1], P[0]); thk = (lam - d.loc[r.k, "lon_rad"]) % (2 * math.pi); th.append(thk)
        lat_err.append(abs(math.asin(math.sin(r.w) * math.sin(r.i)) - d.loc[r.k, "lat_rad"]))
    df["theta"] = th
    for tag, g in df.groupby("tag"):
        t = g.theta.values; dt = g.dt.values[:-1]
        if len(t) > 1: dth_err += list(np.abs(((t[1:] - t[:-1] - W_ROT * dt) + math.pi) % (2 * math.pi) - math.pi))
    frame = dict(max_lat_err_rad=float(np.max(lat_err)), median_theta_step_err_rad=float(np.median(dth_err)), p95_theta_step_err_rad=float(np.percentile(dth_err, 95)))
    print("frame check", frame)
    df["dr_T"] = [tess_kick(cs, Ct, St, r.a, r.e, r.i, r.raan, r.w, r.theta) for r in df.itertuples()]
    df["dr_ZT"] = df.dr_Z + df.dr_T
    dev = df[df.split == "dev"]; x = dev.da - dev.da_Z; y = dev.dr - dev.dr_ZT; beta2 = float((x * y).sum() / (x * x).sum())
    df["dr_ZTD"] = df.dr_ZT + beta2 * (df.da - df.da_Z); dev = df[df.split == "dev"]
    ev = lambda s, c: 1 - ((s.dr - s[c]) ** 2).sum() / (s.dr ** 2).sum()
    res = dict(frame=frame, beta_with_tesseral=beta2)
    for nm, s in (("heldout", df[df.split == "heldout"]), ("dev", dev), ("all", df)):
        res[f"EV_{nm}"] = {c: float(ev(s, c)) for c in ("dr_Z", "dr_T", "dr_ZT", "dr_ZTD")}
        m = np.abs(s.dr) >= 10; res[f"sign_{nm}_10m"] = float((np.sign(s.dr_ZTD[m]) == np.sign(s.dr[m])).mean())
        res[f"rms_resid_{nm}"] = float(np.sqrt(((s.dr - s.dr_ZTD) ** 2).mean())); res[f"rms_obs_{nm}"] = float(np.sqrt((s.dr ** 2).mean()))
    # 10-pass net change (the accumulating part)
    rows = []
    for tag, g in df.groupby("tag"):
        g = g.sort_values("k")
        for s0 in range(0, len(g) - 10 + 1, 10):
            h = g.iloc[s0:s0 + 10]; rows.append((h.dr.sum(), h.dr_Z.sum(), h.dr_T.sum(), h.dr_ZTD.sum(), h.split.iloc[0]))
    o = pd.DataFrame(rows, columns=["obs", "Z", "T", "ZTD", "split"])
    res["net10"] = {nm: dict(n=int(len(s)), EV_Z=float(1 - ((s.obs - s.Z) ** 2).sum() / (s.obs ** 2).sum()), EV_ZTD=float(1 - ((s.obs - s.ZTD) ** 2).sum() / (s.obs ** 2).sum()),
                             rms_obs=float(np.sqrt((s.obs ** 2).mean()))) for nm, s in (("heldout", o[o.split == "heldout"]), ("all", o))}
    df.to_csv(out / "pairs_full_mars_tess.csv", index=False); (out / "k_tess_results.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
