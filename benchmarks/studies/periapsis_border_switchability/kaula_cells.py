"""Per-cell Kaula parameters (Pi_n, rate ratios, horizon effects, Pi_K) and the Cowell numerical check, for one planet's cell set.

Usage: python3 -I kaula_cells.py --planet Earth|Mars --stage analytic|cowell [--results-root DIR] [--gravity-dir DIR] [--workers W] [--nsum 200]
Reads <root>/mars/{pi_cells.csv, planet_constants.csv, mars_cell_labels.csv, j3_exposure.csv, truth/}, <root>/oracle/oracle_cell_labels.csv,
the Earth truth under <root>/border{,2}/grid and the Earth v2 predictions <root>/estimator_v2/pred. Writes <root>/kaula/cells_<planet>_{analytic,cowell}.csv.
Definitions are in project_specs/kaula_harmonic_parameter.md. K is the W_HI geometric-mean fit of kaula_calibrate.py.
"""
import argparse, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import kaula_lib as kl

HORIZON_D = 60.0
NMAX_CMP = 12            # degrees reported individually
CIRC = 12                # Cowell: zonals up to this degree for the "all zonals" run


def paths(root, planet):
    m = root / "mars"
    if planet == "Mars":
        L = pd.read_csv(m / "mars_cell_labels.csv"); L = L[L.status != "missing"]
        hist = lambda t, inc: m / "truth" / t / "passes_3b.csv"; pred = lambda t: m / "pred" / f"{t}.csv"
    else:
        L = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv")
        hist = lambda t, inc: (root / "border" / "grid" / t[6:] if inc == 50.0 else root / "border2" / "grid" / t) / "passes_3b.csv"
        pred = lambda t: root / "estimator_v2" / "pred" / f"{t}.csv"
    P = pd.read_csv(m / "pi_cells.csv"); P = P[P.planet == planet].set_index("tag")
    pc = pd.read_csv(m / "planet_constants.csv").set_index("planet").loc[planet]
    return L, P, pc, hist, pred


def fitK(gdir, planet):
    df, R, mu, _ = kl.read_gravity(planet, gdir); n, A = kl.degree_amplitude(df)
    nz = n[np.isfinite(A) & (A > 0)]; return kl.fit_K(n, A, 10, int(nz.max()))[0], df, R, mu


class GammaCache:
    def __init__(self, nmax): self.d = {}; self.nmax = nmax
    def get(self, e, inc):
        key = (round(inc, 4), round(e, 3))
        if key not in self.d:
            self.d[key] = np.array([0.0, 0.0, 0.0] + [kl.gamma_max(n, key[1], key[0], nw=64, nq=4096 if n < 120 else 8192) for n in range(3, self.nmax + 1)])
        return self.d[key]


def analytic(args):
    root = pathlib.Path(results_root(args.results_root)); out = root / "kaula"; out.mkdir(exist_ok=True)
    L, P, pc, hist, predf = paths(root, args.planet)
    K, gdf, R, mu = fitK(args.gravity_dir, args.planet); Jfile = kl.zonal_J(gdf); nfile = max(Jfile)
    ns = np.arange(3, args.nsum + 1)
    Jk = np.sqrt(2 * ns + 1) * K / ns ** 2                       # Kaula rms J_n (unnormalized)
    Jh = np.array([abs(Jfile[n]) if n <= nfile else Jk[n - 3] for n in ns])   # measured to the file degree, Kaula beyond
    gc = GammaCache(args.nsum); rows = []
    for c in L.itertuples():
        pi = P.loc[c.tag]; a, e = float(pi.a_m), float(pi.e); inc = np.radians(c.inc_deg); p = a * (1 - e * e); rp = a * (1 - e); N = np.sqrt(mu / a ** 3)
        Pi2 = abs(Jfile[2]) * (R / p) ** 2; _, wdot2, _ = kl.zonal_rates(2, e, inc, np.radians(c.omega_deg), Jfile[2] * (R / p) ** 2, N)
        Th = min(HORIZON_D, float(c.t_end_days)) * 86400.0; Teff = Th if abs(wdot2) * Th < 2 else 2.0 / abs(wdot2)
        G = gc.get(e, inc)[3:]; pfac = (R / p) ** ns
        PiN_h, PiN_k = Jh * pfac, Jk * pfac
        rate_h, rate_k = a * N * PiN_h * G, a * N * PiN_k * G    # m/s, envelope rate of the periapsis radius (max over omega)
        Dh, Dk = rate_h * Teff / 1e3, rate_k * Teff / 1e3          # km over the horizon
        PiJ3B = float(pi.Pi_3B); r3b = 1.875 * e * np.sqrt(1 - e * e) * PiJ3B                # dimensionless 3B rate coefficient (max over angles)
        D3B = a * N * r3b * Th / 1e3
        row = dict(tag=c.tag, planet=args.planet, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, epoch=c.epoch, status=c.status, t_end_days=c.t_end_days,
                   a_m=a, e=e, p_m=p, R_over_p=R / p, R_over_rp=R / rp, Pi_J2=Pi2, Pi_3B=PiJ3B, wdot_J2_rad_s=wdot2, Teff_days=Teff / 86400, D_3B_km=D3B)
        for n in range(3, NMAX_CMP + 1):
            j = n - 3
            row.update({f"Pi_{n}": PiN_h[j], f"Gmax_{n}": G[j], f"rho_{n}": PiN_h[j] * G[j] / r3b, f"D_{n}_km": Dh[j], f"Dk_{n}_km": Dk[j], f"Pi_{n}_over_Pi3B": PiN_h[j] / PiJ3B})
        for Nin in (2, 3, 4, 5, 6, 8):
            s = ns > Nin
            row[f"PiK_wc_N{Nin}_km"] = Dh[s].sum(); row[f"PiK_rss_N{Nin}_km"] = np.sqrt((Dh[s] ** 2).sum())
            row[f"PiKk_wc_N{Nin}_km"] = Dk[s].sum(); row[f"PiKk_rss_N{Nin}_km"] = np.sqrt((Dk[s] ** 2).sum())
            row[f"Pitail_N{Nin}"] = PiN_k[s].sum(); s12 = s & (ns <= CIRC); row[f"PiK_wc_N{Nin}_to{CIRC}_km"] = Dh[s12].sum() if s12.any() else 0.0
            row[f"term_last_over_sum_N{Nin}"] = Dk[-1] / max(Dk[s].sum(), 1e-300)
        # horizon integral along the truth history (cancellation-aware), degrees 3..6, osculating elements per pass
        T = pd.read_csv(hist(c.tag, c.inc_deg), usecols=["time_s", "a_m", "e", "i_rad", "omega_rad"]); T = T[T.time_s <= Th + 1]
        for n in (3, 4, 5, 6):
            Jn = Jfile[n]; ed = np.zeros(len(T))
            for k, t in enumerate(T.itertuples()):
                pk = t.a_m * (1 - t.e ** 2); Nk = np.sqrt(mu / t.a_m ** 3)
                ed[k] = t.a_m * kl.zonal_rates(n, t.e, t.i_rad, t.omega_rad, Jn * (R / pk) ** n, Nk, nq=2048)[0] * -1.0    # d rp / dt = -a de/dt
            tt = T.time_s.values; cum = np.concatenate([[0.0], np.cumsum(0.5 * (ed[1:] + ed[:-1]) * np.diff(tt))]) / 1e3
            row[f"Hrms_{n}_km"] = float(np.sqrt(np.mean(cum ** 2))); row[f"Hmax_{n}_km"] = float(np.max(np.abs(cum))); row[f"Hend_{n}_km"] = float(cum[-1])
        # estimator error against the truth (Earth: v2; Mars: frozen v2) over common passes, km
        pr = predf(c.tag)
        if pr.exists():
            Pp = pd.read_csv(pr, usecols=["alt_m"]); Tt = pd.read_csv(hist(c.tag, c.inc_deg), usecols=["alt_m"]); nn = min(len(Pp), len(Tt))
            row["est_rms_err_km"] = float(np.sqrt(np.mean((Pp.alt_m.values[:nn] - Tt.alt_m.values[:nn]) ** 2)) / 1e3)
        rows.append(row)
    df = pd.DataFrame(rows); df.to_csv(out / f"cells_{args.planet}_analytic.csv", index=False)
    print(args.planet, "K =", K, "cells", len(df), "gamma cache", len(gc.d)); return df


def _cowell_one(job):
    import kaula_cowell as kc
    tag, a, e, inc, raan, om, th, Jf, mu, R = job
    J2 = {2: Jf[2]}; J3 = {2: Jf[2], 3: Jf[3]}; J4 = {n: Jf[n] for n in (2, 3, 4)}; Jall = {n: Jf[n] for n in range(2, CIRC + 1)}
    res = {}
    for name, Js in (("J2", J2), ("J23", J3), ("J24", J4), ("Jall", Jall)):
        res[name] = kc.periapsis_series(Js, a, e, inc, raan, om, th, mu, R)[1]
    def rms(x, y): n = min(len(x), len(y)); return float(np.sqrt(np.mean((x[:n] - y[:n]) ** 2)) / 1e3), float(np.max(np.abs(x[:n] - y[:n])) / 1e3)
    d = dict(tag=tag, n_pass=len(res["J2"]))
    for lab, (x, y) in dict(E3=("J23", "J2"), E4=("J24", "J23"), E5p=("Jall", "J24"), Eall=("Jall", "J23")).items():
        d[f"{lab}_rms_km"], d[f"{lab}_max_km"] = rms(res[x], res[y])
    return d


def cowell(args):
    from multiprocessing import Pool
    root = pathlib.Path(results_root(args.results_root)); out = root / "kaula"
    L, P, pc, hist, _ = paths(root, args.planet); gdf, R, mu, _ = kl.read_gravity(args.planet, args.gravity_dir); Jf = kl.zonal_J(gdf)
    jobs = []
    for c in L.itertuples():
        pi = P.loc[c.tag]; jobs.append((c.tag, float(pi.a_m), float(pi.e), np.radians(c.inc_deg), np.radians(c.raan_deg), np.radians(c.omega_deg), min(HORIZON_D, float(c.t_end_days)) * 86400.0, Jf, mu, R))
    if args.limit: jobs = jobs[::max(1, len(jobs) // args.limit)]
    with Pool(args.workers) as pool: rows = pool.map(_cowell_one, jobs, chunksize=2)
    pd.DataFrame(rows).to_csv(out / f"cells_{args.planet}_cowell.csv", index=False); print("cowell rows", len(rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--planet", required=True, choices=["Earth", "Mars"]); ap.add_argument("--stage", required=True, choices=["analytic", "cowell"])
    ap.add_argument("--gravity-dir", default=None); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--nsum", type=int, default=200); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); analytic(a) if a.stage == "analytic" else cowell(a)
