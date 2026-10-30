"""Tesseral parameter Pi_T / D_T per cell (Mars 128, Earth 336), the deterministic per-pass map at Mars, and the pre-registered criteria.

Usage: python3 -I tesseral_param.py [--results-root DIR] [--pi-cells FILE] [--nf 96] [--limit N]
Definitions and criteria: tasks/2026-10-15_periapsis-region-border-targeting/tesseral_preregistration.md (written before this ran).
Reads (results root of the 20x20 truth): mars/{cells.csv, truth_variant_differences.csv, truth/<tag>/passes_{3b,zonal}.csv},
earth_cells_336.csv, earth_truth_old_vs_new.csv, border/grid, border2/grid; --pi-cells: the Stage C pi_cells.csv (Earth a, e; it lives
under the 3x0 results root). Writes <root>/tesseral/{param_cells_Mars.csv, param_cells_Earth.csv, map_passes_Mars.csv, map_bins.csv,
criteria.csv, sigma_n.csv}.
"""
import argparse, datetime as dt, math, pathlib, sys, tomllib
import numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import tesseral_lib as tl
import kaula_lib as kl

L = 20
SPIN = {"Mars": 7.08823596e-5, "Earth": 7.2921158553e-5}   # SOURCED: planets.jl Mars struct; planet_constants.csv (Earth)
HORIZON_D = 60.0
# SOURCED: pck00011.tpc (GRAM Suite SPICE/pck) BODY499_PM and the Mars nutation-precession angles 21..26 with their PM amplitudes (deg, deg/century)
PM0, PM1 = 176.049863, 350.891982443297
PM_TERMS = [(0.000145, 129.071773, 19140.0328244), (0.000157, 36.352167, 38281.0473591), (0.000040, 56.668646, 57420.9295360),
            (0.000001, 67.364003, 76560.2552215), (0.000001, 104.792680, 95700.4387578), (0.584542, 95.391654, 0.5042615)]
TDB_MINUS_UTC = 32 + 32.184    # SOURCED: naif0012.tls DELTA_AT = 32 s from 1999-JAN-1 to 2006-JAN-1, plus TT - TAI = 32.184 s (TDB - TT < 2 ms ignored)


def mars_W(epoch_utc, t_s):
    """Mars prime-meridian angle W (rad) from the IAU_MARS node, at epoch + t_s (pck00011)."""
    j2000 = dt.datetime(2000, 1, 1, 12, 0, 0)
    d = ((epoch_utc - j2000).total_seconds() + TDB_MINUS_UTC + np.asarray(t_s)) / 86400.0; T = d / 36525.0
    W = PM0 + PM1 * d + sum(c * np.sin(np.radians(a + b * T)) for c, a, b in PM_TERMS)
    return np.radians(W % 360.0)


def sigma_n(planet, nmax=50):
    cs = tl.coefficients(planet, nmax); Cb, Sb = cs["Cbar"], cs["Sbar"]
    s = np.zeros(nmax + 1)
    for n in range(2, nmax + 1): s[n] = math.sqrt(((Cb[n, 1:n + 1] ** 2 + Sb[n, 1:n + 1] ** 2).sum()) / (2 * n))
    return s, cs


def gain(m, nu0, nuend, Kh):
    """Accumulation gain G_m^2 (pre-registered closed form)."""
    g = min(1.0 / (2 * math.sin(math.pi * m * nu0) ** 2), Kh ** 2 / 3.0)
    nup = (nuend - nu0) / max(Kh - 1, 1)
    if abs(nup) > 0:
        lo, hi = sorted((m * nu0, m * nuend))
        for j in range(int(math.floor(lo)) + 1, int(math.ceil(hi))):
            Kc = abs(m * nu0 - j) / (m * abs(nup)); g += (1 / (m * abs(nup))) * max(0.0, 1 - Kc / Kh)
    return g


def cell_param(hist, a1, e1, i1, w1, sig, mu, R, spin, nf, theta=None, map_coef=None):
    """Returns dict with D_T, D_T_hist, parts; if theta given (absolute rotation angle per pass) also the deterministic map arrays."""
    t = hist.time_s.values; Kh = len(hist)
    rp1 = a1 * (1 - e1)
    Y1 = tl.term_kicks(a1, e1, i1, 0.0, w1, 0.0, spin, mu, R, L, nf=nf)
    nu = spin * np.diff(t) / (2 * math.pi); nu0, nuend = nu[0], nu[-1]
    G2 = np.array([0.0] + [gain(m, nu0, nuend, Kh) for m in range(1, L + 1)])
    kap = sig[:L + 1, None] * np.abs(Y1)                               # rms per-pass kick, m
    DT = math.sqrt((kap ** 2 * G2[None, :]).sum())
    DT_nores = math.sqrt((kap ** 2 * np.minimum(1 / (2 * np.sin(np.pi * np.arange(L + 1) * nu0) ** 2 + 1e-300), Kh ** 2 / 3)[None, :] * (np.arange(L + 1) > 0)).sum())
    th = theta if theta is not None else spin * t
    S = np.zeros((L + 1, L + 1), complex); acc = []; mp_full = [0.0]; mp_half = []
    for k, r in enumerate(hist.itertuples()):
        Yk = tl.term_kicks(r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, th[k], spin, mu, R, L, nf=nf)
        if map_coef is not None:
            Yh = tl.term_kicks(r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, th[k], spin, mu, R, L, nf=nf, half=True)
            mp_half.append(float(np.real((map_coef * Yh).sum()))); mp_full.append(mp_full[-1] + float(np.real((map_coef * Yk).sum())))
        acc.append((sig[:L + 1, None] ** 2 * np.abs(S) ** 2).sum()); S += Yk
    out = dict(rp1_m=rp1, R_over_rp=R / rp1, nu0=nu0, nu_end=nuend, K_h=Kh, D_T_km=DT / 1e3, D_T_nores_km=DT_nores / 1e3,
               D_T_hist_km=math.sqrt(np.mean(acc)) / 1e3, kick_rms_m=math.sqrt((kap ** 2).sum()),
               D_T_m1to4_km=math.sqrt((kap[:, 1:5] ** 2 * G2[None, 1:5]).sum()) / 1e3,
               **{f"Pi_T_{n}": sig[n] * (R / rp1) ** n for n in (2, 3, 4, 6, 10, 20)},
               **{f"kappa_m{m}_m": math.sqrt((kap[:, m] ** 2).sum()) for m in range(1, 9)},
               **{f"G_m{m}": math.sqrt(G2[m]) for m in range(1, 9)})
    if map_coef is not None: out["map"] = np.array(mp_full[:-1]) + np.array(mp_half)    # prediction at pass k (1-based order)
    return out, Y1


def mars(root, nf, limit, sig, cs50):
    m = root / "mars"; cells = pd.read_csv(m / "cells.csv"); lab = pd.read_csv(m / "truth_variant_differences.csv").set_index("tag")
    cfg = tomllib.load(open(HERE / "mars_generality_config_g20.toml", "rb")); epochs = [dt.datetime.fromisoformat(x) for x in cfg["grid"]["epochs_utc"]]
    cs = tl.coefficients("Mars", L); mu, R = cs["mu"], cs["R"]; coef = (cs["Cbar"] - 1j * cs["Sbar"]); coef[:, 0] = 0
    rows, mp = [], []
    for c in (cells.itertuples() if not limit else list(cells.itertuples())[:limit]):
        z = pd.read_csv(m / "truth" / c.tag / "passes_zonal.csv"); f = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv")
        z = z[z.time_s <= HORIZON_D * 86400 + 1].reset_index(drop=True)
        ep = epochs[int(c.epoch_idx) - 1]; th = mars_W(ep, z.time_s.values)
        o, Y1 = cell_param(z, z.a_m[0], z.e[0], z.i_rad[0], z.omega_rad[0], sig, mu, R, SPIN["Mars"], nf, theta=th, map_coef=coef)
        # Kaula/file tail beyond the truth's degree (real-planet statement, not tested): file to 50 at pass 1, constant-nu gain only
        n = min(len(z), len(f)); d = (f.alt_m.values[:n] - z.alt_m.values[:n])
        for k in range(n):
            mp.append(dict(tag=c.tag, inc=c.inc, k=k + 1, t_d=z.time_s[k] / 86400, truth_diff_m=d[k], map_m=o["map"][k]))
        o.pop("map")
        rows.append(dict(tag=c.tag, inc=c.inc, raan=c.raan, omega=c.omega, epoch_idx=c.epoch_idx, y_km=lab.loc[c.tag, "tesseral_rms_km"],
                         t_end_full=lab.loc[c.tag, "t_end_full"], t_end_zonal=lab.loc[c.tag, "t_end_zonal"], **o))
        print(c.tag, f"D_T {o['D_T_km']:.2f} hist {o['D_T_hist_km']:.2f} y {rows[-1]['y_km']:.2f}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(mp)


def earth(root, nf, limit, sig, pic):
    cells = pd.read_csv(root / "earth_cells_336.csv"); lab = pd.read_csv(root / "earth_truth_old_vs_new.csv").set_index("tag")
    P = pd.read_csv(pic); P = P[P.planet == "Earth"].set_index("tag")
    cs = tl.coefficients("Earth", L); mu, R = cs["mu"], cs["R"]; rows = []
    for c in (cells.itertuples() if not limit else list(cells.itertuples())[:limit]):
        g = (root / "border" / "grid" / c.tag[6:]) if c.inc_deg == 50.0 else (root / "border2" / "grid" / c.tag)
        z = pd.read_csv(g / "passes_3b.csv"); z = z[z.time_s <= HORIZON_D * 86400 + 1].reset_index(drop=True)
        o, _ = cell_param(z, float(P.loc[c.tag, "a_m"]), float(P.loc[c.tag, "e"]), math.radians(c.inc_deg), math.radians(c.omega_deg), sig, mu, R, SPIN["Earth"], nf)
        rows.append(dict(tag=c.tag, inc=c.inc_deg, raan=c.raan_deg, omega=c.omega_deg, y_km=lab.loc[c.tag, "rms_alt_diff_m"] / 1e3, **o))
    return pd.DataFrame(rows)


def criteria(M, E, MP):
    res = []
    def add(i, s, v, ok, thr): res.append(dict(id=i, statement=s, value=v, threshold=thr, verdict="PASS" if ok else "FAIL"))
    for col, tagc in (("D_T_km", ""), ("D_T_hist_km", " (secondary, D_T_hist)")):
        r = M.y_km / M[col]; add("T1" + tagc, "Mars median y / D_T", float(r.median()), 0.5 <= r.median() <= 2, "[0.5, 2]")
        fr = float(((M.y_km >= M[col] / 3) & (M.y_km <= 3 * M[col])).mean()); add("T2" + tagc, "Mars fraction within x3", fr, fr >= 0.8, ">= 0.80")
        rho = float(spearmanr(M[col], M.y_km)[0]); add("T3" + tagc, "Mars Spearman", rho, rho >= 0.3, ">= 0.30")
        pos = M.y_km > 1; rec = float(((M[col] > 1) & pos).sum() / max(pos.sum(), 1)); add("T4" + tagc, f"Mars recall of y > 1 km ({int(pos.sum())} positives; flagged {int((M[col] > 1).sum())})", rec, rec >= 0.95, ">= 0.95")
        fe = float((E[col] < 1).mean()); add("T5" + tagc, "Earth fraction D_T < 1 km", fe, fe >= 0.95, ">= 0.95")
        add("T6" + tagc, "Earth median D_T / median y^E", float(E[col].median() / E.y_km.median()), E[col].median() <= 3 * E.y_km.median(), "<= 3")
        add("T7" + tagc, "median D_T Mars / Earth", float(M[col].median() / E[col].median()), M[col].median() / E[col].median() >= 10, ">= 10")
    b = MP[(MP.t_d >= 0) & (MP.t_d < 5)]; err = math.sqrt(np.mean((b.map_m - b.truth_diff_m) ** 2)) / 1e3; eff = math.sqrt(np.mean(b.truth_diff_m ** 2)) / 1e3
    add("M1", f"map RMS error 0-5 d (km); truth-difference RMS {eff:.3f} km", err, err <= 1 and err <= 0.5 * eff, "<= 1 km and <= 0.5 x effect")
    return pd.DataFrame(res)


def tail(root):
    """Real-planet statement (not a test; the truth stops at 20x20): D_T contribution of degrees 21..50 (Mars50c file values, orders
    1..n) against 2..20, at pass-1 elements and the closed-form gain of each cell; Earth: GGM05C to 50. One cell per inclination."""
    out = []
    for planet, csvn in (("Mars", "param_cells_Mars.csv"), ("Earth", "param_cells_Earth.csv")):
        P = pd.read_csv(root / "tesseral" / csvn); sig, _ = sigma_n(planet, 50); cs = tl.coefficients(planet, 50)
        for inc, g in P.groupby("inc"):
            c = g.iloc[0]
            if planet == "Mars":
                z = pd.read_csv(root / "mars" / "truth" / c.tag / "passes_zonal.csv"); a1, e1, i1, w1 = z.a_m[0], z.e[0], z.i_rad[0], z.omega_rad[0]
            else:
                gdir = (root / "border" / "grid" / c.tag[6:]) if inc == 50.0 else (root / "border2" / "grid" / c.tag)
                z = pd.read_csv(gdir / "passes_3b.csv"); a1, e1, i1, w1 = z.a_m[0], z.e[0], z.i_rad[0], z.omega_rad[0]
            Y = tl.term_kicks(a1, e1, i1, 0.0, w1, 0.0, SPIN[planet], cs["mu"], cs["R"], 50, nf=256)
            G2 = np.array([0.0] + [gain(m, c.nu0, c.nu_end, int(c.K_h)) for m in range(1, 51)])
            k2 = (sig[:, None] * np.abs(Y)) ** 2 * G2[None, :]
            out.append(dict(planet=planet, inc=inc, tag=c.tag, D_T_2to20_km=math.sqrt(k2[:21].sum()) / 1e3, D_T_21to50_km=math.sqrt(k2[21:].sum()) / 1e3,
                            kick_rms_21to50_m=math.sqrt(((sig[21:, None] * np.abs(Y[21:])) ** 2).sum()),
                            **{f"D_T_n2to{N}_km": math.sqrt(k2[:N + 1].sum()) / 1e3 for N in (3, 4, 6, 8, 12)},
                            **{f"D_T_m1to{M_}_km": math.sqrt(k2[:21, :M_ + 1].sum()) / 1e3 for M_ in (2, 3, 4, 6)}))
    T = pd.DataFrame(out); T.to_csv(root / "tesseral" / "tail_21to50.csv", index=False); print(T.round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--pi-cells", default=None); ap.add_argument("--nf", type=int, default=96); ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tail", action="store_true", help="only the degree 21..50 statement (needs the param_cells files)")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); out = root / "tesseral"; out.mkdir(parents=True, exist_ok=True)
    if a.tail: return tail(root)
    sM, csM = sigma_n("Mars"); sE, csE = sigma_n("Earth")
    pd.DataFrame(dict(n=np.arange(2, 51), sigma_Mars=sM[2:], sigma_Earth=sE[2:])).to_csv(out / "sigma_n.csv", index=False)
    M, MP = mars(root, a.nf, a.limit, sM, csM); E = earth(root, a.nf, a.limit, sE, a.pi_cells)
    M.to_csv(out / "param_cells_Mars.csv", index=False); E.to_csv(out / "param_cells_Earth.csv", index=False); MP.to_csv(out / "map_passes_Mars.csv", index=False)
    edges = [0, 5, 10, 20, 30, 40, 50, 60]; bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = MP[(MP.t_d >= lo) & (MP.t_d < hi)]
        if len(b): bins.append(dict(bin=f"{lo}-{hi}", n=len(b), truth_rms_km=math.sqrt(np.mean(b.truth_diff_m ** 2)) / 1e3,
                                    map_err_rms_km=math.sqrt(np.mean((b.map_m - b.truth_diff_m) ** 2)) / 1e3, zero_err_rms_km=math.sqrt(np.mean(b.truth_diff_m ** 2)) / 1e3))
    pd.DataFrame(bins).to_csv(out / "map_bins.csv", index=False)
    C = criteria(M, E, MP); C.to_csv(out / "criteria.csv", index=False)
    pd.set_option("display.width", 250); print(pd.DataFrame(bins).round(3).to_string(index=False)); print(C.to_string(index=False))


if __name__ == "__main__":
    main()
