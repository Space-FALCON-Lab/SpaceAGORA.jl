"""Pass-pair tables for the analytic validation (WP-A3): periapsis-to-periapsis pairs with the per-pair closed-form predictions attached.

Stage O dev flights (Mars, Odyssey parameters, 48 runs; SOURCED tasks/.../stageo_data/dev/runs/<run>/{passes,decisions}.csv). Burn pairs (a nonzero commanded dv
between pass k and k+1, decisions.csv target_pass == k+1) are flagged and excluded from natural-drift tests.
Mars mean-field truth cells with gravity ablations (3b = full, zonal, deg3, j2only; SOURCED analytic_data/snapshot/mars_truth).
Python only; the expensive per-pair tesseral quadrature is cached as CSV under analytic_data/cache/.
"""
import math, os, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analytic_lib as L, stagek_lib as SK, tesseral_lib as TL, apoapsis_pi_lib as AP


def tasks_dir():
    env = os.environ.get("SPACEAGORA_TASKS_DIR")
    if env: return pathlib.Path(env)
    for q in [HERE, *HERE.parents]:
        if (q / "tasks" / "2026-10-15_periapsis-region-border-targeting").is_dir(): return q / "tasks" / "2026-10-15_periapsis-region-border-targeting"
    import subprocess   # worktree: tasks/ lives in the main checkout
    root = pathlib.Path(subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--path-format=absolute", "--git-common-dir"], text=True).strip()).parent
    return root / "tasks" / "2026-10-15_periapsis-region-border-targeting"


TASKS = tasks_dir(); CACHE = TASKS / "analytic_data" / "cache"
CD, AREA = 2.2, 2.2 * 1.7 + 2 * 1.945 * 1.7     # SOURCED stageo_config.toml [vehicle]: cd 2.2, bus 2.2 x 1.7, two panels 1.945 x 1.7 (drag area)
MU_M = 42828370371000.0                          # SOURCED Mars50c.csv header GM (the truth gravity); the Stage O planet record uses the same file
W_MARS = L.PLANETS["mars"]["w_rot"]


def stageo_runs():
    root = TASKS / "stageo_data" / "dev" / "runs"
    for d in sorted(root.iterdir()):
        arm, seed = d.name.split("_pert_ODY_s"); yield d, arm.split("_")[0], int(seed)


def stageo_pairs(use_cache=True, nf=128):
    f = CACHE / "stageo_pairs.csv"
    if use_cache and f.exists(): return pd.read_csv(f)
    cs = TL.coefficients("Mars", 20); Ct, St = TL.select(cs["C"], cs["S"], lambda n, m: m >= 1); P = SK.planet("Mars")
    rows = []
    for d, arm, seed in stageo_runs():
        ps = pd.read_csv(d / "passes.csv"); dec = pd.read_csv(d / "decisions.csv")
        burn = set(dec.loc[dec.delta_v_mps.abs() > 0, "target_pass"].astype(int)); dv = dict(zip(dec.target_pass.astype(int), dec.delta_v_mps))
        om_u = np.unwrap(ps.omega_rad.values); raan_u = np.unwrap(ps.raan_rad.values)
        for k in range(len(ps) - 1):
            r, r1 = ps.iloc[k], ps.iloc[k + 1]
            th = L.theta_p_from_pass(r.raan_rad, r.omega_rad, r.i_rad, r.lon_rad)
            tk = L.tess_arc(cs, Ct, St, r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, th, W_MARS, nf)
            ze, zw, zO, zrp = L.zonal_rates_per_pass(P, r.a_m, r.e, r.i_rad, r.omega_rad)
            beta = r.mass_kg / (CD * AREA); frot = AP.f_rot(r.radius_m, r.a_m, r.vrel_mps, r.lat_rad, W_MARS, MU_M)
            vp = math.sqrt(MU_M * (2 / r.radius_m - 1 / r.a_m))
            rows.append(dict(arm=arm, seed=seed, run=d.name, k=int(r["pass"]), burn=bool((int(r["pass"]) + 1) in burn), dv=dv.get(int(r["pass"]) + 1, 0.0),
                             t=r.time_s, dt=r1.time_s - r.time_s, a=r.a_m, e=r.e, inc=r.i_rad, raan=r.raan_rad, om=r.omega_rad, lat=r.lat_rad, lon=r.lon_rad,
                             rp=r.radius_m, ra=r.ra_m, alt=r.alt_m, mass=r.mass_kg, beta=beta, frot=frot, vp=vp, vrel=r.vrel_mps, h_local=r.h_local_m,
                             rho_mean=r.rho, rho_peak=r.rho_peak_kg_m3, rho_eff=r.rho_eff, q_mean=r.q_w_cm2, q_eff=r.q_eff_w_cm2, r_pert=r.r_pert,
                             o_dr=r1.radius_m - r.radius_m, o_da=r1.a_m - r.a_m, o_dra=r1.ra_m - r.ra_m, o_de=r1.e - r.e, o_dom=om_u[k + 1] - om_u[k],
                             o_dOm=raan_u[k + 1] - raan_u[k], o_di=r1.i_rad - r.i_rad, o_dalt=r1.alt_m - r.alt_m, o_dg=-((r1.radius_m - r1.alt_m) - (r.radius_m - r.alt_m)),
                             z_de=ze, z_dom=zw, z_dOm=zO, z_drp=zrp,
                             t_da=tk["a"], t_de=tk["e"], t_di=tk["i"], t_dOm=tk["raan"], t_dom=tk["omega"], t_drp=tk["rp"], t_dE=tk["dE"], t_dhz=tk["dhz"]))
    df = pd.DataFrame(rows); CACHE.mkdir(parents=True, exist_ok=True); df.to_csv(f, index=False); return df


SNAP = TASKS / "analytic_data" / "snapshot"
MARS_EPOCH_UTC = {1: (2001, 11, 6, 19, 0, 32), 2: (2001, 12, 21, 19, 0, 32), 3: (2002, 2, 4, 19, 0, 32)}   # SOURCED stagef_config.toml [mars] epochs_utc (E1; E2 = +45 d; E3 = +90 d)
MARS_PM = (176.630, 350.89198226)    # SOURCED pck00011.tpc BODY499_PM (deg, deg/day): W = W0 + Wdot d (d days past J2000 TDB)


def mars_theta(et):
    """Planet rotation angle (inertial x of the body-equator frame to the body-fixed x) = prime-meridian angle W; checked against the lon_rad-based value of the F20 tables
    to 2e-5 rad (constant offset removed: W is measured from the node, which is the frame x axis)."""
    return math.radians(MARS_PM[0] + MARS_PM[1] * et / 86400.0) % (2 * math.pi)


def mars_truth_pairs(arms=("3b", "zonal"), use_cache=True, nf=64, tags=None):
    """Mean-field Mars truth (128 cells x 3 gravity arms): per pair (k-1 -> k, as apoapsis_pi_check.pairs) the observed change of r_a + r_p, the Laplace drag prediction from the
    pass's own mean-field heat rate, and the tesseral first-order kicks of the (k-1) state: r_a over the apoapsis-to-apoapsis arc of row k (DERIVED from the data: regressing obs on the drag of passes k-1, k, k+1 gives coefficients -0.05, 1.05, -0.00 on the 3b arm, so apo_alt_m of row k is the
    apoapsis AFTER periapsis k and obs = [apo(k) - apo(k-1)] + [alt(k) - alt(k-1)] = Delta r_a (apoapsis to apoapsis, centered on periapsis k: elements of row k) + Delta r_p (periapsis k-1 to k)) and r_p over the periapsis-to-periapsis arc. Cached."""
    import analytic_ephem as E
    f = CACHE / "mars_truth_pairs.csv"
    if use_cache and f.exists(): return pd.read_csv(f)
    cs = TL.coefficients("Mars", 20); Ct, St = TL.select(cs["C"], cs["S"], lambda n, m: m >= 1); H = L.PLANETS["mars"]["H"]; rows = []
    root = SNAP / "mars_truth"
    for td in sorted(root.iterdir()):
        tag = td.name; ep = int(tag[-1]); et0 = E.utc_to_et(*MARS_EPOCH_UTC[ep])
        if tags is not None and tag not in tags: continue
        for arm in arms:
            d = pd.read_csv(td / f"passes_{arm}.csv")
            for k in range(1, len(d) - 1):
                r0, r = d.iloc[k - 1], d.iloc[k]
                th = mars_theta(et0 + r0.time_s)
                tk = L.tess_arc(cs, Ct, St, r0.a_m, r0.e, r0.i_rad, r0.raan_rad, r0.omega_rad, th, W_MARS, nf)                                   # periapsis to periapsis
                th1 = mars_theta(et0 + r.time_s)
                ta = L.tess_arc(cs, Ct, St, r.a_m, r.e, r.i_rad, r.raan_rad, r.omega_rad, th1, W_MARS, nf, f_lo=-math.pi, f_hi=math.pi)           # apoapsis to apoapsis around periapsis k
                beta = r.mass_kg / (CD * AREA); fr = AP.f_rot(r.radius_m, r.a_m, r.vrel_mps, r.lat_rad, W_MARS, MU_M)
                pred = 2 * AP.delta_a_pass(r.a_m, r.e, r.radius_m, r.q_w_cm2 * 1e4, H, beta, MU_M, fr)          # apoapsis_pi_check convention: the drag of pass k (row k) predicts the pair (k-1, k)
                rows.append(dict(tag=tag, arm=arm, k=int(r["pass"]), a=r0.a_m, e=r0.e, c=r0.a_m * r0.e / H, alt=r.alt_m, pred=pred,
                                 obs=(r.apo_alt_m + r.alt_m) - (r0.apo_alt_m + r0.alt_m), obs_ra=r.apo_alt_m - r0.apo_alt_m,
                                 t_dra=ta["ra"], t_drp=tk["rp"], t_da_apo=ta["a"], t_dE_apo=ta["dE"], t_dhz_apo=ta["dhz"], t_da_per=tk["a"], t_dE_per=tk["dE"], t_dhz_per=tk["dhz"]))
    df = pd.DataFrame(rows); CACHE.mkdir(parents=True, exist_ok=True); df.to_csv(f, index=False); return df


EPOCHS = {   # SOURCED stagef_config.toml epochs_utc per planet (E1..E3)
    "mars": MARS_EPOCH_UTC,
    "earth": {1: (2020, 1, 1, 0, 0, 0), 2: (2020, 2, 15, 0, 0, 0), 3: (2020, 4, 1, 0, 0, 0)},
    "venus": {1: (2014, 5, 19, 4, 0, 0), 2: (2014, 7, 3, 4, 0, 0), 3: (2014, 8, 17, 4, 0, 0)}}
BODIES = {"mars": (10,), "earth": (10, 301), "venus": (10,)}


def stagek_pairs(planet):
    """Stage K pair table (snapshot) of a planet; Mars comes with the tesseral columns (pairs_full_mars_tess.csv)."""
    return pd.read_csv(SNAP / "stagek" / ("pairs_full_mars_tess.csv" if planet == "mars" else f"pairs_full_{planet}.csv"))


def tide_pairs(planet, use_cache=True, nf=48, nodes=7):
    """T8 on the Stage K pair tables: per pair, the direct first-order quadrature of the exact third-body tide along the Kepler arc (perturbers from DE440s, frozen at
    nothing: positions interpolated over the arc by a cubic spline through `nodes` ephemeris points) and the singly averaged closed form at the arc midpoint."""
    import analytic_ephem as E
    from scipy.interpolate import CubicSpline
    f = CACHE / f"tide_pairs_{planet}.csv"
    if use_cache and f.exists(): return pd.read_csv(f)
    df = stagek_pairs(planet); P = SK.planet(planet.capitalize()); mu = P["mu"]; rows = []
    et0 = {ep: E.utc_to_et(*EPOCHS[planet][ep]) for ep in EPOCHS[planet]}
    for r in df.itertuples():
        et = et0[int(r.epoch)] + r.t0; T = 2 * math.pi * math.sqrt(r.a ** 3 / mu); ts = np.linspace(0, T, nodes)
        dr_dir = dr_avg = 0.0
        for b in BODIES[planet]:
            pos = np.array([E.third_body_in_frame(planet, b, et + t) for t in ts]).T          # (3, nodes)
            cs_ = CubicSpline(ts, pos, axis=1); mub = E.GM[b]
            k = L.tide_arc_direct(mu, mub, r.a, r.e, r.i, r.raan, r.w, lambda tau: cs_(tau), nf)
            dr_dir += k["rp"]
            Q = TL.pqw_matrix(r.i, r.raan, r.w); Ph, Qh = Q[:, 0], Q[:, 1]
            dr_avg += L.tide_rp_averaged(r.a, r.e, mu, mub, cs_(0.5 * T), Ph, Qh)
        rows.append((dr_dir, dr_avg))
    df["dr_T8"], df["dr_T8avg"] = np.array(rows).T; CACHE.mkdir(parents=True, exist_ok=True); df.to_csv(f, index=False); return df


PM = {"mars": MARS_PM, "earth": (190.147, 360.9856235), "venus": (160.20, -1.4813688)}   # SOURCED pck00010/00011 BODY<id>_PM (deg, deg/day)


def planet_theta(planet, et):
    """Prime-meridian angle W(et) of the body-equator frame (see mars_theta)."""
    w0, wd = PM[planet]; return math.radians(w0 + wd * et / 86400.0) % (2 * math.pi)


def tess_pairs_other(planet, use_cache=True, nf=64):
    """Tesseral (m >= 1) first-order r_p kick of the Stage K pairs of Earth and Venus (ASSUMED small there; plan 3.3): returns the pair table with dr_T added.
    Gravity: EarthGGM05C (L = 20), MGNP180U (L = 20), from the repo files. Frame: body-equator frame of analytic_ephem, theta = W(t)."""
    import analytic_ephem as E
    f = CACHE / f"tess_pairs_{planet}.csv"
    if use_cache and f.exists(): return pd.read_csv(f)
    df = stagek_pairs(planet); cs = TL.coefficients(planet.capitalize(), 20); Ct, St = TL.select(cs["C"], cs["S"], lambda n, m: m >= 1)
    w_rot = L.PLANETS[planet]["w_rot"]; et0 = {ep: E.utc_to_et(*EPOCHS[planet][ep]) for ep in EPOCHS[planet]}; out = []
    for r in df.itertuples():
        th = planet_theta(planet, et0[int(r.epoch)] + r.t0)
        out.append(L.tess_arc(cs, Ct, St, r.a, r.e, r.i, r.raan, r.w, th, w_rot, nf)["rp"])
    df["dr_T"] = out; CACHE.mkdir(parents=True, exist_ok=True); df.to_csv(f, index=False); return df
