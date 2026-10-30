"""Independent numerical check of the Sun-tide path used by the Venus fixed-epoch study (offline Python, no flights, no GRAM).

The ranking model sums, pass by pass, a first-order Sun kick at frozen osculating elements (omega, Omega, e are not re-evolved by the Sun). This script integrates the real equations of motion
(Venus point mass + J2 from the ranking config + Sun differential acceleration, Sun position from the engine-sourced table), applies the nominal drag decay as an impulsive tangential
delta-v at each periapsis (apoapsis drop from the ranking model's nominal schedule), and records the osculating periapsis radius at every apoapsis. No Sun-frozen-element assumption, no
short-period closed form. It checks (a) the closed-form Sun path, (b) the DP floor computed on the numerical path, for a list of candidates.

Usage: python3 -I venus_fixed_epoch_cowell.py OUT_DIR [--jobs N] [--one]      (OUT_DIR = the directory holding truth.csv and analysis/)
Sun-only variant for the drag impulse is the nominal one; the Cowell path has no zonal J3+, tesseral, or surface channel (not in the equations here; the model's versions of those are 0.005 to 0.1 ln).
"""
import json, math, pathlib, sys, time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V
import venus_fixed_epoch as E

_W = {}


def sun_spline(ctx, off, T_days=62.0, step_h=1.0):
    days = ctx.anchor_days + off + np.arange(0.0, T_days, step_h / 24.0)
    R = np.array([G.sun_in_mars_frame(d) for d in days]); return CubicSpline((days - days[0]) * 86400.0, R)


def propagate(ctx, nom, cand, sun, rtol=1e-11):
    mu, J2, Rg, ms = ctx.mu, ctx.J2, ctx.Rg, G.MU_SUN; o = ctx.cfg["orbit"]
    inc, om, Om = (math.radians(cand[k]) for k in ("inc", "omega", "raan")); ra, rp = o["ra0_m"], ctx.Req + o["h_ref_m"]; a = 0.5 * (ra + rp); va = math.sqrt(mu * (2 / ra - 1 / a))
    P, Q, W = E.pq_vectors(np.array(inc), np.array(Om), np.array(om)); y0 = np.concatenate([-ra * P, -va * Q])

    def rhs(t, y):
        r = y[:3]; rn = math.sqrt(r @ r); acc = -mu * r / rn ** 3
        z2 = (r[2] / rn) ** 2; k = 1.5 * J2 * mu * Rg ** 2 / rn ** 5; acc = acc + k * np.array([r[0] * (5 * z2 - 1), r[1] * (5 * z2 - 1), r[2] * (5 * z2 - 3)])
        Rs = sun(t); d = Rs - r; acc = acc + ms * (d / np.linalg.norm(d) ** 3 - Rs / np.linalg.norm(Rs) ** 3); return np.concatenate([y[3:], acc])

    def per(t, y): return y[:3] @ y[3:]
    per.terminal = True; per.direction = 1                                         # periapsis: r.v from - to +
    def apo(t, y): return y[:3] @ y[3:]
    apo.terminal = True; apo.direction = -1
    def osc(y):
        r, v = y[:3], y[3:]; rn = np.linalg.norm(r); aa = 1 / (2 / rn - v @ v / mu); h = np.cross(r, v); e = math.sqrt(max(0.0, 1 - (h @ h) / (mu * aa))); return aa * (1 - e), aa, e
    t = 0.0; y = y0.copy(); rps = [osc(y)[0]]; n = nom.n; tp = []
    for k in range(n - 1):
        s = solve_ivp(rhs, (t, t + 3 * 86400.0), y, method="DOP853", rtol=rtol, atol=1e-6, events=per); t, y = s.t[-1], s.y[:, -1]; tp.append(t)
        r, v = y[:3], y[3:]; rn = np.linalg.norm(r); aa = 1 / (2 / rn - v @ v / mu); an = aa - nom.S["drop"][k] / 2; vn = math.sqrt(mu * (2 / rn - 1 / an)); y = np.concatenate([r, v / np.linalg.norm(v) * vn])
        s = solve_ivp(rhs, (t, t + 3 * 86400.0), y, method="DOP853", rtol=rtol, atol=1e-6, events=apo); t, y = s.t[-1], s.y[:, -1]; rps.append(osc(y)[0])
    rps = np.array(rps); path = -(rps - rps[0]) / nom.H                                                   # ln q (periapsis descent in scale heights is positive); pass k+1 state after arc k
    return path, np.array(tp)


def one(args):
    c = args; ctx, nom = _W["ctx"], _W["nom"]; sun = sun_spline(ctx, c["epoch_off"]); path, tp = propagate(ctx, nom, c, sun)
    cf = E.closed_form(ctx, nom, c["epoch_off"], np.radians([c["inc"]]), np.radians([c["omega"]]), np.radians([c["raan"]]))["path3"][0]
    tot_num = path + nom.other; tot_cf = cf + nom.other; sn = E.window(tot_num); sc = E.window(tot_cf)
    return dict(id=c["id"], set=c["set"], epoch_off=c["epoch_off"], inc=c["inc"], omega0=c["omega"], raan0=c["raan"], R_num=float(np.ptp(sn)), R_cf=float(np.ptp(sc)), N_num=float(sn[-1] - sn[0]), N_cf=float(sc[-1] - sc[0]),
                floor_num=E.dp_floor(sn, nom.w, nom.W_eff), floor_cf=E.dp_floor(sc, nom.w, nom.W_eff), sun_path_rms_err=float(np.sqrt(np.mean((path - cf) ** 2))), sun_path_end_num=float(path[-1]), sun_path_end_cf=float(cf[-1]),
                sun_path_max_abs_num=float(np.abs(path).max()), t_first_peri_h=float(tp[0] / 3600), t_last_peri_d=float(tp[-1] / 86400))


def main(out_dir, jobs, only_one=False):
    import multiprocessing as mp
    out = pathlib.Path(out_dir); ctx, info = V.make_ctx(); nom = E.Nominal(ctx); _W.update(ctx=ctx, nom=nom); cands = []
    for ep in E.HEAD:
        for om in range(0, 360, 30):
            for ra in (0, 90, 180, 270): cands.append(dict(set="polar", inc=89.876, epoch_off=ep, omega=float(om), raan=float(ra)))
        for inc in (39.2, 63.4, 116.6):
            for om in (0, 90):
                for ra in (0, 180): cands.append(dict(set="other_inc", inc=inc, epoch_off=ep, omega=float(om), raan=float(ra)))
    pk = out / "analysis" / "planner_fine_pool_picks.csv"
    if pk.exists():
        p = pd.read_csv(pk); p = p[p["rank"] < 3]
        for r in p.itertuples(): cands.append(dict(set="planner_pick", inc=float(r.inc), epoch_off=float(r.epoch_off), omega=float(r.omega0), raan=float(r.raan0)))
    for i, c in enumerate(cands): c["id"] = "C%04d" % i
    if only_one: cands = cands[:1]
    t0 = time.time(); print(len(cands), "candidates", flush=True)
    with mp.get_context("fork").Pool(jobs) as pool: rows = pool.map(one, cands, chunksize=2)
    R = pd.DataFrame(rows); R.to_csv(out / "analysis" / "cowell_check.csv", index=False); print("done %.0f s" % (time.time() - t0)); print(R.round(3).head(12).to_string())


if __name__ == "__main__":
    a = sys.argv[1:]; jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 6; main(a[0], jobs, "--one" in a)
