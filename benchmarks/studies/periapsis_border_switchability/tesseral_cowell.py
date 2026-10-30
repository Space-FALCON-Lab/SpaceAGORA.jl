"""Vectorized Cowell runs for the tesseral derivation (Mars50c, body-equator frame, uniform rotation about the pole).

Usage: python3 -I tesseral_cowell.py --set mech|drag [--passes 40] [--results-root DIR]
  mech: 10 Mars orbits (the repo's Odyssey-like a, e at pass 1), no drag, no Sun; per orbit the variants
        zonal (20x0), full (20x20), zonal + order m (m = 1..20), zonal + the tesserals of degree n (n = 2..20).  -> 410 trajectories
  drag: 6 orbits with the single-layer exponential atmosphere of the Mars() struct and constant beta (the Stage C truth drag model,
        no Sun), variants zonal and full, up to 60 d (decay cut 80 km): the orbit period shrinks, so nu = T/T_rot sweeps through commensurabilities.
All trajectories of a set are integrated together (one state vector, shared steps: differences between variants carry no
step-selection noise). Writes <root>/tesseral/cowell_<set>.csv: one row per trajectory and apoapsis (osculating elements) with the
time, rotation angle and position of the following periapsis.
"""
import argparse, math, pathlib, sys, time
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import tesseral_lib as tl

# SOURCED: src/models/environment/ephemerides/planets.jl (Mars struct): spin, exponential atmosphere, equatorial radius, flattening
W_ROT = 7.08823596e-5; RHO_REF = 8.7489231e-7; H_REF = 90.0e3; SCALE_H = 6.308278108e3; RE = 3.396190e6; FLAT = 1 - 3.3762e6 / 3.396190e6
BETA = 20.240070950010093      # SOURCED: <root>/mars/planet_constants.csv (Stage C derived beta, kg/m^2)
A0, E0 = 1.6065203772101851e7, 0.781859753293487   # SOURCED: Mars truth pass 1 of I0252_R000_W000_E1 (osculating a, e at periapsis)

MECH = [(25.2, 0, 0), (25.2, 90, 90), (45.0, 180, 45), (63.4, 270, 0), (63.4, 0, 270), (80.0, 90, 135),
        (93.522, 180, 0), (93.522, 270, 180), (116.6, 0, 90), (116.6, 90, 225)]          # (inc, raan, omega) deg, ASSUMED spread
DRAG = [(25.2, 0, 0), (45.0, 180, 45), (63.4, 270, 0), (93.522, 180, 0), (116.6, 0, 90), (80.0, 90, 135)]


def variants(cs, kind):
    C, S = cs["C"], cs["S"]; v = {"zonal": tl.select(C, S, lambda n, m: m == 0), "full": tl.select(C, S, lambda n, m: True)}
    if kind == "mech":
        for mm in range(1, cs["L"] + 1): v[f"m{mm}"] = tl.select(C, S, lambda n, m, mm=mm: m == 0 or m == mm)
        for nn in range(2, cs["L"] + 1): v[f"n{nn}"] = tl.select(C, S, lambda n, m, nn=nn: m == 0 or n == nn)
    return v


def run(kind, passes, root, rtol=1e-11):
    cs = tl.coefficients("Mars", 20); mu, R = cs["mu"], cs["R"]
    orbs = MECH if kind == "mech" else DRAG; V = variants(cs, kind); vn = list(V)
    traj = [(o, name) for o in range(len(orbs)) for name in vn]; K = len(traj)
    Cst = np.stack([V[name][0] for _, name in traj], axis=-1); Sst = np.stack([V[name][1] for _, name in traj], axis=-1)
    LU, TU = R, math.sqrt(R ** 3 / mu); VU = LU / TU
    y0 = np.zeros((6, K))
    for k, (o, _) in enumerate(traj):
        inc, raan, om = orbs[o]; s = tl.pqw_matrix(*np.radians([inc, raan, om]))
        p = A0 * (1 - E0 ** 2); r = p / (1 - E0)
        y0[:3, k] = s @ np.array([-r, 0.0, 0.0]) / LU; y0[3:, k] = s @ (math.sqrt(mu / p) * np.array([0.0, -(1 - E0), 0.0])) / VU
    drag = kind == "drag"; dead = np.zeros(K, bool)

    def rhs(t, yv):
        y = yv.reshape(6, K); r = y[:3] * LU; th = W_ROT * t * TU
        rb = tl.rot_z(-th * np.ones(K), r); ab = tl.accel_body(rb[0], rb[1], rb[2], cs, Cst, Sst); ai = tl.rot_z(th * np.ones(K), ab)
        rn = np.linalg.norm(r, axis=0); acc = -mu * r / rn ** 3 + ai
        if drag:
            v = y[3:] * VU; vrel = v - np.array([-W_ROT * r[1], W_ROT * r[0], 0 * r[0]]); vr = np.linalg.norm(vrel, axis=0)
            s2 = (r[2] / rn) ** 2; h = rn - RE * (1 - FLAT * s2)           # first-order ellipsoidal altitude (ASSUMED approximation)
            rho = RHO_REF * np.exp(-(h - H_REF) / SCALE_H); acc = acc - 0.5 * rho * vr * vrel / BETA
        out = np.concatenate([y[3:], acc * TU ** 2 / LU]); out[:, dead] = 0.0     # a decayed trajectory is frozen
        return out.ravel()

    T0 = 2 * math.pi * math.sqrt(A0 ** 3 / mu) / TU
    rows = {}; t = 0.0; y = y0.ravel(); seg = 0; nexts = 0; tstart = time.time(); count = np.zeros(K, int)
    for k in range(K): rows[(k, 0)] = _apo_row(traj[k], 0, 0.0, y0[:, k] * np.r_[[LU] * 3, [VU] * 3], mu)
    pend = {}
    while count[~dead].min(initial=passes) < passes and t * TU < 60.0 * 86400:
        def low(tt, yv):          # terminal: an alive trajectory reaches 90 km spherical altitude (ASSUMED decay cut; truth impact 50 km)
            r_ = np.linalg.norm(yv.reshape(6, K)[:3], axis=0) * LU - RE - 90e3
            return float(r_[~dead].min()) if (~dead).any() else 1.0
        low.terminal = True; low.direction = -1
        sol = solve_ivp(rhs, (t, t + 0.5 * T0), y, method="DOP853", rtol=rtol, atol=1e-13, dense_output=True, events=low if drag else None)
        if drag and sol.status == 1:
            r_ = np.linalg.norm(sol.y[:, -1].reshape(6, K)[:3], axis=0); r_[dead] = np.inf; dead[int(np.argmin(r_))] = True
        nexts += sol.nfev
        ts = np.linspace(t, sol.t[-1], 400); Y = sol.sol(ts).reshape(6, K, -1)
        rv = (Y[:3] * Y[3:]).sum(0)                    # (K, nt)
        for k in range(K):
            sgn = np.sign(rv[k]); idx = np.nonzero(sgn[:-1] != sgn[1:])[0]
            for i0 in idx:
                fk = lambda tt, k=k: float((sol.sol(tt).reshape(6, K)[:3, k] * sol.sol(tt).reshape(6, K)[3:, k]).sum())
                tc = brentq(fk, ts[i0], ts[i0 + 1], xtol=1e-13)
                st = sol.sol(tc).reshape(6, K)[:, k] * np.r_[[LU] * 3, [VU] * 3]
                if sgn[i0] < 0:            # periapsis
                    th = W_ROT * tc * TU; rb = tl.rot_z(np.array([-th]), st[:3, None])[:, 0]
                    pend[k] = dict(t_peri_s=tc * TU, theta_peri=th, r_peri_m=float(np.linalg.norm(st[:3])),
                                   lat_peri=math.asin(rb[2] / np.linalg.norm(rb)), lon_peri=math.atan2(rb[1], rb[0]))
                elif tc * TU > 1.0 and count[k] < passes:   # apoapsis after a periapsis
                    count[k] += 1; rows[(k, count[k])] = _apo_row(traj[k], count[k], tc * TU, st, mu)
                    rows[(k, count[k] - 1)].update(pend.pop(k, {}))
        t = sol.t[-1]; y = sol.y[:, -1]; seg += 1
        if seg % 10 == 0: print(f"seg {seg} t={t * TU / 86400:.2f} d passes min {count.min()} nfev {nexts} wall {time.time() - tstart:.0f} s", flush=True)
    df = pd.DataFrame(list(rows.values())); out = root / "tesseral"; out.mkdir(parents=True, exist_ok=True)
    meta = pd.DataFrame([dict(orbit=o, inc_deg=v[0], raan_deg=v[1], omega_deg=v[2]) for o, v in enumerate(orbs)])
    df = df.merge(meta, on="orbit"); df.to_csv(out / f"cowell_{kind}.csv", index=False)
    print("wrote", out / f"cowell_{kind}.csv", len(df), "rows; nfev", nexts, "wall", round(time.time() - tstart), "s (not a timing result)")


def _apo_row(tr, k, t_s, st, mu):
    a, e, i, O, w, f = tl.elements(st[:3], st[3:], mu)
    return dict(orbit=tr[0], variant=tr[1], **{"pass": k}, t_apo_s=t_s, a_m=a, e=e, i_rad=i, raan_rad=O, omega_rad=w, rp_m=a * (1 - e))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--set", required=True, choices=["mech", "drag"]); ap.add_argument("--passes", type=int, default=None)
    a = ap.parse_args(); run(a.set, a.passes or (40 if a.set == "mech" else 400), pathlib.Path(results_root(a.results_root)))
