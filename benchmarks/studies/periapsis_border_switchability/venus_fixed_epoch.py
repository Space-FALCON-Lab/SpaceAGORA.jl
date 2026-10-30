"""Venus fixed-epoch planner study (offline; no flights, no GRAM). New file only: imports geometry_rank (Mars, frozen) and geometry_rank_venus (Venus config, frozen) unchanged.

Question (user): with the initial epoch fixed, can the dimensionless perturbation parameters predict how the perturbations move the periapsis, so a planner can choose
(i, omega0, Omega0) that holds the periapsis altitude? Free initial conditions here: i, omega0, Omega0 only. The orbit size (ra0, h_ref), the campaign (ra_target, D_ref) and
the band stay the ranking model's values (the ranking code takes them from the config; they were not varied).

Commands (python3 -I venus_fixed_epoch.py ...):
  --selftest
  --truth OUT_DIR [--jobs N]     full ranking-model evaluation of the grid x epochs, the ridge scan and the VEX geometry at every epoch (the "truth": drift paths and DP floors)
  --analyze OUT_DIR              dimensionless parameters, predictors, ablation, out-of-sample tests, planner (reads OUT_DIR/truth.*)
Labels: every constant is the ranking model's (geometry_rank_venus_config.toml: SOURCED/DERIVED/ASSUMED there); nothing new is sourced here. See venus_fixed_epoch_results.md.
"""
import argparse, hashlib, json, math, pathlib, sys, time
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V
import stagek_tess as ST
import stagek_lib as K

EPOCHS = [float(x) for x in range(0, 221, 20)]          # days after the VEX anchor epoch 2014-05-19T04:00 (12 epochs; the Venus year is 224.7 d). ASSUMED set.
HEAD = [0.0, 60.0, 120.0, 180.0]                          # headline epochs (Sun phases 0, 96, 192, 288 deg); 0 = the Venus Express anchor epoch
INCS = [20.0, 39.2, 50.0, 63.4, 70.0, 80.0, 85.0, 89.876, 95.0, 100.0, 110.0, 116.6, 140.0]   # ASSUMED set (VEX/critical/Sun-plane values from the Venus config plus a spread)
OMEGAS = [float(x) for x in range(0, 360, 30)]
RAANS = [float(x) for x in range(0, 360, 45)]
NPAD = 80
H_KM = 3.86


# ------------------------------------------------------------------------------------------------ truth: the full ranking model

_W = {}


def make_grid():
    out = []
    for ep in EPOCHS:
        for inc in INCS:
            for om in OMEGAS:
                for ra in RAANS: out.append(dict(set="grid", inc=inc, epoch_off=ep, omega=om, raan=ra))
    for ep in HEAD:                                            # ridge scan: polar, omega0 every 10 deg, Omega0 every 30 deg (out-of-sample in omega0 and Omega0)
        for om in range(0, 360, 10):
            for ra in range(0, 360, 30): out.append(dict(set="ridge", inc=89.876, epoch_off=ep, omega=float(om), raan=float(ra)))
    for ep in EPOCHS: out.append(dict(set="vex", inc=89.876, epoch_off=ep, omega=75.505, raan=104.115))   # the Venus Express geometry at each epoch
    for i, c in enumerate(out): c["id"] = "F%05d" % i
    return out


def eval_truth(c):
    ctx = _W["ctx"]; r, n, S, d = G.evaluate(ctx, c, None, None)
    path = G.path_of(d); slow, _ = G.split_slow_fast(path, ctx.cfg["band"]["fast_window_passes"])
    nn = S["n"]; w = np.zeros(nn); w[1:] = S["H"][1:] / (ctx.ratio * np.array([G.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, nn)]))
    row = dict(id=c["id"], set=c["set"], epoch_off=c["epoch_off"], inc=c["inc"], omega0=c["omega"], raan0=c["raan"], n_passes=nn, q_req=r["q_req"], W_eff=r["W_eff"], floor=r["floor_mps"],
               floor_driftfree=r["floor_driftfree_mps"], R=r["R"], N=r["N"], regime=r["regime"], days=r["days"])
    for k in V.CH:
        if k in d: cum = np.concatenate([[0.0], np.cumsum(d[k])]); row["net_" + k] = float(cum[-1]); row["rng_" + k] = float(np.ptp(cum))
    pad = lambda x: np.concatenate([np.asarray(x, float), np.full(NPAD - len(x), np.nan)])
    sun = np.concatenate([[0.0], np.cumsum(d["third_body"])])
    arr = np.stack([pad(path), pad(slow), pad(sun), pad(S["a"]), pad(S["e"]), pad(S["t"]), pad(w)])
    return row, arr


def run_truth(out_dir, jobs):
    import multiprocessing as mp
    out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True); ctx, info = V.make_ctx(); _W["ctx"] = ctx; grid = make_grid(); t0 = time.time()
    print("campaign", json.dumps(info), "; candidates", len(grid), flush=True)
    with mp.get_context("fork").Pool(jobs) as pool: res = pool.map(eval_truth, grid, chunksize=8)
    df = pd.DataFrame([r for r, _ in res]); df.to_csv(out / "truth.csv", index=False)
    np.savez_compressed(out / "truth_paths.npz", arr=np.stack([a for _, a in res]), ids=df["id"].values.astype(str))
    (out / "truth_info.json").write_text(json.dumps(dict(campaign=info, n=len(df), wall_s=time.time() - t0, jobs=jobs, epochs=EPOCHS, incs=INCS), indent=1))
    print("done %d candidates in %.0f s" % (len(df), time.time() - t0))


# ------------------------------------------------------------------------------------------------ nominal schedule and Sun geometry

class Nominal:
    """Geometry-independent pass schedule (the closed-form decay law at the VEX anchor geometry and its q_req): a_k, e_k, t_k, H, V, lever weights, band. DERIVED from SOURCED constants
    (ranking config); it replaces the per-candidate schedule, whose geometry dependence is through f_rot and V_rel only (checked in --analyze)."""
    def __init__(self, ctx):
        cfg = ctx.cfg; anc = G.anchor_cand(cfg); self.ctx = ctx; b = cfg["band"]
        q = G.solve_q_req(ctx, anc, b["d_ref_days"]); S = G.schedule(ctx, anc, q); self.S = S; self.q_req = q; n = S["n"]; self.n = n
        self.a, self.e, self.t, self.H, self.V = S["a"], S["e"], S["t"], S["H"], S["V"]; self.rp = S["rp"]
        self.W_eff = b["slack_band_factor"] * math.log(b["q_cap_mean_w_cm2"] / q) - G.fast_width(cfg, ctx.ratio)
        self.w = np.zeros(n); self.w[1:] = S["H"][1:] / (ctx.ratio * np.array([G.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, n)]))
        piE = S["drop"] / (2 * S["a"]); self.lam = G.AP.lam_rp(S["e"], piE)               # drag descent Lambda per pass (scale heights), Pi_E = drop / (2 a)
        self.other = np.concatenate([[0.0], np.cumsum(self.lam[:-1] + 3 * np.diff(np.log(S["V"])))])   # geometry-independent drag + speed channels of the path (ln q)
        self.Ts = np.diff(self.t)


def sun_unit(ctx, nom, off):
    days = ctx.anchor_days + off + nom.t / 86400.0
    R = np.array([G.sun_in_mars_frame(d) for d in days]); rb = np.linalg.norm(R, axis=1); return R / rb[:, None], rb


def pq_vectors(inc, raan, om):
    """Unit vectors toward periapsis P, along the semilatus rectum Q, orbit normal W (arrays broadcast; last axis 3)."""
    ci, si, cO, sO, cw, sw = np.cos(inc), np.sin(inc), np.cos(raan), np.sin(raan), np.cos(om), np.sin(om)
    P = np.stack([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si], -1)
    Q = np.stack([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si], -1)
    W = np.stack([sO * si, -cO * si, ci], -1); return P, Q, W


def j2_rates(ctx, nom, inc):
    """Per-pass secular J2 apse rotation and node regression along the nominal schedule (same maps as G.schedule). inc: array (nc,)."""
    p = nom.a * (1 - nom.e ** 2); s = np.sin(inc)[:, None]
    dom = 3 * math.pi * ctx.J2 * (ctx.Rg / p) ** 2 * (2 - 2.5 * s ** 2); dOm = -3 * math.pi * ctx.J2 * (ctx.Rg / p) ** 2 * np.cos(inc)[:, None]
    return dom, dOm


def closed_form(ctx, nom, ep, inc, om0, Om0):
    """Closed-form Sun path and Sun-geometry parameters for arrays of candidates at one fixed epoch ep (days after the anchor). ALL from the dimensionless Sun tide:
    per-pass change of the periapsis height  dh_k = 15 pi a e sqrt(1-e^2) Pi_3B (u.P)(u.Q)   (singly averaged quadrupole, perturber frozen over an orbit; analytic_results.md section 3),
    Pi_3B,k = (mu_Sun / mu)(a_k / r_b,k)^3 (spec 3.1 item 3). Returns dict of arrays (nc, ...)."""
    u, rb = sun_unit(ctx, nom, ep); n = nom.n; nc = len(inc)
    dom, dOm = j2_rates(ctx, nom, inc)
    om = om0[:, None] + np.concatenate([np.zeros((nc, 1)), np.cumsum(dom[:, :-1], 1)], 1); Om = Om0[:, None] + np.concatenate([np.zeros((nc, 1)), np.cumsum(dOm[:, :-1], 1)], 1)
    P, Q, W = pq_vectors(inc[:, None] * np.ones((1, n)), Om, om)
    bP = (P * u[None]).sum(-1); bQ = (Q * u[None]).sum(-1); bW = (W * u[None]).sum(-1)
    Pi3 = G.MU_SUN / ctx.mu * (nom.a / rb) ** 3                                   # Pi_3B per pass
    scale = 15 * math.pi * nom.a * nom.e * np.sqrt(1 - nom.e ** 2) * Pi3 / nom.H  # A3_k: Sun kick scale in scale heights per orbit
    d3 = -(scale[None] * bP * bQ)[:, :-1]                                          # ln q increment of pass k -> k+1 (= -dh/H)
    path3 = np.concatenate([np.zeros((nc, 1)), np.cumsum(d3, 1)], 1)
    return dict(path3=path3, bP=bP, bQ=bQ, bW=bW, scale=scale, Pi3=Pi3, rb=rb, om=om, Om=Om)


def closed_form_evolving(ctx, nom, ep, inc, om0, Om0):
    """Second-order closed form: the Sun re-evolves the eccentricity vector e and angular-momentum vector j between passes (singly averaged quadrupole, perturber frozen over an orbit).
    DERIVED here from the averaged tidal potential  Phi = -(K a^2 / 4)[1 - 6 e^2 + 15 (e.u)^2 - 3 (j.u)^2],  K = mu_b / r_b^3  (Hamilton equations in (e, j)); per orbit, with 2 pi Pi_3B = 2 pi K / n^2:
        de = -(3/2)(2 pi Pi_3B)[ (j.u)(e x u) + 2 (j x e) - 5 (e.u)(j x u) ],   dj = -(3/2)(2 pi Pi_3B)[ (j.u)(j x u) - 5 (e.u)(e x u) ].
    To first order in the small quantities this is the frozen-element kick of closed_form (checked in --selftest). Midpoint (RK2) step; J2 apse rotation and node regression after each step as in the
    ranking model; drag acts at fixed r_p (a shrinks along the nominal schedule). Returns the Sun path (nc, n) in ln q (positive = periapsis descent)."""
    u, rb = sun_unit(ctx, nom, ep); n = nom.n; nc = len(inc); dom, dOm = j2_rates(ctx, nom, inc)
    ehat, _, What = pq_vectors(inc, Om0, om0); rp = np.full(nc, nom.rp[0]); path = np.zeros((nc, n)); zax = np.tile(np.array([0.0, 0.0, 1.0]), (nc, 1))
    dot = lambda x, y: (x * y).sum(1)[:, None]
    def rates(e, j, uk, Kt):
        eu, ju = dot(e, uk), dot(j, uk)
        return (-1.5 * Kt * (ju * np.cross(e, uk) + 2 * np.cross(j, e) - 5 * eu * np.cross(j, uk)), -1.5 * Kt * (ju * np.cross(j, uk) - 5 * eu * np.cross(e, uk)))
    def rot(v, ax, ang):
        c, s = np.cos(ang)[:, None], np.sin(ang)[:, None]; return v * c + np.cross(ax, v) * s + ax * dot(ax, v) * (1 - c)
    for k in range(n - 1):
        a = nom.a[k]; e = 1 - rp / a; e_ = e[:, None] * ehat; j_ = np.sqrt(1 - e ** 2)[:, None] * What; uk = u[k][None]; Kt = 2 * math.pi * G.MU_SUN / ctx.mu * (a / rb[k]) ** 3
        d1e, d1j = rates(e_, j_, uk, Kt); d2e, d2j = rates(e_ + 0.5 * d1e, j_ + 0.5 * d1j, uk, Kt); e2, j2 = e_ + d2e, j_ + d2j
        en = np.linalg.norm(e2, axis=1); ehat = e2 / en[:, None]; What = j2 / np.linalg.norm(j2, axis=1)[:, None]; What = What - dot(What, ehat) * ehat; What /= np.linalg.norm(What, axis=1)[:, None]
        rp_new = a * (1 - en); path[:, k + 1] = path[:, k] - (rp_new - rp) / nom.H[k]; rp = rp_new
        ehat = rot(ehat, What, dom[:, k]); ehat = rot(ehat, zax, dOm[:, k]); What = rot(What, zax, dOm[:, k])
    return path


def dp_floor(path, w, W):
    return G.tube_floor(path, w, W, 801)


def window(path):
    return uniform_filter1d(path, 15, mode="nearest", axis=-1)


def features(ctx, nom, ep, inc, om0, Om0, with_floor=True, nums=True):
    """The dimensionless parameter set (names follow pi_policy_interface_spec.md / the framework spec) for arrays of candidates at one epoch."""
    cf = closed_form(ctx, nom, ep, inc, om0, Om0); nc = len(inc); n = nom.n
    path = cf["path3"] + nom.other[None]; slow = window(path)
    F = {}
    F["Pi_3B0"] = np.full(nc, cf["Pi3"][0]); F["A3"] = np.full(nc, cf["scale"][0])                                # epoch-only scale
    p0 = nom.a[0] * (1 - nom.e[0] ** 2)
    F["Pi_J2"] = np.full(nc, ctx.J2 * (ctx.Rg / p0) ** 2); F["Gamma_3B_J2"] = F["Pi_3B0"] / F["Pi_J2"]               # Gamma_3B/J2 (item 19)
    F["Pi_H"] = np.full(nc, nom.H[0] / nom.rp[0]); F["Pi_f"] = np.full(nc, ctx.f * ctx.Req / nom.H[0]); F["Lam0"] = np.full(nc, nom.lam[0])   # Pi_f = f R / H = 0 at Venus
    bP0, bQ0, bW0 = cf["bP"][:, 0], cf["bQ"][:, 0], cf["bW"][:, 0]
    F["bP0"], F["bQ0"], F["bW0"] = bP0, bQ0, bW0; F["X3_0"] = cf["scale"][0] * bP0 * bQ0                            # signed Sun kick per orbit at the start, in H (item 21 for the tide)
    F["cos2beta0"] = 1 - bW0 ** 2; psi = np.unwrap(np.arctan2(cf["bQ"], cf["bP"]), axis=1)
    F["sin2psi0"], F["cos2psi0"] = np.sin(2 * psi[:, 0]), np.cos(2 * psi[:, 0])
    F["dpsi"] = psi[:, -1] - psi[:, 0]; F["bW_end"] = cf["bW"][:, -1]; F["sin2psi_mid"] = np.sin(2 * psi[:, n // 2]); F["cos2psi_mid"] = np.cos(2 * psi[:, n // 2]); F["bW_mid"] = cf["bW"][:, n // 2]
    F["N3_cf"] = cf["path3"][:, -1]; F["R3_cf"] = np.ptp(window(cf["path3"]), axis=1)                                   # closed-form Sun net drift and range (ln q)
    F["N_cf"] = slow[:, -1] - slow[:, 0]; F["R_cf"] = np.ptp(slow, axis=1)
    F["Rraw_cf"] = np.ptp(path, axis=1); F["maxabs_cf"] = np.abs(slow - slow[:, :1]).max(1)
    F["Sun_signchg"] = (np.diff(np.sign(np.diff(slow, axis=1)), axis=1) != 0).sum(1).astype(float)
    # J2 apse rotation per orbit and the Sun-tide apse rate ratio is not used: J2 does not move r_p at first order (stagek results); keep the per-orbit rotation as the J2 parameter
    dom, _ = j2_rates(ctx, nom, inc); F["dom_J2"] = dom[:, 0]
    F["_slow"] = slow
    return F, cf


# the numeric start-state parameters (zonal J3+, tesseral) are computed per candidate, with the correct Venus rotation rate (the ranking code's tess_kick uses the Mars W_ROT)


def numeric_params(ctx, nom, ep, inc, om0, Om0):
    T0 = nom.Ts[0]; a0, e0, H0 = nom.a[0], nom.e[0], nom.H[0]; day0 = ctx.anchor_days + ep + nom.t[0] / 86400.0; th0 = ctx.theta(day0)
    old = ST.W_ROT; ST.W_ROT = ctx.w
    try:
        z = np.array([a0 * K.zonal_secular(ctx.P, a0, e0, math.radians(i), math.radians(w), nlo=3, nhi=6, nq=256)[0] * T0 / H0 for i, w in zip(inc, om0)])
        t = np.array([-ST.tess_kick(ctx.cs, ctx.Ct, ctx.St, a0, e0, math.radians(i), math.radians(O), math.radians(w), th0, nf=64) / H0 for i, w, O in zip(inc, om0, Om0)])
    finally: ST.W_ROT = old
    return z, t


def build_feature_table(ctx, nom, df):
    rows = []
    for ep, g in df.groupby("epoch_off"):
        inc = np.radians(g.inc.values); om0 = np.radians(g.omega0.values); Om0 = np.radians(g.raan0.values)
        F, cf = features(ctx, nom, ep, inc, om0, Om0)
        # DP floor on the closed-form path (parameter-only): nominal weights and band
        slow = F.pop("_slow"); F["floor_cf"] = np.array([dp_floor(s, nom.w, nom.W_eff) for s in slow])
        z, t = numeric_params(ctx, nom, ep, g.inc.values, g.omega0.values, g.raan0.values); F["X_zon0"] = z; F["X_T0"] = t
        F["rho_zon"] = np.abs(z) / F["A3"]; F["rho_T"] = np.abs(t) / F["A3"]; F["rho_drag"] = F["Lam0"] / F["A3"]
        F["cosi"], F["sini"] = np.cos(inc), np.sin(inc); F["cos_w"], F["sin_w"] = np.cos(om0), np.sin(om0); F["cos2w"], F["sin2w"] = np.cos(2 * om0), np.sin(2 * om0); F["cos_O"], F["sin_O"] = np.cos(Om0), np.sin(Om0)
        f = pd.DataFrame(F, index=g.index); rows.append(f)
    return pd.concat(rows).loc[df.index]


# ------------------------------------------------------------------------------------------------ small gradient-boosted trees (numpy only; no sklearn on this host)

class GBM:
    def __init__(self, n_trees=150, depth=3, lr=0.1, min_leaf=25, nb=32): self.nt, self.d, self.lr, self.ml, self.nb = n_trees, depth, lr, min_leaf, nb

    def _bins(self, X):
        self.edges = [np.unique(np.quantile(X[:, j], np.linspace(0, 1, self.nb + 1)[1:-1])) for j in range(X.shape[1])]
        return np.stack([np.searchsorted(self.edges[j], X[:, j], side="right") for j in range(X.shape[1])], 1)

    def _tree(self, B, r, idx, depth):
        if depth == 0 or len(idx) < 2 * self.ml: return ("leaf", r[idx].mean())
        best = (0.0, None); rs = r[idx]; tot = rs.sum(); n = len(idx)
        for j in range(B.shape[1]):
            b = B[idx, j]; K_ = len(self.edges[j]) + 1
            if K_ < 2: continue
            cnt = np.bincount(b, minlength=K_).astype(float); sm = np.bincount(b, weights=rs, minlength=K_); cl = np.cumsum(cnt)[:-1]; sl = np.cumsum(sm)[:-1]
            ok = (cl >= self.ml) & (n - cl >= self.ml)
            if not ok.any(): continue
            gain = np.where(ok, sl ** 2 / np.maximum(cl, 1) + (tot - sl) ** 2 / np.maximum(n - cl, 1) - tot ** 2 / n, -1)
            k = int(np.argmax(gain))
            if gain[k] > best[0]: best = (gain[k], (j, k))
        if best[1] is None: return ("leaf", rs.mean())
        j, k = best[1]; m = B[idx, j] <= k
        return ("split", j, k, self._tree(B, r, idx[m], depth - 1), self._tree(B, r, idx[~m], depth - 1))

    def _pred(self, t, B):
        if t[0] == "leaf": return np.full(len(B), t[1])
        _, j, k, l, r = t; m = B[:, j] <= k; out = np.empty(len(B)); out[m] = self._pred(l, B[m]); out[~m] = self._pred(r, B[~m]); return out

    def fit(self, X, y):
        X = np.asarray(X, float); B = self._bins(X); self.f0 = y.mean(); pred = np.full(len(y), self.f0); self.trees = []
        for _ in range(self.nt):
            t = self._tree(B, y - pred, np.arange(len(y)), self.d); self.trees.append(t); pred += self.lr * self._pred(t, B)
        return self

    def predict(self, X):
        X = np.asarray(X, float); B = np.stack([np.searchsorted(self.edges[j], X[:, j], side="right") for j in range(X.shape[1])], 1); p = np.full(len(X), self.f0)
        for t in self.trees: p += self.lr * self._pred(t, B)
        return p


def fy(floor): return np.log10(0.1 + np.asarray(floor))
def fy_inv(y): return 10 ** np.asarray(y) - 0.1


def r2(y, p): return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
def mae(y, p): return float(np.abs(y - p).mean())


# ------------------------------------------------------------------------------------------------ selftest

def selftest():
    ctx, info = V.make_ctx(); nom = Nominal(ctx); ok = True
    def chk(name, c, s=""):
        nonlocal ok; ok &= bool(c); print(("PASS " if c else "FAIL ") + name, s)
    # (1) closed-form per-orbit Sun kick against the model's Gauss-quadrature kick (same state)
    rng = np.random.default_rng(1); ratios = []; a, e = nom.a[0], nom.e[0]
    Rs = G.sun_in_mars_frame(ctx.anchor_days); rb = np.linalg.norm(Rs); u = Rs / rb; Pi3 = G.MU_SUN / ctx.mu * (a / rb) ** 3
    cfk, num = [], []
    for _ in range(60):
        i, O, w = rng.uniform(0.3, 2.8), rng.uniform(0, 6.28), rng.uniform(0, 6.28); P, Q, W = pq_vectors(np.array(i), np.array(O), np.array(w))
        cfk.append(15 * math.pi * a * e * math.sqrt(1 - e * e) * Pi3 * float(P @ u) * float(Q @ u)); num.append(G.sun_kick(ctx, a, e, i, O, w, Rs))
    cfk, num = np.array(cfk), np.array(num); sl = (cfk * num).sum() / (cfk * cfk).sum(); corr = np.corrcoef(cfk, num)[0, 1]
    chk("closed-form Sun kick = model Gauss kick (sign, slope within 20%, corr > 0.98)", sl > 0 and abs(sl - 1) < 0.2 and corr > 0.98, "slope %.3f corr %.4f (number > 1 means the quadrature kick is larger)" % (sl, corr))
    # (2) DP floor of a truth path reproduces the model floor
    c = dict(id="T", inc=89.876, epoch_off=0.0, omega=0.0, raan=0.0); r, n, S, d = G.evaluate(ctx, c, None, None)
    path = G.path_of(d); slow, _ = G.split_slow_fast(path, 15); w = np.zeros(S["n"]); w[1:] = S["H"][1:] / (ctx.ratio * np.array([G.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, S["n"])]))
    chk("DP on the saved path/weights/band reproduces the model floor", abs(dp_floor(slow, w, r["W_eff"]) - r["floor_mps"]) < 1e-9, "%.4f vs %.4f" % (dp_floor(slow, w, r["W_eff"]), r["floor_mps"]))
    # (3) feature path vs model path (same candidate)
    F, cf = features(ctx, nom, 0.0, np.array([math.radians(89.876)]), np.array([0.0]), np.array([0.0]))
    chk("closed-form total path within 25% of the model range (anchor-epoch polar, omega0 = Omega0 = 0)", abs(F["R_cf"][0] / r["R"] - 1) < 0.25 or abs(F["R_cf"][0] - r["R"]) < 0.15, "R_cf %.3f model %.3f" % (F["R_cf"][0], r["R"]))
    chk("Venus flattening zero (Pi_f = 0)", F["Pi_f"][0] == 0.0)
    inc, om, Om = np.array([math.radians(89.876)]), np.array([0.3]), np.array([1.0]); k1 = E_ev = closed_form_evolving(ctx, nom, 0.0, inc, om, Om)[0, 1]; k0 = closed_form(ctx, nom, 0.0, inc, om, Om)["path3"][0, 1]
    chk("evolving (vector) form reproduces the frozen-element kick of the first orbit within 3% (derivation sign and factor)", abs(k1 / k0 - 1) < 0.03, "evolving %.5f frozen %.5f ln q" % (k1, k0))
    # (4) GBM learns a toy function
    X = rng.normal(size=(800, 4)); y = np.sin(X[:, 0]) + 0.5 * (X[:, 1] > 0); g = GBM(80).fit(X[:600], y[:600]); chk("GBM toy fit", r2(y[600:], g.predict(X[600:])) > 0.9, "R2 %.3f" % r2(y[600:], g.predict(X[600:])))
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--truth"); ap.add_argument("--analyze"); ap.add_argument("--jobs", type=int, default=6); a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    if a.truth: run_truth(a.truth, a.jobs)
    if a.analyze:
        import venus_fixed_epoch_analysis as A; A.main(a.analyze)
