"""Venus analogue of geometry_rank.py (WP-G1): rank candidate Venus aerobraking geometries by the dV floor. Python only; no flights, no GRAM sampling.

The Mars script and config are untouched; this file imports geometry_rank as G and replaces only the planet-specific parts: the context (Venus gravity, rotation, vehicle), the engine-sourced
Sun/pole/rotation table (sun_table_venus.jl), the reference campaign (Venus Express flown record) and the ranking driver (no mean-field input: the Venus field term T10 is absent).
Constants and labels: geometry_rank_venus_config.toml.

Usage (python3 -I geometry_rank_venus.py ...):
  --selftest                          table, frame, lever, Sun-kick and campaign checks
  --campaign                          print the derived reference campaign (ra_target, D_ref, anchor q_req)
  --rank OUT_DIR [--fine] [--cut]     rank the grid (default: the 256 cut grid); --fine = the 19-inclination follow-up grid (with --cut: omega0/Omega0 every 90 deg)
  --mars-control OUT_DIR              the same driver on the Mars config without the mean field (control: how much of the Mars result survives without T10)
"""
import argparse, hashlib, json, math, os, pathlib, sys, time
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank as G

CFG_PATH = HERE / "geometry_rank_venus_config.toml"
REPO = HERE.parents[2]
NOTES = "tasks/2026-10-15_periapsis-region-border-targeting"


def notes_dir():
    env = os.environ.get("SPACEAGORA_NOTES_DIR")
    if env: return pathlib.Path(env)
    import subprocess
    common = pathlib.Path(subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--path-format=absolute", "--git-common-dir"], text=True).strip())
    return common.parent / NOTES


def load_cfg(path=CFG_PATH):
    cfg = G.load_cfg(path); cfg["mars"] = cfg["venus"]   # geometry_rank.channels reads cfg["mars"]
    return cfg


# ------------------------------------------------------------------------------------------------ engine table (Sun, pole, rotation angle)

def load_table(cfg):
    """Fill geometry_rank._SUN from the Venus table (sha256 checked against its manifest). theta: rotation angle of the Venus-fixed x axis from the node of the J2000 equator on the Venus equator."""
    env = os.environ.get("SPACEAGORA_VENUS_SUN_TABLE"); d = pathlib.Path(env) if env else notes_dir() / cfg["sun"]["table_dir"]
    if not (d / "sun_table.csv").is_file(): raise FileNotFoundError("Venus sun_table.csv not found in %s (run sun_table_venus.jl)" % d)
    for line in (d / "manifest.sha256").read_text().splitlines():
        h, f = line.split(); assert hashlib.sha256((d / f).read_bytes()).hexdigest() == h, "sha256 mismatch for " + f
    t = pd.read_csv(d / "sun_table.csv"); x = t.days_j2000_utc.values
    z = t[["pole_x", "pole_y", "pole_z"]].values; p = t[["pm_x", "pm_y", "pm_z"]].values
    n = np.cross(np.array([0, 0, 1.0]), z); n /= np.linalg.norm(n, axis=1)[:, None]          # node = (-z_y, z_x, 0)
    th = np.unwrap(np.arctan2((np.cross(n, p) * z).sum(1), (n * p).sum(1)))
    G._SUN.clear()
    G._SUN.update(spl=CubicSpline(x, t[["sun_x_m", "sun_y_m", "sun_z_m"]].values), pole=CubicSpline(x, z), theta=CubicSpline(x, th), lo=x[0], hi=x[-1], dir=d, check=d / "sun_check.csv")
    return G._SUN


class VCtx(G.Ctx):
    """Venus constants of one run (same attributes the Mars code reads from Ctx)."""
    def __init__(self, cfg, ratio=None):
        m = cfg["venus"]; v = cfg["vehicle"]; self.cfg = cfg
        self.ratio = float(cfg["floor"]["burn_response_ratio"] if ratio is None else ratio)
        self.P = G.K.planet("Venus", max(20, int(m["j_zonal_max_degree"])))
        self.mu, self.Rg, self.J2 = self.P["mu"], self.P["R"], self.P["J"][2]
        self.Req, self.Rpol = float(m["r_eq_m"]), float(m["r_pol_m"]); self.f = 1 - self.Rpol / self.Req
        self.w = float(m["w_rot_rad_s"])
        area = v["bus_dims_m"][0] * v["bus_dims_m"][2] + 2 * v["panel_dims_m"][1] * v["panel_dims_m"][2]
        self.beta = (v["bus_mass_kg"] + 2 * v["panel_mass_each_kg"] + v["prop_mass_kg"]) / (v["cd"] * area)
        cs = G.TL.coefficients("Venus", int(m["tesseral_degree"])); self.cs = cs; self.Ct, self.St = G.TL.select(cs["C"], cs["S"], lambda n, mm: mm >= 1)
        a = cfg["anchor"]; self.anchor_days = (pd.Timestamp(a["epoch_utc"]) - G.J2000_UTC).total_seconds() / 86400.0
        self.s_t_ody = None; self.sun = bool(cfg["floor"].get("include_sun", False))

    def theta(self, days):
        return float(G._SUN["theta"](days)) % (2 * math.pi)


# ------------------------------------------------------------------------------------------------ reference campaign

def flown_campaign(ctx):
    """Model integration of the Venus Express flown record (passes 2950 to 3005) at the anchor: returns (ra_target, D_ref days, table of q per pass)."""
    cfg = ctx.cfg; c = cfg["campaign"]; o = cfg["orbit"]; tbl = pd.read_csv(REPO / c["pass_table"]).set_index("pass")
    passes = np.arange(int(tbl.index.min()), int(tbl.index.max()) + 1); inc = math.radians(cfg["anchor"]["inclination_deg"]); lat = 0.0
    h = np.interp(passes, tbl.index.values, tbl.h_peri_km.values) * 1e3; pdyn = np.interp(passes, tbl.index.values, tbl.pdyn_Nm2.values)   # missing pass: linear in the neighbors (ASSUMED)
    ra = o["ra0_m"]; t = None; ts = []; qs = []
    for k in range(len(passes)):
        rp = ctx.Req + h[k]; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
        if t is None: t = 0.5 * 2 * math.pi * math.sqrt(a ** 3 / ctx.mu)
        V, frot = G.vrel_parts(ctx, rp, e, inc, lat); q = pdyn[k] * V
        drop = cfg["decay"]["decay_scale"] * G.M.drop_per_pass(a, e, rp, q, G.M.k_theory(cfg["decay"]["h_fallback_m"], ctx.beta, frot), ctx.mu)
        ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); T = 2 * math.pi * math.sqrt(amid ** 3 / ctx.mu)
        ts.append(t); qs.append(q / 1e4); ra = ra_next; t += T
    S = dict(t=np.array(ts), n=len(ts))
    return ra, G.campaign_days(S), np.array(qs)


def derive_campaign(ctx):
    cfg = ctx.cfg; c = cfg["campaign"]
    if c["mode"] == "flown":
        ra_t, d_ref, q = flown_campaign(ctx); info = dict(mode="flown", q_mean_w_cm2=float(q.mean()), q_max_w_cm2=float(q.max()), n_passes=len(q))
    elif c["mode"] == "constant_q":
        cfg["orbit"]["ra_target_m"] = 1.0; cfg["band"]["d_ref_days"] = 1.0
        S = G.schedule(ctx, G.anchor_cand(cfg), c["q_ref_w_cm2"]); k = int(np.searchsorted(S["t"] / 86400.0, c["horizon_days"]))
        ra_t = S["ra"][k] - S["drop"][k]; d_ref = c["horizon_days"]; info = dict(mode="constant_q", q_ref=c["q_ref_w_cm2"])
    else: raise ValueError(c["mode"])
    cfg["orbit"]["ra_target_m"] = float(ra_t); cfg["band"]["d_ref_days"] = float(d_ref)
    info.update(ra_target_m=float(ra_t), d_ref_days=float(d_ref)); return info


def make_ctx(cfg_path=CFG_PATH, **over):
    cfg = load_cfg(cfg_path)
    for sec, kv in over.items():
        for k, v in kv.items(): cfg[sec][k] = v
    cfg["mars"] = cfg["venus"]; load_table(cfg); ctx = VCtx(cfg); info = derive_campaign(ctx); return ctx, info


# ------------------------------------------------------------------------------------------------ ranking driver (no mean field)

CH = ["geometry", "drag", "zonal", "tesseral", "third_body", "speed"]
_W = {}


def _eval_i(i):
    ctx, c = _W["ctx"], _W["grid"][i]; r, n, S, d = G.evaluate(ctx, c, None, ctx.s_t_ody)
    for k in CH:
        if k in d: cum = np.concatenate([[0.0], np.cumsum(d[k])]); r["net_" + k] = float(cum[-1]); r["rng_" + k] = float(np.ptp(cum)); r["rms_inc_" + k] = float(np.sqrt(np.mean(d[k] ** 2)))
    return r


def _ens_i(t):
    i, sun = t; return t, G.jitter_floors(_W["ctx"], _W["grid"][i], None, None, sun)


def _pmap(fn, items, jobs):
    import multiprocessing as mp
    if jobs <= 1: return [fn(x) for x in items]
    with mp.get_context("fork").Pool(jobs) as pool: return list(pool.imap(fn, items, chunksize=4))


def rank(ctx, out_dir, fine=False, cut=True, jobs=4, ensemble=True, tag="", scale_tess=False):
    cfg = ctx.cfg; out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    g = cfg["grid"]; grid = G.make_grid(cfg, g["cut256_omega_step_deg"] if cut else None, g["cut256_raan_step_deg"] if cut else None, fine=fine)
    anc = G.anchor_cand(cfg); _W.update(ctx=ctx, grid=grid + [anc])
    r_anc, *_ = G.evaluate(ctx, anc, None, None)
    if scale_tess: ctx.s_t_ody = r_anc["s_tesseral"]   # Mars control only: the tesseral share of the fast width scales from the anchor, as geometry_rank.anchor_ctx
    fl = G.jitter_floors(ctx, anc, None, None)
    r_anc.update(floor_det_mps=r_anc["floor_mps"], floor_ens_mps=float(fl.mean()), floor_ens_sd_mps=float(fl.std())); pd.DataFrame([r_anc]).to_csv(out / ("anchor%s.csv" % tag), index=False)
    print("anchor (VEX): q_req %.4f, W %.3f, W_eff %.3f, floor det %.4f ens %.4f +- %.4f, regime %s, passes %d, days %.2f" % (r_anc["q_req"], r_anc["W"], r_anc["W_eff"], r_anc["floor_mps"], fl.mean(), fl.std(), r_anc["regime"], r_anc["n_passes"], r_anc["days"]), flush=True)
    df = pd.DataFrame(_pmap(_eval_i, range(len(grid)), jobs)); print("stage 1: %d candidates, %.0f s" % (len(df), time.time() - t0), flush=True)
    df["floor_ens_mps"] = np.nan; df["floor_ens_sd_mps"] = np.nan
    if ensemble:
        order = df.sort_values("floor_mps"); sel = set(order.head(int(cfg["floor"]["refine_stat"])).index)
        for inc, gg in df.groupby("inc_deg"):
            gg = gg.sort_values("floor_mps"); sel |= {gg.index[0], gg.index[len(gg) // 2], gg.index[-1]}
        for (i, _), f in _pmap(_ens_i, [(i, None) for i in sorted(sel)], jobs): df.at[i, "floor_ens_mps"] = f.mean(); df.at[i, "floor_ens_sd_mps"] = f.std()
        print("stage 2: %d re-measured, %.0f s" % (len(sel), time.time() - t0), flush=True)
    df.to_csv(out / ("geometry_rank%s.csv" % tag), index=False)
    return df, r_anc


def window(ctx, out_dir, inc, omega, jobs):
    """Epoch (10 d) x Omega0 (15 deg) scan at one inclination and omega0: how wide is the low-floor window (fragility of the best geometries)."""
    out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    grid = [dict(id="W%04d" % i, inc=inc, epoch_off=float(e), omega=omega, raan=float(r)) for i, (e, r) in enumerate((e, r) for e in range(0, 301, 10) for r in range(0, 360, 15))]
    _W.update(ctx=ctx, grid=grid); df = pd.DataFrame(_pmap(_eval_i, range(len(grid)), jobs)); df.to_csv(out / ("window_i%g_w%g.csv" % (inc, omega)), index=False)
    p = df.pivot_table(index="epoch_off_d", columns="raan0_deg", values="floor_mps"); print(p.round(2).to_string()); print("share of the scan with floor < 0.05: %.3f; < 0.3: %.3f" % ((df.floor_mps < .05).mean(), (df.floor_mps < .3).mean()))


def decide(df, r_anc, cfg):
    f_ref = r_anc["floor_ens_mps"]; x = cfg["decision"]["x_floor_ratio"]
    top = df.dropna(subset=["floor_ens_mps"]).sort_values("floor_ens_mps").head(int(cfg["floor"]["refine_top"]))
    feas = top[top.feasible]
    if feas.empty: return dict(branch="geometry_arms_dropped", reason="no feasible candidate", f_ref=f_ref)
    b = feas.iloc[0]; d = dict(f_ref=float(f_ref), f_ref_sd=float(r_anc["floor_ens_sd_mps"]), f_ref_det=float(r_anc["floor_det_mps"]), f_best=float(b.floor_ens_mps), best_id=b["id"], x=x,
                               rule_passes=bool(b.floor_ens_mps <= x * f_ref), ratio=float(b.floor_ens_mps / f_ref) if f_ref > 0 else None)
    return d


# ------------------------------------------------------------------------------------------------ selftest and CLI

def selftest(ctx, info):
    ok = True
    def chk(name, cond, info_=""):
        nonlocal ok; ok &= bool(cond); print(("PASS " if cond else "FAIL ") + name, info_)
    cfg = ctx.cfg; St = G._SUN; m = cfg["venus"]
    ck = pd.read_csv(St["check"]); xs = ck.days_j2000_utc.values; eng = ck[["sun_x_m", "sun_y_m", "sun_z_m"]].values; sp = St["spl"](xs)
    ang = np.degrees(np.arccos(np.clip((eng * sp).sum(1) / np.linalg.norm(eng, axis=1) / np.linalg.norm(sp, axis=1), -1, 1))); dist = np.abs(np.linalg.norm(sp, axis=1) / np.linalg.norm(eng, axis=1) - 1)
    chk("Venus Sun table vs 300 off-grid engine epochs: angle < 1e-4 deg, distance < 1e-7", ang.max() < 1e-4 and dist.max() < 1e-7, "max %.2e deg, %.2e rel" % (ang.max(), dist.max()))
    pole = St["pole"](xs[:5]); ra = math.degrees(math.atan2(pole[0][1], pole[0][0])) % 360; dec = math.degrees(math.asin(pole[0][2] / np.linalg.norm(pole[0])))
    chk("engine pole is the pck00011 Venus pole (RA 272.76, Dec 67.16, 0.1 deg)", abs(ra - m["pole_ra_deg"]) < 0.1 and abs(dec - m["pole_dec_deg"]) < 0.1, "RA %.3f Dec %.3f" % (ra, dec))
    d0, d1 = St["lo"] + 10, St["hi"] - 10; rate = (St["theta"](d1) - St["theta"](d0)) / (d1 - d0)
    chk("engine rotation rate = pck00011 PM rate -1.4813688 deg/day (1e-3 deg/day)", abs(math.degrees(rate) - m["pm_rate_deg_day"]) < 1e-3, "%.5f deg/day" % math.degrees(rate))
    W_pck = (m["pm_w0_deg"] + m["pm_rate_deg_day"] * d0) % 360; th = math.degrees(St["theta"](d0)) % 360
    print("INFO engine theta at day %.0f: %.3f deg; pck00011 W: %.3f deg; difference %.3f deg (the tesseral channel is weak at Venus: 19 m rms in the Part 1 pairs)" % (d0, th, W_pck, ((th - W_pck + 180) % 360) - 180))
    chk("flattening zero (Rp_e = Rp_p), rotation retrograde", ctx.f == 0.0 and ctx.w < 0)
    a = 0.5 * (cfg["orbit"]["ra0_m"] + ctx.Req + cfg["orbit"]["h_ref_m"]); rp = ctx.Req + cfg["orbit"]["h_ref_m"]; e = (cfg["orbit"]["ra0_m"] - rp) / (cfg["orbit"]["ra0_m"] + rp)
    chk("lever at the start about 16.02 km per m/s (venus_gram_results.md section 6)", abs(G.lever(ctx, a, e) / 1e3 - 16.02) < 0.2, "%.3f km/(m/s)" % (G.lever(ctx, a, e) / 1e3))
    # Sun kick maximum over the angles against the closed form 15 pi a e sqrt(1-e^2) Pi_3B / 2 (u.P)(u.Q)max = 1/2 (analytic_results section 3)
    day = ctx.anchor_days; Rs = G.sun_in_mars_frame(day); rb = np.linalg.norm(Rs); Pi3 = G.MU_SUN / ctx.mu * (a / rb) ** 3; closed = 15 * math.pi * a * e * math.sqrt(1 - e * e) * Pi3 * 0.5
    best = max(abs(G.sun_kick(ctx, a, e, math.radians(i), math.radians(O), math.radians(w), Rs)) for i in (30, 60, 90, 120) for O in range(0, 360, 30) for w in range(0, 360, 30))
    chk("numeric Sun kick maximum within 15% of the singly averaged closed form (the form validated in Part 1, EV 0.976; venus_gram_results.md quotes 1.05 H from a formula with prefactor 2 pi x 1.875, half of this one: reported, not enforced)", 0.85 < best / closed < 1.15, "numeric %.0f m, closed %.0f m, %.2f H" % (best, closed, best / ctx.cfg["decay"]["h_fallback_m"]))
    om = G.M.apsidal_rotation_per_pass(a, e, math.radians(63.4), ctx.J2, ctx.Rg, ctx.mu); chk("J2 apse rotation is zero at 63.4 deg to 1e-4 rad per pass", abs(om) < 1e-4, "%.2e rad" % om)
    z6 = G.K.zonal_secular(ctx.P, a, e, math.radians(89.876), math.radians(75.505), nlo=3, nhi=6)[0]; z20 = G.K.zonal_secular(ctx.P, a, e, math.radians(89.876), math.radians(75.505), nlo=3, nhi=20)[0]
    print("INFO zonal de/dt degrees 3..6 vs 3..20 at the anchor: %.3e vs %.3e per s (j_zonal_max_degree check)" % (z6, z20))
    S = G.schedule(ctx, G.anchor_cand(cfg), 0.16); qr = G.solve_q_req(ctx, G.anchor_cand(cfg), cfg["band"]["d_ref_days"]); d = G.campaign_days(G.schedule(ctx, G.anchor_cand(cfg), qr))
    chk("anchor q_req reproduces D_ref to 0.5% (and the flown record's mean q is of the same order)", abs(d / cfg["band"]["d_ref_days"] - 1) < 5e-3, "q_req %.4f, days %.3f vs %.3f; flown mean q %.4f" % (qr, d, cfg["band"]["d_ref_days"], info.get("q_mean_w_cm2", float("nan"))))
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--campaign", action="store_true"); ap.add_argument("--rank"); ap.add_argument("--fine", action="store_true")
    ap.add_argument("--full", action="store_true", help="with --fine: omega0 30 / Omega0 45 instead of the 90/90 cut"); ap.add_argument("--no-ens", action="store_true"); ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--config", default=str(CFG_PATH)); ap.add_argument("--mars-control"); ap.add_argument("--no-sun", action="store_true"); ap.add_argument("--tag", default=""); ap.add_argument("--window", help="OUT_DIR; scan epoch x Omega0 at --inc/--omega"); ap.add_argument("--inc", type=float, default=89.876); ap.add_argument("--omega", type=float, default=0.0)
    ap.add_argument("--set", action="append", default=[], help="section.key=value override of the config (sensitivity runs)"); a = ap.parse_args()
    if a.mars_control:
        cfg = G.load_cfg(); ctx = G.Ctx(cfg); G._SUN.clear(); G.sun_table(); _W.update(ctx=ctx)
        df, r = rank(ctx, a.mars_control, a.fine, not a.full, a.jobs, not a.no_ens, a.tag, True); print(json.dumps(decide(df, r, cfg) if not a.no_ens else {}, indent=1)); return
    over = {}
    for s in a.set:
        k, v = s.split("="); sec, key = k.split("."); over.setdefault(sec, {})[key] = json.loads(v) if v[0] in "[0123456789-.tf\"" else v
    ctx, info = make_ctx(a.config, **over)
    if a.no_sun: ctx.sun = False
    print("campaign:", json.dumps(info))
    if a.selftest: sys.exit(0 if selftest(ctx, info) else 1)
    if a.window: window(ctx, a.window, a.inc, a.omega, a.jobs); return
    if a.rank:
        df, r = rank(ctx, a.rank, a.fine, not a.full, a.jobs, not a.no_ens, a.tag)
        dec = decide(df, r, ctx.cfg) if not a.no_ens else {}; dec["campaign"] = info; dec["overrides"] = over; dec["no_sun"] = bool(a.no_sun)
        (pathlib.Path(a.rank) / ("geometry_decision%s.json" % a.tag)).write_text(json.dumps(dec, indent=1)); print(json.dumps(dec, indent=1))


if __name__ == "__main__":
    main()
