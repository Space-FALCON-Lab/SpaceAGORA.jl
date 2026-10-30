"""WP-G1: rank candidate Mars aerobraking geometries by the dV floor (plan R2.3 grid, R3.1 to R3.6 metric, regimes, pre-check, decision rule). Python only.

Per candidate (omega0, Omega0, epoch offset, i): a burn-free slow path of ln Q over the main phase from closed forms (J2 apse and node rates; J3+ secular
de/dt; the T4a tesseral per-pass map of stagek_tess.tess_kick; drag descent Lambda H; velocity loss 3 dln V_rel; the oblate-surface term from the apse
rotation; and, when sampled, the GRAM mean-field term T10), then the exact weighted-total-variation tube problem (the dV floor) by dynamic programming.
Reuses stagek_lib, stagem_lib, apoapsis_pi_lib, tesseral_lib and stagek_tess (imports, no copies). Constants and labels: geometry_rank_config.toml.

Usage (python3 -I geometry_rank.py ...):
  --selftest                      closed-form, DP, classifier and (if the dev tables are found) dev-data checks
  --dry-run                       counts against plan R2.3, probe/job cost estimate; computes nothing heavy
  --oracle                        reproduce the dev oracle floor (plan R3.2/R3.6) from the Stage O dev tables
  --tracks OUT.csv                predicted periapsis tracks (cand, k, t, lat, lon) for meanfield_sample.jl (fallback scale height)
  --precheck [--meanfield DIR]    R3.3 pre-check: predicted Odyssey floor vs the dev oracle floor (no --meanfield: preview, field term 0)
  --rank OUT_DIR --meanfield DIR  rank all candidates, apply the R3.5 decision rule
"""
import argparse, glob, itertools, json, math, pathlib, subprocess, sys, time
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagek_lib as K
import stagem_lib as M
import apoapsis_pi_lib as AP
import tesseral_lib as TL
from stagek_tess import tess_kick

CFG_PATH = HERE / "geometry_rank_config.toml"
J2000_UTC = pd.Timestamp("2000-01-01T12:00:00")


def load_cfg(path=CFG_PATH):
    return tomllib.loads(pathlib.Path(path).read_text())


# ------------------------------------------------------------------------------------------------ context

class Ctx:
    """Constants of one run: planet, vehicle, gravity coefficients (built once)."""
    def __init__(self, cfg, ratio=None):
        self.cfg = cfg; m = cfg["mars"]; v = cfg["vehicle"]
        self.ratio = float(cfg["floor"]["burn_response_ratio"] if ratio is None else ratio)
        self.P = K.planet("Mars", max(20, int(m["j_zonal_max_degree"])))
        self.mu, self.Rg, self.J2 = self.P["mu"], self.P["R"], self.P["J"][2]
        self.Req, self.Rpol = float(m["r_eq_m"]), float(m["r_pol_m"]); self.f = 1 - self.Rpol / self.Req
        self.w = float(m["w_rot_rad_s"])
        area = v["bus_dims_m"][0] * v["bus_dims_m"][2] + 2 * v["panel_dims_m"][1] * v["panel_dims_m"][2]
        self.beta = (v["bus_mass_kg"] + 2 * v["panel_mass_each_kg"] + v["prop_mass_kg"]) / (v["cd"] * area)
        cs = TL.coefficients("Mars", int(m["tesseral_degree"]))
        self.cs = cs; self.Ct, self.St = TL.select(cs["C"], cs["S"], lambda n, mm: mm >= 1)
        a = cfg["anchor"]; self.anchor_days = (pd.Timestamp(a["epoch_utc"]) - J2000_UTC).total_seconds() / 86400.0
        self.s_t_ody = None   # tesseral fast rms of the anchor, filled by anchor_ctx
        self.sun = bool(cfg["floor"].get("include_sun", False))   # Sun third-body channel (validated against the dev tables, see geometry_rank_results.md addendum)

    def surface_radius(self, lat_gc):
        return self.Req * (1 - self.f * np.sin(lat_gc) ** 2)   # first order in f (stagem_lib.periapsis_geodetic_altitude form)

    def theta(self, t_days_since_j2000):
        """Rotation angle of the body frame from the inertial x axis (node of the J2000 equator on the Mars equator), rad."""
        m = self.cfg["mars"]
        return math.radians((m["pm_w0_deg"] + m["pm_rate_deg_day"] * t_days_since_j2000 + m["theta_offset_deg"]) % 360.0)


def lever(ctx, a, e):
    return K.lever_m_per_mps(ctx.mu, a, e)


# ------------------------------------------------------------------------------------------------ candidates and counts

def make_grid(cfg, omega_step=None, raan_step=None, fine=False):
    g = dict(cfg["grid"])
    if fine: g["inclinations_deg"] = cfg["grid_fine"]["inclinations_deg"]
    os_ = g["omega0_step_deg"] if omega_step is None else omega_step; rs = g["raan0_step_deg"] if raan_step is None else raan_step
    out = []
    for i, (inc, ep, om, ra) in enumerate(itertools.product(sorted(g["inclinations_deg"]), g["epoch_offsets_days"],
                                                             np.arange(0.0, 360.0, os_), np.arange(0.0, 360.0, rs))):
        out.append(dict(inc=float(inc), epoch_off=float(ep), omega=float(om), raan=float(ra)))
    for i, c in enumerate(out): c["id"] = "G%04d" % i
    return out


def anchor_cand(cfg):
    a = cfg["anchor"]
    return dict(id="ODY", inc=a["inclination_deg"], epoch_off=0.0, omega=a["omega_deg"], raan=a["raan_deg"])


def in_cut(cand, cfg, mode):
    g = cfg["grid"]; os_, rs = (g["cut256_omega_step_deg"], g["cut256_raan_step_deg"]) if mode == "256" else (g["cut384_omega_step_deg"], g["cut384_raan_step_deg"])
    return cand["id"] == "ODY" or (abs(cand["omega"] / os_ - round(cand["omega"] / os_)) < 1e-9 and abs(cand["raan"] / rs - round(cand["raan"] / rs)) < 1e-9)


# ------------------------------------------------------------------------------------------------ burn-free schedule (one candidate)

def vrel_parts(ctx, rp, e, inc, lat_gc):
    vp = math.sqrt(ctx.mu * (1 + e) / rp); hh = rp * vp
    v2 = vp * vp - 2 * ctx.w * hh * math.cos(inc) + (ctx.w * rp) ** 2 * (1 - math.sin(lat_gc) ** 2)
    frot = (vp * vp - ctx.w * hh * math.cos(inc)) / v2     # (v . v_rel) / V_rel^2
    return math.sqrt(v2), frot


def schedule(ctx, cand, q_w_cm2, H_of_k=None):
    """Main-phase pass sequence at constant mean-field heat rate q on the held tube (altitude h_ref over the surface). Returns dict of arrays."""
    cfg = ctx.cfg; o = cfg["orbit"]; inc = math.radians(cand["inc"]); om = math.radians(cand["omega"]); Om = math.radians(cand["raan"])
    ra = o["ra0_m"]; q = q_w_cm2 * 1.0e4; rows = []; t = None
    for k in range(int(o["max_passes"])):
        lat = math.asin(math.sin(inc) * math.sin(om)); rp = ctx.surface_radius(lat) + o["h_ref_m"]
        a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
        if t is None: t = 0.5 * 2 * math.pi * math.sqrt(a ** 3 / ctx.mu) + 0.0   # first periapsis half an orbit after the starting apoapsis
        H = cfg["decay"]["h_fallback_m"] if H_of_k is None else float(H_of_k[min(k, len(H_of_k) - 1)])
        V, frot = vrel_parts(ctx, rp, e, inc, lat)
        drop = cfg["decay"]["decay_scale"] * M.drop_per_pass(a, e, rp, q, M.k_theory(H, ctx.beta, frot), ctx.mu)
        rows.append((t, a, e, rp, ra, om, Om, lat, H, V, frot, drop))
        ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); T = 2 * math.pi * math.sqrt(amid ** 3 / ctx.mu)
        p = a * (1 - e * e); om += M.apsidal_rotation_per_pass(a, e, inc, ctx.J2, ctx.Rg, ctx.mu)
        Om += -3 * math.pi * ctx.J2 * (ctx.Rg / p) ** 2 * math.cos(inc)
        ra = ra_next; t += T
        if ra <= o["ra_target_m"]: break
    A = np.array(rows); names = "t a e rp ra omega raan lat H V frot drop".split()
    S = {n: A[:, i] for i, n in enumerate(names)}; S["inc"] = inc; S["n"] = len(A)
    return S


def campaign_days(S):
    return (S["t"][-1] + 0.5 * (S["t"][-1] - S["t"][-2] if S["n"] > 1 else 0.0)) / 86400.0   # last pass plus half a period (to the closing apoapsis)


def solve_q_req(ctx, cand, d_ref_days, H_of_k=None, q0=0.16):
    """Mean-field heat rate (W/cm2) at which the campaign takes d_ref_days (secant on ln q, time ~ 1/q)."""
    q = q0
    for _ in range(12):
        d = campaign_days(schedule(ctx, cand, q, H_of_k))
        if abs(d / d_ref_days - 1) < 1e-3: break
        q *= d / d_ref_days
    return q


def lon_lat(ctx, cand, S):
    """Body-fixed longitude and geodetic latitude of each periapsis (rad); inertial periapsis direction from (omega, Omega, i), minus the rotation angle."""
    inc = S["inc"]; w, O = S["omega"], S["raan"]
    Px = np.cos(O) * np.cos(w) - np.sin(O) * np.sin(w) * math.cos(inc); Py = np.sin(O) * np.cos(w) + np.cos(O) * np.sin(w) * math.cos(inc)
    lam = np.arctan2(Py, Px)
    days = ctx.anchor_days + cand["epoch_off"] + S["t"] / 86400.0
    th = np.array([ctx.theta(d) for d in days])
    lat_gd = np.arctan(np.tan(S["lat"]) / (1 - ctx.f) ** 2)
    return np.mod(lam - th, 2 * math.pi), lat_gd, th



# ------------------------------------------------------------------------------------------------ Sun third-body kick (robustness view)

MU_SUN = 1.32712440018e20      # m3/s2, SOURCED: IAU 2015 nominal solar GM
AU_M = 1.495978707e11
# Mars heliocentric mean elements (J2000 value, rate per Julian century), RECALLED from the JPL approximate-positions table (Standish); a, e, I, L, varpi, Omega in AU / deg
_MARS_EL = dict(a=(1.52371034, 0.00001847), e=(0.09339410, 0.00007882), I=(1.84969142, -0.00813131), L=(-4.55343205, 19140.30268499), w=(-23.94362959, 0.44441088), O=(49.55953891, -0.29257343))


def sun_vector_ecliptic(days_j2000):
    """Mean-element Sun (RECALLED constants): kept only to quantify how much the sourced table changed the result (selftest and addendum); not used by the model."""
    T = days_j2000 / 36525.0; el = {k: v[0] + v[1] * T for k, v in _MARS_EL.items()}
    a, e, I, L, w, O = el["a"], el["e"], math.radians(el["I"]), math.radians(el["L"]), math.radians(el["w"]), math.radians(el["O"])
    M = (L - w + math.pi) % (2 * math.pi) - math.pi; E = M
    for _ in range(30): E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E); om = w - O
    co, so, cO, sO, ci, si = math.cos(om), math.sin(om), math.cos(O), math.sin(O), math.cos(I), math.sin(I)
    r = np.array([(co * cO - so * sO * ci) * x + (-so * cO - co * sO * ci) * y, (co * sO + so * cO * ci) * x + (-so * sO + co * cO * ci) * y, (so * si) * x + (co * si) * y])
    eps = math.radians(23.43928); r = -r * AU_M
    return np.array([r[0], r[1] * math.cos(eps) - r[2] * math.sin(eps), r[1] * math.sin(eps) + r[2] * math.cos(eps)])   # ICRF-equatorial


_SUN = {}


def sun_table(cfg=None):
    """Engine-sourced Sun position (m, from Mars, J2000) and Mars pole (J2000) as cubic splines of UTC days since J2000, from sun_table.jl (SimulationModel's own SPICE call and frame transform).
    The table's sha256 is checked against its manifest; a missing table is an error (no unsourced fallback)."""
    if "spl" in _SUN: return _SUN
    import hashlib
    from scipy.interpolate import CubicSpline
    base = find_dev_dir()
    d = None if base is None else base.parent.parent / "analytic_data" / "sun_table"
    env = __import__("os").environ.get("SPACEAGORA_SUN_TABLE")
    d = pathlib.Path(env) if env else d
    if d is None or not (d / "sun_table.csv").is_file(): raise FileNotFoundError("sun_table.csv not found (run sun_table.jl; set SPACEAGORA_SUN_TABLE)")
    for line in (d / "manifest.sha256").read_text().splitlines():
        h, f = line.split(); assert hashlib.sha256((d / f).read_bytes()).hexdigest() == h, "sha256 mismatch for " + f
    t = pd.read_csv(d / "sun_table.csv"); x = t.days_j2000_utc.values
    _SUN.update(spl=CubicSpline(x, t[["sun_x_m", "sun_y_m", "sun_z_m"]].values), pole=CubicSpline(x, t[["pole_x", "pole_y", "pole_z"]].values), lo=x[0], hi=x[-1], dir=d, check=d / "sun_check.csv")
    return _SUN


def mars_frame(days):
    """Rows (x, y, z) of the Mars-equator inertial frame in J2000: z = the engine's pole (planet_frame_lpi), x = node of the J2000 equator on the Mars equator (as mars_rot)."""
    z = sun_table()["pole"](days); z = z / np.linalg.norm(z); x = np.array([-z[1], z[0], 0.0]); x /= np.linalg.norm(x)
    return np.array([x, np.cross(z, x), z])


def sun_in_mars_frame(days_j2000):
    S = sun_table()
    if not (S["lo"] <= days_j2000 <= S["hi"]): raise ValueError("epoch outside the Sun table")
    return mars_frame(days_j2000) @ S["spl"](days_j2000)


def sun_kick(ctx, a, e, inc, raan, om, Rsun, nf=64):
    """First-order periapsis-radius change over one apoapsis-to-apoapsis arc from the Sun's differential acceleration (Sun fixed over the orbit), Gauss quadrature as tess_kick."""
    f, w, tau, r, v, rr, h = TL.kepler_arc(a, e, inc, raan, om, ctx.mu, nf)
    d = Rsun[:, None] - r; ai = MU_SUN * (d / np.linalg.norm(d, axis=0) ** 3 - (Rsun / np.linalg.norm(Rsun) ** 3)[:, None])
    rhat = r / rr; hv = np.cross(r.T, v.T).T; hhat = hv / np.linalg.norm(hv, axis=0); shat = np.cross(hhat.T, rhat.T).T
    rates = TL.gauss_rates(a, e, inc, om, f, rr, h, (ai * rhat).sum(0), (ai * shat).sum(0), (ai * hhat).sum(0), ctx.mu)
    return float(rates[5] @ (w * rr ** 2 / h))

# ------------------------------------------------------------------------------------------------ natural slow path

def channels(ctx, cand, S, field=None):
    """Per-pass increments of ln Q of the burn-free path (pass k -> k+1), by channel; arrays of length n-1. field: dict(lnrho0) per pass or None."""
    cfg = ctx.cfg; n = S["n"]; a, e, H, inc = S["a"], S["e"], S["H"], S["inc"]
    d = {}
    Rs = ctx.surface_radius(S["lat"]); d["geometry"] = (Rs[1:] - Rs[:-1]) / H[:-1]                     # h_p = r_p - R_s: surface moves, radius does not
    piE = S["drop"] / (2 * a); d["drag"] = np.array([AP.lam_rp(e[k], piE[k]) for k in range(n)])[:-1]    # drag descent Lambda H per pass
    jmax = int(cfg["mars"]["j_zonal_max_degree"]); T = np.diff(S["t"])
    dZ = np.array([K.zonal_secular(ctx.P, a[k], e[k], inc, S["omega"][k], nlo=3, nhi=jmax, nq=256)[0] * T[k] for k in range(n - 1)])
    d["zonal"] = a[:-1] * dZ / H[:-1]                                                                   # d r_p = -a d e
    _, _, th = lon_lat(ctx, cand, S)
    dT = np.array([tess_kick(ctx.cs, ctx.Ct, ctx.St, a[k], e[k], inc, S["raan"][k], S["omega"][k], th[k], nf=int(cfg["mars"]["tesseral_nodes"])) for k in range(n - 1)])
    d["tesseral"] = -dT / H[:-1]
    if ctx.sun:
        dd = ctx.anchor_days + cand["epoch_off"] + S["t"] / 86400.0
        dS = np.array([sun_kick(ctx, a[k], e[k], inc, S["raan"][k], S["omega"][k], sun_in_mars_frame(dd[k])) for k in range(n - 1)]); d["third_body"] = -dS / H[:-1]
    d["speed"] = 3 * np.diff(np.log(S["V"]))
    d["field"] = np.zeros(n - 1) if field is None else np.diff(np.asarray(field["lnrho0"])[np.minimum(np.arange(n), len(field["lnrho0"]) - 1)])
    return d


def path_of(d, keys=None):
    keys = keys or list(d)
    return np.concatenate([[0.0], np.cumsum(sum(d[k] for k in keys))])


def split_slow_fast(path, win):
    slow = uniform_filter1d(path, int(win), mode="nearest"); return slow, path - slow


def width_5_95(x):
    return float(np.percentile(x, 95) - np.percentile(x, 5))


# ------------------------------------------------------------------------------------------------ the floor, corridor time, drift timescale, regime

def tube_floor(path, w, W, levels=801):
    """Exact minimum of sum_{k>=1} w_k |u_k - u_{k-1}| (m/s; w in m/s per ln unit) with path_k + u_k inside a band of width W (ln); u_0 free (the start altitude is a design variable).
    Dynamic programming on a level grid with the min-plus (distance transform) step; the band is widened by half a grid step (error <= w dl / 2 per burn)."""
    if W <= 0: return float("inf")
    lo = -np.asarray(path, float); hi = lo + W; g = np.linspace(lo.min(), hi.max(), levels); dl = g[1] - g[0]; j = np.arange(levels)
    ok = lambda k: (g >= lo[k] - dl / 2) & (g <= hi[k] + dl / 2)
    f = np.where(ok(0), 0.0, np.inf)
    for k in range(1, len(path)):
        s = w[k] * dl
        fw = np.minimum.accumulate(f - s * j) + s * j                      # best predecessor at or below j
        bw = (np.minimum.accumulate((f + s * j)[::-1])[::-1]) - s * j      # best predecessor at or above j
        f = np.where(ok(k), np.minimum(fw, bw), np.inf)
    return float(f.min())


def zero_burn_fraction(path, W):
    """Largest share of passes whose burn-free path lies inside a band of width W for the best band position (plan R3.2: raw path, raw band, pass-weighted)."""
    p = np.sort(path); cnt = np.searchsorted(p, p + W, side="right") - np.arange(len(p))
    return float(cnt.max() / len(p))


def drift_timescale(slow, t_s):
    """Characteristic period (days) of the linearly detrended slow path and its ratio to the campaign duration: the power-weighted mean frequency of the Hann-windowed
    spectrum over frequencies of at least one cycle per campaign (an argmax would always land on the lowest allowed bin for a 1/f-like path)."""
    t = np.asarray(t_s, float); T = t[-1] - t[0]; N = max(64, 4 * len(t)); tu = np.linspace(t[0], t[-1], N); y = np.interp(tu, t, slow)
    y = y - np.polyval(np.polyfit(tu, y, 1), tu); y = y * np.hanning(N); pad = 8 * N
    P = np.abs(np.fft.rfft(y, pad)) ** 2; fr = np.fft.rfftfreq(pad, (tu[1] - tu[0])); m = fr >= 1.0 / T
    if P[m].sum() <= 0: return float("nan"), float("nan")
    per = 1.0 / (np.sum(fr[m] * P[m]) / np.sum(P[m]))
    return per / 86400.0, per / T


def classify(F, F_loss, D_need, W_eff, R, tol):
    """Regime (plan R3.4) from the floors: helps (below the drift-free floor), interferes (above), neutral (equal within tol)."""
    if W_eff <= 0: return "infeasible"
    if F < F_loss - tol: return "helps"
    if F > F_loss + tol: return "interferes"
    return "neutral"


# ------------------------------------------------------------------------------------------------ one candidate, end to end

def fast_width(cfg, ratio):
    b = cfg["band"]; a, c = b["fast_width_ln_at_ratio_054"], b["fast_width_ln_at_ratio_1"]
    return a + (c - a) * (ratio - 0.54) / (1.0 - 0.54)     # linear between the two derived points (ASSUMED interpolation)


def evaluate(ctx, cand, field=None, s_t_ody=None, jitter=0.0):
    cfg = ctx.cfg; b = cfg["band"]; fl = cfg["floor"]
    H_of_k = None if field is None else field["H"]
    q_req = solve_q_req(ctx, cand, b["d_ref_days"], H_of_k)
    S = schedule(ctx, cand, q_req * (1.0 + jitter), H_of_k)   # jitter: time-scale perturbation of the pass epochs (tesseral phase), W and q_req unchanged
    d = channels(ctx, cand, S, field)
    n = path_of(d); slow, fast = split_slow_fast(n, b["fast_window_passes"])
    tess = path_of(d, ["tesseral"]); s_t = float(np.std(tess - split_slow_fast(tess, b["fast_window_passes"])[0]))
    fw0 = fast_width(cfg, ctx.ratio)
    if s_t_ody is None: fw = fw0
    else:
        lam = min(1.0, (s_t_ody / b["fast_sd_ln_ody"]) ** 2); fw = fw0 * math.sqrt((1 - lam) + lam * (s_t / s_t_ody) ** 2)   # ASSUMED structure: only the tesseral share of the fast variance scales
    W = b["slack_band_factor"] * math.log(b["q_cap_mean_w_cm2"] / q_req); W_eff = W - fw
    # lever weights: burn between pass k-1 and k acts through the state of pass k-1; m/s per ln unit = H / (ratio * lever)
    w = np.zeros(S["n"]); w[1:] = S["H"][1:] / (ctx.ratio * np.array([lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, S["n"])]))
    F = tube_floor(slow, w, W_eff, fl["n_levels"])
    W_ts = math.log(b["q_cap_mean_w_cm2"] / q_req) - fw      # time-safe band [q_req, q_cap]: the lowest admissible q still meets D_ref, whatever the path does inside it
    F_ts = tube_floor(slow, w, W_ts, fl["n_levels"]) if W_ts > 0 else float("inf")
    d_cap = campaign_days(schedule(ctx, cand, b["q_cap_mean_w_cm2"], H_of_k))   # fastest possible campaign: constant q at the cap
    loss_slow = split_slow_fast(path_of(d, ["speed"]), b["fast_window_passes"])[0]; F_loss = tube_floor(loss_slow, w, W_eff, fl["n_levels"])
    D_need = float(-d["speed"].sum()); R = float(slow.max() - slow.min()); N = float(slow[-1] - slow[0])
    per_d, per_ratio = drift_timescale(slow, S["t"])
    return dict(id=cand["id"], omega0_deg=cand["omega"], raan0_deg=cand["raan"], epoch_off_d=cand["epoch_off"], inc_deg=cand["inc"], n_passes=S["n"],
                days=campaign_days(S), q_req=q_req, W=W, fast_w=fw, W_eff=W_eff, floor_mps=F, floor_driftfree_mps=F_loss, D_need=D_need, R=R, N=N,
                zero_burn_slack=zero_burn_fraction(n, W), zero_burn_odyssey_band=zero_burn_fraction(n, b["odyssey_band_ln"]),
                floor_timesafe_mps=F_ts, days_at_qcap=d_cap, frot_mean=float(S["frot"].mean()), vrel_start_mps=float(S["V"][0]), vrel_end_mps=float(S["V"][-1]),
                net_third_body_ln=float(d["third_body"].sum()) if "third_body" in d else 0.0, tau_drift_days=per_d, tau_drift_over_campaign=per_ratio, s_tesseral=s_t,
                regime=classify(F, F_loss, D_need, W_eff, R, cfg["floor"]["regime_tol_mps"]),
                feasible=bool(W_eff > 0 and q_req < b["q_cap_mean_w_cm2"] and b["q_cap_mean_w_cm2"] * math.exp(b["ln_qs_over_qmean_max"]) <= b["q_never_exceed_qs_w_cm2"]),
                field_mode="meanfield" if field is not None else "preview_no_field"), n, S, d


def jitter_floors(ctx, cand, field=None, frac=None, sun=None):
    """Floors over a phase-jitter ensemble. The tesseral kicks depend on the pass epochs (rotation phase), which the model reproduces only to a few percent of elapsed time
    (model/flown 1.04 to 1.05 in the selftest), so one realization is one draw of a noisy quantity; the dev oracle is a 12-seed mean."""
    jc = ctx.cfg["floor"]; frac = jc["phase_jitter_frac"] if frac is None else frac; old = ctx.sun
    if sun is not None: ctx.sun = sun
    try:
        return np.array([evaluate(ctx, cand, field, ctx.s_t_ody, j)[0]["floor_mps"] for j in np.linspace(-frac, frac, int(jc["phase_jitter_n"]))])
    finally:
        ctx.sun = old


def anchor_ctx(ctx, field=None):
    """Anchor evaluation: the deterministic run (sets the tesseral fast-rms reference) and its jitter ensemble (floor_mps is the ensemble mean)."""
    r, *_ = evaluate(ctx, anchor_cand(ctx.cfg), field); ctx.s_t_ody = r["s_tesseral"]; fl = jitter_floors(ctx, anchor_cand(ctx.cfg), field)
    r["floor_det_mps"] = r["floor_mps"]; r["floor_mps"] = float(np.mean(fl)); r["floor_jitter_sd_mps"] = float(np.std(fl)); return r


# ------------------------------------------------------------------------------------------------ dev oracle (plan R3.2 / R3.6)

def find_dev_dir():
    import os
    env = os.environ.get("SPACEAGORA_STAGEO_DEV")
    if env and pathlib.Path(env).is_dir(): return pathlib.Path(env)
    try:
        common = pathlib.Path(subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--path-format=absolute", "--git-common-dir"], text=True).strip())
    except Exception:
        return None
    p = common.parent / "tasks" / "2026-10-15_periapsis-region-border-targeting" / "stageo_data" / "dev"
    return p if p.is_dir() else None


def dev_natural(ctx, run, ratio):
    """Natural (burn-removed) ln q_mean path of one dev run and the burn cost weights (m/s per ln unit). Burn k acts between passes k-1 and k."""
    P = pd.read_csv(run + "/passes.csv").dropna(subset=["ra_m"]); D = pd.read_csv(run + "/decisions.csv"); n = len(P)
    lnq = np.log(P.q_w_cm2.values); H = P.h_local_m.values; rp = P.rp_m.values; ra = P.ra_m.values; burn = np.zeros(n); w = np.zeros(n)
    for i in range(1, n):
        a = 0.5 * (ra[i - 1] + rp[i - 1]); e = (ra[i - 1] - rp[i - 1]) / (ra[i - 1] + rp[i - 1]); w[i] = H[i] / (ratio * lever(ctx, a, e))
    for r in D.itertuples():
        i = int(r.target_pass) - 1
        if r.delta_v_mps != 0 and 1 <= i < n: burn[i] = -r.delta_v_mps / w[i]      # w = H / (ratio lever): ratio lever dv / H = dv / w
    return P, lnq - np.cumsum(burn), w


def oracle_floors(ctx, dev, arm, ratio, slack=0.72, levels=301):
    out = []
    for run in sorted(glob.glob(str(dev / "runs" / (arm + "*")))):
        P, nat, w = dev_natural(ctx, run, ratio); slow, fast = split_slow_fast(nat, ctx.cfg["band"]["fast_window_passes"]); fw = width_5_95(fast)
        out.append(dict(floor=tube_floor(slow, w, slack - fw, levels), fast_w=fw, zb_slack=zero_burn_fraction(nat, slack), zb_ody=zero_burn_fraction(nat, ctx.cfg["band"]["odyssey_band_ln"]),
                        net=float(nat[-1] - nat[0]), rng=float(nat.max() - nat.min())))
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ mean field (T10) input

def load_meanfield(directory, cfg):
    """Per candidate: H (local scale height from the three sampled altitudes) and ln rho at the center altitude, by pass index."""
    files = sorted(glob.glob(str(pathlib.Path(directory) / "mf_*.csv")))
    if not files: raise FileNotFoundError("no mf_*.csv in %s" % directory)
    df = pd.concat([pd.read_csv(f) for f in files]); off = cfg["sampling"]["alt_offsets_m"]; dz = 0.5 * (off[2] - off[0])
    out = {}
    for cid, g in df.groupby("cand"):
        piv = g.pivot_table(index="k", columns="alt_idx", values="rho_kg_m3").sort_index()
        l = np.log(piv.values); out[cid] = dict(H=-2 * dz / (l[:, 2] - l[:, 0]), lnrho0=l[:, 1], k=piv.index.values)
    return out


# ------------------------------------------------------------------------------------------------ decision rule (plan R3.5)

def decision(rank_df, f_ody, x):
    feas = rank_df[(rank_df.feasible) & np.isfinite(rank_df.floor_mps)]
    if feas.empty: return dict(branch="geometry_arms_dropped", reason="no feasible candidate", f_ody=f_ody, f_best=None)
    best = feas.loc[feas.floor_mps.idxmin()]; keep = best.floor_mps <= x * f_ody
    return dict(branch="geometry_arms_mandatory" if keep else "geometry_arms_dropped", f_ody=f_ody, f_best=float(best.floor_mps), ratio=float(best.floor_mps / f_ody),
                x=x, best_id=best["id"])


def precheck(cfg, ratio, f_pred):
    ref = cfg["precheck"]["oracle_floor_ratio_1_mps"] if abs(ratio - 1.0) < 1e-9 else cfg["precheck"]["oracle_floor_ratio_054_mps"]
    tol = cfg["precheck"]["tolerance"]; rel = f_pred / ref - 1
    return dict(f_pred=f_pred, f_ref=ref, rel_err=rel, tolerance=tol, passed=bool(abs(rel) <= tol))


# ------------------------------------------------------------------------------------------------ counts and tracks

def counts(cfg):
    full = make_grid(cfg); n_full = len(full); s = cfg["sampling"]; na = len(s["alt_offsets_m"])
    c256 = make_grid(cfg, cfg["grid"]["cut256_omega_step_deg"], cfg["grid"]["cut256_raan_step_deg"]); c384 = make_grid(cfg, cfg["grid"]["cut384_omega_step_deg"], cfg["grid"]["cut384_raan_step_deg"])
    return dict(full=n_full, cut256=len(c256), cut384=len(c384), calls_full=n_full * s["passes_nominal"] * na, calls_256=len(c256) * s["passes_nominal"] * na,
                calls_384=len(c384) * s["passes_nominal"] * na, calls_anchor=s["passes_nominal"] * na)


def write_tracks(ctx, cfg, out, fine=False):
    rows = []; g = make_grid(cfg, fine=fine) + [anchor_cand(cfg)]
    for c in g:
        q = solve_q_req(ctx, c, cfg["band"]["d_ref_days"]); S = schedule(ctx, c, q); lon, lat, _ = lon_lat(ctx, c, S)
        ep = (pd.Timestamp(cfg["anchor"]["epoch_utc"]) + pd.Timedelta(days=c["epoch_off"])).isoformat()
        for k in range(S["n"]):
            rows.append((c["id"], ep, c["epoch_off"], c["omega"], c["raan"], c["inc"], k, S["t"][k], lat[k], lon[k], cfg["orbit"]["h_ref_m"],
                         int(in_cut(c, cfg, "256")), int(in_cut(c, cfg, "384"))))
    df = pd.DataFrame(rows, columns="cand epoch_utc epoch_off_d omega0_deg raan0_deg inc_deg k t_s lat_gd_rad lon_rad h_center_m in_cut256 in_cut384".split())
    df.to_csv(out, index=False); return df


# ------------------------------------------------------------------------------------------------ CLI

def cmd_dry_run(cfg):
    c = counts(cfg); s = cfg["sampling"]
    print("candidate grid (R2.3): omega0 x Omega0 x epoch x i = %d x %d x %d x %d = %d" % (len(np.arange(0, 360, cfg["grid"]["omega0_step_deg"])), len(np.arange(0, 360, cfg["grid"]["raan0_step_deg"])),
          len(cfg["grid"]["epoch_offsets_days"]), len(cfg["grid"]["inclinations_deg"]), c["full"]), "(plan: 1,536) ->", "MATCH" if c["full"] == 1536 else "MISMATCH")
    print("plus 1 Odyssey anchor candidate (pre-check, F_ODY), not in the grid count")
    print("sampling calls (candidates x %d passes x %d altitudes): full %d (plan: about 0.92 million) -> %s; anchor %d" % (s["passes_nominal"], len(s["alt_offsets_m"]), c["calls_full"],
          "MATCH" if round(c["calls_full"] / 1e6, 2) == 0.92 else "MISMATCH", c["calls_anchor"]))
    print("cut grid if the probe exceeds %.1f ms/call: omega0 90 / Omega0 90 -> %d candidates (%d calls), plan states 256; omega0 60 / Omega0 90 (the plan's parenthesis) -> %d candidates (%d calls)" % (
        s["probe_limit_ms"], c["cut256"], c["calls_256"], c["cut384"], c["calls_384"]))
    print("  FLAG: the plan's parenthesis (omega0 every 60 deg, Omega0 every 90 deg) gives 384, not 256; 256 needs omega0 every 90 deg. cut_mode=%s in the config" % cfg["grid"]["cut_mode"])
    lo, hi = s["est_ms_per_call_low"], s["est_ms_per_call_high"]
    for nm, calls in (("full", c["calls_full"]), ("cut 256", c["calls_256"])):
        print("job estimate %-7s %d calls: %.1f to %.1f min single process at the ASSUMED %.1f to %.1f ms/call (%.2f to %.2f worker-h); %d workers: %.1f to %.1f min of coordinator wall (contended, not a performance number)" % (
            nm, calls, calls * lo / 6e4, calls * hi / 6e4, lo, hi, calls * lo / 3.6e6, calls * hi / 3.6e6, s["workers"], calls * lo / 6e4 / s["workers"], calls * hi / 6e4 / s["workers"]))
    print("probe: %d calls as %d cases of %d (route :process asserted); cut rule: if median ms/call > %.1f the full job is refused without --cut" % (s["probe_calls"], s["probe_workers"], s["probe_calls"] // s["probe_workers"], s["probe_limit_ms"]))
    ctx = Ctx(cfg); S = schedule(ctx, anchor_cand(cfg), 0.163); t0 = time.time()
    for k in range(10): tess_kick(ctx.cs, ctx.Ct, ctx.St, S["a"][0], S["e"][0], S["inc"], S["raan"][0], S["omega"][0], 0.3)
    ms = (time.time() - t0) * 100
    print("Python: tess_kick measured %.2f ms/call (this machine, idle-ish); full rank %d candidates x ~%d passes = about %.1f min tesseral (plan: about 12 min); not a timing result" % (
        ms, c["full"] + 1, S["n"], (c["full"] + 1) * S["n"] * ms / 6e4))
    print("anchor schedule (q = 0.163): %d passes, %.1f days, r_a end %.0f km" % (S["n"], campaign_days(S), (S["ra"][-1] - S["drop"][-1]) / 1e3))
    print("dry run: nothing flown, no Julia, no GRAM")


def cmd_oracle(cfg, ctx):
    dev = find_dev_dir()
    if dev is None: print("dev tables not found (set SPACEAGORA_STAGEO_DEV)"); return None
    res = {}
    for arm, refs in (("N2I", {1.0: 1.54, 0.54: 1.61}), ("FO", {1.0: 3.11, 0.54: 4.23})):
        for ratio, ref in refs.items():
            d = oracle_floors(ctx, dev, arm, ratio, levels=cfg["floor"]["n_levels"])
            print("%-3s ratio %.2f: floor %.3f (sd %.3f, n=%d) vs plan %.2f (%+.1f%%); fast width %.3f; zero-burn share slack/Odyssey band %.3f/%.3f; net %.2f ln, range %.2f ln" % (
                arm, ratio, d.floor.mean(), d.floor.std(), len(d), ref, 100 * (d.floor.mean() / ref - 1), d.fast_w.mean(), d.zb_slack.mean(), d.zb_ody.mean(), d.net.mean(), d.rng.mean()))
            res[(arm, ratio)] = d
    return res


def cmd_precheck(cfg, ctx, meanfield):
    field = None
    if meanfield:
        mf = load_meanfield(meanfield, cfg); field = mf["ODY"]
    r = anchor_ctx(ctx, field)
    pc = precheck(cfg, ctx.ratio, r["floor_mps"]); mode = "meanfield" if field else "PREVIEW (no mean-field term; not the R3.3 pre-check)"
    print("Odyssey anchor [%s]: %d passes, %.1f days, q_req %.3f (plan R1.3: about 0.163), W %.3f, fast width %.3f, W_eff %.3f, floor %.3f m/s, drift-free floor %.3f, regime %s" % (
        mode, r["n_passes"], r["days"], r["q_req"], r["W"], r["fast_w"], r["W_eff"], r["floor_mps"], r["floor_driftfree_mps"], r["regime"]))
    print("  path: net %.2f ln, range %.2f ln, D_need %.2f ln, zero-burn share %.2f (slack) / %.2f (Odyssey band), drift period %.1f d = %.2f of the campaign" % (
        r["N"], r["R"], r["D_need"], r["zero_burn_slack"], r["zero_burn_odyssey_band"], r["tau_drift_days"], r["tau_drift_over_campaign"]))
    print("  floor = mean over the phase-jitter ensemble (deterministic run %.3f, jitter sd %.3f)\n  pre-check vs dev oracle %.2f m/s: predicted %.3f, error %+.1f%%, tolerance +-%.0f%% -> %s" % (r["floor_det_mps"], r["floor_jitter_sd_mps"], pc["f_ref"], pc["f_pred"], 100 * pc["rel_err"], 100 * pc["tolerance"], "PASS" if pc["passed"] else "FAIL"))
    return r, pc


_G = {}


def _eval_i(i):
    c = _G["grid"][i]; return evaluate(_G["ctx"], c, _G["mf"][c["id"]], _G["ctx"].s_t_ody)[0]


def _ens_task(t):
    i, frac, sun = t; c = _G["grid"][i]; return t, jitter_floors(_G["ctx"], c, _G["mf"][c["id"]], frac, sun)


def _pmap(fn, items, jobs):
    import multiprocessing as mp
    if jobs <= 1: return [fn(x) for x in items]
    with mp.get_context("fork").Pool(jobs) as pool: return list(pool.imap(fn, items, chunksize=4))


def indistinguishable(df, best_id):
    """Candidates whose ensemble mean is within 2 x the combined jitter sd of the best one (the sd of one realization is the relevant noise, not the se of the mean)."""
    e = df.dropna(subset=["floor_ens_mps"]); b = e[e.id == best_id].iloc[0]
    return e[(e.floor_ens_mps - b.floor_ens_mps) <= 2 * np.sqrt(e.floor_ens_sd_mps ** 2 + b.floor_ens_sd_mps ** 2)]


def cmd_rank(cfg, ctx, out_dir, meanfield, fine=False, jobs=1):
    mf = load_meanfield(meanfield, cfg); out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True); fl = cfg["floor"]
    a, pc = cmd_precheck(cfg, ctx, meanfield)
    grid = [c for c in make_grid(cfg, fine=fine) if c["id"] in mf]; byid = {c["id"]: c for c in grid}
    _G.update(ctx=ctx, grid=grid, mf=mf); t0 = time.time()
    print("stage 1: %d candidates, %d jobs" % (len(grid), jobs), flush=True)
    df = pd.DataFrame(_pmap(_eval_i, range(len(grid)), jobs)); print("  stage 1 done, %.0f s" % (time.time() - t0), flush=True)
    df["floor_ens_mps"] = np.nan; df["floor_ens_sd_mps"] = np.nan
    order = df.sort_values("floor_mps"); sel = set(order.head(int(fl["refine_stat"])).index)
    for inc, g in df.groupby("inc_deg"):
        g = g.sort_values("floor_mps"); sel |= {g.index[0], g.index[len(g) // 2], g.index[-1]}
    # winner's curse and class noise: ensemble re-measure of the top refine_stat candidates and of the best / median / worst of every inclination
    res = _pmap(_ens_task, [(i, None, None) for i in sorted(sel)], jobs)
    for (i, _, _), f in res: df.at[i, "floor_ens_mps"] = f.mean(); df.at[i, "floor_ens_sd_mps"] = f.std()
    print("  stage 2 done (%d re-measured), %.0f s" % (len(sel), time.time() - t0), flush=True)
    df.to_csv(out / "geometry_rank.csv", index=False)
    top20 = df.loc[order.head(int(fl["refine_top"])).index].assign(floor_mps=lambda x: x.floor_ens_mps)
    dec = decision(top20, a["floor_mps"], cfg["decision"]["x_floor_ratio"]); dec["precheck"] = pc; dec["anchor_jitter_sd_mps"] = a["floor_jitter_sd_mps"]
    dec["rule_uses"] = "jitter-ensemble means of the top %d deterministic candidates and of the anchor" % fl["refine_top"]
    ts = df[np.isfinite(df.floor_timesafe_mps) & df.feasible]; anc = evaluate(ctx, anchor_cand(cfg), mf["ODY"], ctx.s_t_ody)[0]
    dec["timesafe_band_sensitivity"] = dict(n_timesafe_feasible=int(len(ts)), f_best_det=float(ts.floor_timesafe_mps.min()) if len(ts) else None, f_ody_det=float(anc["floor_timesafe_mps"]),
                                            ratio=float(ts.floor_timesafe_mps.min() / anc["floor_timesafe_mps"]) if len(ts) and np.isfinite(anc["floor_timesafe_mps"]) else None)
    best_id = top20.sort_values("floor_ens_mps").iloc[0]["id"]; ind = indistinguishable(df, best_id)
    dec["best_id"] = best_id; dec["n_indistinguishable_from_best"] = int(len(ind)); dec["n_remeasured"] = int(df.floor_ens_mps.notna().sum()); dec["n_candidates"] = int(len(df))
    (out / "geometry_decision.json").write_text(json.dumps(dec, indent=1)); print(json.dumps(dec, indent=1))
    # robustness for the top 10 (ensemble ranking): without the Sun term, and with larger tesseral timing error
    top10 = df.dropna(subset=["floor_ens_mps"]).sort_values("floor_ens_mps").head(10); idx = [list(df.index).index(i) for i in top10.index]
    tasks = [(i, f, sun) for i in idx for (f, sun) in ((None, False), (0.02, True), (0.10, True), (0.20, True))]
    rob = {t: f for t, f in _pmap(_ens_task, tasks, jobs)}; rows = []
    for i in idx:
        r = df.iloc[i]; row = dict(id=r["id"], inc=r.inc_deg, omega0=r.omega0_deg, raan0=r.raan0_deg, epoch=r.epoch_off_d, floor_ens=r.floor_ens_mps, ens_sd=r.floor_ens_sd_mps, net_sun_ln=r.net_third_body_ln)
        for (f, sun), nm in (((None, False), "no_sun"), ((0.02, True), "jit2"), ((0.10, True), "jit10"), ((0.20, True), "jit20")):
            v = rob[(i, f, sun)]; row[nm + "_mean"] = v.mean(); row[nm + "_sd"] = v.std()
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "robustness_top10.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(pd.DataFrame(rows).round(3).to_string())
    if not pc["passed"]: print("PRE-CHECK FAILED: the ranking is not trusted; only G1 flown, exploratory (plan R3.3)")
    print("stage 3 done, %.0f s" % (time.time() - t0))


def cmd_sun_compare(cfg, ctx, meanfield, old_dir, out_dir, jobs):
    """Re-evaluate the 50 best (old ensemble floor) candidates plus the anchor with the current Sun table; compare with the stored ranking (deterministic floor, ensemble floor, net Sun term)."""
    mf = load_meanfield(meanfield, cfg); old = pd.read_csv(pathlib.Path(old_dir) / "geometry_rank.csv"); out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    anchor_ctx(ctx, mf["ODY"]); top = old.dropna(subset=["floor_ens_mps"]).sort_values(["floor_ens_mps", "floor_mps"]).head(50)
    byid = {c["id"]: c for c in make_grid(cfg, fine=True)}; grid = [byid[i] for i in top["id"]]; _G.update(ctx=ctx, grid=grid, mf=mf)
    new = pd.DataFrame(_pmap(_eval_i, range(len(grid)), jobs)); ens = _pmap(_ens_task, [(i, None, None) for i in range(len(grid))], jobs)
    new["floor_ens_mps"] = [f.mean() for _, f in ens]; new["floor_ens_sd_mps"] = [f.std() for _, f in ens]
    m = top.reset_index(drop=True)[["id", "floor_mps", "floor_ens_mps", "floor_ens_sd_mps", "net_third_body_ln"]].merge(new[["id", "floor_mps", "floor_ens_mps", "floor_ens_sd_mps", "net_third_body_ln"]], on="id", suffixes=("_old", "_new"))
    m["d_net_sun_ln"] = m.net_third_body_ln_new - m.net_third_body_ln_old; m["d_floor_det"] = m.floor_mps_new - m.floor_mps_old; m["d_floor_ens"] = m.floor_ens_mps_new - m.floor_ens_mps_old
    m.to_csv(out / "sun_compare_top50.csv", index=False); pd.set_option("display.width", 250)
    print(m[["d_net_sun_ln", "d_floor_det", "d_floor_ens"]].abs().describe().round(5).to_string()); print("mean ens sd (old):", m.floor_ens_sd_mps_old.mean().round(4))
    print("rank by new ens floor, top 10:", list(m.sort_values(["floor_ens_mps_new", "floor_mps_new"]).id.head(10)), "old:", list(m.id.head(10)))
    return m


def cmd_class_best(cfg, ctx, meanfield, rank_dir, out_dir, incs, det_max, jobs):
    """Jitter-ensemble floors of every candidate of the given inclination classes whose deterministic floor is below det_max (so ties at zero are all measured); writes class_remeasure.csv."""
    mf = load_meanfield(meanfield, cfg); d = pd.read_csv(pathlib.Path(rank_dir) / "geometry_rank.csv"); out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    anchor_ctx(ctx, mf["ODY"]); sel = d[d.inc_deg.isin(incs) & (d.floor_mps < det_max)].copy(); byid = {c["id"]: c for c in make_grid(cfg, fine=True)}
    _G.update(ctx=ctx, grid=[byid[i] for i in sel["id"]], mf=mf); res = _pmap(_ens_task, [(i, None, None) for i in range(len(sel))], jobs)
    sel["floor_ens_mps"] = [f.mean() for _, f in res]; sel["floor_ens_sd_mps"] = [f.std() for _, f in res]; sel.to_csv(out / "class_remeasure.csv", index=False)
    for inc, g in sel.groupby("inc_deg"): print(inc, len(g), "\n", g.sort_values(["floor_ens_mps", "floor_ens_sd_mps"])[["id", "omega0_deg", "raan0_deg", "epoch_off_d", "floor_mps", "floor_ens_mps", "floor_ens_sd_mps"]].head(5).round(4).to_string())


# ------------------------------------------------------------------------------------------------ selftest

def selftest(cfg):
    from scipy.optimize import linprog
    ok = True
    def chk(name, cond, info=""):
        nonlocal ok; ok &= bool(cond); print(("PASS " if cond else "FAIL ") + name, info)
    rng = np.random.default_rng(1)
    # 1. DP against an independent LP (min sum w|du|, band constraints, free start)
    for trial in range(3):
        n = 30; path = np.cumsum(rng.normal(0, 0.15, n)); w = np.concatenate([[0], rng.uniform(0.1, 0.5, n - 1)]); W = 0.5
        c = np.concatenate([np.zeros(n), w[1:]]); A = []; b = []
        for k in range(1, n):
            r = np.zeros(2 * n - 1); r[k] = 1; r[k - 1] = -1; r[n + k - 1] = -1; A.append(r); b.append(0)
            r = np.zeros(2 * n - 1); r[k] = -1; r[k - 1] = 1; r[n + k - 1] = -1; A.append(r); b.append(0)
        bounds = [(-path[k], -path[k] + W) for k in range(n)] + [(0, None)] * (n - 1)
        lp = linprog(c, A_ub=np.array(A), b_ub=b, bounds=bounds, method="highs").fun
        dp = tube_floor(path, w, W, 4001)
        chk("tube_floor vs LP (trial %d)" % trial, abs(dp - lp) <= 0.01 * max(lp, 0.05) + 1e-3, "dp %.4f lp %.4f" % (dp, lp))
    chk("tube_floor flat path inside band is zero", tube_floor(np.zeros(20), np.ones(20), 0.3) == 0.0)
    chk("tube_floor step beyond band: (step - W) x weight", abs(tube_floor(np.r_[np.zeros(10), np.ones(10) * 1.0], np.ones(20) * 2.0, 0.4, 2001) - 2.0 * 0.6) < 0.01)
    # 2. zero-burn share against brute force
    p = rng.normal(0, 1, 60); W = 1.1; brute = max(np.mean((p >= b0) & (p <= b0 + W)) for b0 in p)
    chk("zero_burn_fraction vs brute force", abs(zero_burn_fraction(p, W) - brute) < 1e-12)
    # 3. drift timescale: known period in time on a nonuniform grid
    t = np.cumsum(np.linspace(6e4, 2.2e4, 200)); T = t[-1]; per = T / 3.0
    pd_, pr = drift_timescale(np.sin(2 * math.pi * t / per) + 0.2 * (t / T), t)
    chk("drift_timescale recovers a period of campaign/3", abs(pr - 1 / 3) < 0.05, "ratio %.3f" % pr)
    # 4. classifier on synthetic paths (loss ramp, drift cancelling it, large oscillation)
    k = np.arange(200); w = np.ones(200) * 1.0; ramp = -np.linspace(0, 1.2, 200)
    Wf = 0.5; Fl = tube_floor(ramp, w, Wf); Fh = tube_floor(ramp + np.linspace(0, 1.2, 200), w, Wf); Fi = tube_floor(ramp + 1.5 * np.sin(k / 10.0), w, Wf)
    chk("classify helps / interferes / neutral", classify(Fh, Fl, 1.2, Wf, 0.0, 0.05) == "helps" and classify(Fi, Fl, 1.2, Wf, 3.0, 0.05) == "interferes"
        and classify(Fl, Fl, 1.2, Wf, 1.2, 0.05) == "neutral" and classify(0, 0, 0, -0.1, 0, 0.05) == "infeasible", "floors %.2f %.2f %.2f" % (Fl, Fh, Fi))
    # 5. counts against plan R2.3
    c = counts(cfg); chk("grid counts 1536 / 256 / 384; 921,600 sampling calls", (c["full"], c["cut256"], c["cut384"], c["calls_full"]) == (1536, 256, 384, 921600), str(c))
    chk("decision rule thresholds", decision(pd.DataFrame(dict(id=["a", "b"], feasible=[True, True], floor_mps=[1.0, 0.7])), 1.0, 0.7)["branch"] == "geometry_arms_mandatory"
        and decision(pd.DataFrame(dict(id=["a"], feasible=[True], floor_mps=[0.71])), 1.0, 0.7)["branch"] == "geometry_arms_dropped")
    ctx = Ctx(cfg)
    # 6. physics reuse: J2 apse rate sign and size, lever
    om_dot = M.apsidal_rotation_per_pass(1.55e7, 0.776, math.radians(93.522), ctx.J2, ctx.Rg, ctx.mu)
    chk("J2 apsidal rotation per pass is negative for i = 93.5 deg, a few mrad", -0.01 < om_dot < 0, "%.2e rad" % om_dot)
    try:
        St = sun_table(); ck = pd.read_csv(St["check"]); xs = ck.days_j2000_utc.values; eng = ck[["sun_x_m", "sun_y_m", "sun_z_m"]].values; sp = St["spl"](xs)
        ang = np.degrees(np.arccos(np.clip((eng * sp).sum(1) / np.linalg.norm(eng, axis=1) / np.linalg.norm(sp, axis=1), -1, 1))); dist = np.abs(np.linalg.norm(sp, axis=1) / np.linalg.norm(eng, axis=1) - 1)
        chk("Sun table (engine SPICE call) vs 300 off-grid engine epochs 2001-10 to 2002-12: angle < 1e-4 deg, distance < 1e-8", ang.max() < 1e-4 and dist.max() < 1e-8, "max %.2e deg, %.2e rel" % (ang.max(), dist.max()))
        pole = St["pole"](xs[:5]); chk("engine pole is the IAU 2000 pole (RA 317.68, Dec 52.89 deg, 0.1 deg)", abs(math.degrees(math.atan2(pole[0][1], pole[0][0])) % 360 - 317.68) < 0.1 and abs(math.degrees(math.asin(pole[0][2])) - 52.89) < 0.1)
    except FileNotFoundError as ex: chk("Sun table present", False, str(ex))
    dev = find_dev_dir()
    if dev is None: print("SKIP dev-data checks (SPACEAGORA_STAGEO_DEV not set and the main checkout tables not found)"); print("SELFTEST", "PASS" if ok else "FAIL"); return ok
    run = sorted(glob.glob(str(dev / "runs" / "N2I*ODY_s1")))[0]; P = pd.read_csv(run + "/passes.csv")
    # 7. frame: body rotation angle and geodetic latitude against the pass table
    th_err = []
    for r in P.iloc[[0, 1, 2, 50, 100, 150]].itertuples():
        dd = ctx.anchor_days + r.time_s / 86400.0
        Px = math.cos(r.raan_rad) * math.cos(r.omega_rad) - math.sin(r.raan_rad) * math.sin(r.omega_rad) * math.cos(r.i_rad)
        Py = math.sin(r.raan_rad) * math.cos(r.omega_rad) + math.cos(r.raan_rad) * math.sin(r.omega_rad) * math.cos(r.i_rad)
        th_err.append(((math.atan2(Py, Px) - r.lon_rad - ctx.theta(dd) + math.pi) % (2 * math.pi)) - math.pi)
    lat_err = [abs(math.atan(math.tan(math.asin(math.sin(r.i_rad) * math.sin(r.omega_rad))) / (1 - ctx.f) ** 2) - r.lat_rad) for r in P.itertuples()]
    chk("rotation angle W(t) + offset reproduces lam - lon_body (max error < 0.02 deg)", np.max(np.abs(th_err)) < math.radians(0.02), "max %.4f deg" % math.degrees(np.max(np.abs(th_err))))
    chk("geodetic latitude of the periapsis vs table (median < 0.05 deg)", np.median(lat_err) < math.radians(0.05), "median %.4f deg" % math.degrees(np.median(lat_err)))
    # 8. decay scale: observed apoapsis drop over the K-theory map at q_mean
    ratios = []
    for rn in sorted(glob.glob(str(dev / "runs" / "N2I*"))):
        T = pd.read_csv(rn + "/passes.csv").dropna(subset=["ra_m"])
        for i in range(len(T) - 1):
            a_, e_, rp_ = T.a_m[i], T.e[i], T.rp_m[i]; _, fr = vrel_parts(ctx, rp_, e_, T.i_rad[i], math.asin(math.sin(T.i_rad[i]) * math.sin(T.omega_rad[i])))
            obs = T.ra_m[i] - T.ra_m[i + 1]; pred = M.drop_per_pass(a_, e_, rp_, T.q_w_cm2[i] * 1e4, M.k_theory(T.h_local_m[i], ctx.beta, fr), ctx.mu)
            if obs > 0 and pred > 0: ratios.append(obs / pred)
    med = float(np.median(ratios)); chk("decay scale: median observed/K-theory drop in [1.08, 1.18] (plan R1.1: 1.131)", 1.08 <= med <= 1.18, "%.3f over %d passes" % (med, len(ratios)))
    # 9. dev oracle reproduction (plan R3.2): fast width, floor within the plan's +-30%, zero-burn share conventions
    d1 = oracle_floors(ctx, dev, "N2I", 1.0, levels=301); d5 = oracle_floors(ctx, dev, "N2I", 0.54, levels=301)
    chk("fast width N2-I at ratio 1 / 0.54 reproduces 0.257 / 0.221 (+-0.01)", abs(d1.fast_w.mean() - 0.257) < 0.01 and abs(d5.fast_w.mean() - 0.221) < 0.01, "%.3f %.3f" % (d1.fast_w.mean(), d5.fast_w.mean()))
    chk("oracle floor N2-I within 30% of 1.54 (ratio 1) and 1.61 (0.54)", abs(d1.floor.mean() / 1.54 - 1) < 0.3 and abs(d5.floor.mean() / 1.61 - 1) < 0.3, "%.3f %.3f" % (d1.floor.mean(), d5.floor.mean()))
    chk("zero-burn share N2-I Odyssey/slack band reproduces 0.62 / 0.58 (+-0.03)", abs(d1.zb_ody.mean() - 0.62) < 0.03 and abs(d1.zb_slack.mean() - 0.58) < 0.03, "%.3f %.3f" % (d1.zb_ody.mean(), d1.zb_slack.mean()))
    # 10. schedule: campaign length at the dev flights' time-weighted harmonic q against the flown days (plan R1.2: the constant-Q model is 2 to 7% optimistic)
    ratios_days = []
    for rn in sorted(glob.glob(str(dev / "runs" / "N2I*")))[:4]:
        T = pd.read_csv(rn + "/passes.csv").dropna(subset=["ra_m"]); dt = np.diff(np.r_[0, T.time_s.values]); dt[0] = dt[1]; qh = dt.sum() / (dt / T.q_w_cm2.values).sum()
        S = schedule(ctx, anchor_cand(cfg), qh); ratios_days.append(campaign_days(S) / (T.time_s.values[-1] / 86400.0))
    chk("schedule days at the flown harmonic-mean q within 15% of the flown days (4 N2-I runs)", np.max(np.abs(np.array(ratios_days) - 1)) < 0.15, "model/flown %s" % np.round(ratios_days, 3))
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--oracle", action="store_true"); ap.add_argument("--tracks"); ap.add_argument("--precheck", action="store_true"); ap.add_argument("--rank"); ap.add_argument("--meanfield")
    ap.add_argument("--burn-response-ratio", type=float); ap.add_argument("--no-sun", action="store_true"); ap.add_argument("--sun-compare"); ap.add_argument("--class-best"); ap.add_argument("--incs"); ap.add_argument("--det-max", type=float, default=0.3); ap.add_argument("--old-rank"); ap.add_argument("--fine", action="store_true"); ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--config", default=str(CFG_PATH)); a = ap.parse_args(); cfg = load_cfg(a.config)
    if a.dry_run: return cmd_dry_run(cfg)
    if a.selftest: sys.exit(0 if selftest(cfg) else 1)
    ctx = Ctx(cfg, a.burn_response_ratio)
    if a.no_sun: ctx.sun = False
    if a.class_best: cmd_class_best(cfg, ctx, a.meanfield, a.old_rank, a.class_best, [float(x) for x in a.incs.split(',')], a.det_max, a.jobs)
    elif a.sun_compare: cmd_sun_compare(cfg, ctx, a.meanfield, a.old_rank, a.sun_compare, a.jobs)
    elif a.oracle: cmd_oracle(cfg, ctx)
    elif a.tracks: df = write_tracks(ctx, cfg, a.tracks, a.fine); print("tracks: %d rows, %d candidates -> %s" % (len(df), df.cand.nunique(), a.tracks))
    elif a.precheck: cmd_precheck(cfg, ctx, a.meanfield)
    elif a.rank: cmd_rank(cfg, ctx, a.rank, a.meanfield or sys.exit("--rank needs --meanfield DIR"), a.fine, a.jobs)
    else: ap.print_help()


if __name__ == "__main__":
    main()
