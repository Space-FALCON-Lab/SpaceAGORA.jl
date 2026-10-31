"""Revision 21 WP21-B: rendezvous POSITION scoring (analysis only: Python, no Julia, no flights, reads flight directories).

Target. The apoapsis direction at the rendezvous time t*: initial elements of the geometry (omega0, Omega0, i0; SOURCED planner_config / geometry_gate_table) propagated with the secular J2 rates evaluated along the
model-world NOMINAL plan's (a, e) history (the open-loop plan from the calibrated first-pass level, planner_rev20_replay.nominal_arrival_d recipe, whose arrival is the guard's t*) up to t*; i stays i0 (the plan has no inclination state).
Flight. The apoapsis direction at the arrival t_arr = t(closing periapsis) + T/2 from the flight's final osculating (omega, Omega, i) (passes.csv columns omega_rad, raan_rad, i_rad at the closing pass) advanced by
half a J2 orbit (periapsis -> apoapsis).
Error = great-circle angle between the two unit vectors (apoapsis direction = minus the perifocal P axis); arc = angle x RA_T (7,408 km, SOURCED).

J2 secular rates (SOURCED form, as planner_rev18_flights.PLB_apsidal; Kaula 1966 / Vallado ch. 9): per orbit d(omega) = 3 pi J2 (R/p)^2 (2 - 5/2 sin^2 i), d(Omega) = -3 pi J2 (R/p)^2 cos i, d(i) = 0, p = a (1 - e^2).
Drag (DERIVED, first order): the drag acceleration is in the orbit plane and opposite to the relative velocity, which is the velocity up to the (small) atmosphere rotation. The tangential part gives da, de but d(omega) from a purely
tangential force is (2 sin f / e) F_t / (n a ...) with f the true anomaly of the force: it is odd in f for a force symmetric about periapsis, so its integral over a pass is zero at first order; the residual is second order in the
rotation of the atmosphere (Omega_planet r / v ~ 0.5 percent at periapsis) and acts through di and dOmega, which the predictor residual (b) measures, not assumes.
Tesseral (bound): tesseral_bound() integrates the first-order per-pass kicks of (i, Omega, omega) from the 20x20 field (tesseral_lib.gauss_kick, tesseral-only coefficients, Mars50c) along the nominal history for several
planet-rotation phases (the phase at epoch is not known here: ASSUMED sampled), reporting sum |kick| (worst case, DERIVED), rms (random phase) and the phase-sampled accumulated angle.

Decomposition (angles in the tangent plane at the target; vectors L_b(u) = log map of u at base b projected to the target tangent plane):
  (a) timing = J2 nominal at t_arr minus J2 nominal at t*      (= J2 rate x (t_arr - t*))
  (c) history = J2 from the flight's own (a, e, i) history at t_arr minus J2 nominal at t_arr
  (b) predictor residual = actual flight apoapsis direction minus the J2 prediction from its own history
  total = a + b + c to second order in the angles (they telescope exactly on the sphere; the numerical closure residual is the tangent-plane composition error, reported per flight).
Independent check (not by construction): (a) against the analytic rate: |d(apo)/dt| from the partial derivatives of the apoapsis direction at the nominal final rates times |t_arr - t*|.

Usage: python3 -I planner_rev21_position.py selftest | run [--out-dir DIR] | tesseral [--out-dir DIR]
Env: SPACEAGORA_PLANNER_ROOT (flight root, default ~/Documents/SpaceAGORA-planner-data/planner), SPACEAGORA_BORDER_TASK_DIR (task directory, planner_lib.task_dir)."""
import argparse, glob, json, math, os, pathlib, re, sys, tomllib
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev18_offline as O
import planner_rev20_window as W

TD = pathlib.Path(PLB.task_dir()); OUT = TD / "rev21_data" / "position"
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
CFG = tomllib.load(open(HERE / "planner_config.toml", "rb"))["planner"]
MU = O.MU; RA_T = 7.408e6        # SOURCED target apoapsis radius (default; the per-geometry target radius is target_radius())
LOW_TARGET_ID = {"G0486": "CRIT63_2000km", "G1230": "CRIT117_1000km", "G1038": "ODY_500km"}      # Revision 25 low-altitude block: token -> target_id of target_orbits.json (ASSUMED mapping from the Revision 25 text: 63.4 deg/2,000 km, 116.6 deg/1,000 km, 93.5 deg/500 km)
TARGET_OVERRIDE = {}      # token -> target radius (m); set by the scorer CLI (--target-radius TOK=m) or by a caller
J2 = float(W.J2); RJ = float(W.RJ)       # SOURCED: Kaula Mars table (planner_rev20_window)
TOL_DEG = (0.5, 1.0, 2.0)        # tol_pos_deg = 1.0 ASSUMED (named), 0.5 and 2.0 shown
TOL_BASE = 1.0
W_ROT = 7.08823596e-5            # SOURCED planner_config / planets.jl Mars spin (rad/s)
GEOM_ID = {"ODY": "ODY", "G5863": "G1", "G0975": "G0975", "G0063": "G0063", "G0390": "G0390", "G1216": "G1216", "G0528": "G0528", "G0890": "G0890", "G0975M": "G0975M", "G0063M": "G0063M", "G0486": "G0486", "G1230": "G1230", "G1038": "G1038", "G0140_AODY": "G0140_AODY", "G0042_A2000": "G0042_A2000", "G3257_AODY": "G3257_AODY", "G3161_A2000": "G3161_A2000"}      # run-id prefix -> planner token
GATE_ID = {"G1": "G5863"}


# ------------------------------------------------------------------------------------------------ geometry
def deg(x): return math.degrees(x)


def apo_hat(om, Om, inc):
    """Unit vector to the apoapsis (minus the perifocal P axis) in the inertial frame of the elements."""
    P = np.array([math.cos(Om) * math.cos(om) - math.sin(Om) * math.sin(om) * math.cos(inc),
                  math.sin(Om) * math.cos(om) + math.cos(Om) * math.sin(om) * math.cos(inc),
                  math.sin(om) * math.sin(inc)])
    return -P


def ang(u, v): return math.acos(max(-1.0, min(1.0, float(np.dot(u, v)))))


def logmap(base, u):
    """Tangent vector at base pointing to u with length = great-circle angle."""
    c = float(np.dot(u, base)); th = math.acos(max(-1.0, min(1.0, c))); v = u - c * base; n = float(np.linalg.norm(v))
    return np.zeros(3) if n < 1e-15 else th * v / n


def tproj(base, v): return v - float(np.dot(v, base)) * base


def j2_orbit(a, e, inc):
    """Secular J2 change per orbit (d omega, d Omega) in rad (SOURCED form)."""
    p = a * (1 - e * e); k = 3 * math.pi * J2 * (RJ / p) ** 2
    return k * (2 - 2.5 * math.sin(inc) ** 2), -k * math.cos(inc)


def _table_files():
    """Ground-table files of the task directory: the generic and Revision 23 files first, then any other planner_ground_tables_*.json (e.g. the Revision 25 low-altitude pins), SOURCED."""
    fs = [TD / "planner_ground_tables_generic.json", TD / "planner_ground_tables_rev23.json"]
    rest = [f for f in sorted(TD.glob("planner_ground_tables_*.json")) if f not in fs and f.name not in ("planner_ground_tables.json",)]
    return fs + sorted(rest, key=lambda f: (0 if f.name.endswith("mf.json") else 1, f.name))      # mean-field tables (…mf.json, e.g. rev25mf) before the preview pins of the same tokens


_TF = {}


def table_entry(tok):
    """Raw ground-table dict of a token from the table files (None if absent)."""
    for f in _table_files():
        try:
            if f not in _TF: _TF[f] = json.load(open(f))["geometries"]
            g = _TF[f].get(tok)
        except Exception: continue
        if g: return g
    return None


def target_radius(tok):
    """Target apoapsis radius (m) of a geometry token: explicit override, then the table's own field (ra_target_m / target_radius_m), then target_orbits.json (low-altitude tokens, circle_r_km), else the default 7,408,000 m (SOURCED stageo_config.toml)."""
    if tok in TARGET_OVERRIDE: return float(TARGET_OVERRIDE[tok])
    g = table_entry(tok) if tok not in ("ODY", "G1") else None
    if g:
        for k in ("ra_target_m", "target_radius_m"):
            if g.get(k): return float(g[k])
    if tok in LOW_TARGET_ID:
        try:
            for t in json.load(open(TD / "target_orbits.json"))["targets"]:
                if t["target_id"] == LOW_TARGET_ID[tok]: return float(t["circle_r_km"]) * 1e3
        except Exception: pass
    return RA_T


def geometry_elements(tok):
    """(i0, omega0, Omega0) in rad of a planner token. SOURCED: planner_config [planner.geometry.*] and geometry_gate_table.json."""
    if tok in CFG["geometry"]:
        g = CFG["geometry"][tok]; return math.radians(g["inc_deg"]), math.radians(g["omega_deg"]), math.radians(g["raan_deg"])
    J = json.load(open(TD / "geometry_gate_table.json"))
    for g in J["geometries"]:
        if g["id"] == GATE_ID.get(tok, tok):
            return math.radians(g["inclination_deg"]), math.radians(g["omega0_deg"]), math.radians(g["raan0_deg"])
    g = table_entry(tok)       # Revision 23 and later tokens: elements stored with the table (SOURCED)
    if g: return math.radians(g["inc_deg"]), math.radians(g["omega_deg"]), math.radians(g["raan_deg"])
    raise KeyError(tok)


class J2Track:
    """Piecewise J2 secular propagation from an initial (omega0, Omega0, i0) along a pass history. Orbit k carries the per-orbit increment of its own (a_k, e_k, i_k) spread over
    [lo_k, hi_k] (lo_1 = 0, hi_k = (t_k + t_{k+1})/2, hi_last = t_last + T_last/2); beyond the last orbit the last rate continues. Times in seconds from the epoch (the starting apoapsis)."""

    def __init__(self, t, a, e, inc, om0, Om0, i0):
        n = len(t); T = 2 * np.pi * np.sqrt(np.asarray(a, float) ** 3 / MU); hi = np.empty(n); hi[:-1] = 0.5 * (np.asarray(t[:-1]) + np.asarray(t[1:])); hi[-1] = t[-1] + 0.5 * T[-1]
        lo = np.concatenate(([0.0], hi[:-1])); self.lo, self.hi, self.T = lo, hi, T
        d = np.array([j2_orbit(a[k], e[k], inc[k]) for k in range(n)]); self.dom, self.dOm = d[:, 0], d[:, 1]
        self.om0, self.Om0, self.i0 = om0, Om0, i0; self.n = n

    def at(self, t):
        f = np.clip((t - self.lo) / (self.hi - self.lo), 0.0, 1.0)
        extra = max(0.0, t - self.hi[-1]) / self.T[-1]
        return self.om0 + float(np.sum(f * self.dom)) + extra * self.dom[-1], self.Om0 + float(np.sum(f * self.dOm)) + extra * self.dOm[-1]

    def apo(self, t, inc=None):
        om, Om = self.at(t); return apo_hat(om, Om, self.i0 if inc is None else inc)

    def rate_vec(self, t):
        """Time derivative (rad/s) of the apoapsis direction at t from the local orbit's rates (analytic, central difference in the elements, step 1e-6 rad)."""
        k = int(np.clip(np.searchsorted(self.hi, t), 0, self.n - 1)); om, Om = self.at(t); r_om = self.dom[k] / self.T[k]; r_Om = self.dOm[k] / self.T[k]; h = 1e-6
        d_om = (apo_hat(om + h, Om, self.i0) - apo_hat(om - h, Om, self.i0)) / (2 * h); d_Om = (apo_hat(om, Om + h, self.i0) - apo_hat(om, Om - h, self.i0)) / (2 * h)
        return d_om * r_om + d_Om * r_Om


# ------------------------------------------------------------------------------------------------ nominal plan history
def _rp():
    import planner_rev20_replay as RP
    return RP


GENERIC = ("G0975", "G0063", "G0390", "G1216")      # tokens with a ground table in planner_ground_tables_generic.json (SOURCED)
KAPPA_CAL = TD / "rev21_data" / "position" / "kappa_generic.json"


REV23 = ("G0528", "G0890", "G0975M", "G0063M")      # tables in planner_ground_tables_rev23.json (with the GRAM mean-field term; kappa 1 ASSUMED)


def get_tab(tok):
    """Ground table of a token: ODY, G1, G0975 via planner_rev20_replay.get_tab; every other token from the task-directory table files (generic, Revision 23, pins), with the same derived fields as planner_rev20_replay.load_tab."""
    if tok in ("ODY", "G1", "G0975"): return _rp().get_tab(tok)
    g = table_entry(tok)
    if g is None: raise KeyError(f"no ground table for {tok}")
    tab = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in g.items()}; tab["token"] = tok
    n = tab["n"]; ra = tab["ra"]; tab["d0"] = np.concatenate((np.diff(-ra), [float(ra[-2] - ra[-1])])); tab["_ra_inc"] = ra[::-1]; tab["_d0_inc"] = tab["d0"][::-1]
    tab["_k_inc"] = np.arange(n, 0, -1, dtype=float)[::-1] * 1.0; tab["_kk"] = np.arange(1, n + 1, dtype=float); return tab


KAPPA_ALIAS = {"G0975M": "G0975"}      # the flights' config maps the mean-field G0975M to the G0975 calibration (kappa 0.9518, MEASURED from the gen1 flights); logged guard t* of G0975M flights reproduces with it


def kappa_for(tok):
    """Drop-per-q factor: MEASURED rev19 calibration (ODY, G1); a generic token's value from rev21_data/position/kappa_generic.json when its calibration was accepted (IQR tight), else 1.0 (ASSUMED)."""
    RP = _rp()
    if tok in ("ODY", "G1"): return RP.KAPPA[tok]
    try:
        c = json.load(open(KAPPA_CAL)).get(KAPPA_ALIAS.get(tok, tok))
        if c and c.get("accepted"): return float(c["kappa_med"])
    except Exception: pass
    return 1.0


def calibrate_kappa(tok, dirs):
    """planner_rev18_calibrate method: per flight kappa = sum(flown apoapsis drops) / sum(model drops at the current apoapsis radius, d0 interpolated on the table, x q / q_req) over the counted passes; geometry value = median over flights
    (IQR reported). Accepted when (q75 - q25) / median <= KAPPA_IQR_MAX (ASSUMED 0.05)."""
    tab = get_tab(tok); ks = []
    for d in dirs:
        P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra_all = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra_all) & (ra_all <= target_radius(tok)))[0]
        if not len(hit): continue
        n = int(hit[0] + 1); ra = ra_all[:n]; q = P.q_w_cm2.values.astype(float)[:n]; ra_prev = np.concatenate(([float(tab["ra"][0])], ra[:-1])); drop = ra_prev - ra
        d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); m_st = np.interp(ra_prev, tab["ra"][::-1], d0[::-1]) * q / tab["q_req"]; ok = np.isfinite(drop) & (m_st > 0)
        ks.append(float(np.sum(drop[ok]) / np.sum(m_st[ok])))
    k = np.array(ks); q25, q75 = (float(np.percentile(k, 25)), float(np.percentile(k, 75))) if len(k) else (float("nan"),) * 2; med = float(np.median(k)) if len(k) else float("nan")
    return dict(geom=tok, n=len(k), kappa_med=med, q25=q25, q75=q75, min=float(k.min()) if len(k) else float("nan"), max=float(k.max()) if len(k) else float("nan"), iqr_rel=(q75 - q25) / med if len(k) else float("nan"),
                accepted=bool(len(k) >= 4 and (q75 - q25) / med <= KAPPA_IQR_MAX), iqr_max_assumed=KAPPA_IQR_MAX, method="planner_rev18_calibrate: sum(flown drops)/sum(model drops), median over flights")


KAPPA_IQR_MAX = 0.05      # ASSUMED acceptance: relative IQR of the per-flight kappa


def headline_ic(tok):
    """The flown S-fixed initial-condition plan of a geometry (same selection as planner_rev20_score.ic_dv_pre): the configured headline file first, then the other headline_*.json pins of the task directory (mean-field pins before preview pins), SOURCED."""
    cand = [TD / CFG["ic"]["path"]] + sorted((f for f in (TD / "rev20_data" / "ic").glob("headline_*.json") if f != TD / CFG["ic"]["path"]), key=lambda f: (0 if f.stem.endswith("mf") else 1, f.name))
    for f in cand:
        try: J = json.load(open(f))
        except Exception: continue
        c = [p for p in J["plans"] if p["geom"] == tok and p["scenario"] == "S-fixed" and p["fb_term"] is False and p["horizon"] == "arrival" and abs(p["delta_ln"]) < 1e-12 and abs(p["h_des_m"] - 250000.0) < 1e-6 and p["t_star_d"] is None]
        if len(c) == 1: return c[0]
    return None


def arrival_history(tab, x, q_first_ln, kappa, ra_target=RA_T):
    """Copy of planner_rev18_offline.arrival that also returns the pass history: t_k (periapsis time, s), a_k, e_k (orbit k before the drop), the arrival day and pass count."""
    n = tab["n"]; q_req = tab["q_req"]; ra = float(tab["ra"][0]); t = float(tab["t"][0]); path = tab["path"]; kmax = int(n * 3)
    d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])])); ra_t = tab["ra"][::-1]; d0_r = d0[::-1]; T_, A_, E_ = [], [], []
    for k in range(1, kmax + 1):
        kk = min(k, n); xk = float(x[min(k, len(x)) - 1]); nu = xk - (q_first_ln + float(path[kk - 1] - path[0]))
        rp = float(tab["rp"][kk - 1]) - float(tab["H"][kk - 1]) * nu; a = 0.5 * (ra + rp)
        drop = kappa * float(np.interp(ra, ra_t, d0_r)) * math.exp(xk) / q_req; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2 * math.pi * math.sqrt(amid ** 3 / MU)
        T_.append(t); A_.append(a); E_.append((ra - rp) / (ra + rp))
        if ra_next <= ra_target:
            return dict(t=np.array(T_), a=np.array(A_), e=np.array(E_), days=(t + 0.5 * Tk) / 86400.0, N=k, T_last=Tk)
        ra = ra_next; t += Tk
    return dict(t=np.array(T_), a=np.array(A_), e=np.array(E_), days=float("nan"), N=kmax, T_last=Tk)


_NOM = {}


def nominal(tok, ic_x0=None):
    """Model-world nominal for a geometry: dict(track=J2Track on the plan history, days=t* (the nominal arrival of the S-fixed IC plan, guard recipe), hist=history). SOURCED recipe: planner_rev20_replay.tstar_for.
    ic_x0 (a flight's logged first-pass ln q anchor, e.g. the S-free arm) replaces the headline S-fixed shift: the nominal of the plan that flight flies."""
    key = (tok, None if ic_x0 is None else round(float(ic_x0), 9))
    if key in _NOM: return _NOM[key]
    RP = _rp(); ic = headline_ic(tok) if ic_x0 is None else dict(x0=float(ic_x0)); x0 = math.log(RP.QFIRST_CAL.get(tok, O.Q_FIRST)) + RP.ic_shift(ic)      # Q_FIRST 0.212 (MEASURED calibration target) where no per-geometry q_first exists
    tab = get_tab(tok); sw, P = RP.arm_def("OLR"); inp = PLB.ol_inputs(sw, tab, P, x0); pl = PLB.dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"])
    x = np.concatenate(([x0], np.asarray(inp["level"]) + pl["nu"])); h = arrival_history(tab, x, x0, kappa_for(tok), target_radius(tok)); i0, om0, Om0 = geometry_elements(tok)
    tr = J2Track(h["t"], h["a"], h["e"], np.full(len(h["t"]), i0), om0, Om0, i0)
    _NOM[key] = dict(track=tr, days=h["days"], hist=h, tok=tok, N=h["N"]); return _NOM[key]


# ------------------------------------------------------------------------------------------------ one flight
def flight_position(d, tok, t_star_d=None, nom=None, ra_target=None, stop_extra_m=0.0):
    """Position error of one finished flight directory against the nominal of its geometry. Returns None when the target apoapsis was never reached."""
    nom = nom or nominal(tok); tstar = nom["days"] * 86400.0 if t_star_d is None else t_star_d * 86400.0
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); rt = target_radius(tok) if ra_target is None else float(ra_target); hit = np.where(np.isfinite(ra) & (ra <= rt + stop_extra_m))[0]
    if not len(hit): return None
    n = int(hit[0] + 1); Pn = P.iloc[:n]; t = Pn.time_s.values.astype(float); a = Pn.a_m.values.astype(float); e = Pn.e.values.astype(float); inc = Pn.i_rad.values.astype(float)
    om = Pn.omega_rad.values.astype(float); Om = Pn.raan_rad.values.astype(float); i0, om0, Om0 = geometry_elements(tok)
    fl = J2Track(t, a, e, inc, om0, Om0, i0); t_arr = float(t[-1] + 0.5 * fl.T[-1]); nt = nom["track"]
    dom, dOm = j2_orbit(a[-1], e[-1], inc[-1]); u_fl = apo_hat(om[-1] + 0.5 * dom, Om[-1] + 0.5 * dOm, inc[-1])      # actual: final osculating elements advanced half an orbit to the apoapsis
    u_tgt = nt.apo(tstar); u_nA = nt.apo(t_arr); u_fp = fl.apo(t_arr)
    tot = logmap(u_tgt, u_fl); A = logmap(u_tgt, u_nA); C = tproj(u_tgt, logmap(u_nA, u_fp)); B = tproj(u_tgt, logmap(u_fp, u_fl)); res = tot - (A + B + C)
    om_p, Om_p = fl.at(t_arr); ang_tot = ang(u_fl, u_tgt); dt_d = (t_arr - tstar) / 86400.0
    a_an = float(np.linalg.norm(nt.rate_vec(0.5 * (t_arr + tstar)))) * abs(t_arr - tstar)          # independent analytic J2-rate x dt
    wrap = lambda x: (x + math.pi) % (2 * math.pi) - math.pi
    return dict(n=n, t_arr_d=t_arr / 86400.0, tstar_d=tstar / 86400.0, dt_d=dt_d, err_deg=deg(ang_tot), arc_km=ang_tot * rt / 1e3,
                a_timing_deg=deg(float(np.linalg.norm(A))), b_pred_deg=deg(float(np.linalg.norm(B))), c_hist_deg=deg(float(np.linalg.norm(C))),
                sum_abs_deg=deg(float(np.linalg.norm(A) + np.linalg.norm(B) + np.linalg.norm(C))), closure_deg=deg(float(np.linalg.norm(res))), closure_rel=float(np.linalg.norm(res) / ang_tot) if ang_tot > 0 else 0.0,
                a_analytic_deg=deg(a_an), a_ratio=float(np.linalg.norm(A) / a_an) if a_an > 0 else float("nan"),
                d_omega_deg=deg(wrap(om[-1] + 0.5 * dom - om_p)), d_raan_deg=deg(wrap(Om[-1] + 0.5 * dOm - Om_p)), d_inc_deg=deg(inc[-1] - fl.i0),
                j2_apse_deg_per_d=deg(float(np.linalg.norm(nt.rate_vec(tstar)))) * 86400.0)


# ------------------------------------------------------------------------------------------------ tesseral bound
def tesseral_bound(tok, n_phase=12, nf=64, stride=1):
    """First-order tesseral bound on the apoapsis-direction error accumulated along the NOMINAL history. Per pass: kicks (d i, d Omega, d omega) of tesseral-only Mars50c 20x20 (tesseral_lib.gauss_kick, apoapsis-to-apoapsis
    Gauss quadrature, elements of the nominal J2 track, planet rotation phase theta0 + W_ROT t); angle of the apoapsis-direction displacement per pass. sum_abs = worst case (phases aligned, no feedback), rms = quadrature sum, and
    the phase-sampled coherent accumulation |sum of vectors| (n_phase values of theta0 in [0, 2 pi), ASSUMED sampling). ASSUMED: no drag feedback on the pass times; osculating elements at periapsis stand in for those at the apoapsis."""
    import tesseral_lib as tl
    cs = tl.coefficients("Mars", 20); C, S = tl.select(cs["C"], cs["S"], lambda n, m: m >= 1); nm = nominal(tok); h = nm["hist"]; tr = nm["track"]; i0 = tr.i0
    idx = np.arange(0, len(h["t"]), stride); disp = np.zeros((n_phase, len(idx), 3)); per = np.zeros((n_phase, len(idx)))
    for ip in range(n_phase):
        th0 = 2 * math.pi * ip / n_phase
        for jj, k in enumerate(idx):
            om, Om = tr.at(float(h["t"][k])); kk = tl.gauss_kick(float(h["a"][k]), float(h["e"][k]), i0, Om, om, th0 + W_ROT * float(h["t"][k]), W_ROT, cs, C, S, nf=nf)
            u0 = apo_hat(om, Om, i0); u1 = apo_hat(om + kk["omega"], Om + kk["raan"], i0 + kk["i"]); disp[ip, jj] = logmap(u0, u1); per[ip, jj] = ang(u0, u1)
    sc = float(stride)
    return dict(geom=tok, passes_used=int(len(idx)), stride=stride, phases=n_phase, nf=nf, sum_abs_deg=[deg(float(x) * sc) for x in per.sum(1)], rms_deg=[deg(float(math.sqrt((x ** 2).sum())) * math.sqrt(sc)) for x in per],
                coherent_deg=[deg(float(np.linalg.norm(disp[ip].sum(0))) * sc) for ip in range(n_phase)], per_pass_rms_deg=deg(float(np.sqrt((per ** 2).mean()))))


# ------------------------------------------------------------------------------------------------ discovery and statistics
def parse_id(name):
    m = re.match(r"^(\w+?)_PL_(.+)_pert_[A-Za-z0-9]+_s(\d+)(?:_(L\d\w*))?$", name)
    if not m: return None
    return dict(gid=m.group(1), tok=GEOM_ID.get(m.group(1)), arm=m.group(2), seed=int(m.group(3)), level=m.group(4) or "")


def finished(d): return all((d / f).exists() for f in ("passes.csv", "decisions.csv", "measured.csv"))


def block_dirs(root, block):
    """Flight directories of a block (a directory under the root, or a path relative to it, e.g. rev20val, rev21c, ../trx_pull4/planner/devR)."""
    return [p for p in sorted((pathlib.Path(root) / block).glob("*")) if p.is_dir() and parse_id(p.name) and parse_id(p.name)["tok"]]


def score_block(root, block, tstar_override=None):
    """Position rows of every finished, reached flight of a block. tstar_override: dict tok -> nominal t* (days) (planner_rev20_score passes the logged one)."""
    rows = []
    for d in block_dirs(root, block):
        k = parse_id(d.name)
        if not finished(d): continue
        r = flight_position(d, k["tok"], None if not tstar_override else tstar_override.get(k["tok"]))
        if r is None: continue
        r.update(block=block, run=d.name, arm=k["arm"], geom=k["tok"], level=k["level"], seed=k["seed"]); rows.append(r)
    return rows


def pct(x, q): return float(np.percentile(x, q)) if len(x) else float("nan")


def cell_stats(rows, keyf):
    """Per-cell median, 90th percentile, fraction within tol (0.5, 1.0, 2.0 deg) of err_deg, plus medians of the decomposition."""
    cells = {}
    for r in rows: cells.setdefault(keyf(r), []).append(r)
    out = {}
    for c, rs in cells.items():
        e = np.array([r["err_deg"] for r in rs]); g = lambda k: float(np.median([r[k] for r in rs]))
        out[c] = dict(n=len(rs), med_deg=float(np.median(e)), p90_deg=pct(e, 90), max_deg=float(e.max()), within={f"{t:g}": float((e <= t).mean()) for t in TOL_DEG}, med_arc_km=g("arc_km"), med_a_deg=g("a_timing_deg"), med_b_deg=g("b_pred_deg"), med_c_deg=g("c_hist_deg"),
                      med_abs_dt_d=float(np.median([abs(r["dt_d"]) for r in rs])), med_closure_deg=g("closure_deg"), max_closure_deg=float(max(r["closure_deg"] for r in rs)))
    return out


def fmt_cell_table(stats, title_key="cell"):
    h = f"| {title_key} | n | median err (deg) | p90 (deg) | within 0.5 | within 1.0 | within 2.0 | median arc (km) | med (a) timing | med (b) predictor | med (c) history | med abs dt (d) |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
    return h + "\n".join(f"| {c} | {s['n']} | {s['med_deg']:.2f} | {s['p90_deg']:.2f} | {s['within']['0.5']:.2f} | {s['within']['1']:.2f} | {s['within']['2']:.2f} | {s['med_arc_km']:.0f} | {s['med_a_deg']:.2f} | {s['med_b_deg']:.2f} | {s['med_c_deg']:.2f} | {s['med_abs_dt_d']:.2f} |" for c, s in stats.items())


# ------------------------------------------------------------------------------------------------ selftest
def selftest():
    # 1 geometry
    assert np.allclose(apo_hat(0.0, 0.0, 1.0), [-1, 0, 0]); assert np.allclose(apo_hat(math.pi / 2, 0.0, math.pi / 2), [0, 0, -1])
    for om, Om, i in [(0.3, 1.1, 1.6), (4.0, 5.5, 2.0)]: assert abs(np.linalg.norm(apo_hat(om, Om, i)) - 1) < 1e-14
    u = apo_hat(0.3, 0.2, 1.0); v = apo_hat(0.31, 0.2, 1.0); assert abs(ang(u, v) - np.linalg.norm(logmap(u, v))) < 1e-12 and np.linalg.norm(logmap(u, u)) == 0.0
    # 2 J2 rates: critical inclination zero, retrograde Omega sign, matches rev18 form
    ic = math.asin(math.sqrt(0.8)); assert abs(j2_orbit(1.0e7, 0.5, ic)[0]) < 1e-15
    assert j2_orbit(1.6e7, 0.78, math.radians(93.5))[1] > 0 and j2_orbit(1.6e7, 0.78, math.radians(60))[1] < 0
    import planner_rev18_flights as FL
    a, e, i = 1.6e7, 0.78, 1.6; assert abs(j2_orbit(a, e, i)[0] - FL.PLB_apsidal(a, e, i)) < 1e-14
    # 3 track: constant history reproduces rate x t, continuation beyond the last orbit
    n = 40; T = 2 * math.pi * math.sqrt(1.5e7 ** 3 / MU); t = (np.arange(n) + 0.5) * T; tr = J2Track(t, np.full(n, 1.5e7), np.full(n, 0.7), np.full(n, 1.5), 0.2, 0.3, 1.5)
    dom, dOm = j2_orbit(1.5e7, 0.7, 1.5); om, Om = tr.at(10 * T); assert abs(om - (0.2 + 10 * dom)) < 1e-12 and abs(Om - (0.3 + 10 * dOm)) < 1e-12
    om, Om = tr.at(n * T + 3 * T); assert abs(om - (0.2 + (n + 3) * dom)) < 1e-12
    # analytic rate matches the track difference
    r = tr.rate_vec(5 * T); d = tr.apo(5 * T + 100.0) - tr.apo(5 * T - 100.0); assert np.linalg.norm(r * 200.0 - d) < 1e-9
    # 4 nominal arrival_history == planner_rev18_offline.arrival and t* of the guard recipe
    RP = _rp(); tab = RP.get_tab("ODY"); x0 = math.log(RP.QFIRST_CAL["ODY"]); sw, P = RP.arm_def("OLR"); inp = PLB.ol_inputs(sw, tab, P, x0); pl = PLB.dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"])
    x = np.concatenate(([x0], np.asarray(inp["level"]) + pl["nu"])); h = arrival_history(tab, x, x0, RP.KAPPA["ODY"]); o = O.arrival(tab, x, x0, kappa=RP.KAPPA["ODY"])
    assert h["N"] == o["N"] and abs(h["days"] - o["days"]) < 1e-12 and abs(h["days"] - RP.nominal_arrival_d("ODY", x0)) < 1e-12
    # 5 decomposition closure on synthetic geometry: perturb timing + history + residual
    ub = apo_hat(1.9, 0.5, 1.63); v1 = apo_hat(1.9 + 0.004, 0.5 - 0.001, 1.63); v2 = apo_hat(1.9 + 0.01, 0.5, 1.63 + 0.002); v3 = apo_hat(1.9 + 0.012, 0.5 + 0.003, 1.63 + 0.004)
    tot = logmap(ub, v3); s = logmap(ub, v1) + tproj(ub, logmap(v1, v2)) + tproj(ub, logmap(v2, v3)); assert np.linalg.norm(tot - s) < 5e-5 * 3 and np.linalg.norm(tot) > 5e-3
    # 6 zero error for a flight equal to its nominal is checked on data when present
    for t in ("G0063", "G0390", "G1216", "G0975"):
        tb = get_tab(t); nm = nominal(t); assert nm["N"] > 10 and np.isfinite(nm["days"]) and tb["token"] == t and kappa_for(t) > 0
    k = calibrate_kappa("G0975", [d for d in block_dirs(ROOT, "gen1")][:0]); assert k["n"] == 0
    print("selftest ok")
    return True


# ------------------------------------------------------------------------------------------------ CLI
BLOCKS = ["rev20val", "rev20olt", "rev17", "gen1", "dev", "devT", "devR", "../trx_pull4/planner/devR"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "run", "tesseral", "kappa"]); ap.add_argument("--out-dir"); ap.add_argument("--root"); ap.add_argument("--blocks", nargs="*"); ap.add_argument("--n-phase", type=int, default=12); ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args(); out = pathlib.Path(a.out_dir) if a.out_dir else OUT; root = pathlib.Path(a.root) if a.root else ROOT
    if a.cmd == "selftest": selftest(); return
    out.mkdir(parents=True, exist_ok=True)
    if a.cmd == "kappa":
        ds = [d for d in block_dirs(root, "gen1") if parse_id(d.name)["tok"] == "G0975" and finished(d)]; c = calibrate_kappa("G0975", ds); KAPPA_CAL.parent.mkdir(parents=True, exist_ok=True)
        old = json.load(open(KAPPA_CAL)) if KAPPA_CAL.exists() else {}; old["G0975"] = c; (out / "kappa_generic.json").write_text(json.dumps(old, indent=1)); print(c); return
    if a.cmd == "tesseral":
        res = {tok: tesseral_bound(tok, a.n_phase, stride=a.stride) for tok in ("ODY", "G1", "G0975")}
        (out / "tesseral_bound.json").write_text(json.dumps(res, indent=1)); [print(k, "worst-case sum %.3f deg, rms %.3f, coherent max %.3f" % (max(v["sum_abs_deg"]), max(v["rms_deg"]), max(v["coherent_deg"]))) for k, v in res.items()]; return
    rows = []
    for b in (a.blocks or BLOCKS): rows += score_block(root, b)
    (out / "position_flights.json").write_text(json.dumps(rows, indent=1))
    noms = {t: dict(tstar_d=nominal(t)["days"], N=nominal(t)["N"], kappa=kappa_for(t)) for t in ("ODY", "G1") + GENERIC}
    summ = dict(tol_deg_assumed=TOL_BASE, nominal=noms, by_block_arm_geom_level={f"{c[0]}|{c[1]}|{c[2]}|{c[3]}": v for c, v in cell_stats(rows, lambda r: (r["block"], r["arm"], r["geom"], r["level"])).items()})
    (out / "position_summary.json").write_text(json.dumps(summ, indent=1)); print(len(rows), "flights ->", out)


if __name__ == "__main__":
    main()
