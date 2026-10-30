"""Analyzer of halflever_twin_run.jl (GRAM-on burn / no-burn twin, half-lever test I1).
    python3 -I halflever_twin_analyze.py --selftest
    python3 -I halflever_twin_analyze.py --dir <halflever_twin_data> [--dev <stageo_data/dev/runs>]
Every number is DERIVED from the run files (states in the Mars body-equator inertial frame, mu = the harmonics GM, SOURCED from halflever_vacuum_results.md). Writes <dir>/twin_summary.csv,
<dir>/twin_passes.csv (per pass differences), <dir>/twin_crossings.csv (altitude-crossing log of every run) and <dir>/twin_crossdiff.csv (burn minus no-burn at every crossing).
Method notes. Exact lever = osculating a(1-e) change of an impulsive dv along the velocity at the NO-BURN flight's first apoapsis event (as halflever_vacuum_analyze.py). Drag Delta v and energy
are trapezoid integrals of the logged drag acceleration (CorridorBallisticDrag formula with the GRAM mean density times the held perturbation factor). The drag prediction of the periapsis change is
first order: sum over the no-burn flight's steps of (d rp_osc / d v) . (a_drag_burn - a_drag_noburn) dt, the two accelerations aligned by time from each flight's own periapsis event."""
import csv, glob, math, os, sys
import numpy as np

MU = 4.2828370371e13        # SOURCED: harmonics GM, halflever_vacuum_results.md
LEVELS_KM = [600, 400, 300, 250, 200, 150, 120]
GAP_S = 30.0                # a gap larger than this between logged steps starts a new in-atmosphere segment
DEPART_TOL_M = 15.0         # ASSUMED: departure tolerance in osculating rp; the vacuum run's inbound time-aligned envelope is 3420.2 to 3430.3 m (halflever_vacuum_data, DERIVED)
EXACT_VAC_P25 = 3420.54     # SOURCED: halflever_vacuum_results.md, +0.25 m/s


def rd(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fcol(rows, k):
    return np.array([float(r[k]) for r in rows])


def elements_v(r, v, mu=MU):
    """Vectorized a, e, omega (rad, [0, 2 pi)), i for arrays r, v of shape (n, 3)."""
    rn = np.linalg.norm(r, axis=1); v2 = (v * v).sum(1); h = np.cross(r, v); hn = np.linalg.norm(h, axis=1)
    ev = np.cross(v, h) / mu - r / rn[:, None]; e = np.linalg.norm(ev, axis=1); a = 1.0 / (2.0 / rn - v2 / mu)
    nv = np.stack([-h[:, 1], h[:, 0], np.zeros(len(rn))], 1); nn = np.linalg.norm(nv, axis=1)
    cw = np.clip((nv * ev).sum(1) / (nn * np.maximum(e, 1e-15)), -1, 1); om = np.arccos(cw); om = np.where(ev[:, 2] >= 0, om, 2 * np.pi - om)
    return a, e, om, np.arccos(np.clip(h[:, 2] / hn, -1, 1))


def rp_of(r, v):
    a, e, _, _ = elements_v(np.atleast_2d(r), np.atleast_2d(v)); return (a * (1 - e))[0] if np.ndim(r) == 1 else a * (1 - e)


def exact_prediction(r, v, dv):
    """Change of osculating a(1-e) for an impulsive dv (m/s, signed) along the velocity at state (r, v)."""
    r = np.asarray(r, float); v = np.asarray(v, float); vn = np.linalg.norm(v)
    return rp_of(r, v * (1 + dv / vn)) - rp_of(r, v)


class Run:
    def __init__(self, path):
        self.path = path; self.name = os.path.basename(path)
        self.events = rd(os.path.join(path, "events.csv")); self.passes = rd(os.path.join(path, "passes.csv")); self.burn = rd(os.path.join(path, "burn.csv"))
        S = rd(os.path.join(path, "steps.csv")); self.cols = list(S[0].keys()) if S else []
        self.S = {k: fcol(S, k) for k in self.cols}
        self.t = self.S["t_s"]; self.r = np.stack([self.S["rx"], self.S["ry"], self.S["rz"]], 1); self.v = np.stack([self.S["vx"], self.S["vy"], self.S["vz"]], 1)
        self.ad = np.stack([self.S["adx"], self.S["ady"], self.S["adz"]], 1)
        self.a, self.e, self.om, self.inc = elements_v(self.r, self.v); self.rp = self.a * (1 - self.e); self.speed = np.linalg.norm(self.v, axis=1)
        self.rho_eff = self.S["rho_mean"] * self.S["r_pert"]
        # segments and pass labels (a segment is one atmosphere pass; label = the periapsis event inside it)
        brk = np.where(np.diff(self.t) > GAP_S)[0] + 1; self.seg_start = np.r_[0, brk]; self.seg_end = np.r_[brk, len(self.t)]
        self.per = [(int(e["pass"]), float(e["t_s"])) for e in self.events if e["kind"] == "periapsis"]
        self.apo = [(int(e["pass"]), float(e["t_s"])) for e in self.events if e["kind"] == "apoapsis"]
        self.seg_pass = []
        for s, e_ in zip(self.seg_start, self.seg_end):
            k = [p for p, tp in self.per if self.t[s] <= tp <= self.t[e_ - 1]]; self.seg_pass.append(k[0] if k else -1)
        # cumulative drag integrals (along-track dv, scalar magnitude, energy), from the first logged step, trapezoid within segments
        f_t = (self.ad * self.v).sum(1) / self.speed; f_m = np.linalg.norm(self.ad, axis=1); f_e = (self.ad * self.v).sum(1)
        self.cum_t = np.zeros(len(self.t)); self.cum_m = np.zeros(len(self.t)); self.cum_e = np.zeros(len(self.t)); c = np.zeros(3)
        for s, e_ in zip(self.seg_start, self.seg_end):
            self.cum_t[s], self.cum_m[s], self.cum_e[s] = c
            for i in range(s + 1, e_):
                dt = self.t[i] - self.t[i - 1]
                c = c + 0.5 * dt * np.array([f_t[i] + f_t[i - 1], f_m[i] + f_m[i - 1], f_e[i] + f_e[i - 1]]); self.cum_t[i], self.cum_m[i], self.cum_e[i] = c
        self.f_t = f_t

    def per_t(self, k):
        return dict(self.per).get(k)

    def crossings(self):
        """Rows of the altitude-crossing log: linear interpolation of every quantity between the two logged steps that bracket the level."""
        out = []; alt = self.S["alt_m"]
        for sg, (s, e_) in enumerate(zip(self.seg_start, self.seg_end)):
            for i in range(s, e_ - 1):
                for L in LEVELS_KM:
                    lv = L * 1e3
                    if (alt[i] - lv) * (alt[i + 1] - lv) > 0 or alt[i] == alt[i + 1]:
                        continue
                    w = (lv - alt[i]) / (alt[i + 1] - alt[i]); lin = lambda x: (1 - w) * x[i] + w * x[i + 1]
                    inbound = alt[i + 1] < alt[i]
                    out.append(dict(run=self.name, pass_=self.seg_pass[sg], level_km=L, dir="in" if inbound else "out", t_s=lin(self.t), radius_m=lin(np.linalg.norm(self.r, axis=1)),
                        speed_mps=lin(self.speed), a_m=lin(self.a), e=lin(self.e), rp_osc_m=lin(self.rp), omega_rad=lin(self.om), rho_mean=lin(self.S["rho_mean"]), r_pert=lin(self.S["r_pert"]),
                        rho_eff=lin(self.rho_eff), cum_drag_dv_alongtrack_mps=lin(self.cum_t), cum_drag_dv_scalar_mps=lin(self.cum_m), cum_drag_energy_j_per_kg=lin(self.cum_e),
                        step_gap_s=float(self.t[i + 1] - self.t[i])))
        return out

    def peak(self, k):
        """Pass peak density (mean and effective = mean x perturbation) over the logged steps of pass k."""
        sg = [j for j, p in enumerate(self.seg_pass) if p == k]
        if not sg:
            return float("nan"), float("nan")
        s, e_ = self.seg_start[sg[0]], self.seg_end[sg[0]]; return float(self.S["rho_mean"][s:e_].max()), float(self.rho_eff[s:e_].max())


def tau_series(run, k):
    """(tau = t - t_periapsis_event, index range) of the logged steps of pass k."""
    sg = [j for j, p in enumerate(run.seg_pass) if p == k]
    if not sg:
        return None
    s, e_ = run.seg_start[sg[0]], run.seg_end[sg[0]]; return run.t[s:e_] - run.per_t(k), s, e_


def drag_prediction(RA, RB, k):
    """First-order effect of the drag-acceleration difference (burn minus no-burn) in pass k on osculating rp: (up to the periapsis event, to the end of the logged pass). Sensitivity from the no-burn state."""
    tb = tau_series(RB, k); ta = tau_series(RA, k)
    if tb is None or ta is None:
        return float("nan"), float("nan")
    tau_b, sb, eb = tb; tau_a, sa, ea = ta
    dA = np.stack([np.interp(tau_b, tau_a, RA.ad[sa:ea, j]) for j in range(3)], 1); d = dA - RB.ad[sb:eb]
    r = RB.r[sb:eb]; v = RB.v[sb:eb]; base = RB.rp[sb:eb]; eps = 1e-3; S = np.zeros_like(v)
    for j in range(3):
        dv = np.zeros(3); dv[j] = eps; S[:, j] = (rp_of(r, v + dv) - base) / eps
    dtw = np.gradient(RB.t[sb:eb]); c = (S * d).sum(1) * dtw
    return float(c[tau_b <= 0].sum()), float(c.sum())


def pair(A, B, dv_cmd, exact_B, summary, passes_rows, cross_rows, seed_label):
    """Compare burn run A with its no-burn twin B."""
    bj = max(int(r["j"]) for r in A.burn); bA = [r for r in A.burn if int(r["j"]) == bj][0]; bBl = [r for r in B.burn if int(r["j"]) == bj]; bB = bBl[0] if bBl else None; treq = float(bA["request_t_s"])
    out = dict(seed=seed_label, run=A.name, dv_cmd_mps=dv_cmd, burn_apoapsis=bj, first_post_burn_pass=bj + 1)
    # bit-identity before the burn request
    mask = B.t < treq; nA = int((A.t < treq).sum()); nB = int(mask.sum())
    same = nA == nB and all(np.array_equal(A.S[k][:nA], B.S[k][:nB]) for k in A.cols)
    out.update(request_t_s=treq, n_steps_before_request=nB, steps_identical_before_request=bool(same))
    # pass rows whose apoapsis closes before the request are bit-identical; the row of the pass just before the burn differs only in apo_alt_m (the apoapsis event falls inside the burn window)
    npre = bj; diffcols = sorted({k for i in range(min(npre, len(A.passes), len(B.passes))) for k in A.passes[i] if A.passes[i][k] != B.passes[i][k]})
    out["pre_burn_pass_rows_differing_columns"] = "+".join(diffcols) if diffcols else "none"; out["pre_burn_pass_rows_identical_but_apo_alt"] = bool(set(diffcols) <= {"apo_alt_m", "ra_m"})
    out["request_t_equal"] = None if bB is None else bool(float(bA["request_t_s"]) == float(bB["request_t_s"])); out["request_nu_equal"] = None if bB is None else bool(float(bA["nu_request_rad"]) == float(bB["nu_request_rad"]))
    evb = [e for e in B.events if float(e["t_s"]) < treq]; eva = [e for e in A.events if float(e["t_s"]) < treq]
    out["events_identical_before_request"] = bool(evb == eva)
    # exact lever and burn window
    ws, we = float(bA["win_start_s"]), float(bA["win_stop_s"]); out.update(win_start_s=ws, win_stop_s=we, burn_len_s=we - ws, dv_used_mps=float("nan"))
    exact = exact_B(dv_cmd, bj); out["exact_lever_m"] = exact
    # per-pass event differences
    rpA = {int(e["pass"]): np.array([float(e["rx"]), float(e["ry"]), float(e["rz"])]) for e in A.events if e["kind"] == "periapsis"}
    rpB = {int(e["pass"]): np.array([float(e["rx"]), float(e["ry"]), float(e["rz"])]) for e in B.events if e["kind"] == "periapsis"}
    pA = {int(r["pass"]): r for r in A.passes}; pB = {int(r["pass"]): r for r in B.passes}
    apA = {int(e["pass"]): float(np.linalg.norm([float(e["rx"]), float(e["ry"]), float(e["rz"])])) for e in A.events if e["kind"] == "apoapsis"}
    apB = {int(e["pass"]): float(np.linalg.norm([float(e["rx"]), float(e["ry"]), float(e["rz"])])) for e in B.events if e["kind"] == "apoapsis"}
    k0 = bj + 1   # first post-burn pass; the pair estimator of the earlier analyses uses rp(k0) - rp(k0 - 1) of the burn flight against lever x dv
    if k0 in rpA and bj in rpA and k0 in rpB and bj in rpB:
        nA = float(np.linalg.norm(rpA[k0]) - np.linalg.norm(rpA[bj])); nB = float(np.linalg.norm(rpB[k0]) - np.linalg.norm(rpB[bj]))
        out.update(pair_change_burn_flight_m=nA, pair_change_noburn_flight_m=nB, pair_naive_ratio=nA / exact if exact else float("nan"),
                   pair_ratio_with_minus113_baseline=(nA + 113.0) / exact if exact else float("nan"))   # -113 m: F-O clean-pass median, halflever_attribution_results.md (SOURCED)
    for k in sorted(set(rpA) & set(rpB)):
        d_rp = float(np.linalg.norm(rpA[k]) - np.linalg.norm(rpB[k])); a, b = pA.get(k), pB.get(k)
        row = dict(seed=seed_label, run=A.name, dv_cmd_mps=dv_cmd, pass_=k, lag=k - (bj + 1), d_rp_event_m=d_rp, ratio_to_exact=d_rp / exact if exact else float("nan"))
        if a and b:
            for key, nm in (("a_m", "d_a_m"), ("e", "d_e"), ("omega_rad", "d_omega_rad"), ("alt_m", "d_alt_m"), ("rho", "d_rho_mean"), ("r_pert", "d_r_pert"), ("rho_eff", "d_rho_eff"), ("time_s", "d_t_peri_s")):
                row[nm] = float(a[key]) - float(b[key])
            row["rho_mean_ratio_AB"] = float(a["rho"]) / float(b["rho"]); row["rho_eff_ratio_AB"] = float(a["rho_eff"]) / float(b["rho_eff"]); row["h_local_B_m"] = float(b["h_local_m"])
            row["pred_rho_ratio_from_dalt"] = math.exp(-row["d_alt_m"] / float(b["h_local_m"]))
        row["d_ra_after_event_m"] = (apA.get(k, float("nan")) - apB.get(k, float("nan")))
        row["pk_rho_mean_A"], row["pk_rho_eff_A"] = A.peak(k); row["pk_rho_mean_B"], row["pk_rho_eff_B"] = B.peak(k)
        up, full = drag_prediction(A, B, k); row["drag_pred_to_peri_m"] = up; row["drag_pred_full_pass_m"] = full
        passes_rows.append(row)
    # crossing differences
    cA = A.crossings(); cB = B.crossings(); keyB = {(c["pass_"], c["level_km"], c["dir"]): c for c in cB}; cr = []
    for c in cA:
        b = keyB.get((c["pass_"], c["level_km"], c["dir"]))
        if b is None:
            continue
        d = dict(seed=seed_label, run=A.name, dv_cmd_mps=dv_cmd, pass_=c["pass_"], level_km=c["level_km"], dir=c["dir"])
        for k in ("t_s", "radius_m", "speed_mps", "a_m", "e", "rp_osc_m", "omega_rad", "rho_mean", "r_pert", "rho_eff", "cum_drag_dv_alongtrack_mps", "cum_drag_dv_scalar_mps", "cum_drag_energy_j_per_kg"):
            d["d_" + k] = c[k] - b[k]
        d["rp_osc_B_m"] = b["rp_osc_m"]; d["rho_eff_B"] = b["rho_eff"]; d["max_step_gap_s"] = max(c["step_gap_s"], b["step_gap_s"]); cr.append(d)
    cross_rows.extend(cr)
    return out


def run_all(D, dev=None):
    names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(D, "s*_*")) if os.path.isfile(os.path.join(p, "steps.csv")))
    runs = {n: Run(os.path.join(D, n)) for n in names}; summary = []; passes_rows = []; cross_rows = []; allcross = []
    for n, R in runs.items():
        allcross.extend(R.crossings())
    seeds = sorted({n.split("_")[0] for n in runs}); dvs = {"none": 0.0, "p25": 0.25, "m25": -0.25, "p14": 0.14, "f136": 0.136, "late25": 0.25}
    for s in seeds:
        B = runs.get(s + "_none")
        if B is None:
            continue
        evs = [e for e in B.events if e["kind"] == "apoapsis"]
        def exact_B(dv, j, evs=evs):
            ev = evs[j - 1]; return float(exact_prediction(np.array([float(ev[k]) for k in ("rx", "ry", "rz")]), np.array([float(ev[k]) for k in ("vx", "vy", "vz")]), dv))
        for tag, dv in dvs.items():
            A = runs.get(s + "_" + tag)
            if A is None or tag == "none":
                continue
            summary.append(pair(A, B, dv, exact_B, summary, passes_rows, cross_rows, s))
    return runs, summary, passes_rows, cross_rows, allcross


def dev_twins(R):
    """Existing Stage O dev flights as twins: the N2 arms burn at apoapsis 1 (decisions.csv target_pass 2), F-O does not (dv = 0 there), same seed, so bit-identical up to the burn.
    Ratio = (pass-2 radius of the N2 flight minus that of the F-O flight of the same seed) / (lever at the pass-1 (a, e) x the first N2 dv). SOURCED files, DERIVED ratios."""
    out = []
    for arm, pat in (("N2T", "N2T_cc01_a010_t152_w04000"), ("N2I", "N2I_cc12_a010_t137_w04000"), ("N21", "N21_cc01_a010_t178_w04000")):
        for k in range(1, 13):
            try:
                fo = rd(os.path.join(R, "FO_fixed1000_pert_ODY_s%d" % k, "passes.csv")); n2 = rd(os.path.join(R, "%s_pert_ODY_s%d" % (pat, k), "passes.csv")); dec = rd(os.path.join(R, "%s_pert_ODY_s%d" % (pat, k), "decisions.csv"))
            except FileNotFoundError:
                continue
            if int(dec[0]["target_pass"]) != 2 or float(dec[0]["delta_v_mps"]) == 0:
                continue
            dv = float(dec[0]["delta_v_mps"]); a = float(fo[0]["a_m"]); e = float(fo[0]["e"]); lev = 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)
            d = float(n2[1]["radius_m"]) - float(fo[1]["radius_m"]); out.append(dict(arm=arm, seed=k, dv_mps=dv, d_rp_pass2_m=d, lever_formula_m=lev * dv, ratio=d / (lev * dv)))
    return out


def vacuum_compare(vac_dir, rp_e=3396190.0):
    """Altitude-aligned osculating-rp difference (burn minus no-burn) of the VACUUM flights (halflever_vacuum_data nom / none) at the crossing levels, inbound and outbound, radius level = rp_e + L
    (equatorial radius ASSUMED for the level; geodetic and spherical altitude differ by at most a few km, which moves the crossing by about 1 s and the osculating rp by about 0.01 m)."""
    def load(n):
        rows = rd(os.path.join(vac_dir, n, "grid.csv")); G = np.array([[float(r[k]) for k in ("t_s", "rx", "ry", "rz", "vx", "vy", "vz")] for r in rows])
        a, e, _, _ = elements_v(G[:, 1:4], G[:, 4:7]); return G[:, 0], np.linalg.norm(G[:, 1:4], axis=1), a * (1 - e)
    tA, rA, pA = load("nom"); tB, rB, pB = load("none"); out = []
    def cross(t, r, p, R, inb):
        for i in range(len(t) - 1):
            if t[i] > 80000 and (r[i] - R) * (r[i + 1] - R) <= 0 and ((r[i + 1] < r[i]) == inb):
                w = (R - r[i]) / (r[i + 1] - r[i]); return p[i] + w * (p[i + 1] - p[i])
    for L in LEVELS_KM:
        for inb in (True, False):
            a = cross(tA, rA, pA, rp_e + L * 1e3, inb); b = cross(tB, rB, pB, rp_e + L * 1e3, inb); out.append(dict(level_km=L, dir="in" if inb else "out", vac_d_rp_osc_m=a - b, vac_minus_exact_m=a - b - EXACT_VAC_P25))
    return out


def dev_pair_estimator(R):
    """The pair estimator of the earlier analyses applied to the FIRST burn of every dev flight: (rp at the target pass minus rp one pass earlier) / (formula lever at the earlier pass's (a, e) x dv).
    F-O first burn: apoapsis 2 (target pass 3); N2 arms: apoapsis 1 (target pass 2). The causal ratio of the same burns is the twin ratio (about 1.000)."""
    pats = {"FO": "FO_fixed1000", "N2T": "N2T_cc01_a010_t152_w04000", "N2I": "N2I_cc12_a010_t137_w04000", "N21": "N21_cc01_a010_t178_w04000"}; out = []
    for arm, pat in pats.items():
        for k in range(1, 13):
            try:
                P = rd(os.path.join(R, "%s_pert_ODY_s%d" % (pat, k), "passes.csv")); dec = rd(os.path.join(R, "%s_pert_ODY_s%d" % (pat, k), "decisions.csv"))
            except FileNotFoundError:
                continue
            d = [x for x in dec if abs(float(x["delta_v_mps"])) > 0]
            if not d:
                continue
            tp = int(d[0]["target_pass"]); dv = float(d[0]["delta_v_mps"]); a = float(P[tp - 2]["a_m"]); e = float(P[tp - 2]["e"]); lev = 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)
            ch = float(P[tp - 1]["radius_m"]) - float(P[tp - 2]["radius_m"]); out.append(dict(arm=arm, seed=k, target_pass=tp, dv_mps=dv, pair_change_m=ch, pair_ratio=ch / (lev * dv), pair_ratio_minus113_baseline=(ch + 113.0) / (lev * dv)))
    return out


def write(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def selftest():
    # exact lever: tangential impulse at the apoapsis of a known orbit equals 4 a sqrt((1-e)/(1+e))/sqrt(mu/a) dv to the linearization
    ra, rp = 28.533e6, 3.49e6; a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp); va = math.sqrt(MU * (2 / ra - 1 / a))
    ex = exact_prediction([ra, 0, 0], [0, va, 0], 0.25); lin = 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a) * 0.25
    assert abs(ex / lin - 1) < 5e-3, (ex, lin)
    assert abs(exact_prediction([ra, 0, 0], [0, va, 0], -0.25) / ex + 1) < 5e-3
    # elements_v on a perifocal orbit with omega = 1 rad
    om, nd, inc = 1.0, 0.3, 0.8
    Rz = lambda x: np.array([[math.cos(x), -math.sin(x), 0], [math.sin(x), math.cos(x), 0], [0, 0, 1]]); Rx = lambda x: np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
    R = Rz(nd) @ Rx(inc) @ Rz(om); rpv = R @ np.array([rp, 0, 0]); vp = R @ np.array([0, math.sqrt(MU * (2 / rp - 1 / a)), 0])
    aa, ee, ww, ii = elements_v(rpv[None], vp[None]); assert abs(aa[0] / a - 1) < 1e-9 and abs(ee[0] - e) < 1e-9 and abs(ww[0] - om) < 1e-9 and abs(ii[0] - inc) < 1e-9, (aa, ee, ww, ii)
    print("halflever_twin_analyze selftest ok: exact lever %.2f m for +0.25 m/s (formula %.2f), elements ok" % (ex, lin))


def main():
    if "--selftest" in sys.argv:
        return selftest()
    D = sys.argv[sys.argv.index("--dir") + 1]
    runs, summary, passes_rows, cross_rows, allcross = run_all(D)
    write(os.path.join(D, "twin_summary.csv"), summary); write(os.path.join(D, "twin_passes.csv"), passes_rows)
    write(os.path.join(D, "twin_crossings.csv"), allcross); write(os.path.join(D, "twin_crossdiff.csv"), cross_rows)
    if "--dev" in sys.argv:
        dt = dev_twins(sys.argv[sys.argv.index("--dev") + 1]); write(os.path.join(D, "twin_devpairs.csv"), dt); rr = np.array([r["ratio"] for r in dt])
        print("dev N2-vs-F-O first-burn pairs: n=%d ratio mean %.4f min %.4f max %.4f" % (len(rr), rr.mean(), rr.min(), rr.max()))
    if "--dev" in sys.argv:
        dp = dev_pair_estimator(sys.argv[sys.argv.index("--dev") + 1]); write(os.path.join(D, "twin_devpair_estimator.csv"), dp)
        for arm in ("FO", "N2T", "N2I", "N21"):
            rr = np.array([r["pair_ratio"] for r in dp if r["arm"] == arm]); rr.size and print("dev first-burn pair estimator %s: n=%d mean %.3f min %.3f max %.3f" % (arm, rr.size, rr.mean(), rr.min(), rr.max()))
    if "--vac" in sys.argv:
        vc = vacuum_compare(sys.argv[sys.argv.index("--vac") + 1]); write(os.path.join(D, "twin_vacuum_compare.csv"), vc)
    for s in summary:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()})
    for r in passes_rows:
        if r["pass_"] in (1, 2, 3):
            print(r["run"], "pass", r["pass_"], "d_rp_event %.2f m ratio %.4f d_omega %.3e d_a %.2f d_alt %.1f rho_eff_ratio %.4f drag_pred_to_peri %.2f full %.2f" % (r["d_rp_event_m"], r["ratio_to_exact"],
                r.get("d_omega_rad", float("nan")), r.get("d_a_m", float("nan")), r.get("d_alt_m", float("nan")), r.get("rho_eff_ratio_AB", float("nan")), r["drag_pred_to_peri_m"], r["drag_pred_full_pass_m"]))


if __name__ == "__main__":
    main()
