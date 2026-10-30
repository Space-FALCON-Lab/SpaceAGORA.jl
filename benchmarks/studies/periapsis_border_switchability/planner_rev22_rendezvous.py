"""Revision 22 (WP22): rendezvous at a general point of the final aerobraking ellipse (analysis only: Python, no Julia, no flights).

Pure functions (SI, radians unless named _deg):
  state on an ellipse at true anomaly nu (perifocal and inertial), Kepler time between two anomalies, secular J2 advance of (Omega, omega),
  relative plane (inclination, line of nodes) of two orbits, radius crossing, relative-velocity cost (including the out-of-plane part), optimal nu_r.
Event model (Revision 22 item 1): the spacecraft stays on the final ellipse F after the closing apoapsis (nu = 180 deg) and meets a target at the true anomaly nu_r in [180, 360) deg,
i.e. before the next periapsis (the next drag pass changes the orbit); states below H_MIN_ALT are excluded (ASSUMED, below). Rendezvous burn = |v_target - v_F| at the meeting point (position
continuity: the two orbits must contain the point; the burn changes velocity only). Coplanar: any radius in both radius ranges is a crossing (the target's apse direction is a free
family parameter, so each crossing radius is realizable); non-coplanar: the point must lie on the line of nodes (its true anomaly is the node argument u_n, one of the two nodes is always on [180, 360)).
The target PHASE is free in the sweeps (ASSUMED: a target on its orbit can be made to be at the point; the phase-locked case is scored for the nominal ellipse in the flight scoring).

Usage: python3 -I planner_rev22_rendezvous.py selftest | run [--out-dir DIR]
Env: SPACEAGORA_PLANNER_ROOT (flight root), SPACEAGORA_BORDER_TASK_DIR (task dir, planner_lib.task_dir). Imports planner_rev21_position read-only."""
import json, math, os, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev21_position as POS

MU = POS.MU; RJ = POS.RJ; RA_T = POS.RA_T            # SOURCED: planner config mu, Kaula reference radius (3394.2 km), target apoapsis radius 7,408 km
TWO_PI = 2 * math.pi
H_MIN_ALT = 200e3        # ASSUMED lowest rendezvous altitude (outside the drag pass); sensitivity 150 and 300 km in the run
F_RP, F_RA = 3.504e6, RA_T      # ASSUMED representative final ellipse (r_p about 3,500 km, Revision 22 prediction header); flown ellipses are scored with their own values
# ASSUMED sweep ranges (named, Revision 22 item 4)
SWEEP = dict(ra_t_fracs=[0.8, 0.9, 1.0, 1.1, 1.2], rp_t_n=10, di_deg=[0.1, 0.5, 1.0, 2.0, 5.0], un_step_deg=5.0, circ_r_max_frac=1.3, circ_n=121, nr=4001)


# ------------------------------------------------------------------------------------------------ two-body primitives
def deg(x): return math.degrees(x)


def ell(rp, ra):
    """(a, e) from periapsis and apoapsis radii."""
    return 0.5 * (rp + ra), (ra - rp) / (ra + rp)


def radius(a, e, nu): return a * (1 - e * e) / (1 + e * np.cos(nu))


def vel_rt(a, e, nu):
    """(radial, transverse) velocity components on an ellipse at true anomaly nu."""
    h = math.sqrt(MU * a * (1 - e * e)); return (MU / h) * e * np.sin(nu), (MU / h) * (1 + e * np.cos(nu))


def mean_anomaly(e, nu):
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2)); return (E - e * np.sin(E)) % TWO_PI


def kepler_time(a, e, nu1, nu2):
    """Forward time (s) from true anomaly nu1 to nu2 (both rad; result in [0, T))."""
    n = math.sqrt(MU / a ** 3); return float(((mean_anomaly(e, nu2) - mean_anomaly(e, nu1)) % TWO_PI) / n)


def period(a): return TWO_PI * math.sqrt(a ** 3 / MU)


def rot(i, Om, om):
    """Perifocal -> inertial rotation matrix (columns P, Q, W)."""
    co, so, ci, si, cO, sO = math.cos(om), math.sin(om), math.cos(i), math.sin(i), math.cos(Om), math.sin(Om)
    P = np.array([cO * co - sO * so * ci, sO * co + cO * so * ci, so * si]); Q = np.array([-cO * so - sO * co * ci, -sO * so + cO * co * ci, co * si]); return np.column_stack((P, Q, np.cross(P, Q)))


def state_inertial(el, nu):
    """Inertial (r, v) of an element dict (a, e, i, Om, om) at true anomaly nu (scalar)."""
    a, e = el["a"], el["e"]; vr, vt = vel_rt(a, e, nu); r = radius(a, e, nu); R = rot(el["i"], el["Om"], el["om"])
    c, s = math.cos(nu), math.sin(nu); rh = R @ np.array([c, s, 0.0]); th = R @ np.array([-s, c, 0.0]); return r * rh, vr * rh + vt * th


def j2_rates(a, e, i):
    """Secular J2 rates (d omega/dt, d Omega/dt) in rad/s (SOURCED form, planner_rev21_position.j2_orbit per orbit divided by the period)."""
    d = POS.j2_orbit(a, e, i); T = period(a); return d[0] / T, d[1] / T


def advance(el, dt):
    """Element dict advanced by the secular J2 drift over dt seconds (a, e, i constant)."""
    wo, wO = j2_rates(el["a"], el["e"], el["i"]); out = dict(el); out["om"] = el["om"] + wo * dt; out["Om"] = el["Om"] + wO * dt; return out


# ------------------------------------------------------------------------------------------------ plane geometry
def relative_plane(F, T):
    """Relative geometry of two orbit planes. Returns (di, u_n): di = angle between the orbit normals (rad); u_n = argument, in F's perifocal frame measured from F's periapsis, of the
    line-of-nodes direction L = n_F x n_T (the other node is u_n + pi). Coplanar (di < 1e-12): (0.0, None)."""
    RF = rot(F["i"], F["Om"], F["om"]); RT = rot(T["i"], T["Om"], T["om"]); nf, nt = RF[:, 2], RT[:, 2]; c = float(np.dot(nf, nt)); di = math.acos(max(-1.0, min(1.0, c)))
    L = np.cross(nf, nt); nL = float(np.linalg.norm(L))
    if di < 1e-12 or nL < 1e-15: return 0.0, None
    L = L / nL; return di, math.atan2(float(np.dot(L, RF[:, 1])), float(np.dot(L, RF[:, 0]))) % TWO_PI


def nu_range(a, e, h_min=H_MIN_ALT, n=7201):
    """Admissible true anomalies on the descending leg: nu in [pi, 2pi) with r(nu) >= RJ + h_min. Returns (nu_lo, nu_hi) = (pi, last admissible)."""
    nu = np.linspace(math.pi, TWO_PI, n); ok = radius(a, e, nu) >= RJ + h_min; return math.pi, float(nu[ok][-1])


def nu_at_radius(a, e, r):
    """True anomaly on the descending leg (pi..2pi) at radius r, or None when r is outside [r_p, r_a]."""
    c = (a * (1 - e * e) / r - 1) / e if e > 0 else None
    if c is None or c < -1 - 1e-12 or c > 1 + 1e-12: return None
    return TWO_PI - math.acos(max(-1.0, min(1.0, c)))


# ------------------------------------------------------------------------------------------------ cost
def cost_at(Fae, nu, Tae, di=0.0, branch=-1):
    """|v_target - v_F| (m/s) at F's true anomaly nu, the target (a_t, e_t) crossing the same radius with the radial-velocity sign `branch` (-1 descending, +1 ascending), its plane tilted by di about
    the radial direction (the line of nodes at the node; di = 0 coplanar). NaN when the radius is outside the target's range. A circular target has vr = 0 and the branch is irrelevant."""
    (a, e), (at, et) = Fae, Tae; r = radius(a, e, nu); vr, vf = vel_rt(a, e, nu); pt = at * (1 - et * et)
    if et > 0:
        c = (pt / r - 1) / et
        if abs(c) > 1 + 1e-12: return float("nan")
        nut = branch * math.acos(max(-1.0, min(1.0, c))); vrt, vtt = vel_rt(at, et, nut)
    else:
        if abs(r - at) > 1e-6 * at: return float("nan")
        vrt, vtt = 0.0, math.sqrt(MU / at)
    return math.sqrt((vrt - vr) ** 2 + (vtt * math.cos(di) - vf) ** 2 + (vtt * math.sin(di)) ** 2)


def _best_over_radius(Fae, Tae, r_lo, r_hi, nr, di=0.0):
    """Minimum cost over crossing radii in [r_lo, r_hi] and both branches (grid then three zoom passes). Returns (cost, nu_r_rad, branch) or None."""
    a, e = Fae; best = None; lo, hi = r_lo, r_hi
    if hi < lo: return None
    for zoom in range(4):
        rs = np.linspace(lo, hi, nr if zoom == 0 else 201)
        for br in (-1, 1):
            for r in rs:
                nu = nu_at_radius(a, e, float(r))
                if nu is None: continue
                c = cost_at(Fae, nu, Tae, di, br)
                if c == c and (best is None or c < best[0] - 1e-12): best = (c, nu, br, float(r))
        if best is None: return None
        step = (hi - lo) / (nr if zoom == 0 else 200); lo, hi = max(r_lo, best[3] - 2 * step), min(r_hi, best[3] + 2 * step)
    return best[0], best[1], best[2]


def _crossing_limits(Fae, Tae, h_min=H_MIN_ALT):
    """Radius interval on F's descending leg that lies in the target's radius range and above the altitude floor."""
    (a, e), (at, et) = Fae, Tae; rp_f, ra_f = a * (1 - e), a * (1 + e); rp_t, ra_t = at * (1 - et), at * (1 + et)
    return max(rp_f, rp_t, RJ + h_min), min(ra_f, ra_t)


def optimal_coplanar(Fae, Tae, nr=2001, h_min=H_MIN_ALT):
    """Coplanar rendezvous: best crossing over the admissible arc. Returns dict(cost, nu_r_deg, branch, t_s) or None (no crossing: infeasible without a transfer)."""
    lo, hi = _crossing_limits(Fae, Tae, h_min); res = _best_over_radius(Fae, Tae, lo, hi, nr)
    if res is None: return None
    c, nu, br = res; a, e = Fae; return dict(cost=c, nu_r_deg=deg(nu), branch=br, t_s=kepler_time(a, e, math.pi, nu))


def apoapsis_cost(Fae, Tae, di=0.0):
    """Cost of meeting exactly at F's apoapsis (nu = 180 deg); NaN when the target's radius range does not contain r_a."""
    return cost_at(Fae, math.pi, Tae, di, -1)


def node_cost(Fae, Tae, di, u_n, h_min=H_MIN_ALT):
    """Non-coplanar: the meeting point is the node of F on the descending leg (u_n or u_n + pi, the one in [pi, 2pi)). Returns dict(cost, nu_r_deg, t_s, node_radius) or None (node radius outside the
    target's range, or below the altitude floor)."""
    a, e = Fae; nu = u_n if u_n >= math.pi else u_n + math.pi
    if nu >= TWO_PI: nu -= TWO_PI
    r = float(radius(a, e, nu))
    if r < RJ + h_min: return None
    cs = [cost_at(Fae, nu, Tae, di, b) for b in (-1, 1)]; cs = [c for c in cs if c == c]
    if not cs: return None
    return dict(cost=min(cs), nu_r_deg=deg(nu), t_s=kepler_time(a, e, math.pi, nu), node_radius=r)


def coaxial_crossing(Fae, Tae, h_min=H_MIN_ALT):
    """Coaxial coplanar orbits (same apse line, Delta omega = 0): crossing of r_F(nu) = r_T(nu) on the descending leg. Returns dict(cost, nu_r_deg) or None."""
    (a, e), (at, et) = Fae, Tae; p, pt = a * (1 - e * e), at * (1 - et * et); den = p * et - pt * e
    if abs(den) < 1e-12: return None
    c = (pt - p) / den
    if abs(c) > 1: return None
    nu = TWO_PI - math.acos(c); r = float(radius(a, e, nu))
    if r < RJ + h_min: return None
    vr, vt = vel_rt(a, e, nu); vrt, vtt = vel_rt(at, et, nu); return dict(cost=math.hypot(vrt - vr, vtt - vt), nu_r_deg=deg(nu))


# ------------------------------------------------------------------------------------------------ selftest
def selftest():
    ok = []
    def chk(name, cond, info=""): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + (" | " + info if info else ""))
    F = ell(F_RP, F_RA); a, e = F
    # 1. tangential apoapsis: circular target at r_a, cost = v_circ - v_a; optimum at nu = 180 (the only crossing)
    va = vel_rt(a, e, math.pi)[1]; vc = math.sqrt(MU / F_RA); o = optimal_coplanar(F, (F_RA, 0.0))
    chk("circular target at r_a: cost = v_circ - v_a", o and abs(o["cost"] - (vc - va)) < 1e-6, f"cost {o['cost']:.4f} v_c {vc:.3f} v_a {va:.3f} diff {vc - va:.3f}")
    chk("circular target at r_a: optimum nu_r = 180 deg", o and abs(o["nu_r_deg"] - 180) < 1e-6)
    # 2. identical coplanar orbits: cost 0 at every admissible nu on the descending leg
    lo, hi = nu_range(a, e); cs = [cost_at(F, nu, F, 0.0, -1) for nu in np.linspace(lo, hi, 9)]
    chk("identical orbits: cost 0", max(cs) < 1e-6, f"max {max(cs):.2e}"); o = optimal_coplanar(F, F)
    chk("identical orbits: optimal cost 0", o and o["cost"] < 1e-6)
    # 3. plane change: 2 v_theta sin(di/2) at any node (2 v sin(di/2) at apoapsis where v_r = 0)
    for dd in (0.5, 2.0, 5.0):
        di = math.radians(dd); c = cost_at(F, math.pi, F, di, -1); chk(f"plane change {dd} deg at apoapsis = 2 v_a sin(di/2)", abs(c - 2 * va * math.sin(di / 2)) < 1e-9, f"{c:.5f}")
    nu = math.radians(250.0); vr, vt = vel_rt(a, e, nu); di = math.radians(1.0); c = cost_at(F, nu, F, di, -1)
    chk("plane change at nu=250: 2 v_theta sin(di/2)", abs(c - 2 * vt * math.sin(di / 2)) < 1e-9, f"{c:.5f}")
    # 4. relative plane recovers di and the node lies in both planes; J2 drift consistent with POS
    Fe = dict(a=a, e=e, i=math.radians(93.5), Om=math.radians(10.0), om=math.radians(40.0)); Te = dict(a=a, e=e, i=math.radians(93.5), Om=math.radians(10.0 + 2.0), om=math.radians(100.0))
    di, un = relative_plane(Fe, Te); di_exact = 2 * math.asin(math.sin(math.radians(93.5)) * math.sin(math.radians(1.0)))
    chk("relative_plane di = 2 asin(sin i sin(dOm/2))", abs(di - di_exact) < 1e-12, f"{deg(di):.6f} deg")
    r_n, _ = state_inertial(Fe, un); RT = rot(Te["i"], Te["Om"], Te["om"]); chk("F's node point lies in the target plane", abs(float(np.dot(r_n, RT[:, 2]))) < 1e-6 * np.linalg.norm(r_n), f"{float(np.dot(r_n, RT[:, 2])):.2e} m")
    d1 = advance(Fe, 86400.0); wo, wO = j2_rates(a, e, Fe["i"]); chk("J2 advance equals rate x dt and POS per-orbit value", abs(d1["Om"] - Fe["Om"] - wO * 86400.0) < 1e-15 and abs(wO * period(a) - POS.j2_orbit(a, e, Fe["i"])[1]) < 1e-15)
    # 5. Kepler time apoapsis -> periapsis = T/2; monotone
    chk("Kepler time apoapsis to periapsis = T/2", abs(kepler_time(a, e, math.pi, TWO_PI) - period(a) / 2) < 1e-6 * period(a), f"{kepler_time(a, e, math.pi, TWO_PI):.3f} s")
    ts = [kepler_time(a, e, math.pi, nu) for nu in np.linspace(math.pi + 0.01, TWO_PI - 0.01, 50)]; chk("Kepler time monotone", all(np.diff(ts) > 0))
    # 6. state check: inertial energy and angular momentum equal the elements'
    r, v = state_inertial(Fe, 4.0); chk("inertial state: energy", abs(0.5 * np.dot(v, v) - MU / np.linalg.norm(r) + MU / (2 * a)) < 1e-6 * MU / a); chk("inertial state: |h|", abs(np.linalg.norm(np.cross(r, v)) - math.sqrt(MU * a * (1 - e * e))) < 1e-3)
    # 7. coaxial crossing equals an optimal search restricted to equal radii (target (rp+200 km, ra+400 km) at Delta omega = 0)
    T = ell(F_RP + 2e5, F_RA + 4e5); cx = coaxial_crossing(F, T); chk("coaxial crossing has equal radii", cx is None or abs(radius(a, e, math.radians(cx["nu_r_deg"])) - radius(T[0], T[1], math.radians(cx["nu_r_deg"]))) < 1e-3)
    # 8. optimal search vs brute force over nu
    best = min((cost_at(F, nu, T, 0.0, b), nu) for nu in np.linspace(lo, hi, 20001) for b in (-1, 1) if cost_at(F, nu, T, 0.0, b) == cost_at(F, nu, T, 0.0, b)); o = optimal_coplanar(F, T)
    chk("optimal_coplanar <= brute force over nu (1e-3 m/s)", o["cost"] <= best[0] + 1e-3, f"{o['cost']:.5f} vs {best[0]:.5f}")
    print("SELFTEST", "PASS" if all(ok) else "FAIL", f"({sum(ok)}/{len(ok)})"); return all(ok)


# ------------------------------------------------------------------------------------------------ item 2: target-family sweep and region map
def sweep(F=None, h_min=H_MIN_ALT):
    """Target-family sweep (ranges in SWEEP, all ASSUMED). Returns dict with circular curve, coplanar map and non-coplanar map."""
    F = F or ell(F_RP, F_RA); a, e = F; rp, ra = F[0] * (1 - F[1]), F[0] * (1 + F[1]); S = SWEEP; va = vel_rt(a, e, math.pi)[1]
    # circular targets, coplanar
    circ = []
    for rt in np.linspace(rp, S["circ_r_max_frac"] * ra, S["circ_n"]):
        o = optimal_coplanar(F, (float(rt), 0.0), h_min=h_min); cA = apoapsis_cost(F, (float(rt), 0.0))
        circ.append(dict(r_t_km=float(rt) / 1e3, cost=None if o is None else o["cost"], nu_r_deg=None if o is None else o["nu_r_deg"], t_s=None if o is None else o["t_s"], cost_apo=None if cA != cA else cA))
    # circular targets with a plane mismatch (Revision 22a): the meeting point is the node, so its radius must equal r_t; the node argument is then fixed by r_t.
    # For each (di, r_t): required node true anomaly, cost; plus the node-offset table (u_n scan: radial gap |r_F(node) - r_t| for r_t = r_a and r_t = mean radius).
    non = []
    for dd in S["di_deg"]:
        di = math.radians(dd)
        for rt in np.linspace(rp, S["circ_r_max_frac"] * ra, S["circ_n"]):
            nu = nu_at_radius(a, e, float(rt))
            if nu is None or radius(a, e, nu) < RJ + h_min: non.append(dict(di_deg=dd, r_t_km=float(rt) / 1e3, feasible=False)); continue
            c = cost_at(F, nu, (float(rt), 0.0), di, -1); non.append(dict(di_deg=dd, r_t_km=float(rt) / 1e3, feasible=c == c, cost=c, u_n_required_deg=deg(nu), cost_per_deg=c / dd))
    offs = []
    for lab, rt in (("r_a", ra), ("mean_radius", a)):
        for un in np.arange(180.0, 360.0, S["un_step_deg"]):
            r = float(radius(a, e, math.radians(un))); offs.append(dict(target=lab, r_t_km=rt / 1e3, u_n_deg=float(un), node_radius_km=r / 1e3, radial_gap_km=(r - rt) / 1e3))
    cop = [dict(c, cls=None) for c in circ]
    return dict(F=dict(rp_km=rp / 1e3, ra_km=ra / 1e3, a_km=a / 1e3, e=e, v_apo=va, v_circ_ra=math.sqrt(MU / ra), h_min_m=h_min), sweep=S, circular=circ, coplanar=cop, noncoplanar=non, node_offsets=offs)


def summarize_sweep(SW):
    """Region statistics, P22-1..P22-3 numbers."""
    F = SW["F"]; out = {}; cir = SW["circular"]; ra_km = F["ra_km"]
    p1 = [c for c in cir if abs(c["r_t_km"] - ra_km) < 1e-6]
    if p1: out["P22_1_cost_circ_at_ra"] = p1[0]["cost"]; out["P22_1_nu_deg"] = p1[0]["nu_r_deg"]; out["P22_1_vc_minus_va"] = F["v_circ_ra"] - F["v_apo"]
    inner = [c for c in cir if c["cost"] is not None and c["r_t_km"] < ra_km - 1e-6]; ref = out.get("P22_1_cost_circ_at_ra")
    out["P22_2_inner_n"] = len(inner); out["P22_2_apoapsis_optimal_count"] = sum(1 for c in inner if abs(c["nu_r_deg"] - 180) < 0.5)
    out["P22_2_cheaper_than_ra_match_count"] = sum(1 for c in inner if ref is not None and c["cost"] < ref); out["P22_2_min_cost"] = min(c["cost"] for c in inner) if inner else None
    out["P22_2_min_cost_r_t_km"] = min(inner, key=lambda c: c["cost"])["r_t_km"] if inner else None
    non = SW["noncoplanar"]; out["noncoplanar_by_di"] = {}
    for dd in SWEEP["di_deg"]:
        fe = [r for r in non if r["di_deg"] == dd and r["feasible"]]; ap = [r for r in fe if abs(r["r_t_km"] - ra_km) < 1e-6]
        out["noncoplanar_by_di"][str(dd)] = dict(n=sum(1 for r in non if r["di_deg"] == dd), feasible=len(fe), cost_min=min(r["cost"] for r in fe), min_cost_per_deg=min(r["cost_per_deg"] for r in fe),
                                                 cost_at_ra=ap[0]["cost"] if ap else None, cost_per_deg_at_ra=ap[0]["cost_per_deg"] if ap else None, cost_med=float(np.median([r["cost"] for r in fe])), cost_max=max(r["cost"] for r in fe))
    return out



# ------------------------------------------------------------------------------------------------ item 3: sensitivity to the apoapsis target
def apoapsis_sensitivity(rp=F_RP, ra_list_km=(6000, 6500, 7000, 7383, 7408, 8000, 9000), h_min=H_MIN_ALT):
    """Best rendezvous cost of fixed reference targets as the final apoapsis radius of F moves (r_p fixed). Targets are fixed in space (circular only, Revision 22a): circular at 7,408 km; circular at the nominal mean radius. `gap_km` = r_t - r_a(F) for the circular target at 7,408 km (a crossing needs r_a(F) >= r_t: below it no intersection exists). Also the tracking case: circular target at F's own r_a."""
    tg = dict(circ_7408=(RA_T, 0.0), circ_mean=(0.5 * (rp + RA_T), 0.0)); rows = []
    for ra_km in ra_list_km:
        F = ell(rp, ra_km * 1e3); row = dict(ra_F_km=ra_km, v_apo=vel_rt(F[0], F[1], math.pi)[1])
        for k, T in tg.items():
            o = optimal_coplanar(F, T, h_min=h_min); row[k] = None if o is None else dict(cost=o["cost"], nu_r_deg=o["nu_r_deg"])
        row["gap_km"] = (RA_T - ra_km * 1e3) / 1e3; trk = optimal_coplanar(F, (ra_km * 1e3, 0.0), h_min=h_min); row["circ_at_own_ra"] = None if trk is None else trk["cost"]; rows.append(row)
    d = (optimal_coplanar(ell(rp, 7408e3 + 1e3), (7409e3, 0.0))["cost"] - optimal_coplanar(ell(rp, 7408e3 - 1e3), (7407e3, 0.0))["cost"]) / 2.0
    return dict(rp_km=rp / 1e3, rows=rows, d_cost_per_km_tracking_circ=d)


# ------------------------------------------------------------------------------------------------ item 4: node controllability by arrival time
def node_controllability(tok, range_d=None, h_des_rp=F_RP):
    """Node drift of the final ellipse over an arrival-time range. Rates from planner_rev21_position (nominal J2 track); final-orbit rate (upper bound per day: the drift accelerates as the orbit shrinks) and
    campaign-mean rate (lower bound per day). Returns degrees. rel_to_circ: node rate of the final ellipse minus that of a circular orbit at 7,408 km with the same inclination (a target that also drifts)."""
    nom = POS.nominal(tok); tr = nom["track"]; i0 = tr.i0; ts = nom["days"] * 86400.0
    om_f, Om_f = tr.at(ts); Om0 = tr.Om0
    a, e = ell(h_des_rp, RA_T); wo, wO = j2_rates(a, e, i0); wO_c = j2_rates(RA_T, 0.0, i0)[1]; T = period(a)
    mean_rate = (Om_f - Om0) / ts           # rad/s, campaign mean (negative for prograde)
    out = dict(geom=tok, inc_deg=deg(i0), nominal_arrival_d=nom["days"], T_final_d=T / 86400.0, dOm_final_deg_per_d=deg(wO) * 86400.0, dOm_mean_campaign_deg_per_d=deg(mean_rate) * 86400.0,
               dom_final_deg_per_d=deg(wo) * 86400.0, dOm_per_final_orbit_deg=deg(wO) * T, dOm_circ_ra_deg_per_d=deg(wO_c) * 86400.0, rel_to_circ_deg_per_d=deg(wO - wO_c) * 86400.0)
    for lab, dt_d in (("window_3_orbits", 3 * T / 86400.0), ("one_day", 1.0)):
        out[lab] = dict(dt_d=dt_d, dOm_lo_deg=abs(deg(mean_rate)) * dt_d * 86400.0 / 86400.0 if False else abs(deg(mean_rate) * 86400.0) * dt_d, dOm_hi_deg=abs(deg(wO) * 86400.0) * dt_d)
    if range_d:
        lo, hi = range_d; w = hi - lo
        out["range_d"] = [lo, hi]; out["range_width_d"] = w; out["range_orbits"] = w / (T / 86400.0)
        out["range_dOm_lo_deg"] = abs(deg(mean_rate) * 86400.0) * w; out["range_dOm_hi_deg"] = abs(deg(wO) * 86400.0) * w
        s = math.sin(i0)
        for k in ("lo", "hi"):
            dO = math.radians(out[f"range_dOm_{k}_deg"]); out[f"range_di_{k}_deg"] = deg(2 * math.asin(max(-1, min(1, s * math.sin(dO / 2)))))      # relative inclination created by a node shift dOm at equal i
    return out


def tstar_ranges():
    """Feasible t* range of WP1 per geometry (rev20_data/ic/curves_*.json rows: min and max t_arr_d over feasible plans; SOURCED from those files)."""
    out = {}; d = pathlib.Path(POS.TD) / "rev20_data" / "ic"
    for g in ("ODY", "G1", "G0975"):
        try:
            J = json.load(open(d / f"curves_{g}.json")); ts = [r["t_arr_d"] for r in J["rows"] if r.get("feasible", True) and r.get("t_arr_d") == r.get("t_arr_d") and not r.get("beyond_table")]
            out[g] = (min(ts), max(ts), len(ts))
        except Exception as ex: out[g] = None
    return out


# ------------------------------------------------------------------------------------------------ item 5: flight scoring
def final_ellipse(d):
    """Final ellipse of a finished flight directory at its closing apoapsis: dict(a, e, i, Om, om, rp, ra, t_arr (s), tok). Closing pass = first pass with r_a <= 7,408 km (as planner_rev21_position); periapsis radius
    from that pass, apoapsis radius = its post-pass ra_m (SOURCED columns of passes.csv); i, Omega, omega osculating at that periapsis advanced by half a J2 orbit (same convention as WP21-B)."""
    import pandas as pd
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= RA_T))[0]
    if not len(hit): return None
    n = int(hit[0]); row = P.iloc[n]; rp = float(row.rp_m); rar = float(row.ra_m); a, e = ell(rp, rar); a0, e0, i0 = float(row.a_m), float(row.e), float(row.i_rad)
    dom, dOm = POS.j2_orbit(a0, e0, i0); T0 = period(a0)
    return dict(a=a, e=e, i=i0, Om=float(row.raan_rad) + 0.5 * dOm, om=float(row.omega_rad) + 0.5 * dom, rp=rp, ra=rar, t_arr=float(row.time_s) + 0.5 * T0, n_pass=n + 1)


def score_flight(d, tok, h_min=H_MIN_ALT, n_nu=721):
    """Generalized-event score of one flight against the three reference targets (all DERIVED). Targets lie in the NOMINAL plan's plane at the rendezvous time (inclination i0, Omega and omega from the nominal J2
    track, extended linearly past its end): R1 circular at 7,408 km; R3 circular at the nominal mean radius (a). Rendezvous point: nu_r on F in [180 deg, nu_max]; miss = 3-D distance to the target curve
    (phase free) or to the target vehicle (phase-locked); the best nu_r minimizes the miss; cost = |v_target - v_F| there (in the inertial frame). Also reported at the apoapsis (nu_r = 180 deg)."""
    Fe = final_ellipse(d)
    if Fe is None: return None
    nom = POS.nominal(tok); tr = nom["track"]; h = nom["hist"]; i0 = tr.i0; ts = nom["days"] * 86400.0
    a_n = float(h["a"][-1]); e_n = float(h["e"][-1]); rp_n = a_n * (1 - e_n); aN, eN = ell(rp_n, RA_T)
    lo, hi = nu_range(Fe["a"], Fe["e"], h_min); nus = np.linspace(lo, hi, n_nu); dts = np.array([kepler_time(Fe["a"], Fe["e"], math.pi, float(x)) for x in nus])
    om_rate_n, Om_rate_n = j2_rates(aN, eN, i0); om_n0, Om_n0 = tr.at(Fe["t_arr"])
    # F positions/velocities at t_arr + dts (J2-advanced frame)
    PF = np.zeros((n_nu, 3)); VF = np.zeros((n_nu, 3))
    for k in range(n_nu):
        el = advance(Fe, float(dts[k])); PF[k], VF[k] = state_inertial(el, float(nus[k]))
    out = dict(geom=tok, run=d.name, rp_km=Fe["rp"] / 1e3, ra_km=Fe["ra"] / 1e3, i_deg=deg(Fe["i"]), t_arr_d=Fe["t_arr"] / 86400.0, dt_d=(Fe["t_arr"] - ts) / 86400.0, n_pass=Fe["n_pass"])
    nt = np.linspace(0, TWO_PI, 1441)[:-1]
    def target_el(dt): return dict(a=aN, e=eN, i=i0, Om=om_n0 * 0 + Om_n0 + Om_rate_n * dt, om=om_n0 + om_rate_n * dt)
    def curve_miss(a_t, e_t, circ):
        """Phase-free nearest point on the target curve for each F point, with a local refinement (two zoom passes)."""
        best = np.full(n_nu, np.inf); bnu = np.zeros(n_nu); bv = np.zeros((n_nu, 3))
        for k in range(n_nu):
            el = dict(target_el(float(dts[k])), a=a_t, e=e_t); R = rot(el["i"], el["Om"], el["om"]); lo_, hi_ = 0.0, TWO_PI; grid = nt
            for z in range(3):
                r = radius(a_t, e_t, grid); pos = (R[:, :2] @ np.vstack((r * np.cos(grid), r * np.sin(grid)))).T; dd = np.linalg.norm(pos - PF[k], axis=1); j = int(np.argmin(dd)); c = float(grid[j]); step = float(grid[1] - grid[0]) if len(grid) > 1 else 0.0
                grid = np.linspace(c - 2 * step, c + 2 * step, 81)
            best[k] = dd[j]; bnu[k] = c % TWO_PI
            vr, vt = vel_rt(a_t, e_t, bnu[k]); bv[k] = R[:, :2] @ np.array([vr * math.cos(bnu[k]) - vt * math.sin(bnu[k]), vr * math.sin(bnu[k]) + vt * math.cos(bnu[k])])
        return best, bv
    def pack(miss, V, tag):
        j = int(np.argmin(miss)); return {f"{tag}_miss_best_km": float(miss[j]) / 1e3, f"{tag}_nu_r_deg": deg(float(nus[j])), f"{tag}_t_r_s": float(dts[j]), f"{tag}_cost_best_mps": float(np.linalg.norm(V[j] - VF[j])),
                                          f"{tag}_miss_apo_km": float(miss[0]) / 1e3, f"{tag}_cost_apo_mps": float(np.linalg.norm(V[0] - VF[0]))}
    for tag, a_t in dict(R1=RA_T, R3=0.5 * (rp_n + RA_T)).items():
        miss, V = curve_miss(a_t, 0.0, True); out.update(pack(miss, V, tag))
    # P22-5 split for R1: miss with the target plane frozen at t* (the plan's plane) vs timing part = arc between the nominal plane normals at t* and at the flight's arrival
    Rs = rot(i0, *(lambda o, O: (O, o))(*tr.at(ts))); Ra = rot(i0, Om_n0, om_n0); out["R1_timing_part_km"] = RA_T * math.acos(max(-1.0, min(1.0, float(np.dot(Rs[:, 2], Ra[:, 2]))))) / 1e3
    zs = np.array([abs(float(np.dot(PF[k], Rs[:, 2]))) for k in range(n_nu)]); rr = np.array([abs(float(np.linalg.norm(PF[k])) - RA_T) for k in range(n_nu)]); mf = np.sqrt(zs ** 2 + rr ** 2); out["R1_miss_frozen_plane_km"] = float(mf.min()) / 1e3
    # WP21-B position error for the same flight (apoapsis direction, deg and arc) for cross-reference
    try:
        pe = POS.flight_position(d, tok, None); out["wp21b_err_deg"] = pe["err_deg"]; out["wp21b_arc_km"] = pe["arc_km"]
    except Exception: out["wp21b_err_deg"] = None; out["wp21b_arc_km"] = None
    return out


def score_block(root, block, h_min=H_MIN_ALT):
    rows = []
    for dd in POS.block_dirs(root, block):
        k = POS.parse_id(dd.name)
        if not POS.finished(dd): continue
        try: r = score_flight(dd, k["tok"], h_min)
        except Exception as ex: rows.append(dict(run=dd.name, error=f"{type(ex).__name__}: {ex}")); continue
        if r is None: continue
        r.update(block=block, arm=k["arm"], level=k["level"], seed=k["seed"]); rows.append(r)
    return rows


# ------------------------------------------------------------------------------------------------ region map text and the run
def region_text(SW):
    """Region map of circular coplanar targets: r_t from r_p to 1.3 r_a. Class A = apoapsis (only r_t = r_a), B = crossing at another nu_r, C = no crossing (needs a transfer)."""
    ra_km = SW["F"]["ra_km"]; lines = []
    for c in SW["circular"]:
        cls = "C" if c["cost"] is None else ("A" if abs(c["nu_r_deg"] - 180) < 0.5 else "B")
        lines.append(f"r_t {c['r_t_km']:8.1f} km  {cls}  " + ("-" if c["cost"] is None else f"cost {c['cost']:7.1f} m/s  nu_r {c['nu_r_deg']:6.2f} deg  t {c['t_s']:7.0f} s"))
    return "\n".join(lines)


def cell_stats(rows, keys):
    out = {}
    for k in keys:
        v = np.array([r[k] for r in rows if r.get(k) is not None], float)
        out[k] = None if not len(v) else dict(n=int(len(v)), median=float(np.median(v)), p90=float(np.percentile(v, 90)), max=float(v.max()), min=float(v.min()))
    return out


def run(out_dir=None):
    out = pathlib.Path(out_dir) if out_dir else pathlib.Path(POS.TD) / "rev22_data"; out.mkdir(parents=True, exist_ok=True); J = lambda n, o: json.dump(o, open(out / n, "w"), indent=1)
    SW = sweep(); SU = summarize_sweep(SW)
    F0 = ell(F_RP, F_RA); dn = [(optimal_coplanar(F0, (float(r), 0.0), nr=51)["cost"], float(r)) for r in np.linspace(RJ + H_MIN_ALT + 1e3, F_RA * (1 - 1e-6), 1500)]
    SU["P22_2_dense_inner_min_cost"] = min(dn)[0]; SU["P22_2_dense_inner_min_r_t_km"] = min(dn)[1] / 1e3; SU["P22_2_dense_inner_n"] = len(dn)
    J("sweep.json", SW); J("sweep_summary.json", SU); open(out / "region_map_coplanar.txt", "w").write(region_text(SW) + "\n")
    hs = {}
    for hm in (150e3, 200e3, 300e3):
        s2 = summarize_sweep(sweep(h_min=hm)); hs[str(int(hm / 1e3))] = dict(P22_2_min_cost=s2["P22_2_min_cost"], P22_2_cheaper_count=s2["P22_2_cheaper_than_ra_match_count"], P22_2_apo_optimal=s2["P22_2_apoapsis_optimal_count"])
    J("sweep_hmin_sensitivity.json", hs)
    J("apoapsis_sensitivity.json", apoapsis_sensitivity()); rng = tstar_ranges(); J("tstar_ranges.json", rng)
    NC = {g: node_controllability(g, None if not rng.get(g) else (rng[g][0], rng[g][1])) for g in ("ODY", "G1", "G0975", "G0063", "G0390", "G1216")}; J("node_controllability.json", NC)
    allrows = {}
    for blk in ("rev20val", "rev21c", "rev21v"):
        if (POS.ROOT / blk).exists():
            rows = score_block(POS.ROOT, blk); allrows[blk] = rows; J(f"flights_{blk}.json", rows)
    keys = ["R1_miss_best_km", "R1_cost_best_mps", "R1_cost_apo_mps", "R1_miss_frozen_plane_km", "R1_timing_part_km", "R3_miss_best_km", "R3_cost_best_mps", "wp21b_arc_km", "wp21b_err_deg"]
    cells = {}
    for blk, rows in allrows.items():
        by = {}
        for r in rows:
            if "error" in r: continue
            by.setdefault((r["geom"], r["arm"], r["level"]), []).append(r)
        for k, v in by.items(): cells[f"{blk}|{k[0]}|{k[1]}|{k[2]}"] = dict(n=len(v), **cell_stats(v, keys))
    J("flight_cells.json", cells)
    print("written", out)
    return SW, SU, hs, NC, allrows, cells


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest": sys.exit(0 if selftest() else 1)
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        od = sys.argv[sys.argv.index("--out-dir") + 1] if "--out-dir" in sys.argv else None; run(od)
    else: print(__doc__)
