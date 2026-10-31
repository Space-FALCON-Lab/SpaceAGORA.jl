"""Revision 26, WP26-A (analysis only, Python, no Julia, no flights): arrival POSITION error of the guarded OLRT3 arm by counterfactual replay, and the levers L1-L5 of Revision 26.

Position model. A replayed path gives the pass history (periapsis times t_k, a_k = (ra_(k-1) + rp_k)/2, e_k) and the arrival time. The apoapsis direction at arrival is J2-propagated along that replayed history from the
geometry's initial elements (planner_rev21_position.J2Track, i = the source flight's own i history), plus the predictor residual (b) of the SOURCE flight (everything J2 does not explain: tesseral, drag on i/Omega, short period;
a tangent-plane vector at the target, ASSUMED independent of the policy). Target: the nominal apoapsis direction at t* (planner_rev21_position.nominal). Error = |A + C + B| in the tangent plane at the target (A timing,
C history, B residual; Revision 21 WP21-B definition). Validation: own-path replay against planner_rev21_position.flight_position.
Reuses planner_rev20_replay (guard machinery, GUARD_V3) read-only.
"""
import json, math, pathlib, sys, time
import numpy as np
import pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_rev20_replay as R
import planner_rev21_position as PS
import planner_rev21_guard as G

OUT = PS.TD / "rev26_data"
R.DEV_EXT = True


def src_residual(rid, tok, nom):
    """(b) of the SOURCE flight as a tangent vector at the target of the nominal: B = tproj(u_tgt, logmap(u_fp, u_fl)), plus the source's own error decomposition for the validation."""
    d = R.run_dir(rid); P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float)
    n = int(np.where(ra <= PS.target_radius(tok))[0][0] + 1); Pn = P.iloc[:n]
    t = Pn.time_s.values.astype(float); a = Pn.a_m.values.astype(float); e = Pn.e.values.astype(float); inc = Pn.i_rad.values.astype(float); om = Pn.omega_rad.values.astype(float); Om = Pn.raan_rad.values.astype(float)
    i0, om0, Om0 = PS.geometry_elements(tok); fl = PS.J2Track(t, a, e, inc, om0, Om0, i0); t_arr = float(t[-1] + 0.5 * fl.T[-1])
    dom, dOm = PS.j2_orbit(a[-1], e[-1], inc[-1]); u_fl = PS.apo_hat(om[-1] + 0.5 * dom, Om[-1] + 0.5 * dOm, inc[-1]); u_fp = fl.apo(t_arr)
    tstar = nom["days"] * 86400.0; u_tgt = nom["track"].apo(tstar)
    rp = Pn.rp_m.values.astype(float); ra_n = ra[:n]; ra_prev = np.concatenate(([float(R.get_tab(tok)["ra"][0])], ra_n[:-1])); a_form = 0.5 * (ra_prev + rp); e_form = (ra_prev - rp) / (ra_prev + rp)
    return dict(B=PS.tproj(u_tgt, PS.logmap(u_fp, u_fl)), Braw=PS.logmap(u_fp, u_fl), inc=inc, n=n, err_src=PS.deg(PS.ang(u_fl, u_tgt)), ka=a / a_form, de=e - e_form)


def replay_position(F, rid, tok, nom, ic, guard, tstar_d, noise="pass", seed=0, sr=None, arm="OLRT3"):
    """Replay OLRT3 (guard None or a guard dict) with the trace and return the replay result plus the position decomposition (degrees)."""
    r = R.replay(F, tok, arm, ic=ic, noise=noise, guard=guard, tstar_d=tstar_d, seed=seed, trace=True)
    if not r["reached"]:
        return r, None
    sr = sr or src_residual(rid, tok, nom)
    tr = r["trace"]; p = np.array([x[0] for x in tr]); ra = np.array([x[3] for x in tr]); t = np.array([x[4] for x in tr]); rp = np.array([x[5] for x in tr])
    ra_prev = np.concatenate(([float(R.get_tab(tok)["ra"][0])], ra[:-1])); ext = lambda v: np.concatenate((v, np.full(max(0, len(t) - len(v)), v[-1])))[:len(t)]
    a = 0.5 * (ra_prev + rp) * ext(sr["ka"]); e = (ra_prev - rp) / (ra_prev + rp) + ext(sr["de"])      # anchored to the source flight's osculating (a, e) at the periapsis event (own path: exact)
    inc = np.concatenate((sr["inc"], np.full(max(0, len(t) - len(sr["inc"])), sr["inc"][-1])))[:len(t)]
    i0, om0, Om0 = PS.geometry_elements(tok); fl = PS.J2Track(t, a, e, inc, om0, Om0, i0)
    t_arr = r["t_arr_d"] * 86400.0; tstar = tstar_d * 86400.0; nt = nom["track"]
    u_tgt = nt.apo(tstar); u_nA = nt.apo(t_arr); u_fp = fl.apo(t_arr)
    A = PS.logmap(u_tgt, u_nA); C = PS.tproj(u_tgt, PS.logmap(u_nA, u_fp)); B = sr["B"]
    Bt = PS.tproj(u_fp, sr["Braw"]); nb = float(np.linalg.norm(Bt)); u_fl = u_fp if nb < 1e-15 else math.cos(nb) * u_fp + math.sin(nb) * Bt / nb       # actual direction = own-history J2 direction moved by the source residual
    tot = PS.logmap(u_tgt, u_fl)
    s_hat = PS.tproj(u_tgt, nt.rate_vec(tstar)); s_hat = s_hat / np.linalg.norm(s_hat); c_hat = np.cross(u_tgt, s_hat)
    rate = float(np.linalg.norm(nt.rate_vec(tstar)))          # rad/s of the target apse along-track
    dOm = PS.deg((fl.at(t_arr)[1] - nt.at(tstar)[1] + math.pi) % (2 * math.pi) - math.pi)
    pos = dict(err_deg=PS.deg(PS.ang(u_fl, u_tgt)), along_deg=PS.deg(float(np.dot(tot, s_hat))), cross_deg=PS.deg(float(np.dot(tot, c_hat))), a_deg=PS.deg(float(np.linalg.norm(A))), b_deg=PS.deg(float(np.linalg.norm(B))),
               c_deg=PS.deg(float(np.linalg.norm(C))), dOm_deg=dOm, rate_deg_per_orbit=PS.deg(rate * 2 * math.pi * math.sqrt((0.5 * (7.408e6 + 3.5e6)) ** 3 / PS.MU)),
               tot_vec=tot.tolist(), s_hat=s_hat.tolist(), c_hat=c_hat.tolist(), rate=rate)
    return r, pos


def validate_own(geom, rids):
    nom = PS.nominal(geom); rows = []
    for rid in rids:
        F = R.get_flight(rid); sr = src_residual(rid, geom, nom)
        r, pos = replay_position(F, rid, geom, nom, None, None, nom["days"], sr=sr, arm=R.split_id(rid)["arm"])
        # flown position by the WP21-B code on the same flight (t* = nominal)
        fp = PS.flight_position(R.run_dir(rid), geom, nom["days"])
        rows.append((rid, fp["err_deg"], pos["err_deg"]))
    return rows


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        for g in ("ODY", "G1"):
            ids = R.dev_paths_ext(g)
            rows = validate_own(g, ids)
            d = np.array([x[2] - x[1] for x in rows]); print(g, len(rows), "own-path replay position minus flown (deg): median %.3f max abs %.3f" % (np.median(d), np.abs(d).max()), [(x[0][-9:], round(x[1], 2), round(x[2], 2)) for x in rows[:5]])


# ------------------------------------------------------------------------------------------------ L3: second guard objective on the accumulated apse rotation (position stage)
POS_DEFAULT = dict(pos_band_deg=0.3,     # ASSUMED: the position stage acts when the predicted along-track HISTORY term (c) at the predicted arrival exceeds this plus the funnel (deg); (b, tesseral ~0.3 deg) is not predictable
                   kC_deg_per_d=0.03,   # ASSUMED: funnel slope of the history-term band per day of predicted time to go
                   tol_t_d=0.10,         # ASSUMED: arrival-time residual (d) below which the time objective counts as met
                   M_frac=0.5,           # ASSUMED: the plan horizon is split in two segments at this fraction of the predicted passes to arrival; each segment gets its own tube tightening
                   fd_step=0.1,          # ASSUMED: finite-difference step of the Newton iteration (in units of shift_scale)
                   shift_scale=0.1,      # ASSUMED: level shift (ln q) per unit of the two shape parameters
                   hw=0.08,              # ASSUMED: half width (ln q) of the band the DP follows around the shifted level
                   nit=3,                # ASSUMED: Newton iterations
                   min_steps=12,         # ASSUMED: no position stage with fewer predicted passes to arrival (no authority left)
                   latch=True,           # ASSUMED: keep the last shape parameters while the history term is inside its funnel (against re-engaging)
                   warm=True,            # ASSUMED: start the Newton iteration from the previous shape parameters
                   max_solves=14)        # ASSUMED solves per apoapsis (v3: 12)


def predict_hist(tab, kh, ra0, rp0, t0, x_plan, rp_plan, term):
    """planner_rev20_replay.predict_arrival that also returns the per-pass history (periapsis time, apoapsis radius before the pass, periapsis radius) of the predicted passes."""
    ra = ra0; rp_prev = rp0; t = t0; nx = len(x_plan); T_, RA_, RP_ = [], [], []
    for m in range(400):
        xm = x_plan[min(m, nx - 1)]; rpm = rp_plan[min(m, nx - 1)]
        t += 0.5 * R.Tk(ra, rp_prev) + 0.5 * R.Tk(ra, rpm); T_.append(t); RA_.append(ra); RP_.append(rpm)
        dr = kh * R.d0_at(tab, ra) * math.exp(xm) / tab["q_req"]
        if term is not None:
            Reff = ra - R.RA_T + term["term_aim_below_m"]
            if 0 < Reff <= term["term_activate_steps"] * term["term_step_m"]:
                dr = Reff / max(1, math.ceil(Reff / term["term_step_m"]))
        ra = ra - dr
        if ra <= R.RA_T:
            return t + 0.5 * R.Tk(ra, rpm), len(T_), np.array(T_), np.array(RA_), np.array(RP_)
        rp_prev = rpm
    return float("nan"), 400, np.array(T_), np.array(RA_), np.array(RP_)


def make_pos_fn(tok, nom, sr, pp=None):
    """Closure bound to one replayed path (source residual sr supplies the osculating-(a, e) anchors and the inclination history)."""
    pp = {**POS_DEFAULT, **(pp or {})}; i0, om0, Om0 = PS.geometry_elements(tok); nt = nom["track"]; tab0 = R.get_tab(tok)

    def stage(c):
        gd = c["gd"]; P = c["P"]; tab = c["tab"]; lvl = c["lvl"]; A_ = c["A_"]; B_ = c["B_"]; L_ = c["L_"]; j = c["j"]; info = c["info"]
        n = len(lvl); cap = int(P["update_cap_passes"]); gain = float(P["update_gain"]); K2 = max(int(gd["bias_steps"]), cap); m_ = np.arange(n)
        corr = gain * c["bias"] * (np.minimum(m_, K2) - np.minimum(m_, cap)); lvl0 = lvl if c["lvl_alt"] is None else c["lvl_alt"]; lvl_g = lvl0 + corr
        A2 = A_ + lvl - lvl_g; B2 = B_ + lvl - lvl_g; kh = c["kap"] if c["kap_meas"] is None else c["kap_meas"]; term = P if gd["term_model"] else None
        tstar_pos = c["tstar_s"] - (gd["t_offset_d"].get(tok, 0.0) if isinstance(gd["t_offset_d"], dict) else float(gd["t_offset_d"])) * 86400.0 - info.get("aim_shift_s", 0.0); tt = c["tstar_s"]
        u_tgt = nt.apo(tstar_pos); s_hat = PS.tproj(u_tgt, nt.rate_vec(tstar_pos)); s_hat = s_hat / np.linalg.norm(s_hat); rate_dd = PS.deg(float(np.linalg.norm(nt.rate_vec(tstar_pos)))) * 86400.0
        ex = lambda v, k: np.concatenate((v, np.full(max(0, k - len(v)), v[-1])))[:k] if k > 0 else v[:0]
        pen = P["dp_burn_penalty_mps"] if gd.get("pen_g") is None else gd["pen_g"]
        width = float(np.max(B2 - A2)); calls = [0]
        ras = np.array(c["ras"]); rps = np.array(c["rps"]); ts = np.array(c["ts"]); jp = len(ts)          # passes flown so far
        ra_prev_past = np.concatenate(([float(tab0["ra"][0])], ras[:-1])) if jp else np.zeros(0)

        def evaluate(nu):
            x = lvl_g + nu; k0 = min(j - 1, tab["n"])
            rp = c["rp_now"] + np.array([R.tfun(tab, "rp", min(j + m, tab["n"])) - R.tfun(tab, "rp", k0) for m in range(n)]) - c["H"] * nu
            ta, N, Tn, RAn, RPn = predict_hist(tab, kh, c["ra_now"], c["rp_now"], ts[-1], x, rp, term)
            if not math.isfinite(ta):
                return None
            t_all = np.concatenate((ts, Tn)); ra_prev = np.concatenate((ra_prev_past, RAn)); rp_all = np.concatenate((rps, RPn)); K = len(t_all)
            ka = ex(sr["ka"], K); de = ex(sr["de"], K); inc = ex(sr["inc"], K)
            a = 0.5 * (ra_prev + rp_all) * ka; e = (ra_prev - rp_all) / (ra_prev + rp_all) + de
            tr = PS.J2Track(t_all, a, e, inc, om0, Om0, i0); u_fp = tr.apo(ta); E = PS.logmap(u_tgt, u_fp); Cv = PS.tproj(u_tgt, PS.logmap(nt.apo(ta), u_fp))
            return ta, N, PS.deg(float(np.dot(E, s_hat))), PS.deg(float(np.linalg.norm(E - np.dot(E, s_hat) * s_hat))), PS.deg(float(np.dot(Cv, s_hat)))

        ev_t = evaluate(c["nu_t"]) if c["nu_t"] is not None else None
        if ev_t is None or ev_t[1] < pp["min_steps"]:
            return c["nu_t"], c["eng_t"]
        ta_t, N0, al_t, cr_t, Cal_t = ev_t
        info.setdefault("pos_C", []).append((j, N0, round(Cal_t, 3), round(al_t, 3)))
        tau_d = max((ta_t - c["ts"][-1]) / 86400.0, 0.0); band_c = pp["pos_band_deg"] + pp["kC_deg_per_d"] * tau_d          # funnel: the prediction of the history term is uncertain far from the end
        e_prev = info.get("pos_e")
        if abs(Cal_t) <= band_c and not (pp["latch"] and e_prev is not None):
            return c["nu_t"], c["eng_t"]                                              # the history term is inside its funnel: the time guard (v3) decides
        M = max(1, int(pp["M_frac"] * N0)); seg = (np.arange(n) >= M).astype(int)

        x_ref = lvl_g + (c["nu_t"] if c["nu_t"] is not None else 0.0)                    # level path of the plan the time guard chose

        def plan(e):
            # shift actuator: the DP follows a narrow band (half width hw) around x_ref + e_seg inside the tube, so the burns realize the level shift at least dV
            cen = x_ref - lvl_g + np.where(seg == 0, e[0], e[1]) * pp["shift_scale"]
            A_t = np.clip(cen - pp["hw"], A2, B2); B_t = np.clip(cen + pp["hw"], A2, B2)
            B_t = np.maximum(B_t, A_t)
            calls[0] += 1
            return dpf_(A_t, B_t)
        dpf_ = lambda At, Bt: R.dpf(At, Bt, L_, pen)

        def resid(e):
            nu = plan(e)
            if nu is None:
                return None, None
            ev = evaluate(nu)
            if ev is None:
                return None, None
            ta, N, al, cr, Cal = ev
            return np.array([(ta - tt) / 86400.0 * rate_dd, Cal]), nu
        e = np.zeros(2) if (e_prev is None or not pp['warm']) else np.array(e_prev, float); r, nu = resid(e)
        if r is None:
            return c["nu_t"], c["eng_t"]
        best = (float(np.hypot(*r)), nu, e.copy()); r_t = np.array([(ta_t - tt) / 86400.0 * rate_dd, Cal_t])
        for it in range(int(pp["nit"])):
            if calls[0] + 2 > pp["max_solves"]:
                break
            J = np.zeros((2, 2)); ok = True
            for kx in range(2):
                e2 = e.copy(); e2[kx] += pp["fd_step"]; r2, _ = resid(e2)
                if r2 is None:
                    ok = False; break
                J[:, kx] = (r2 - r) / pp["fd_step"]
            if not ok or abs(np.linalg.det(J)) < 1e-9:
                break
            d = np.linalg.solve(J, -r); e = np.clip(e + d, -1.0, 1.0)
            if gd.get("pos_debug"):
                info.setdefault("pos_dbg", []).append((j, [round(float(v), 2) for v in r], np.round(J, 1).tolist(), [round(float(v), 3) for v in d]))
            if calls[0] >= pp["max_solves"]:
                break
            r, nu = resid(e)
            if r is None:
                break
            if float(np.hypot(*r)) < best[0]:
                best = (float(np.hypot(*r)), nu, e.copy())
            if float(np.hypot(*r)) < 0.15:
                break
        info["pos_solves"] = info.get("pos_solves", 0) + calls[0]
        info["pos_e"] = None if abs(Cal_t) <= 0.5 * band_c and abs(best[2][0]) + abs(best[2][1]) < 1e-9 else [float(v) for v in best[2]]
        info.setdefault("pos_log", []).append((j, N0, round(Cal_t, 2), round(cr_t, 2), round((ta_t - tt) / 86400.0, 2), round(best[0], 2), [round(float(v), 2) for v in best[2]]))
        if best[0] < float(np.hypot(*r_t)) - 1e-9:
            info["pos_engaged"] = info.get("pos_engaged", 0) + 1
            return best[1], True
        return c["nu_t"], c["eng_t"]
    return stage


# Guard v4 (Revision 26): FROZEN after the dev-only tuning (p_stage*.jsonl in rev26_data, 24 position-stage configurations on 27 dev paths, pass indexing, S-fixed start): guard v3 plus the position stage with the constants below.
# Selection rule (fixed before the freeze): 0 hard violations, window 27/27, then the highest within-1-deg fraction at ODY with L1, then the smallest median extra campaign dV over v3.
POS_V4 = dict(shift_scale=0.3, fd_step=0.3, hw=0.15, latch=False, warm=False)
GUARD_V4 = dict(G.GUARD_V3)          # the guard parameters are v3's; the position stage is attached per path with make_pos_fn(..., POS_V4)


def held_eval(geoms=("ODY", "G1"), noise="pass"):
    """The single frozen held test (seeds 13-28, natural paths of the held OLR flights): OLRT3 alone, guard v3, guard v4 (v3 + position stage), S-fixed start, t* = the nominal. Run once."""
    out = {}
    for g in geoms:
        ids = R.held_paths(g); ic, nom, ts = ctx(g); rows = {k: [] for k in ("alone", "v3", "v4")}
        for rid in ids:
            F = R.get_flight(rid); m = R.split_id(rid); sd = m["seed"] * 100; sr = src_residual(rid, g, nom)
            for k, gd in (("alone", None), ("v3", G.GUARD_V3), ("v4", dict(G.GUARD_V3, pos_fn=make_pos_fn(g, nom, sr, POS_V4)))):
                r, ps = replay_position(F, rid, g, nom, ic, gd, ts, noise=noise, seed=sd, sr=sr)
                if ps is None:
                    continue
                e1, tau = l1_error(ps, r["t_arr_d"] - ts, r["w_max_d"])
                rows[k].append(dict(rid=rid, dv=r["dv_campaign"], dt=r["t_arr_d"] - ts, win=bool(abs(r["t_arr_d"] - ts) <= r["w_max_d"]), hard=bool(r["hard_violation"]), apo=bool(r["ok_apo"]), err=ps["err_deg"], err_L1=e1, along=ps["along_deg"], cross=ps["cross_deg"], burns=r["burns"], geom=g))
        out[g] = rows
    return out


# ------------------------------------------------------------------------------------------------ node check v2 (Revision 26c): variant v4n = v4 + node check at 1.0 deg; v4 itself is not changed
NODE_TOL_DEG = 1.0      # the user's node tolerance


def make_aim_fn(tok, nom, sr, tol_deg=NODE_TOL_DEG):
    """Mirror of planner_node_v2_aim: forecast the node error at the predicted arrival of the planner's plan from the J2 track (flown history + predicted passes); inside the tolerance the aim stays at t*, else it moves toward
    the zero-error time ta - err/rate, clamped to [t* - W, t* + W] (W = 3 final periods = the window)."""
    i0, om0, Om0 = PS.geometry_elements(tok); nt = nom["track"]; tab0 = R.get_tab(tok); wrap = lambda x: (x + 180.0) % 360.0 - 180.0

    def aim(c):
        gd = c["gd"]; P = c["P"]; tab = c["tab"]; lvl = c["lvl"]; j = c["j"]; n = len(lvl); tstar = c["tstar_s"]
        W = R.W_MAX_PERIODS * R.Tk(R.RA_T, R.GUARD_W_RP_M); term = P if gd["term_model"] else None; kh = c["kap"] if c["kap_meas"] is None else c["kap_meas"]
        x = lvl + c["nu_ref"]; k0 = min(j - 1, tab["n"])
        rp = c["rp_now"] + np.array([R.tfun(tab, "rp", min(j + m, tab["n"])) - R.tfun(tab, "rp", k0) for m in range(n)]) - c["H"] * c["nu_ref"]
        ts = np.array(c["ts"]); ras = np.array(c["ras"]); rps = np.array(c["rps"]); ta, N, Tn, RAn, RPn = predict_hist(tab, kh, c["ra_now"], c["rp_now"], ts[-1], x, rp, term)
        if not math.isfinite(ta):
            return tstar
        ra_prev_past = np.concatenate(([float(tab0["ra"][0])], ras[:-1])); t_all = np.concatenate((ts, Tn)); ra_prev = np.concatenate((ra_prev_past, RAn)); rp_all = np.concatenate((rps, RPn)); K = len(t_all)
        ex = lambda v: np.concatenate((v, np.full(max(0, K - len(v)), v[-1])))[:K]
        a = 0.5 * (ra_prev + rp_all) * ex(sr["ka"]); e = (ra_prev - rp_all) / (ra_prev + rp_all) + ex(sr["de"]); tr = PS.J2Track(t_all, a, e, ex(sr["inc"]), om0, Om0, i0)
        tstar_nom = tstar - (gd["t_offset_d"].get(tok, 0.0) if isinstance(gd["t_offset_d"], dict) else float(gd["t_offset_d"])) * 86400.0
        Om_t = PS.deg(nt.at(tstar_nom)[1]); err = wrap(PS.deg(tr.at(ta)[1]) - Om_t); kk = int(np.clip(np.searchsorted(tr.hi, ta), 0, tr.n - 1)); rate = PS.deg(tr.dOm[kk] / tr.T[kk])
        info = c["info"]; info.setdefault("node_log", []).append((j, round(err, 3), round((ta - tstar) / 86400.0, 2)))
        if abs(err) <= tol_deg:
            return tstar
        if rate == 0.0:
            return tstar
        info["node_engaged"] = info.get("node_engaged", 0) + 1
        return float(min(max(ta - err / rate, tstar - W), tstar + W))
    return aim


# ------------------------------------------------------------------------------------------------ dev evaluation of a guard variant (S-fixed start, t* = the nominal of the S-fixed plan)
_CTX = {}


def ctx(tok):
    if tok not in _CTX:
        ic = R.get_ic(tok, "S-fixed"); nom = PS.nominal(tok)
        _CTX[tok] = (ic, nom, R.tstar_for(tok, ic))
    return _CTX[tok]


def l1_error(pos, dt_d, w_d):
    """L1: the rendezvous time t*' is chosen in [t* - W, t* + W] and within W of the arrival: minimize |tot - rate tau s_hat| over tau (days). Returns (error deg, tau d)."""
    tot = np.array(pos["tot_vec"]); s = np.array(pos["s_hat"]); rate_day = pos["rate"] * 86400.0
    lo = max(-w_d, dt_d - w_d); hi = min(w_d, dt_d + w_d)
    if lo > hi:
        return pos["err_deg"], 0.0
    al = float(np.dot(tot, s)); cr = float(np.linalg.norm(tot - al * s))
    # target moves along +s by rate*tau: error vector = tot - rate*tau*s ; along component al - rate*tau
    tau = float(np.clip(al / rate_day, lo, hi)) if rate_day > 0 else 0.0
    return PS.deg(math.hypot(al - rate_day * tau, cr)), tau


def eval_dev(guard, geoms=("ODY", "G1"), noise="pass", paths=None, base=None, pos=None):
    rows = []
    for g in geoms:
        ic, nom, ts = ctx(g)
        for rid in (paths or R.dev_paths_ext(g)):
            F = R.get_flight(rid); m = R.split_id(rid); sd = m["seed"] * 100
            sr = src_residual(rid, g, nom)
            gdp = guard if pos is None else dict(guard, pos_fn=make_pos_fn(g, nom, sr, pos))
            r, ps = replay_position(F, rid, g, nom, ic, gdp, ts, noise=noise, seed=sd, sr=sr)
            if ps is None:
                continue
            e1, tau = l1_error(ps, r["t_arr_d"] - ts, r["w_max_d"])
            rows.append(dict(geom=g, rid=rid, dv=r["dv_campaign"], t_arr=r["t_arr_d"], dt=r["t_arr_d"] - ts, win=bool(abs(r["t_arr_d"] - ts) <= r["w_max_d"]), hard=bool(r["hard_violation"]), apo=bool(r["ok_apo"]),
                             err=ps["err_deg"], along=ps["along_deg"], cross=ps["cross_deg"], a=ps["a_deg"], b=ps["b_deg"], c=ps["c_deg"], err_L1=e1, tau_L1=tau, rate_deg_d=ps["rate"] * 86400 * 57.29578, pos_eng=r["info"].get("pos_engaged", 0), burns=r["burns"]))
    return rows


def summ(rows, key="err"):
    out = {}
    for g in sorted(set(r["geom"] for r in rows)) + ["ALL"]:
        rr = [r for r in rows if g == "ALL" or r["geom"] == g]; e = np.array([r[key] for r in rr])
        out[g] = dict(n=len(rr), med=float(np.median(e)), p75=float(np.percentile(e, 75)), p90=float(np.percentile(e, 90)), within1=float((e < 1.0).mean()), win=float(np.mean([r["win"] for r in rr])), dv=float(np.median([r["dv"] for r in rr])))
    return out


if __name__ == "__main__" and sys.argv[1] == "base":
    rows = eval_dev(G.GUARD_V3)
    for k in ("err", "err_L1"):
        print(k, {g: {a: round(b, 2) for a, b in v.items()} for g, v in summ(rows, k).items()})
    for g in ("ODY", "G1"):
        rr = [r for r in rows if r["geom"] == g]
        print(g, "median |along| %.2f |cross| %.2f  a %.2f b %.2f c %.2f  rate %.2f deg/d  W %.2f d" % (np.median([abs(r["along"]) for r in rr]), np.median([abs(r["cross"]) for r in rr]), np.median([r["a"] for r in rr]), np.median([r["b"] for r in rr]), np.median([r["c"] for r in rr]), rr[0]["rate_deg_d"], 0.42))
    json.dump(rows, open(str(OUT / "dev_base_v3.json"), "w"), indent=1, default=float) if OUT.exists() else None


# ------------------------------------------------------------------------------------------------ dev tuning driver
_BASE = {}


def base_rows(kind):
    if kind not in _BASE:
        _BASE[kind] = eval_dev(None if kind == "alone" else G.GUARD_V3) if kind != "alone" else eval_dev_alone()
    return _BASE[kind]


def eval_dev_alone():
    rows = []
    for g in ("ODY", "G1"):
        ic, nom, ts = ctx(g)
        for rid in R.dev_paths_ext(g):
            F = R.get_flight(rid); m = R.split_id(rid); sr = src_residual(rid, g, nom)
            r, ps = replay_position(F, rid, g, nom, ic, None, ts, seed=m["seed"] * 100, sr=sr)
            if ps is None:
                continue
            e1, tau = l1_error(ps, r["t_arr_d"] - ts, r["w_max_d"])
            rows.append(dict(geom=g, rid=rid, dv=r["dv_campaign"], dt=r["t_arr_d"] - ts, win=bool(abs(r["t_arr_d"] - ts) <= r["w_max_d"]), hard=bool(r["hard_violation"]), apo=bool(r["ok_apo"]), err=ps["err_deg"], err_L1=e1, along=ps["along_deg"], cross=ps["cross_deg"], pos_eng=0, burns=r["burns"]))
    return rows


def tune_pos(configs, name, guard=None):
    OUT.mkdir(parents=True, exist_ok=True)
    v3 = base_rows("v3"); alone = base_rows("alone")
    for pp in configs:
        t0 = time.time(); rows = eval_dev(guard or G.GUARD_V3, pos=pp)
        d3 = [a["dv"] - b["dv"] for a, b in zip(rows, v3)]; d0 = [a["dv"] - b["dv"] for a, b in zip(rows, alone)]
        row = dict(pos=pp, guard=None if guard is None else {k: v for k, v in guard.items() if k not in G.GUARD_V3 or guard[k] != G.GUARD_V3[k]}, n=len(rows), win=sum(r["win"] for r in rows), hard=sum(r["hard"] for r in rows), apo=sum(r["apo"] for r in rows),
                   err=summ(rows, "err"), err_L1=summ(rows, "err_L1"), extra_vs_v3=float(np.median(d3)), extra_vs_alone=float(np.median(d0)), pos_eng=float(np.mean([r["pos_eng"] for r in rows])), wall=time.time() - t0)
        with open(OUT / f"{name}.jsonl", "a") as fh:
            fh.write(json.dumps(row, default=float) + "\n")
        e, l = row["err"], row["err_L1"]
        print("L3", pp, "| win %d/%d hard %d apo %d | err med ODY %.2f G1 %.2f within1 ODY %.2f G1 %.2f ALL %.2f | L1+L3 med ODY %.2f G1 %.2f within1 ALL %.2f | extra dV vs v3 %.2f vs alone %.2f | %.0fs" % (
            row["win"], row["n"], row["hard"], row["apo"], e["ODY"]["med"], e["G1"]["med"], e["ODY"]["within1"], e["G1"]["within1"], e["ALL"]["within1"], l["ODY"]["med"], l["G1"]["med"], l["ALL"]["within1"], row["extra_vs_v3"], row["extra_vs_alone"], row["wall"]), flush=True)


if __name__ == "__main__" and sys.argv[1] == "held":
    res = held_eval(); OUT.mkdir(parents=True, exist_ok=True); json.dump(res, open(str(OUT / "held_v4.json"), "w"), indent=1, default=float)
    for g, rows in res.items():
        for k, rr in rows.items():
            e = np.array([r["err"] for r in rr]); e1 = np.array([r["err_L1"] for r in rr])
            print(g, k, "n", len(rr), "win %d hard %d apo %d | err median %.2f within1 %.2f | L1 median %.2f within1 %.2f | median campaign dV %.2f" % (sum(r["win"] for r in rr), sum(r["hard"] for r in rr), sum(r["apo"] for r in rr), np.median(e), (e < 1).mean(), np.median(e1), (e1 < 1).mean(), np.median([r["dv"] for r in rr])))
    sys.exit(0)
if __name__ == "__main__" and sys.argv[1] == "tune":
    tune_pos(eval(sys.argv[3]), sys.argv[2])


# ------------------------------------------------------------------------------------------------ flown guarded flights: L1 and the node check evaluated on the measured flights (no replay)
def flown_vectors(d, tok, nom=None):
    nom = nom or PS.nominal(tok); tstar = nom["days"] * 86400.0; rt = PS.target_radius(tok)
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= rt))[0]
    if not len(hit):
        return None
    n = int(hit[0] + 1); Pn = P.iloc[:n]; t = Pn.time_s.values.astype(float); a = Pn.a_m.values.astype(float); e = Pn.e.values.astype(float); inc = Pn.i_rad.values.astype(float); om = Pn.omega_rad.values.astype(float); Om = Pn.raan_rad.values.astype(float)
    i0, om0, Om0 = PS.geometry_elements(tok); fl = PS.J2Track(t, a, e, inc, om0, Om0, i0); t_arr = float(t[-1] + 0.5 * fl.T[-1]); nt = nom["track"]
    dom, dOm = PS.j2_orbit(a[-1], e[-1], inc[-1]); u_fl = PS.apo_hat(om[-1] + 0.5 * dom, Om[-1] + 0.5 * dOm, inc[-1]); u_tgt = nt.apo(tstar); tot = PS.logmap(u_tgt, u_fl)
    rv = PS.tproj(u_tgt, nt.rate_vec(tstar)); rate = float(np.linalg.norm(rv)); w_d = 3 * fl.T[-1] / 86400.0
    wrap = lambda x: (x + math.pi) % (2 * math.pi) - math.pi; om_t, Om_t = nt.at(tstar)
    out = dict(err=PS.deg(PS.ang(u_fl, u_tgt)), dt=(t_arr - tstar) / 86400.0, w_d=w_d, rate_deg_d=PS.deg(rate) * 86400.0, dOm_deg=PS.deg(wrap(Om[-1] + 0.5 * dOm - Om_t)), dom_deg=PS.deg(wrap(om[-1] + 0.5 * dom - om_t)), inc=PS.deg(i0))
    if rate > 1e-12:
        s = rv / rate; al = float(np.dot(tot, s)); cr = float(np.linalg.norm(tot - al * s)); out.update(along=PS.deg(al), cross=PS.deg(cr))
        lo = max(-w_d, out["dt"] - w_d); hi = min(w_d, out["dt"] + w_d); rd = PS.deg(rate) * 86400.0
        tau = float(np.clip(PS.deg(al) / rd, lo, hi)) if lo <= hi else 0.0
        out["err_L1"] = math.hypot(PS.deg(al) - rd * tau, PS.deg(cr)); out["tau"] = tau
    else:
        out.update(along=float("nan"), cross=out["err"], err_L1=out["err"], tau=0.0)
    return out


def flown_table(blocks=("rev20val", "rev21v", "rev23geo", "rev25low")):
    rows = []
    for b in blocks:
        for d in PS.block_dirs(PS.ROOT, b):
            k = PS.parse_id(d.name)
            if "W" not in k["arm"] or not PS.finished(d):
                continue
            try:
                v = flown_vectors(d, k["tok"])
            except Exception as ex:
                continue
            if v:
                v.update(block=b, run=d.name, geom=k["tok"], arm=k["arm"], level=k["level"], seed=k["seed"]); rows.append(v)
    return rows


if __name__ == "__main__" and sys.argv[1] == "flown":
    rows = flown_table(); OUT.mkdir(parents=True, exist_ok=True); json.dump(rows, open(str(OUT / "flown_guarded_position.json"), "w"), indent=1, default=float)
    cells = {}
    for r in rows:
        cells.setdefault((r["geom"], r["arm"], round(r["inc"], 1)), []).append(r)
    for k, rs in sorted(cells.items()):
        e = np.array([r["err"] for r in rs]); e1 = np.array([r["err_L1"] for r in rs]); dO = np.array([abs(r["dOm_deg"]) for r in rs])
        print(k, "n", len(rs), "median err %.2f -> L1 %.2f ; within1 %d -> %d ; median |along| %.2f |cross| %.2f |dOmega| %.2f ; rate %.2f deg/d" % (np.median(e), np.median(e1), (e < 1).sum(), (e1 < 1).sum(), np.nanmedian([abs(r["along"]) for r in rs]), np.median([r["cross"] for r in rs]), np.median(dO), np.median([r["rate_deg_d"] for r in rs])))


# ------------------------------------------------------------------------------------------------ maneuver cost of what the levers cannot remove (Revision 22 / Revision 20 tools)
def maneuver_costs(rows, target_deg=0.9):
    """Per flown flight: (i) node: relative inclination di = 2 asin(sin i sin(dOmega/2)) needed to bring the cross-track error to target_deg, cost = 2 v_a sin(di/2) at the apoapsis (the planner_rev22_rendezvous plane-change term, exact
    at an apoapsis node); (ii) along-track: the residual along-track error after L1 expressed as a time shift D = residual/rate and costed with planner_rev20_window.phasing (shift of the arrival by D, extra over the first burn;
    side analysis, outside the Revision 20d scope)."""
    import planner_rev20_window as W
    ra, rp = 7.408e6, 3.5e6; a = 0.5 * (ra + rp); va = math.sqrt(PS.MU * (2 / ra - 1 / a)); rp_des = rp + 145e3
    tab = W.phasing_table(ra, rp, rp_des, rp + 70e3, rp + 1000e3 - 3.5e6 + rp - 0.0) if False else W.phasing_table(ra, rp, rp_des, 3.5e6 + 70e3, 3.5e6 + 1000e3)
    out = []
    for r in rows:
        cross_need = max(0.0, r["cross"] - target_deg)
        di = math.radians(cross_need); cost_node = 2 * va * math.sin(di / 2) if cross_need > 0 else 0.0                    # plane change that removes the cross-track excess (cross-track angle ~ relative inclination at the apoapsis)
        resid_al = max(0.0, abs(r["along"]) - target_deg - r["rate_deg_d"] * r["w_d"]); D = resid_al / r["rate_deg_d"] * 86400.0 if r["rate_deg_d"] > 0 else 0.0
        cph = 0.0 if D == 0 else min(W.phasing_cost(tab, D)[0], W.phasing_cost(tab, -D)[0])
        out.append(dict(run=r["run"], geom=r["geom"], cost_node=cost_node, cost_phasing=cph if math.isfinite(cph) else float("inf"), D_s=D))
    return out


# ------------------------------------------------------------------------------------------------ v4n: dev, then the single held run
def run_variants(paths_fn, geoms=("ODY", "G1"), noise="pass"):
    out = {}
    for g in geoms:
        ic, nom, ts = ctx(g); rows = {k: [] for k in ("alone", "v4", "v4n")}
        for rid in paths_fn(g):
            F = R.get_flight(rid); m = R.split_id(rid); sd = m["seed"] * 100; sr = src_residual(rid, g, nom)
            for k, gd in (("alone", None), ("v4", dict(G.GUARD_V3, pos_fn=make_pos_fn(g, nom, sr, POS_V4))), ("v4n", dict(G.GUARD_V3, pos_fn=make_pos_fn(g, nom, sr, POS_V4), aim_fn=make_aim_fn(g, nom, sr)))):
                r, ps = replay_position(F, rid, g, nom, ic, gd, ts, noise=noise, seed=sd, sr=sr)
                if ps is None:
                    continue
                e1, tau = l1_error(ps, r["t_arr_d"] - ts, r["w_max_d"])
                rows[k].append(dict(rid=rid, dv=r["dv_campaign"], dt=r["t_arr_d"] - ts, w=r["w_max_d"], win=bool(abs(r["t_arr_d"] - ts) <= r["w_max_d"]), hard=bool(r["hard_violation"]), apo=bool(r["ok_apo"]), err=ps["err_deg"], err_L1=e1, along=ps["along_deg"], cross=ps["cross_deg"], dOm=ps["dOm_deg"],
                                    node_eng=r["info"].get("node_engaged", 0), burns=r["burns"], infeas=r["info"].get("infeasible", 0)))
        out[g] = rows
    return out


def report_variants(res):
    for g, rows in res.items():
        for k, rr in rows.items():
            e = np.array([x["err"] for x in rr]); e1 = np.array([x["err_L1"] for x in rr])
            print(g, k, "n", len(rr), "win %d hard %d apo %d | within1 raw %d L1 %d | err med %.2f L1 med %.2f | median |node| %.2f |along| %.2f cross %.2f | dV %.2f | node engaged paths %d | infeasible paths %d" % (
                sum(x["win"] for x in rr), sum(x["hard"] for x in rr), sum(x["apo"] for x in rr), (e < 1).sum(), (e1 < 1).sum(), np.median(e), np.median(e1), np.median([abs(x["dOm"]) for x in rr]), np.median([abs(x["along"]) for x in rr]),
                np.median([abs(x["cross"]) for x in rr]), np.median([x["dv"] for x in rr]), sum(1 for x in rr if x["node_eng"] > 0), sum(1 for x in rr if x["infeas"] > 0)))
        a, b = rows["v4"], rows["v4n"]
        print(g, "window lost by v4n relative to v4:", [(x["rid"][-9:], round(y["dt"], 2), round(y["w"], 2)) for x, y in zip(a, b) if x["win"] and not y["win"]], "| extra dV v4n-v4 median %.2f" % np.median([y["dv"] - x["dv"] for x, y in zip(a, b)]))


if __name__ == "__main__" and sys.argv[1] in ("v4n_dev", "v4n_held"):
    res = run_variants(R.dev_paths_ext if sys.argv[1] == "v4n_dev" else R.held_paths); OUT.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(str(OUT / (sys.argv[1] + ".json")), "w"), indent=1, default=float); report_variants(res)
