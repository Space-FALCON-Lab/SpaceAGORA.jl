"""Revision 20, WP1 (analysis only, Python, no Julia, no flights): fuel-optimal initial conditions for an aerobraking campaign that ends at a rendezvous apoapsis.

Decision variables (see tasks/.../results_rev20_ic.md section 0 for the amendments to the plan's provisional list):
  x0     = ln q of the first counted pass (pass 1, the start level); the initial periapsis altitude follows: h_p0 = H_REF - H_1 (x0 - x_ref) (first order).
  nu_pre = x0 - x_ref, a tangential burn at the STARTING apoapsis before pass 1, priced w_pre |nu_pre| with w_pre = H_1 / (ratio lever(a0, e0)) (DERIVED, planner_lib.burn_weight at table row 0).
  The campaign DP of ol_schedule starts AFTER pass 1 (steps are passes 2..): nu_{-1} = 0 relative to the calibrated pass-1 level; there is no walk-in. Here the initial offset s = nu_pre becomes the state before step 0.
Scenarios: "ref" (reference start x_ref, no pre-burn: planner_lib.ol_schedule), "S-fixed" (capture periapsis given at x_ref; any change pays w_pre |s|), "S-free" (x0 free at no campaign cost).

free_start_dp is dp_exact (planner_lib) / dp_lin (planner_rev18_offline) with a priced initial step. Exactness: the objective w_pre |s| + lin0 s + sum_k [L_k |nu_k - nu_{k-1}| + lin_k nu_k] with box constraints
A_k <= nu_k <= B_k, A0 <= s <= B0 is a convex piecewise-linear function of (s, nu_0, ...). The value function f_k(nu) of the dynamic program is convex piecewise linear; the infimal convolution with L_k |.| only clips slopes to [-L_k, L_k],
so its breakpoints are a subset of those of f_{k-1}; the box constraints add only the points A_k, B_k; the initial function has breakpoints {0, A0, B0}. The minimum of a convex piecewise-linear function over an interval is attained at a breakpoint or an
endpoint, therefore the state set {0} U {A_k} U {B_k} U {A0, B0} suffices (same argument as dp_exact). With a fixed charge per burn (pen > 0) the problem is not convex; the vertex set is then a heuristic exactly as in dp_exact (all headline results use pen = 0, the
model-world convention of Revisions 18 and 19).

Usage: python3 -I planner_rev20_ic.py --selftest | --run [--out-dir DIR]"""
import argparse, hashlib, json, math, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev18_offline as O
import planner_analysis as AN

INF = math.inf
MU = O.MU
RA_T = AN.RA_TARGET                        # SOURCED: target apoapsis radius (7,408 km)
P = O.P
Q_FIRST = O.Q_FIRST                        # MEASURED calibration target of the first pass (0.212 W/cm2)
X_REF = math.log(Q_FIRST)
H_REF = 105319.6                           # SOURCED geometry_rank_config.toml h_ref_m (Odyssey pass-1 altitude); ASSUMED to transfer to the other geometries
HORIZON_D = 130.0                          # SOURCED planner_config.toml horizon_days
QUANT_PERIODS = 1.0                        # ASSUMED: arrival is quantized to the final-orbit period, so a window s is evaluated as max(s, QUANT_PERIODS x T_final)
W_MAX_PERIODS = 3.0                        # ASSUMED (Revision 18): W_max = 3 final-orbit periods
S_GRID_D = (0.0, "WMAX", 1.0, 2.0)         # ASSUMED slack set of the task (days)
KAPPA_G0975 = 1.0                          # ASSUMED: no flight calibration exists for the generic table (preview_no_field), so the uncalibrated table is used
GEOMS = ("ODY", "G1", "G0975")
TD = pathlib.Path(PLB.task_dir())
LAM_MAX = 50.0                             # m/s per orbit; bracket limit before the pinned plan is tried (NUMERICAL)
SHRUNK_LO = None


# ----------------------------------------------------------------------------------------------------- inputs

def load_tab(token):
    if token in ("ODY", "G1"):
        tab = O.load_tab(token)
    elif token in ("G0975", "G0063", "G0390", "G1216"):                  # generic tables (Revision 21 adds the other three)
        g = json.load(open(TD / "planner_ground_tables_generic.json"))["geometries"][token]
        tab = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in g.items()}
    else:
        raise ValueError(token)
    tab["token"] = token
    return tab


def load_kappa():
    """SOURCED: rev19_data/rev18_calibrate.json kappa_med per geometry (MEASURED from the flights); G0975 ASSUMED 1."""
    d = json.load(open(TD / "rev19_data" / "rev18_calibrate.json"))["geometry"]
    k = {g["geom"]: float(g["kappa_med"]) for g in d}
    k["G0975"] = KAPPA_G0975
    for g in ("G0063", "G0390", "G1216"):
        k[g] = 1.0                         # ASSUMED: no flight calibration (Revision 21)
    return k


def first_burn_dv(ra, rp, rp_des):
    """Tangential periapsis-raise burn at apoapsis: dv = v_a(a') - v_a(a) (same as planner_rev18_flights.first_burn_dv)."""
    va = lambda r_p: math.sqrt(MU * (2 / ra - 2 / (ra + r_p)))
    return va(rp_des) - va(rp)


def w_pre_of(tab):
    """DERIVED: m/s per ln unit of a burn at the starting apoapsis (table row 0 is the start orbit)."""
    return PLB.burn_weight(float(tab["H"][0]), P["burn_response_ratio"], MU, float(tab["a"][0]), float(tab["e"][0]))


# ----------------------------------------------------------------------------------------------------- the DP with a priced initial step

def free_start_dp(A, B, L, w_pre, lin=None, A0=-INF, B0=INF, lin0=0.0, tol=1.0e-9, pen=0.0):
    """dp_exact with a priced initial offset s (the pre-burn nu_pre before step 0). Minimizes w_pre |s| + lin0 s + sum_k [L_k |nu_k - nu_{k-1}| + lin_k nu_k] (nu_{-1} = s), A_k <= nu_k <= B_k, A0 <= s <= B0.
    w_pre = inf: s = 0 forced (A0, B0, lin0 ignored): identical arithmetic to planner_lib.dp_exact (lin None) / planner_rev18_offline.dp_lin (lin given). w_pre = 0: free start. pen: fixed cost per burn (steps; the
    pre-burn pays it too when s != 0) as in dp_exact. Returns dict(nu, s, obj) or None if infeasible."""
    A = [float(x) for x in A]; B = [float(x) for x in B]; L = [float(x) for x in L]; n = len(A)
    finite = math.isfinite(w_pre)
    extra = [x for x in (A0, B0) if math.isfinite(x)] if finite else []
    val = sorted(set([0.0] + A + B + extra)); ns = len(val); i0 = val.index(0.0)
    if not finite:
        prev = [INF] * ns; prev[i0] = 0.0
    else:
        prev = [(w_pre * abs(v) + lin0 * v + (pen if v != 0.0 else 0.0)) if (A0 - tol <= v <= B0 + tol) else INF for v in val]
    parent = [None] * n
    for k in range(n):
        Lk = L[k]; cur = prev[:]; arg = list(range(ns))
        for j in range(1, ns):
            c = cur[j - 1] + Lk * (val[j] - val[j - 1])
            if c < cur[j]:
                cur[j] = c; arg[j] = arg[j - 1]
        for j in range(ns - 2, -1, -1):
            c = cur[j + 1] + Lk * (val[j + 1] - val[j])
            if c < cur[j]:
                cur[j] = c; arg[j] = arg[j + 1]
        if pen > 0.0:
            for j in range(ns):
                if prev[j] <= pen + cur[j]:
                    cur[j] = prev[j]; arg[j] = j
                else:
                    cur[j] = pen + cur[j]
        if lin is not None:
            lk = float(lin[k])
            for j in range(ns):
                cur[j] = cur[j] + lk * val[j]
        for j in range(ns):
            if val[j] < A[k] - tol or val[j] > B[k] + tol:
                cur[j] = INF
        parent[k] = arg; prev = cur
    j = min(range(ns), key=lambda x: prev[x])
    if not math.isfinite(prev[j]):
        return None
    obj = float(prev[j]); nu = [0.0] * n
    for k in range(n - 1, -1, -1):
        nu[k] = val[j]; j = parent[k][j]
    return dict(nu=np.array(nu), s=float(val[j]), obj=obj)


# ----------------------------------------------------------------------------------------------------- the case: one geometry, scenario, options

class Case:
    """One (geometry, scenario) planning problem. scenario in {'ref', 'S-fixed', 'S-free'}; horizon 'arrival' (DP over passes 2..N, the flown cost; fixed point in N) or 'full' (DP over the whole table plus tail, as
    ol_schedule; dV counted through the arrival pass). delta_ln: x0 must stay at least delta below the upper (reserve-shrunk) tube edge. pen: fixed cost per burn."""

    def __init__(self, token, scenario, h_des_m=250e3, fb_term=False, kappa=None, horizon="arrival", delta_ln=0.0, pen=0.0, ref_dk=None):
        self.token = token; self.tab = load_tab(token); self.scenario = scenario; self.h_des = float(h_des_m); self.fb_term = bool(fb_term)
        self.kappa = float(load_kappa()[token] if kappa is None else kappa); self.horizon = horizon; self.delta = float(delta_ln); self.pen = float(pen)
        tab = self.tab; inp = PLB.ol_inputs(dict(), tab, P, X_REF); self.passes = inp["passes"]; self.level = np.array(inp["level"]); self.lo = np.array(inp["lo"]); self.hi = np.array(inp["hi"]); self.w = np.array(inp["w"])
        self.tw = PLB.tie_weights(self.w, P["dp_tie_eps"]); self.A = self.lo - self.level; self.B = self.hi - self.level
        self.lo_s = float(self.lo[0]); self.hi_s = float(self.hi[0]); self.hi0 = math.log(P["q_cap_mean_w_cm2"])        # shrunk edges (constants over passes); hi0 = physical upper edge
        self.w_pre_full = w_pre_of(tab)
        self.w_pre = {"ref": INF, "S-fixed": self.w_pre_full, "S-free": 0.0}[scenario]
        self.A0 = self.lo_s - X_REF; self.B0 = self.hi_s - self.delta - X_REF
        self.d0 = np.concatenate((np.diff(-tab["ra"]), [float(tab["ra"][-2] - tab["ra"][-1])]))
        self.dk = None                    # reference drop weights (set after the lambda = 0 plan)
        self.memo = {}
        self.N_guess = None
        p0 = self.solve(0.0); self._set_dk(p0)
        self.nominal = p0

    # --- helpers
    def _set_dk(self, plan):
        x = plan["x"]; n = self.tab["n"]
        dk = np.array([self.d0[min(k, n) - 1] * math.exp(x[min(k, len(x)) - 1]) / self.tab["q_req"] for k in range(1, len(self.passes) + 2)])
        N = plan["N"]; self.dk = dk / dk[min(len(dk) - 1, N - 1)]

    def _fb_slope(self, N):
        """DERIVED: dFB/d(nu_N) in m/s per ln unit: H_N dFB/drp at the table periapsis, desired periapsis radius fixed."""
        kk = min(N, self.tab["n"]); rp_t = float(self.tab["rp"][kk - 1]); H = float(self.tab["H"][kk - 1]); rp_des = rp_t - H_REF + self.h_des
        d = (first_burn_dv(RA_T, rp_t - 250.0, rp_des) - first_burn_dv(RA_T, rp_t + 250.0, rp_des)) / 500.0
        return H * d

    def _once(self, m, lam, pin, fb_N):
        """DP over passes 2..m+1 (m steps). Returns (dp result, x_full, arrival)."""
        A = self.A[:m].copy(); B = self.B[:m].copy(); A0 = self.A0; B0 = self.B0
        if pin == "hi":
            A = B.copy(); A0 = B0
        elif pin == "lo":
            B = A.copy(); B0 = A0
        lin = np.zeros(m); lin0 = 0.0
        if lam != 0.0:
            lin = lin - lam * self.dk[1:1 + m]; lin0 = -lam * float(self.dk[0])
        if fb_N is not None and self.fb_term:
            lin[min(fb_N - 2, m - 1)] += self._fb_slope(fb_N)
        r = free_start_dp(A, B, self.tw[:m], self.w_pre, lin if (lam != 0.0 or (fb_N is not None and self.fb_term)) else None, A0, B0, lin0, pen=self.pen)
        if r is None:
            return None, None, None
        x = np.concatenate(([X_REF + r["s"]], self.level[:m] + r["nu"]))
        return r, x, O.arrival(self.tab, x, X_REF, kappa=self.kappa)

    def _account(self, r, x, arr, m):
        N = arr["N"]; nu = r["nu"]; s = r["s"]
        du = np.diff(np.concatenate(([s], nu))); steps = np.arange(2, 2 + m)
        mk = steps <= N
        dv_c = float(np.sum((self.w[:m] * np.abs(du))[mk]))
        burns = int(((np.abs(du) * self.w[:m] >= P["dp_min_burn_mps"]) & mk).sum())
        nuN = float(nu[N - 2]) if 2 <= N <= m + 1 else float(s)
        kk = min(N, self.tab["n"]); H = float(self.tab["H"][kk - 1]); rp_t = float(self.tab["rp"][kk - 1])
        rp_f = rp_t - H * nuN; alt_f = H_REF - H * nuN; fb_signed = first_burn_dv(RA_T, rp_f, rp_t - H_REF + self.h_des); fb = abs(fb_signed)      # a lowering burn (final periapsis already above h_des) also costs its magnitude
        xs = float(s); dv_pre_charged = self.w_pre_full * abs(xs) if self.scenario == "S-fixed" else 0.0
        a_f = 0.5 * (RA_T + rp_f); T_final = 2 * math.pi * math.sqrt(a_f ** 3 / MU)
        upto = x[:N]; hi_path = np.concatenate(([self.hi_s], self.hi))[:N]
        return dict(N=N, days=arr["days"], x=x, nu=nu, s=xs, x0=float(x[0]), dv_campaign=dv_c, burns=burns, dv_pre=dv_pre_charged, dv_pre_equiv=self.w_pre_full * abs(xs), fb=fb, fb_signed=fb_signed, total=dv_c + dv_pre_charged + fb,
                    beyond_table=bool(N > len(self.passes) + 1), rp_final=rp_f, alt_final_m=alt_f, ln_q_final=float(x[min(N, len(x)) - 1]), T_final_s=T_final, margin_hi_x0=self.hi_s - float(x[0]), margin_hi_phys_x0=self.hi0 - float(x[0]),
                    path_margin_min=float(np.min(hi_path - upto)), obj=r["obj"], ra_end=arr["ra_end"])

    def solve(self, lam=0.0, pin=None):
        key = (round(lam, 15), pin)
        if key in self.memo:
            return self.memo[key]
        nfull = len(self.passes)
        if self.horizon == "full":
            N = self.N_guess or nfull; res = None
            for _ in range(5):
                r, x, arr = self._once(nfull, lam, pin, N)
                if r is None:
                    break
                res = (r, x, arr, nfull)
                if arr["N"] == N or not self.fb_term:
                    break
                N = arr["N"]
        else:
            Nh = self.N_guess or 133; tried = {}; res = None
            for _ in range(10):
                m = min(max(Nh - 1, 1), nfull)
                r, x, arr = self._once(m, lam, pin, m + 1)
                if r is None:
                    break
                tried[m + 1] = (r, x, arr, m)
                if arr["N"] == m + 1:
                    res = tried[m + 1]; break
                Nh = arr["N"] if math.isfinite(arr["days"]) else Nh + 5
                if Nh in tried:
                    break
            if res is None and tried:                         # no self-consistent horizon: take the cheapest plan whose horizon covers its own arrival pass
                cand = [v for k, v in tried.items() if v[2]["N"] <= k]
                res = min(cand, key=lambda v: v[0]["obj"]) if cand else None
        dropped = False
        if res is None and self.fb_term:                      # no self-consistent arrival pass with the FB term (it pulls the final level down, which delays the arrival past its own horizon): drop the term, flag it
            self.fb_term = False
            try:
                saveN = self.N_guess; res = None
                Nh = saveN or 133
                for _ in range(10):
                    m = min(max(Nh - 1, 1), nfull); r, x, arr = self._once(m, lam, pin, None)
                    if r is None:
                        break
                    res = (r, x, arr, m)
                    if arr["N"] == m + 1 or self.horizon == "full":
                        break
                    Nh = arr["N"] if math.isfinite(arr["days"]) else Nh + 5
            finally:
                self.fb_term = True
            dropped = True
        out = None
        if res is not None:
            r, x, arr, m = res; out = self._account(r, x, arr, m); out["lam"] = lam; out["pin"] = pin; out["fb_term_dropped"] = dropped
            if lam == 0.0 and pin is None:
                self.N_guess = out["N"]
        self.memo[key] = out
        return out

    # --- the arrival window by the Revision 18 price (bracket + bisection, memoized); a jump of the vertex plan over the window is closed by a convex mixture of the two bracketing plans
    def _mix(self, spec_a, spec_b, plan_a, plan_b, t_lo, t_hi, need_early):
        """Plans a (arrival outside the window, cheaper) and b (arrival beyond the window) at a common horizon; convex combination (1 - th) a + th b of the offsets (feasible: the constraints are boxes); smallest th that lands in the window.
        The mixture of two Lagrangian-optimal plans at the critical price is Lagrangian-optimal when the arrival is linear in the offsets (approximately true; the arrival and the dV of the mixture are evaluated exactly)."""
        m = min(max(len(plan_a["nu"]), len(plan_b["nu"])), len(self.passes))
        pad = lambda nu: np.concatenate((nu, np.clip(np.full(m - len(nu), nu[-1]), self.A[len(nu):m], self.B[len(nu):m]))) if len(nu) < m else nu[:m]      # steps beyond a plan's own horizon are filled by its last offset clipped into the tube
        ra = dict(nu=pad(plan_a["nu"]), s=plan_a["s"]); rb = dict(nu=pad(plan_b["nu"]), s=plan_b["s"])
        def at(th):
            nu = (1 - th) * ra["nu"] + th * rb["nu"]; s = (1 - th) * ra["s"] + th * rb["s"]; r = dict(nu=nu, s=s, obj=float("nan"))
            x = np.concatenate(([X_REF + s], self.level[:m] + nu)); arr = O.arrival(self.tab, x, X_REF, kappa=self.kappa)
            return r, x, arr
        ok = (lambda arr: arr["days"] <= t_hi) if need_early else (lambda arr: arr["days"] >= t_lo)
        lo_t, hi_t = 0.0, 1.0
        if not ok(at(1.0)[2]):
            return None
        for _ in range(50):
            mid = 0.5 * (lo_t + hi_t)
            if ok(at(mid)[2]):
                hi_t = mid
            else:
                lo_t = mid
        r, x, arr = at(hi_t); out = self._account(r, x, arr, m); out["lam"] = spec_b[0]; out["pin"] = spec_b[1]; out["fb_term_dropped"] = False; out["mix_theta"] = hi_t
        return out

    def window_plan(self, t_star_d, s_d):
        s_eff = max(float(s_d), QUANT_PERIODS * self.nominal["T_final_s"] / 86400.0)
        t_lo, t_hi = t_star_d - s_eff, t_star_d + s_eff
        p0 = self.nominal
        within = lambda p: p is not None and math.isfinite(p["days"]) and t_lo <= p["days"] <= t_hi and not p["beyond_table"]
        if within(p0):
            return dict(plan=p0, feasible=True, s_eff_d=s_eff, binding=False)
        need_early = p0["days"] > t_hi                         # lam > 0 rewards an earlier arrival
        sg = 1.0 if need_early else -1.0
        past = (lambda p: p["days"] <= t_hi and not p["beyond_table"]) if need_early else (lambda p: p["days"] >= t_lo and not p["beyond_table"])      # reached the window (possibly overshooting it)
        prev = 0.0; lam = 1.0e-3; hit = None
        while lam <= LAM_MAX:
            p = self.solve(sg * lam)
            if p is not None and math.isfinite(p["days"]) and past(p):
                hit = lam; break
            prev = lam; lam *= 2.0
        if hit is None:
            pin = self.solve(0.0, "hi" if need_early else "lo")
            if pin is None or not (math.isfinite(pin["days"]) and past(pin)):
                return dict(plan=pin, feasible=False, s_eff_d=s_eff, binding=True)
            pa = self.solve(sg * LAM_MAX); spec_a = (sg * LAM_MAX, None); spec_b = (0.0, "hi" if need_early else "lo"); pb = pin; bind = True
        else:
            lo_l, hi_l = prev, hit
            for _ in range(40):
                mid = 0.5 * (lo_l + hi_l); p = self.solve(sg * mid)
                if p is not None and math.isfinite(p["days"]) and past(p):
                    hi_l = mid
                else:
                    lo_l = mid
                if hi_l - lo_l < 1.0e-9 * max(hi_l, 1.0e-3):
                    break
            pa = self.solve(sg * lo_l) if lo_l > 0 else p0; spec_a = (sg * lo_l, None); pb = self.solve(sg * hi_l); spec_b = (sg * hi_l, None); bind = False
        mx = self._mix(spec_a, spec_b, pa, pb, t_lo, t_hi, need_early) if pa is not None else None        # the mixture lands on the window edge nearest to the cheaper plan a: never dearer than the plan b that overshoots into the window
        cand = [q for q in (pb, mx) if within(q)]
        if cand:
            return dict(plan=min(cand, key=lambda q: q["total"]), feasible=True, s_eff_d=s_eff, binding=bind)
        return dict(plan=mx if mx is not None else pb, feasible=False, s_eff_d=s_eff, binding=bind)


# ----------------------------------------------------------------------------------------------------- public interface

LABELS = dict(scenario="ASSUMED (scenario definition)", x0="DERIVED (DP optimum over the start level; reference x_ref MEASURED)", h_p0_m="DERIVED (H_REF SOURCED for Odyssey, ASSUMED to transfer; first order in x0 - x_ref)",
              nu_pre="DERIVED (x0 - x_ref)", dv_pre_mps="DERIVED (w_pre |nu_pre|; S-free: not charged by definition)", campaign_dv_mps="DERIVED (model world, noise free, exact weights, through the arrival pass)",
              fb_mps="DERIVED (exact vis-viva periapsis raise to h_des; h_des ASSUMED)", total_dv_mps="DERIVED (campaign + charged pre-burn + FB)", N="DERIVED (arrival pass)", t_arr_d="DERIVED (state-dependent arrival x kappa; kappa MEASURED per geometry, G0975 ASSUMED 1)")


def select_ic(geom_token, scenario, t_star_d=None, s_d=None, h_des_m=250e3, fb_term=False, kappa=None, *, horizon="arrival", delta_ln=0.0, pen=0.0, case=None):
    """Fuel-optimal initial conditions. scenario 'ref' | 'S-fixed' | 'S-free'. With t_star_d (and s_d, default 0) the arrival window |t_arr - t*| <= max(s, QUANT_PERIODS x T_final) is enforced by the Revision 18 price
    (bracket + bisection); reported infeasible (feasible False, tube binding) when no plan reaches it. Returns the ic_plan dict."""
    c = case or Case(geom_token, scenario, h_des_m, fb_term, kappa, horizon, delta_ln, pen)
    if t_star_d is None:
        w = dict(plan=c.nominal, feasible=True, s_eff_d=None, binding=False)
    else:
        w = c.window_plan(float(t_star_d), 0.0 if s_d is None else float(s_d))
    p = w["plan"]
    if p is None:
        return dict(scenario=scenario, geom=geom_token, feasible=False, tube_binding=True, t_star_d=t_star_d, s_d=s_d, fb_term=bool(fb_term), horizon=horizon, h_des_m=float(h_des_m), delta_ln=float(delta_ln), labels=LABELS)      # no plan reaches the window (not even the pinned one)
    d = dict(scenario=scenario, geom=geom_token, x0=p["x0"], h_p0_m=H_REF - float(c.tab["H"][0]) * p["s"], nu_pre=p["s"], dv_pre_mps=p["dv_pre"], dv_pre_equivalent_mps=p["dv_pre_equiv"],
             campaign_dv_mps=p["dv_campaign"], fb_mps=p["fb"], fb_signed_mps=p["fb_signed"], total_dv_mps=p["total"], N=int(p["N"]), t_arr_d=p["days"], burns=p["burns"], T_final_s=p["T_final_s"], rp_final_m=p["rp_final"],
             alt_final_km=p["alt_final_m"] / 1e3, margin_to_shrunk_upper_edge_x0=p["margin_hi_x0"], margin_to_physical_upper_edge_x0=p["margin_hi_phys_x0"], path_margin_min=p["path_margin_min"],
             lam=p["lam"], fb_term_dropped=bool(p.get("fb_term_dropped", False)), beyond_table=bool(p["beyond_table"]), mix_theta=p.get("mix_theta"), feasible=bool(w["feasible"]), tube_binding=bool(w["binding"]), horizon=horizon, h_des_m=float(h_des_m), fb_term=bool(fb_term), kappa=c.kappa, delta_ln=float(delta_ln), w_pre=c.w_pre_full,
             t_star_d=t_star_d, s_d=s_d, s_eff_d=w["s_eff_d"], dt_d=None if t_star_d is None else p["days"] - t_star_d, labels=LABELS)
    return d


# ----------------------------------------------------------------------------------------------------- reference burn dict (exactness check ii)

def burn_dict_from_plan(case, nu, s=0.0):
    du = np.diff(np.concatenate(([s], nu))); sched = {}
    for p_, d_, w_ in zip(case.passes, du, case.w):
        if abs(d_) * w_ >= P["dp_min_burn_mps"]:
            sched[p_] = -float(d_) * float(w_)
    return sched


# ----------------------------------------------------------------------------------------------------- selftest

def selftest():
    rng = np.random.default_rng(20)
    # 1. w_pre = inf equals dp_exact and dp_lin on random instances (exact equality of nu and obj)
    for n in (4, 12, 30):
        for _ in range(15):
            lev, lo, hi, w = PLB._random_instance(rng, n); A = lo - lev; B = hi - lev
            r0 = PLB.dp_exact(A, B, w); r1 = free_start_dp(A, B, w, INF)
            assert np.array_equal(r0["nu"], r1["nu"]) and r0["obj"] == r1["obj"] and r1["s"] == 0.0, (n,)
            lin = rng.normal(0, 0.05, n); nl = O.dp_lin(A, B, w, lin); r2 = free_start_dp(A, B, w, INF, lin)
            assert np.array_equal(nl, r2["nu"])
            pen = 0.07; r3 = PLB.dp_exact(A, B, w, pen=pen); r4 = free_start_dp(A, B, w, INF, pen=pen); assert np.array_equal(r3["nu"], r4["nu"]) and r3["obj"] == r4["obj"]
    # 2. priced start against brute force on random instances: min over s on a fine grid of w_pre |s| + dp_exact(A - s, B - s) (s grid contains the vertices)
    for w_pre in (0.0, 0.3, 2.5):
        for _ in range(10):
            n = 10; lev, lo, hi, w = PLB._random_instance(rng, n); A = lo - lev; B = hi - lev; A0, B0 = -0.6, 0.5
            r = free_start_dp(A, B, w, w_pre, None, A0, B0)
            grid = np.unique(np.concatenate((np.linspace(A0, B0, 2001), [0.0, A0, B0])))
            best = INF
            for s in grid:
                q = PLB.dp_exact(A - s, B - s, w)
                if q is not None:
                    best = min(best, w_pre * abs(s) + q["obj"])
            assert r["obj"] <= best + 1e-9 and best - r["obj"] < 2.0e-3 * max(w_pre, float(np.max(w)) , 1.0), (w_pre, r["obj"], best)
    print("selftest 1-2 ok: w_pre = inf equals dp_exact / dp_lin / dp_exact(pen) exactly; priced start equals brute force on random instances")
    # 3. AC1 on the three geometries
    rows = []
    for tok in GEOMS:
        for hz in ("full", "arrival"):
            # (ii) nu_pre forbidden, x0 = x_ref reproduces planner_lib.ol_schedule exactly (pen 0 and the config default pen)
            for pen in (0.0, P["dp_burn_penalty_mps"]):
                c = Case(tok, "ref", horizon="full", pen=pen); Pp = dict(P, dp_burn_penalty_mps=pen)
                sched, _ = PLB.ol_schedule(dict(), c.tab, Pp, X_REF, P["dp_tie_eps"])
                r = free_start_dp(c.A, c.B, c.tw, INF, pen=pen); mine = burn_dict_from_plan(c, r["nu"])
                assert mine == sched, (tok, pen, len(mine), len(sched))
            # (i) free-start optimum = brute-force grid over x0 (>= 201 points across the tube), pre-burn priced: S-free (w_pre 0) and S-fixed
            for sc in ("S-free", "S-fixed"):
                c = Case(tok, sc, horizon=hz); p = c.nominal
                m = len(p["nu"]); wp = c.w_pre
                grid = np.unique(np.concatenate((np.linspace(c.lo_s, c.hi_s, 401), [X_REF, c.lo_s, c.hi_s]))) - X_REF
                best = INF; bs = None
                for s in grid:
                    q = PLB.dp_exact(c.A[:m] - s, c.B[:m] - s, c.tw[:m])
                    v = (wp * abs(s) if math.isfinite(wp) else 0.0) + q["obj"]
                    if v < best:
                        best, bs = v, s
                # compare objective of the free DP on the same horizon (horizon fixed at the converged m)
                rr = free_start_dp(c.A[:m], c.B[:m], c.tw[:m], wp, None, c.A0, c.B0)
                assert rr["obj"] <= best + 1e-9 and best - rr["obj"] < 0.01, (tok, hz, sc, rr["obj"], best)
                rows.append((tok, hz, sc, len(grid), rr["obj"], best, rr["s"], bs))
                # (iii) S-fixed never worse than the reference start
                if sc == "S-fixed":
                    cr = Case(tok, "ref", horizon=hz); assert p["obj"] <= cr.nominal["obj"] + 1e-9 and p["total"] <= cr.nominal["total"] + 1.0e-6, (tok, hz, p["total"], cr.nominal["total"])   # 1e-6: the DP minimizes with tie-break weights (relative 1e-6)
    print("selftest 3 ok: (ii) ol_schedule burn dict reproduced exactly (pen 0 and 0.01) on ODY, G1, G0975; (i) free optimum vs brute-force grid (401+ points) within 0.01 m/s; (iii) S-fixed not worse than reference")
    for r in rows:
        print("  %-6s %-8s %-8s grid %d  DP %.6f  brute %.6f  (diff %.2e)  s* %.4f  s_grid* %.4f" % (r[0], r[1], r[2], r[3], r[4], r[5], r[5] - r[4], r[6], r[7]))
    # 4. interface and window sanity
    ic = select_ic("ODY", "S-free"); assert set(("scenario", "x0", "h_p0_m", "nu_pre", "dv_pre_mps", "campaign_dv_mps", "fb_mps", "fb_signed_mps", "total_dv_mps", "N", "t_arr_d", "labels")) <= set(ic)
    w = select_ic("ODY", "S-free", t_star_d=ic["t_arr_d"] - 5.0, s_d=0.0); assert w["feasible"] and abs(w["dt_d"]) <= w["s_eff_d"] + 1e-9, w["dt_d"]
    print("selftest 4 ok: ic_plan keys; window enforced (ODY S-free, t* = nominal - 5 d: dt %.3f d, total %.4f m/s)" % (w["dt_d"], w["total_dv_mps"]))
    print("planner_rev20_ic selftest ok")


# ----------------------------------------------------------------------------------------------------- runs (data for results_rev20_ic.md)

OUT = TD / "rev20_data" / "ic"
SCENARIOS = ("ref", "S-fixed", "S-free")
COLS = ("scenario", "geom", "x0", "h_p0_m", "nu_pre", "dv_pre_mps", "dv_pre_equivalent_mps", "campaign_dv_mps", "fb_mps", "fb_signed_mps", "total_dv_mps", "N", "t_arr_d", "burns", "alt_final_km", "margin_to_shrunk_upper_edge_x0",
        "margin_to_physical_upper_edge_x0", "path_margin_min", "lam", "fb_term_dropped", "beyond_table", "mix_theta", "feasible", "tube_binding", "horizon", "h_des_m", "fb_term", "kappa", "delta_ln", "w_pre", "t_star_d", "s_d", "s_eff_d", "dt_d")


def slim(ic):
    return {k: ic.get(k) for k in COLS if k in ic}


def _dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True); (OUT / name).write_text(json.dumps(obj, indent=1, default=float))


def run_headline():
    res = dict(labels=LABELS, plans=[], hdes=[], delta=[], pen=[], crosscheck={})
    for tok in GEOMS:
        for hz in ("arrival", "full"):
            for sc in SCENARIOS:
                for fbt in (False, True):
                    res["plans"].append(slim(select_ic(tok, sc, fb_term=fbt, horizon=hz)))
        for h in (120e3, 200e3, 250e3, 300e3, 400e3):
            for sc in SCENARIOS:
                for fbt in (False, True):
                    res["hdes"].append(slim(select_ic(tok, sc, h_des_m=h, fb_term=fbt)))
        for dl in (0.0, 0.1, 0.22):
            for sc in ("S-fixed", "S-free"):
                for fbt in (False, True):
                    res["delta"].append(slim(select_ic(tok, sc, fb_term=fbt, delta_ln=dl)))
        for sc in SCENARIOS:
            res["pen"].append(slim(select_ic(tok, sc, pen=P["dp_burn_penalty_mps"])))
    res["p20_2"] = []
    for tok in GEOMS:
        for hz in ("arrival", "full"):
            r = Case(tok, "ref", horizon=hz); f = Case(tok, "S-fixed", horizon=hz); a_, b_ = r.nominal, f.nominal
            res["p20_2"].append(dict(geom=tok, horizon=hz, w2=float(r.w[0]), w_pre=r.w_pre_full, ref_pass2_burn_mps=float(abs(a_["nu"][0]) * r.w[0]), ref_pass2_nu=float(a_["nu"][0]), sfixed_pre_dv_mps=b_["dv_pre"], sfixed_nu_pre=b_["s"],
                                     sfixed_pass2_burn_mps=float(abs(b_["nu"][0] - b_["s"]) * f.w[0]), ref_total=a_["total"], sfixed_total=b_["total"], saving_total=a_["total"] - b_["total"], ref_campaign=a_["dv_campaign"], sfixed_campaign_plus_pre=b_["dv_campaign"] + b_["dv_pre"]))
    try:                                                  # cross-check of the 'full' reference against the Revision 19 model_plan (same inputs; pen 0)
        import planner_rev19_offline as R19
        for tok in ("ODY", "G1"):
            O.KAPPA[tok] = load_kappa()[tok]; tab = O.load_tab(tok)
            for term in (False, True):
                m = R19.model_plan(tab, term); mine = select_ic(tok, "ref", fb_term=term, horizon="full")
                res["crosscheck"]["%s_fbterm_%s" % (tok, term)] = dict(rev19_total=m["total"], rev19_campaign=m["dv_campaign"], rev19_fb=m["first_burn"], rev19_days=m["days"], rev19_N=m["N"],
                                                                       rev20_total=mine["total_dv_mps"], rev20_campaign=mine["campaign_dv_mps"], rev20_fb=mine["fb_mps"], rev20_days=mine["t_arr_d"], rev20_N=mine["N"])
    except Exception as ex:                                # pragma: no cover
        res["crosscheck"]["error"] = repr(ex)
    _dump("headline.json", res)
    return res


def run_curves(tok, ngrid=20):
    cases = {(sc, fbt): Case(tok, sc, fb_term=fbt) for sc in SCENARIOS for fbt in (False, True)}
    base = {sc: cases[(sc, False)] for sc in SCENARIOS}
    dd = lambda sc, pin: (lambda q: q["days"] if (q is not None and math.isfinite(q["days"]) and not q["beyond_table"]) else INF)(base[sc].solve(0.0, pin))      # INF: the pinned plan never reaches the target apoapsis inside the table
    def slowest(sc):                                     # latest arrival reachable inside the table: pinned-lo plan, else the slowest finite plan of the negative-price scan
        d = dd(sc, "lo")
        if math.isfinite(d):
            return d
        best = base[sc].nominal["days"]
        for lam in [1e-3 * 2 ** i for i in range(0, 17)]:
            q = base[sc].solve(-lam)
            if q is not None and math.isfinite(q["days"]) and not q["beyond_table"]:
                best = max(best, q["days"])
        return best
    tmin = min(dd(sc, "hi") for sc in SCENARIOS); tmax = min(HORIZON_D, max(slowest(sc) for sc in SCENARIOS))
    Wmax_d = W_MAX_PERIODS * base["ref"].nominal["T_final_s"] / 86400.0
    grid = sorted(set([round(t, 4) for t in np.linspace(tmin, tmax, ngrid)] + [round(base[sc].nominal["days"], 4) for sc in SCENARIOS]))
    s_list = [Wmax_d if s == "WMAX" else s for s in S_GRID_D]
    rows = []
    for (sc, fbt), c in cases.items():
        for t in grid:
            best = None
            for sd in s_list:                          # ascending slack: a plan feasible for a smaller window is feasible for a larger one, so the dV never rises with s
                ic = slim(select_ic(tok, sc, t_star_d=t, s_d=sd, fb_term=fbt, case=c))
                if best is not None and best.get("feasible") and (not ic.get("feasible") or ic["total_dv_mps"] > best["total_dv_mps"] + 1e-12):
                    ic = dict(best, s_d=sd, carried_from_smaller_window=True)
                best = ic; rows.append(ic)
    tot = lambda sc, pin: (lambda q: q["total"] if q is not None else None)(base[sc].solve(0.0, pin))
    edges = {sc: dict(t_min_d=dd(sc, "hi"), t_max_d=dd(sc, "lo"), dv_at_hi=tot(sc, "hi"), dv_at_lo=tot(sc, "lo"), nominal_d=base[sc].nominal["days"]) for sc in SCENARIOS}
    _dump("curves_%s.json" % tok, dict(geom=tok, W_max_d=Wmax_d, T_final_s=base["ref"].nominal["T_final_s"], grid=grid, s_list=s_list, edges=edges, rows=rows))
    run_post(tok)


def apply_dominance(C):
    """Post-process a curves dict: a plan feasible for a smaller window is feasible for a larger one (dV never rises with s), and the scenarios are nested: a ref plan is an S-fixed and an S-free plan (nu_pre = 0), an S-fixed plan is an S-free plan
    with its pre-burn uncharged. The window search is heuristic (price + mixture), so rows are replaced by the nested better plan and flagged `carried_from`. Idempotent."""
    rows = C["rows"]; key = lambda r: (r["scenario"], bool(r["fb_term"]), round(r["t_star_d"], 6), round(r["s_d"], 9))
    idx = {key(r): r for r in rows if "t_star_d" in r and r.get("t_star_d") is not None}
    good = lambda r: r.get("feasible") and r.get("total_dv_mps") is not None
    for _ in range(2):
        for r in rows:
            if r.get("t_star_d") is None:
                continue
            cands = []
            for s2 in C["s_list"]:
                if s2 <= r["s_d"] + 1e-12:
                    q = idx.get((r["scenario"], bool(r["fb_term"]), round(r["t_star_d"], 6), round(s2, 9)))
                    if q is not None and good(q):
                        cands.append((q["total_dv_mps"], q, "smaller_window"))
                    if r["scenario"] in ("S-fixed", "S-free"):
                        q = idx.get(("ref", bool(r["fb_term"]), round(r["t_star_d"], 6), round(s2, 9)))
                        if q is not None and good(q):
                            cands.append((q["total_dv_mps"] - (0.0), q, "ref"))
                    if r["scenario"] == "S-free":
                        q = idx.get(("S-fixed", bool(r["fb_term"]), round(r["t_star_d"], 6), round(s2, 9)))
                        if q is not None and good(q):
                            cands.append((q["total_dv_mps"] - q["dv_pre_mps"], q, "S-fixed"))
            if not cands:
                continue
            best = min(cands, key=lambda c: c[0])
            if (not good(r)) or best[0] < r["total_dv_mps"] - 1e-9:
                q = best[1]; new = dict(q); new.update(scenario=r["scenario"], s_d=r["s_d"], carried_from=best[2])
                if r["scenario"] == "S-free" and q["scenario"] != "S-free":
                    new["total_dv_mps"] = best[0]; new["dv_pre_mps"] = 0.0
                if r["scenario"] == "S-fixed" and q["scenario"] == "ref":
                    new["dv_pre_equivalent_mps"] = 0.0
                r.clear(); r.update(new)
    C["dominance_applied"] = True
    return C


def run_post(tok):
    C = json.loads((OUT / ("curves_%s.json" % tok)).read_text())
    # rows lacking keys (no plan at all) get their coordinates back from the nested loop order of run_curves
    k = 0
    for sc in SCENARIOS:
        for fbt in (False, True):
            for t in C["grid"]:
                for sd in C["s_list"]:
                    r = C["rows"][k]; r.setdefault("scenario", sc); r.setdefault("fb_term", fbt); r.setdefault("t_star_d", t); r.setdefault("s_d", sd); k += 1
    apply_dominance(C)
    summ = dict(geom=tok, rows=[])
    for fbt in (False, True):
        for sd in C["s_list"]:
            for sc in ("S-fixed", "S-free"):
                sav = []; nfeas = 0; ntot = 0
                for t in C["grid"]:
                    q = [r for r in C["rows"] if r["scenario"] == sc and r["fb_term"] == fbt and abs(r["t_star_d"] - t) < 1e-9 and abs(r["s_d"] - sd) < 1e-9][0]
                    b = [r for r in C["rows"] if r["scenario"] == "ref" and r["fb_term"] == fbt and abs(r["t_star_d"] - t) < 1e-9 and abs(r["s_d"] - sd) < 1e-9][0]
                    ntot += 1
                    if q.get("feasible") and b.get("feasible"):
                        nfeas += 1; sav.append(b["total_dv_mps"] - q["total_dv_mps"])
                sav = np.array(sav) if sav else np.array([float("nan")])
                summ["rows"].append(dict(fb_term=fbt, s_d=sd, scenario=sc, n_both_feasible=nfeas, n_grid=ntot, saving_median=float(np.nanmedian(sav)), saving_max=float(np.nanmax(sav)), saving_min=float(np.nanmin(sav))))
    _dump("curves_%s.json" % tok, C); _dump("summary_%s.json" % tok, summ)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--run", choices=("headline", "curves", "post")); ap.add_argument("--geom"); a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.run == "headline":
        run_headline(); return
    if a.run == "post" and a.geom:
        run_post(a.geom); return
    if a.run == "curves" and a.geom:
        run_curves(a.geom); return
    ap.error("--selftest | --run headline | --run curves --geom ODY|G1|G0975")


if __name__ == "__main__":
    main()
