"""Revision 20 WP3: end-of-campaign window toolbox (pure functions, no I/O, no Julia, usable by any arm). Python only.

Problem. A rendezvous needs the closing apoapsis at clock time t* within +/- W of it (W = W_MAX_PERIODS final-orbit periods, ASSUMED 3). A campaign closes at t_c (flown or planned). The toolbox moves the arrival by
  (a) PAUSE: at an apoapsis DURING the campaign, a tangential burn raises the periapsis by DH_POPUP_M (out of the drag region), the vehicle coasts k whole orbits of the pop-up orbit, a second burn
      (same magnitude: the apoapsis radius is unchanged) lowers the periapsis back. Cost = 2 x exact vis-viva burn. Delay = k x T(pop-up orbit): pause orbits are EXTRA clock time (arrival t_c + k T_pop).
      To first order the number of drag passes is unchanged (the pop-up density is exp(-dh/H) times the flown one; DH_POPUP_M keeps the per-orbit apoapsis drop <= DROP_TOL_M), so the campaign after the
      pause is the flown campaign shifted in time (second order, not modeled: the density depends on the epoch). A delay is positive only: a pause fixes EARLY arrival only.
  (b) PHASING after the closing apoapsis: k_ph orbits of a phasing orbit (apoapsis radius unchanged, periapsis altitude from min(h_des, alt_f + DH_POPUP_M) to H_PH_MAX_M) instead of the standing orbit (the h_des orbit of the
      first burn). Convention (ASSUMED): the arrival moves by D = k_ph (T_ph - T_ref), T_ref = period of the h_des orbit; D < 0 (shorter phasing orbit, lower periapsis) advances the arrival, D > 0 delays it.
      Cost = exact two burns (rp_f -> rp_ph -> rp_des) minus the first burn (rp_f -> rp_des) that every case pays (>= 0 by the triangle inequality). k_ph <= K_PH_MAX (ASSUMED).
      The 8.9 s per m/s lever of Revision 18 is the small-shift limit of this exact relation (checked in the selftest).
  (c) cheapest(): minimum toolbox cost over pause location (any apoapsis before the closing one), whole orbits k, and phasing, that brings |t_arr - t*| <= W; None = window failure.
Apse angle: each extra orbit rotates the apse line by the secular J2 rate of ITS orbit; cheapest() returns rot (rad added by pause orbits and phasing excess) and apse_clock (rad, mismatch against the
final-orbit J2 rate over the same elapsed time). Reported, not controlled.
Usage: python3 -I planner_rev20_window.py --selftest"""
import argparse, math, pathlib, sys
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import geometry_rank as G

MU = A.MU                                  # SOURCED: planner_analysis (Mars GM, as in the flights)
PL = G.K.planet("Mars", 20); J2 = PL["J"][2]; RJ = PL["R"]     # SOURCED: Kaula Mars table (as planner_rev18_flights)

# ---- named parameters; every value is ASSUMED unless marked otherwise ----
W_MAX_PERIODS = 3.0       # ASSUMED (Revision 18): window half width in final-orbit periods
H_DES_M = 250e3           # ASSUMED (Revision 18): periapsis altitude after the first burn; sensitivity 120-400 km
DH_POPUP_M = 70e3         # ASSUMED, derived default: H ln(drop_p99 / DROP_TOL_M) = 6.41 km x ln(331 km / 10 m) = 66.6 km, rounded up to 70 km. drop_p99 = 99th percentile of the flown per-pass
                          #   apoapsis drop (75177 passes of the 366 flights, MEASURED), H = median flown local scale height (MEASURED). planner_rev20_baselines.py --dh-table reproduces both inputs.
DROP_TOL_M = 10.0         # ASSUMED: per-orbit apoapsis drop at the pop-up altitude regarded as negligible (100 pause orbits = 1 km = 4% of the 25 km final-apoapsis tolerance)
                          # Phasing-orbit periapsis altitude floor: min(h_des, alt_f + DH_POPUP_M) (ASSUMED rule: the same drag-safe margin as the pause; alt_f = flown final periapsis altitude)
H_PH_MAX_M = 1000e3       # ASSUMED: highest periapsis altitude of a phasing orbit
K_PH_MAX = 50             # ASSUMED: most phasing orbits (about 7 days at 3.4 h)
N_GRID = 1500             # numerical: grid points per sign of the shift for the phasing cost table
K_SEARCH = 2              # numerical: whole-orbit counts searched around the nearest one for each pause location


def va(ra, rp):
    """Apoapsis speed of the orbit (ra, rp): sqrt(mu (2/ra - 1/a)), a = (ra + rp)/2 (exact vis-viva)."""
    return np.sqrt(MU * (2.0 / ra - 2.0 / (ra + rp)))


def burn(ra, rp1, rp2):
    """Exact tangential apoapsis burn magnitude (m/s) moving the periapsis radius from rp1 to rp2 at fixed ra."""
    return np.abs(va(ra, rp2) - va(ra, rp1))


def period(ra, rp):
    return 2.0 * np.pi * np.sqrt(((ra + rp) / 2.0) ** 3 / MU)


def apsidal_per_orbit(ra, rp, inc):
    """Secular J2 apsidal rotation per orbit (rad): 3 pi J2 (R/p)^2 (2 - 5/2 sin^2 i), p = a (1 - e^2) (same form as planner_rev18_flights.PLB_apsidal)."""
    a = (ra + rp) / 2.0; e = (ra - rp) / (ra + rp); p = a * (1 - e * e)
    return 3 * np.pi * J2 * (RJ / p) ** 2 * (2 - 2.5 * np.sin(inc) ** 2)


def lever_T(ra, rp):
    """Small-shift limit, s of period per m/s of apoapsis burn: 3 T sqrt((1-e)/(1+e)) / V_c (Revision 18)."""
    a = (ra + rp) / 2.0; e = (ra - rp) / (ra + rp)
    return 3 * period(ra, rp) * np.sqrt((1 - e) / (1 + e)) / np.sqrt(MU / a)


def pause_option(ra, rp, inc, dh=None):
    """Pause at apoapses of orbits (ra, rp) (scalars or arrays): dict(cost = 2 burns (m/s), T_pop (s), rot (rad per pop-up orbit))."""
    dh = DH_POPUP_M if dh is None else dh; rp2 = rp + dh
    return dict(cost=2.0 * burn(ra, rp, rp2), T_pop=period(ra, rp2), rot=apsidal_per_orbit(ra, rp2, inc))


def phasing_table(ra_f, rp_f, rp_des, rp_lo, rp_hi, k_max=None, n_grid=None):
    """Cost of shifting the arrival by D seconds with k_ph orbits of a phasing orbit, minimized over k_ph = 1..k_max, on a grid of D (ascending, 0 included). cost = extra m/s over the first burn (inf = unreachable)."""
    k_max = K_PH_MAX if k_max is None else k_max; n_grid = N_GRID if n_grid is None else n_grid
    T_ref = period(ra_f, rp_des); fb = burn(ra_f, rp_f, rp_des)
    D = np.unique(np.concatenate([np.linspace((period(ra_f, rp_lo) - T_ref) * k_max, 0.0, n_grid), np.linspace(0.0, (period(ra_f, rp_hi) - T_ref) * k_max, n_grid)]))
    cost = np.full(D.shape, np.inf); kb = np.zeros(D.shape, int); rpb = np.full(D.shape, np.nan)
    for k in range(1, k_max + 1):
        T_ph = T_ref + D / k; a_ph = (np.maximum(T_ph, 1.0) / (2 * np.pi)) ** (2.0 / 3.0) * MU ** (1.0 / 3.0); rp_ph = 2 * a_ph - ra_f
        ok = (T_ph > 0) & (rp_ph >= rp_lo - 1e-6) & (rp_ph <= rp_hi + 1e-6); rs = np.where(ok, rp_ph, rp_des)
        c = np.where(ok, burn(ra_f, rp_f, rs) + burn(ra_f, rs, rp_des) - fb, np.inf)
        better = c < cost; cost = np.where(better, c, cost); kb = np.where(better, k, kb); rpb = np.where(better, rp_ph, rpb)
    cost[np.abs(D) < 1e-9] = 0.0; kb[np.abs(D) < 1e-9] = 0
    return dict(D=D, cost=cost, k=kb, rp_ph=rpb, T_ref=T_ref, fb=fb, ra_f=ra_f, rp_des=rp_des)


def phasing_cost(tab, D):
    """(extra cost m/s, k_ph, rp_ph) for a shift D (s); inf outside the reachable range. Takes the grid point at or beyond |D| (cost is monotone in |D|), so it never understates."""
    Dg = tab["D"]
    if D < Dg[0] or D > Dg[-1]: return math.inf, 0, float("nan")
    i = min(max(int(np.searchsorted(Dg, D)), 1), len(Dg) - 1); j = i if abs(Dg[i]) >= abs(Dg[i - 1]) else i - 1
    return float(tab["cost"][j]), int(tab["k"][j]), float(tab["rp_ph"][j])


def prepare(ra, rp, inc, ra_f, rp_f, alt_f, inc_f, h_des_m=None):
    """Bundle a campaign for cheapest(): arrays of the apoapses before the closing one (ra, rp, inc of the orbit after each flown periapsis pass), the final state, and its phasing table at h_des_m."""
    h = H_DES_M if h_des_m is None else h_des_m; r_surf = rp_f - alt_f; rp_des = r_surf + h
    tab = phasing_table(ra_f, rp_f, rp_des, r_surf + min(h, alt_f + DH_POPUP_M), r_surf + H_PH_MAX_M)
    return dict(ra=np.asarray(ra, float), rp=np.asarray(rp, float), inc=np.asarray(inc, float), ra_f=ra_f, rp_f=rp_f, inc_f=inc_f, rp_des=rp_des, ptab=tab,
                T_f=period(ra_f, rp_f), rot_f=float(apsidal_per_orbit(ra_f, rp_f, inc_f)))


def cheapest(t_c, t_star, W, st, dh=None):
    """Cheapest toolbox plan bringing |t_arr - t*| <= W (seconds). st = prepare(...). Returns dict(cost m/s, kind none|phasing|pause|pause+phasing, p, k, k_ph, D, t_arr, rot, apse_clock) or None = window failure.
    Pause tried at every stored apoapsis in order of cost with k whole orbits around the nearest count (k >= 1); phasing absorbs what is left beyond the window edge."""
    r0 = t_star - t_c; best = [None]; ptab = st["ptab"]
    def consider(cost, kind, p, k, T_pop, rot_pop, rr):
        D = 0.0 if abs(rr) <= W else rr - math.copysign(W, rr)
        c_ph, kph, rpph = (0.0, 0, float("nan")) if D == 0.0 else phasing_cost(ptab, D)
        tot = cost + c_ph
        if math.isfinite(tot) and (best[0] is None or tot < best[0]["cost"]):
            rot_ph = kph * (float(apsidal_per_orbit(st["ra_f"], rpph, st["inc_f"])) - float(apsidal_per_orbit(st["ra_f"], st["rp_des"], st["inc_f"]))) if kph else 0.0
            rot = k * rot_pop + rot_ph; clock = k * (rot_pop - st["rot_f"] * T_pop / st["T_f"]) + rot_ph       # phasing orbits are inside the shift D: mismatch vs the standing orbit
            best[0] = dict(cost=tot, kind=kind if kph == 0 or kind == "phasing" else kind + "+phasing", p=p, k=k, k_ph=kph, D=D, t_arr=t_c + k * T_pop + D, rot=rot, apse_clock=clock)
    if abs(r0) <= W: return dict(cost=0.0, kind="none", p=-1, k=0, k_ph=0, D=0.0, t_arr=t_c, rot=0.0, apse_clock=0.0)
    consider(0.0, "phasing", -1, 0, 0.0, 0.0, r0)
    if r0 > 0 and len(st["ra"]):
        po = st[("po", dh)] if ("po", dh) in st else st.setdefault(("po", dh), pause_option(st["ra"], st["rp"], st["inc"], dh)); k0 = np.floor(r0 / po["T_pop"]).astype(int)
        for p in np.argsort(po["cost"]):
            if best[0] is not None and po["cost"][p] >= best[0]["cost"]: break          # sorted by pause cost: nothing cheaper remains
            for k in range(max(1, k0[p] - 1), k0[p] + 2 + K_SEARCH):
                consider(float(po["cost"][p]), "pause", int(p), int(k), float(po["T_pop"][p]), float(po["rot"][p]), r0 - k * float(po["T_pop"][p]))
    return best[0]


def selftest():
    ra = 7.408e6; rp = 3.5e6; inc = math.radians(93.5); rp_des = rp + 145e3
    assert abs(burn(ra, rp, rp + 5e4) - burn(ra, rp + 5e4, rp)) < 1e-12
    # small-shift lever: exact period change of an apoapsis burn equals lever_T x dv (1%), and the final-orbit lever is 8-10 s per m/s (Revision 18: 8.896 median)
    a = (ra + rp) / 2; e = (ra - rp) / (ra + rp); dv = 0.05; dr = dv * 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)
    assert abs(burn(ra, rp, rp + dr) / dv - 1) < 0.02 and abs((period(ra, rp + dr) - period(ra, rp)) / dv / lever_T(ra, rp) - 1) < 0.02
    assert 8.0 < lever_T(ra, rp) < 10.0, lever_T(ra, rp)
    p1 = pause_option(ra, rp, inc); assert abs(p1["cost"] - 2 * burn(ra, rp, rp + DH_POPUP_M)) < 1e-12 and p1["T_pop"] > period(ra, rp)
    assert pause_option(2.8e7, rp, inc)["cost"] < p1["cost"]                                    # a bigger (earlier) orbit is cheaper to pop up
    tab = phasing_table(ra, rp, rp_des, rp_des - 75e3, rp_des + 750e3)
    assert tab["cost"][np.argmin(np.abs(tab["D"]))] == 0.0 and np.all(tab["cost"] >= -1e-9)
    for D in (-3000.0, -400.0, 700.0, 5000.0):                                                  # returned (k, rp_ph) reproduce the shift on the grid point at or beyond D
        c, k, rpp = phasing_cost(tab, D); j = np.searchsorted(tab["D"], D); Dg = tab["D"][max(j - 1, 0):j + 1]
        assert math.isfinite(c) and k >= 1 and min(abs(k * (period(ra, rpp) - tab["T_ref"]) - d) for d in Dg) < 1e-6 * abs(D) + 1e-3
    assert phasing_cost(tab, 1e9)[0] == math.inf
    t5 = phasing_table(ra, rp, rp_des, rp_des - 75e3, rp_des + 750e3, k_max=5)
    assert phasing_cost(tab, -3000.0)[0] < phasing_cost(t5, -3000.0)[0]                          # more orbits are cheaper
    c, k, _ = phasing_cost(tab, 500.0); assert abs(c / (2 * 500.0 / (k * lever_T(ra, rp_des))) - 1) < 0.15, (c, k)    # small shift: 2 |D| / (k lever)
    st = prepare([2.8e7, 1.5e7, 9e6], [3.5e6] * 3, [inc] * 3, ra, rp, 105e3, inc); W = W_MAX_PERIODS * period(ra, rp)
    assert cheapest(1000.0, 1000.0 + 0.5 * W, W, st)["kind"] == "none"
    early = cheapest(0.0, 3 * 86400.0, W, st); assert early["kind"].startswith("pause") and abs(early["t_arr"] - 3 * 86400.0) <= W + 1.0 and early["cost"] > 0 and early["p"] == 0
    assert cheapest(3 * 86400.0, 0.0, W, st) is None                                             # a 3 d lag exceeds the phasing authority
    ph = cheapest(0.0, -(W + 1500.0), W, st); assert ph["kind"] == "phasing" and abs(ph["t_arr"] - (-(W + 1500.0))) <= W + 1.0 and abs(ph["D"] + 1500.0) < 1e-6
    print("selftest OK: final-orbit lever %.3f s per m/s; pause cost at (7408 km, 3500 km) %.3f m/s; early plan %s %.3f m/s k=%d; phasing 1500 s advance %.3f m/s k_ph=%d"
          % (lever_T(ra, rp_des), p1["cost"], early["kind"], early["cost"], early["k"], ph["cost"], ph["k_ph"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    selftest() if a.selftest else ap.print_help()
