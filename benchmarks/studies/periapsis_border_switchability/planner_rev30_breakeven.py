"""Revision 30.1 (analysis only: Python, no Julia, no flights): break-even approach-side cost C* per target cell.

C* is the approach-side cost of choosing the insertion node at which insertion-plane matching ties with the in-orbit alternative:
  matching  = aerobraking dV with i and Omega_0 fixed (Revision 28C proxy) + residual plane-change cost + C*
  in orbit  = aerobraking dV of the free geometry (class best) + in-orbit plane change + phasing
so C* = (aero_free + plane_in + phasing) - (aero_fixed + residual plane cost). No value of C is assumed: C* is the number to compare with a mission-design estimate.
Inputs: rev28_data/rev28_plane.json (targets, fixed-node proxies), rev28_data/flown_residuals.csv (MEASURED flown node/plane errors), target_orbits.json (25 targets).
In-orbit figures: SOURCED only for G1 at 7,408 km (Revision 26a, first order, rev22 tools): plane change 17-25 m/s, phasing 15-78 m/s. Transfer to other cells is ASSUMED: the same equivalent plane error di_eq (2 v_a sin(di_eq/2) = 17 / 25 m/s at G1's v_a) corrected at that
cell's apoapsis speed; phasing 0 (lower bound) to 78 m/s (upper) because it is sourced only at G1. Residual plane cost: the flown di of the cell's token(s) (exact cell) or of the class at 7,408 km (transferred, node error in deg ASSUMED altitude independent) evaluated with the exact
apoapsis-node law at the cell's circle radius (r_p 3,504 km ASSUMED); median and p90 over the pooled flights.
Usage: python3 -I planner_rev30_breakeven.py run"""
import json, math, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_select_geometry as SG
import planner_rev22_rendezvous as R22

TD = SG.TD; OUT = TD / "rev30_data"
RP_FINAL = 3.504e6                          # ASSUMED (Revision 22)
IN_ORBIT_PLANE = (17.0, 25.0); IN_ORBIT_PHASING = (15.0, 78.0)     # SOURCED Revision 26a (G1, 7,408 km)
G1_RA = 7.408e6
CLASS_POOL = {"ODY": ["ODY", "G0975", "G0890"], "CRIT63": ["G0528"], "CRIT117": ["G5863"]}     # flown tokens of the class at 7,408 km
EXACT = {("ODY", 500.0): ["G1038"], ("CRIT63", 2000.0): ["G0486"], ("CRIT117", 1000.0): ["G1230"]}
FLOWN_ROWS = [("ODY", "ODY"), ("G0975", "ODY"), ("G0890", "ODY"), ("G1038", "ODY"), ("G0528", "CRIT63"), ("G0486", "CRIT63"), ("G5863", "CRIT117"), ("G1230", "CRIT117")]
TOKEN_CELL = {"ODY": ("ODY", 4011.81), "G0975": ("ODY", 4011.81), "G0890": ("ODY", 4011.81), "G1038": ("ODY", 500.0), "G0528": ("CRIT63", 4011.81), "G0486": ("CRIT63", 2000.0), "G5863": ("CRIT117", 4011.81), "G1230": ("CRIT117", 1000.0)}


def va_of(ra, rp=RP_FINAL):
    return math.sqrt(R22.MU * (2 / ra - 2 / (ra + rp)))


def di_eq():
    v = va_of(G1_RA); return tuple(2 * math.degrees(math.asin(d / (2 * v))) for d in IN_ORBIT_PLANE)


def in_orbit(ra, sourced):
    """(lo, hi) in-orbit plane change + phasing at a cell (m/s)."""
    if sourced: return IN_ORBIT_PLANE[0] + IN_ORBIT_PHASING[0], IN_ORBIT_PLANE[1] + IN_ORBIT_PHASING[1]
    lo, hi = (2 * va_of(ra) * math.sin(math.radians(d) / 2) for d in di_eq()); return lo, hi + IN_ORBIT_PHASING[1]


def costs_at(df, ra):
    return np.array([SG.exact_plane_increment(float(d), ra, RP_FINAL) for d in df.di_deg])


def run():
    P = json.load(open(TD / "rev28_data" / "rev28_plane.json")); F = pd.read_csv(TD / "rev28_data" / "flown_residuals.csv"); T = json.load(open(TD / "target_orbits.json"))["targets"]
    prox = {p["target_id"]: p for p in P["proxy"]}; rows = []
    for e in T:
        cls = e["class_tag"]; alt = e["alt_km"]; ra = e["radius_m"]; pr = prox.get(e["target_id"]); key = (cls, 500.0 if abs(alt - 500) < 1 else 1000.0 if abs(alt - 1000) < 1 else 2000.0 if abs(alt - 2000) < 1 else None)
        exact = EXACT.get(key); pool = exact or (CLASS_POOL.get(cls) if abs(alt - 4011.81) < 1 or True else None)
        r = dict(target_id=e["target_id"], class_tag=cls, alt_km=alt, radius_m=ra, flown_tokens=pool, flown_basis=("exact cell" if exact else ("class at 7,408 km, transferred" if pool else "none")))
        if abs(alt - 4011.81) < 1 and cls in CLASS_POOL: r["flown_basis"] = "exact cell (7,408 km)"
        if pr is None or not pool: r.update(status="undetermined (no flown node-error distribution)"); rows.append(r); continue
        g = F[F.token.isin(pool)]; c = costs_at(g, ra); sourced = (cls == "CRIT117" and abs(alt - 4011.81) < 1)
        lo, hi = in_orbit(ra, sourced); fixed = pr["median_best_interp"]; free = pr["free_best"]
        stats = {"median": float(np.median(c)), "p90": float(np.percentile(c, 90))}
        r.update(n_flights=int(len(g)), residual_node_deg_median=float(g.residual_deg.abs().median()), residual_cost_median=stats["median"], residual_cost_p90=stats["p90"], aero_free=free, aero_fixed=fixed, in_orbit_lo=lo, in_orbit_hi=hi, in_orbit_sourced=sourced)
        for nm in ("median", "p90"):
            base = free + 0.0 - (fixed + stats[nm]); r["Cstar_%s_lo" % nm] = lo + base; r["Cstar_%s_hi" % nm] = hi + base
        r["status"] = "wins at C=0" if r["Cstar_p90_lo"] > 0 else ("never wins (C* <= 0 even for the highest in-orbit figure)" if r["Cstar_median_hi"] <= 0 else "undetermined (sign depends on the in-orbit figure or on the p90)")
        rows.append(r)
    fl = []
    for tok, cls in FLOWN_ROWS:
        g = F[F.token == tok]; cl, alt = TOKEN_CELL[tok]; ra = SG.radius_of(alt) if alt != 4011.81 else 7.408e6; c = costs_at(g, ra); sourced = (tok == "G5863"); lo, hi = in_orbit(ra, sourced)
        m, p9 = float(np.median(c)), float(np.percentile(c, 90)); fl.append(dict(token=tok, class_tag=cl, alt_km=alt, n=int(len(g)), residual_cost_median=m, residual_cost_p90=p9, in_orbit_lo=lo, in_orbit_hi=hi, in_orbit_sourced=sourced,
                                                                                Cstar_median_lo=lo - m, Cstar_median_hi=hi - m, Cstar_p90_lo=lo - p9, Cstar_p90_hi=hi - p9))
    OUT.mkdir(parents=True, exist_ok=True); json.dump(dict(targets=rows, flown=fl, di_eq_deg=di_eq(), params=dict(IN_ORBIT_PLANE=IN_ORBIT_PLANE, IN_ORBIT_PHASING=IN_ORBIT_PHASING, RP_FINAL=RP_FINAL)), open(OUT / "breakeven.json", "w"), indent=1, default=float)
    return rows, fl


if __name__ == "__main__":
    rows, fl = run()
    print("di_eq", di_eq())
    for r in rows: print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k in ("target_id", "flown_basis", "n_flights", "residual_cost_median", "residual_cost_p90", "aero_free", "aero_fixed", "in_orbit_lo", "in_orbit_hi", "Cstar_median_lo", "Cstar_median_hi", "Cstar_p90_lo", "Cstar_p90_hi", "status")})
    for r in fl: print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
