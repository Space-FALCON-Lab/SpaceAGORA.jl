"""Self-tests S-M0..S-M6 of stagem_lib.py (pure Python, no truth data, no Julia).
Criteria: tasks/2026-10-15_periapsis-region-border-targeting/apoapsis_controller_plan.md (local note), section 6 and the WP hand-off table.
Usage: python3 -I stagem_selftest.py        (exit status 0 when every check passes; S-M6 is a Julia check, run by stagem_run.jl)
"""
import math, pathlib, sys
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import stagem_lib as L

MU = dict(mars=4.2828370371e13,        # SOURCED: Mars50c.csv header gm_m3s2 (the Mars truth gravity)
          earth=3.98600436233e14,      # SOURCED: Earth() struct
          venus=3.24858592e14)         # as used in the planning arithmetic (Stage M plan section 4)
# Stage I fixed-arm pass tables, median over the 48 cells per planet (planning arithmetic input): periapsis radius [m], fitted K, initial apoapsis [m]
RP0 = dict(mars=3501385.596638003, earth=6523114.899408773, venus=6187935.558017712)
K_FIT = dict(mars=41.2059727288554, earth=53.630783953686226, venus=23.995849948876433)
RA0_MARS = 28574250.47971137
Q_MEAN = 0.13930491646573473           # Stage I Mars fixed-arm median of the per-cell mean Q [W/cm2] (the plan's "0.139" rows)
RA_T = 3900e3                          # SOURCED: Odyssey Walkout end (plan D1)
H_MARS, BETA_MARS, F_ROT_MARS = 6308.278, 20.24, 1.05   # SOURCED: plan section 2.1 (Mars() H, Stage I beta, f_rot)
R_EQ_MARS, F_MARS = 3396.19e3, 0.005886                  # SOURCED: Mars() struct radii (Addendum B)
J2_MARS, R_J2_MARS = math.sqrt(5) * 0.000875919760273, 3394.2e3   # DERIVED: J2 = -sqrt(5) C20 of Mars50c.csv (fully normalized), its reference radius
E_Q = dict(mars=(-0.5824, 0.5342), earth=(-0.126, 0.124), venus=(-0.277, 0.256))   # one-pass ln-density error 5/95% (Stage I refit adapt_params.csv, P2_10; plan 3.1)
T_FLOOR_M = (-1140.0, 700.0)           # SOURCED: plan section 3.2, a 3.5-5 Mm, 20x20 dev pairs, 5/95% [m]

results = []
def check(tag, ok, msg):
    results.append(ok); print(("PASS " if ok else "FAIL ") + tag + " " + msg)

def cmap(planet, ra0, ra_t, Q_wcm2, rp=None):
    return L.campaign_map(ra0, ra_t, rp or RP0[planet], Q_wcm2 * 1e4, K_FIT[planet], MU[planet])

# ---- S-M0: campaign map reproduces plan section 4 (passes within 1, days within 0.15 d) -----------------
rows = [("full Mars Q.15", RA0_MARS, .15, 334, 86.5), ("full Mars Q.139", RA0_MARS, Q_MEAN, 360, 93.2),
        ("S-P Q.15", 10038e3, .15, 160, 21.0), ("S-P Q.139", 10038e3, Q_MEAN, 172, 22.6),
        ("S-E Q.15", 4904e3, .15, 28, 2.5), ("S-E Q.139", 4904e3, Q_MEAN, 30, 2.6)]
for name, ra0, Q, n_ref, d_ref in rows:
    n, d, _ = cmap("mars", ra0, RA_T, Q)
    check("S-M0", abs(n - n_ref) <= 1 and abs(d - d_ref) <= 0.15, "%s passes %d (ref %d) days %.2f (ref %.1f)" % (name, n, n_ref, d, d_ref))
e_p = (10038e3 - RP0["mars"]) / (10038e3 + RP0["mars"]); e_t = (RA_T - RP0["mars"]) / (RA_T + RP0["mars"])   # Pi-matched G scope: Mars S-P start e, target e_t
for planet, n_ref, d_ref in [("earth", 1002, 109), ("venus", 1830, 204)]:
    rp = RP0[planet]; n, d, _ = cmap(planet, rp * (1 + e_p) / (1 - e_p), rp * (1 + e_t) / (1 - e_t), .15)
    check("S-M0", abs(n - n_ref) <= 1 and abs(d - d_ref) <= 0.6, "G %s passes %d (ref %d) days %.1f (ref %d)" % (planet, n, n_ref, d, d_ref))

# ---- S-M1: theory gain -------------------------------------------------------------------------------
k = L.k_theory(H_MARS, BETA_MARS, F_ROT_MARS)
check("S-M1", abs(k - 41.3) <= 0.1, "k_theory(6308.278, 20.24, 1.05) = %.3f (ref 41.3 +- 0.1)" % k)

# ---- S-M2: Odyssey final raise -> r_p,final --------------------------------------------------------------
rpf = L.rp_final_from_raise(3.9e6, 3.505e6, 19.82, MU["mars"])
check("S-M2", abs(rpf - 3588.0e3) <= 100.0, "rp_final = %.3f km (ref 3588.0 +- 0.1); raise_dv round trip %.5f m/s" % (rpf / 1e3, L.raise_dv(3.9e6, 3.505e6, rpf, MU["mars"])))

# ---- S-M3: trim lever at the target orbit (analytic and by differencing trim_dv) ------------------------
lev = L.trim_lever(RA_T, RP0["mars"], MU["mars"])
fd = 2000.0 / (L.trim_dv(RA_T + 1e3, RA_T, RP0["mars"], MU["mars"]) + L.trim_dv(RA_T - 1e3, RA_T, RP0["mars"], MU["mars"]))
check("S-M3", abs(lev - 4592) <= 1 and abs(fd - lev) <= 5, "trim lever %.1f m per m/s (ref 4592 +- 1), differenced %.1f; %.4f m/s per km" % (lev, fd, 1e3 / lev))

# ---- S-M4: endgame_decision on synthetic sequences ------------------------------------------------------
du = L.drop_per_pass(0.5 * (RA_T + RP0["mars"]), e_t, RP0["mars"], 1e4, K_FIT["mars"], MU["mars"])   # drop per W/cm2 at the target orbit
qm, qe = 0.15, L.q_end(5e3, E_Q["mars"][1], E_Q["mars"][0], du, 0.02, 0.15); qn = qe / 2
check("S-M4", abs(qe - 0.0234) < 0.001 and abs(du - 303.8e3) < 0.2e3, "Mars drop_unit %.1f km, Q_end %.4f (plan: 303.8 km, about 0.023)" % (du / 1e3, qe))
for planet in ("earth", "venus"):
    rp = RP0[planet]; ra = rp * (1 + e_t) / (1 - e_t); a = 0.5 * (ra + rp)
    d1 = L.drop_per_pass(a, e_t, rp, 1e4, K_FIT[planet], MU[planet]); qq = L.q_end(5e3, E_Q[planet][1], E_Q[planet][0], d1, 0.02, 0.15)
    check("S-M4", qq == 0.15, "%s Q_end = %.4f (plan: no walk-down, Q_main)" % (planet, qq))
import random
ok_free = ok_noisy = True; worst = 0.0; qmax = 0.0; npl = 0
for D0 in [8e3, 31e3, 52e3, 100e3, 137.5e3, 400e3, 1.3e6]:
    for noise in (0.0, 1.0):
        rng = random.Random(int(D0)); D = D0; passes = 0; phases = []
        while True:
            act, q, ph = L.endgame_decision(D, du, qm, qe, qn)
            if act == "stop": break
            qmax = max(qmax, q); phases.append(ph); npl += 1
            D -= q * du * (math.exp(rng.gauss(0, 0.3)) if noise else 1.0); passes += 1
            if passes > 1000: ok_noisy = False; break
        if not noise:
            worst = max(worst, abs(D)); ok_free &= abs(D) <= du * qn / 2 + 1e-9
check("S-M4", ok_free and ok_noisy and qmax <= qm, "noise-free end |D| <= drop(Q_min)/2 = %.2f km (worst %.2f km); noisy runs all stop; max planned Q %.3f <= %.2f (%d plans)" % (du * qn / 2e3, worst / 1e3, qmax, qm, npl))
seq = []
D = 100e3
while True:
    act, q, ph = L.endgame_decision(D, du, qm, qe, qn)
    if act == "stop": break
    seq.append(ph); D -= q * du
check("S-M4", seq[-2:] == ["penultimate", "last"] and set(seq[:-2]) <= {"main"}, "phase order for D0 100 km: %s" % ",".join(seq[:2] + ["..."] + seq[-2:]))
check("S-M4", L.stop_nearest(10e3, 45.6e3) and not L.stop_nearest(30e3, 45.6e3), "stop_nearest (S0) rounds to the nearest pass")

# ---- S-M5: precision_mc against the plan 3.4 / A.7 walk-down rows (tolerance 0.01) ------------------------
mc = lambda d, **kw: L.precision_mc(d, *E_Q["mars"], *T_FLOOR_M, 200000, 20261022, **kw)
for name, d, w10, w5, kw in [("walk-down Q.05", 15.2e3, .94, .68, {}), ("walk-down Q.03", 9.1e3, .99, .89, {}),
                             ("M default Q_end", qe * du, .997, .95, {}), ("S0 at Q .15", 45.6e3, .42, .21, dict(s0=True)),
                             ("S0 at Q .112", 34.0e3, .55, .28, dict(s0=True))]:
    r = mc(d, **kw)
    check("S-M5", abs(r["within10"] - w10) <= .01 and abs(r["within5"] - w5) <= .01, "%s d %.1f km: within10 %.3f (ref %.3f) within5 %.3f (ref %.3f) mean|err| %.2f km" % (name, d / 1e3, r["within10"], w10, r["within5"], w5, r["mean_abs"] / 1e3))
r = mc(qe * du); check("S-M5", abs(r["mean_abs"] / 1e3 - 1.98) <= 0.1, "M default mean |err| %.2f km (A.7: 1.98); trim %.2f m/s (A.7: 0.43)" % (r["mean_abs"] / 1e3, r["mean_abs"] / lev))

# ---- k_update: EWMA convergence and guards ---------------------------------------------------------------
a, e = 0.5 * (10038e3 + RP0["mars"]), e_p; q = 0.15e4
obs = L.drop_per_pass(a, e, RP0["mars"], q, 41.2, MU["mars"]); kk = 30.0
for _ in range(60): kk = L.k_update(kk, obs, a, e, RP0["mars"], q, MU["mars"])
check("k_update", abs(kk - 41.2) < 1e-3 and L.k_update(30.0, -1.0, a, e, RP0["mars"], q, MU["mars"]) == 30.0 and L.k_update(30.0, float("nan"), a, e, RP0["mars"], q, MU["mars"]) == 30.0,
      "EWMA converges to 41.2 (%.4f); bad observations leave K unchanged" % kk)

# ---- J2 apsidal rotation over the S-P scope (A.2 item 3: about 48 deg in 161 passes at i 88.6-98.6) -----
for inc_deg in (88.6, 93.5, 98.6):
    ra, dw, n = 6770e3 * 2 - RP0["mars"], 0.0, 0     # ASSUMED start a = 6770 km (plan A.2: "a 6,770 to 3,700 km")
    while ra > RA_T:
        a = 0.5 * (ra + RP0["mars"]); e = (ra - RP0["mars"]) / (ra + RP0["mars"])
        dw += L.apsidal_rotation_per_pass(a, e, math.radians(inc_deg), J2_MARS, R_J2_MARS, MU["mars"])
        ra -= L.drop_per_pass(a, e, RP0["mars"], .15e4, K_FIT["mars"], MU["mars"]); n += 1
    print("INFO J2 apsidal rotation i %.1f: %.1f deg in %d passes" % (inc_deg, math.degrees(abs(dw)), n))
    check("J2", 35 <= math.degrees(abs(dw)) <= 60, "|d omega| %.1f deg (plan: about 48)" % math.degrees(abs(dw)))
h0 = L.periapsis_geodetic_altitude(R_EQ_MARS + 92e3, math.radians(90), math.radians(90), R_EQ_MARS, F_MARS)
check("J2", abs(h0 - 112.0e3) < 50, "geodetic altitude at lat 90 deg N, rp = R_eq + 92 km: %.2f km (Addendum B exact conversion: 112.00)" % (h0 / 1e3))
print("SKIP S-M6 Julia check (B1 edges bitwise equal stagef_edges(:fixed); stagem_run.jl --phase=selftest)")
print("%d/%d checks passed" % (sum(results), len(results))); sys.exit(0 if all(results) else 1)
