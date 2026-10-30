"""Half-lever check (open issue I1): is the realized periapsis response to an apoapsis burn really 0.50 of the vis-viva lever?

Usage: python3 -I halflever_check.py [--selftest] [--dev DIR]   (DIR = tasks/.../stageo_data/dev)
Selftest (no data): the closed-form lever vs a Cartesian two-body burn, exact vs linearized.
Data: per burn (decisions.csv row with delta_v != 0, target pass k), compares the event-radius change rp_k - rp_(k-1) with
lever * dv, against a no-burn counterfactual baseline taken from non-burn pairs of the same run.
"""
import sys, math, pathlib, csv, glob, statistics as st
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagek_lib as K
MU = K.planet("Mars")["mu"]

def lever_formula(a, e, mu=MU):  # m per (m/s), the spec form 4 a sqrt((1-e)/(1+e)) / sqrt(mu/a)
    return 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(mu / a)

def drp_exact(ra, rp, dv, mu=MU):
    """Exact r_p change of a tangential apoapsis burn dv (r_a unchanged), from vis-viva."""
    a = 0.5 * (ra + rp); va = math.sqrt(mu * (2 / ra - 1 / a)); v2 = va + dv
    a2 = 1 / (2 / ra - v2 * v2 / mu)
    return (2 * a2 - ra) - rp

def drp_cartesian(ra, rp, dv, mu=MU):
    """Independent: state at apoapsis, add dv along velocity, recompute r_p from the energy and angular momentum vectors."""
    a = 0.5 * (ra + rp); va = math.sqrt(mu * (2 / ra - 1 / a))
    r = np.array([ra, 0, 0.]); v = np.array([0, va + dv, 0.])
    h = np.cross(r, v); ev = np.cross(v, h) / mu - r / np.linalg.norm(r); e2 = np.linalg.norm(ev)
    a2 = -mu / (2 * (v @ v / 2 - mu / ra))
    return a2 * (1 - e2) - rp

def selftest():
    ok = True
    for e in (0.63, 0.78, 0.85, 0.9):
        a = 1.5e7; ra, rp = a * (1 + e), a * (1 - e)
        for dv in (0.1, 1.0, 5.0):
            ex, ca = drp_exact(ra, rp, dv), drp_cartesian(ra, rp, dv); lin = lever_formula(a, e) * dv
            print(f"e={e} dv={dv}: exact {ex:.3f} cart {ca:.3f} linear {lin:.3f} lin/exact {lin/ex:.5f}")
            ok &= abs(ex - ca) < 1e-3 and abs(lin / ex - 1) < 0.01 * dv + 1e-4
    # the "half" candidates: the spec lever vs 2a(1-e)Pi (framework 5.1) and vs 2*da only
    a, e = 1.5e7, 0.78; L = lever_formula(a, e); print("5.1 form / spec:", 2 * a * (1 - e) / math.sqrt(MU / a) / L)
    # gravity-only vacuum burn (J2..J20), Odyssey-like orbit: the event radius must move by the two-body lever to a few percent
    dr, dosc, L = gravity_only_burn_response(1.5576e7, 0.7764, 1.632, 0.49, 1.91, 0.25)
    print(f"gravity-only: event dr {dr:.1f} osc d(a(1-e)) {dosc:.1f} lever*dv {L:.1f} ratio {dr/L:.3f}")
    ok &= abs(dr / L - 1) < 0.05
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok


def gravity_only_burn_response(a, e, inc, Om, w, dv, nmax=20):
    """Vacuum zonal (J2..J20) orbit: periapsis state -> apoapsis (impulse dv along velocity) -> next periapsis. Returns
    (r_event change vs the no-burn run, osculating a(1-e) change, the two-body lever*dv). Tests whether gravity or the event definition halves the response."""
    from scipy.integrate import solve_ivp
    P = K.planet("Mars", nmax); rhs = K.make_rhs(P, nmax)
    r0, v0 = K.periapsis_state(MU, a, e, inc, Om, w); T = 2 * math.pi * math.sqrt(a ** 3 / MU)
    def run(dvv):
        u = np.concatenate([r0, v0])
        s = solve_ivp(rhs, [0, 1e-3 * T], u, method="DOP853", rtol=1e-11, atol=1e-8).y[:, -1]
        ev = lambda t, x: x[0] * x[3] + x[1] * x[4] + x[2] * x[5]
        ev.direction = -1
        s1 = solve_ivp(rhs, [1e-3 * T, 1.5 * T], s, method="DOP853", rtol=1e-11, atol=1e-8, events=ev)
        t1, x = s1.t_events[0][0], s1.y_events[0][0].copy()
        x[3:] += dvv * x[3:] / np.linalg.norm(x[3:])
        ev.direction = 1
        s2 = solve_ivp(rhs, [t1, t1 + 1.5 * T], x, method="DOP853", rtol=1e-11, atol=1e-8, events=ev)
        y = s2.y_events[0][0]; el = K.elements(MU, np.array(y[:3]), np.array(y[3:]))
        return np.linalg.norm(y[:3]), el[0] * (1 - el[1])
    b, c = run(dv), run(0.0)
    return b[0] - c[0], b[1] - c[1], lever_formula(a, e) * dv

def load(run):
    P = list(csv.DictReader(open(run / "passes.csv"))); D = list(csv.DictReader(open(run / "decisions.csv")))
    return P, D

def analyze(dev, field="rp_m", win=6):
    rows = []
    for run in sorted(glob.glob(str(dev / "runs" / "*"))):
        run = pathlib.Path(run)
        try: P, D = load(run)
        except Exception: continue
        arm = run.name.split("_")[0]
        rp = np.array([float(p[field]) for p in P]); ra = np.array([float(p["ra_m"]) for p in P])
        alt = np.array([float(p["alt_m"]) for p in P])
        dvk = {}
        for d in D:
            if float(d["delta_v_mps"]) != 0: dvk[int(d["target_pass"])] = float(d["delta_v_mps"])
        for k, dv in dvk.items():
            i = k - 1  # index of pass k
            if i < 1 or i >= len(P): continue
            rpb, rpa = rp[i - 1], rp[i]
            a_pre = 0.5 * (ra[i - 1] + rpb)  # pre-burn orbit: apoapsis radius of the pass before, periapsis radius of that pass
            e_pre = (ra[i - 1] - rpb) / (ra[i - 1] + rpb)
            L = lever_formula(a_pre, e_pre) * dv
            Ex = drp_exact(ra[i - 1], rpb, dv)
            # counterfactual baseline: median change of non-burn pairs in the window, and the same for altitude
            nb = [rp[j] - rp[j - 1] for j in range(max(1, i - win), min(len(P), i + win + 1)) if j != i and (j + 1) not in dvk]
            nba = [alt[j] - alt[j - 1] for j in range(max(1, i - win), min(len(P), i + win + 1)) if j != i and (j + 1) not in dvk]
            # 3-pass level shift (burn between i-1 and i), neighbors without burns only
            pre = [rp[j] for j in range(max(0, i - 3), i) if all((jj + 1) not in dvk for jj in range(max(0, i - 3) + 1, i))]
            post = [rp[j] for j in range(i, min(len(P), i + 3)) if all((jj + 1) not in dvk for jj in range(i + 1, min(len(P), i + 3)))]
            rows.append(dict(run=run.name, arm=arm, k=k, dv=dv, L=L, Ex=Ex, naive=rpa - rpb,
                             base=float(np.median(nb)) if nb else float("nan"),
                             dalt=alt[i] - alt[i - 1], basealt=float(np.median(nba)) if nba else float("nan"),
                             shift3=(np.mean(post) - np.mean(pre)) if len(pre) == 3 and len(post) == 3 else float("nan"),
                             # also what apoapsis radius did around the burn (tangential at apoapsis: r_a unchanged)
                             dra=ra[i] - ra[i - 1]))
    return rows

def summarize(rows, label):
    R = lambda f: np.array([r[f] for r in rows])
    L, Ex = R("L"), R("Ex"); out = {}
    def rep(name, resp):
        m = ~np.isnan(resp) & (np.abs(L) > 0)
        x, y = L[m], resp[m]; r = y / x
        beta = (x @ y) / (x @ x); c = np.corrcoef(x, y)[0, 1]
        print(f"{label:6s} {name:26s} n={m.sum():4d} median ratio {np.median(r):.3f} slope(through 0) {beta:.3f} corr {c:.3f}")
        out[name] = (int(m.sum()), float(np.median(r)), float(beta), float(c))
    rep("naive rp change", R("naive"))
    rep("minus no-burn baseline", R("naive") - R("base"))
    rep("3-pass level shift", R("shift3"))
    rep("altitude naive", R("dalt"))
    rep("altitude minus baseline", R("dalt") - R("basealt"))
    print(f"{label:6s} exact/linear lever mean {np.mean(Ex/L):.4f} (nonlinearity of the burn size)")
    return out

def impulse_response(a, e, f, alpha, dv, mu=MU):
    """Exact (d r_p, d e, d omega) of an in-plane impulse dv at true anomaly f, angle alpha from the velocity (+ = toward the left normal)."""
    p = a * (1 - e * e); r = p / (1 + e * math.cos(f))
    rv = np.array([r * math.cos(f), r * math.sin(f), 0.]); vv = math.sqrt(mu / p) * np.array([-math.sin(f), e + math.cos(f), 0.])
    def el(rr, v):
        h = np.cross(rr, v); ev = np.cross(v, h) / mu - rr / np.linalg.norm(rr); aa = -mu / (2 * (v @ v / 2 - mu / np.linalg.norm(rr)))
        return aa, np.linalg.norm(ev), math.atan2(ev[1], ev[0])
    u = vv / np.linalg.norm(vv); n = np.array([-u[1], u[0], 0.])
    a0, e0, o0 = el(rv, vv); a1, e1, o1 = el(rv, vv + dv * (math.cos(alpha) * u + math.sin(alpha) * n))
    return a1 * (1 - e1) - a0 * (1 - e0), e1 - e0, (o1 - o0 + math.pi) % (2 * math.pi) - math.pi

def direction_scan(dev, win=6):
    """FO burn pairs: baseline-corrected change of (r_p, e, omega) across the burn pair, regressed on the exact response of an impulse at
    true anomaly f and angle alpha from the velocity. A tangential apoapsis burn gives zero d omega; the measured d omega is not zero."""
    B = []
    for run in sorted(glob.glob(str(dev / "runs" / "FO_*"))):
        P, D = load(pathlib.Path(run)); col = lambda n: np.array([float(p[n]) for p in P])
        rp, a, e, om = col("rp_m"), col("a_m"), col("e"), col("omega_rad")
        dvk = {int(d["target_pass"]): float(d["delta_v_mps"]) for d in D if float(d["delta_v_mps"]) != 0}
        w = lambda x: (x + math.pi) % (2 * math.pi) - math.pi
        for k, dv in dvk.items():
            j = k - 1
            if j < 1 or j >= len(P): continue
            nb = [x for x in range(max(1, j - win), min(len(P), j + win + 1)) if x != j and (x + 1) not in dvk]
            if not nb: continue
            m = [rp[j] - rp[j - 1] - np.median([rp[x] - rp[x - 1] for x in nb]), e[j] - e[j - 1] - np.median([e[x] - e[x - 1] for x in nb]),
                 w(om[j] - om[j - 1]) - np.median([w(om[x] - om[x - 1]) for x in nb])]
            B.append((a[j - 1], e[j - 1], dv, m))
    meas = np.array([b[3] for b in B]); sl = lambda y, x: (x @ y) / (x @ x) if x @ x > 0 else float("nan")
    x = np.array([b[2] for b in B])
    print("measured slopes per m/s of dv: dr_p %.0f m, de %.2e, domega %.2e rad (corr %.2f %.2f %.2f)" % (
        sl(meas[:, 0], x), sl(meas[:, 1], x), sl(meas[:, 2], x), *[np.corrcoef(x, meas[:, i])[0, 1] for i in range(3)]))
    out = []
    for fd in (180, 175, 170, 165):
        for ad in range(-70, 71, 10):
            pr = np.array([impulse_response(a, e, math.radians(fd), math.radians(ad), dv) for a, e, dv, _ in B])
            out.append((fd, ad, [round(sl(meas[:, i], pr[:, i]), 2) for i in range(3)]))
    print("measured / predicted slope for (rp, e, omega); 1.0 = consistent. tangential at apoapsis (180, 0):", [o for o in out if o[:2] == (180, 0)])
    for o in out:
        if all(abs(v - 1) < 0.35 for v in o[2] if v == v): print("  consistent-ish:", o)

if __name__ == "__main__":
    if "--selftest" in sys.argv: sys.exit(0 if selftest() else 1)
    dev = pathlib.Path(sys.argv[sys.argv.index("--dev") + 1])
    direction_scan(dev)
    rows = analyze(dev)
    for arm in ("FO", "N2T", "N2I", "N21"):
        sub = [r for r in rows if r["arm"] == arm]
        if sub: summarize(sub, arm)
    summarize(rows, "ALL")
    ra = np.array([r["dra"] for r in rows if r["arm"] == "FO"])
    print("FO: r_a change across the burn pair, median", np.median(ra))
