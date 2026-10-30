"""Analyzer of halflever_vacuum.jl output (half-lever test I1): realized next-periapsis change against the exact vis-viva impulse and the formula lever.
    python3 -I halflever_vacuum_analyze.py --selftest
    python3 -I halflever_vacuum_analyze.py --dir <halflever_vacuum output dir>
Every number is DERIVED from the run files (states in the Mars body-equator inertial frame, mu = meta.csv mu_harm). Prints a table and writes analysis.csv."""
import csv, math, os, sys
import numpy as np

G0 = 9.80665
ISP = 220.0   # SOURCED: stageo_config.toml [planets.mars.vehicle] (also in burn.csv)


def elements(r, v, mu):
    """a, e, i, node, omega, nu (rad; omega, nu in [0, 2pi)) from the state."""
    r = np.asarray(r, float); v = np.asarray(v, float); rn = np.linalg.norm(r)
    h = np.cross(r, v); hn = np.linalg.norm(h); ev = np.cross(v, h) / mu - r / rn; e = np.linalg.norm(ev)
    a = 1.0 / (2.0 / rn - v @ v / mu); inc = math.acos(h[2] / hn)
    nvec = np.array([-h[1], h[0], 0.0]); nn = np.linalg.norm(nvec)
    node = math.atan2(nvec[1], nvec[0]) % (2 * math.pi)
    om = math.acos(np.clip(nvec @ ev / (nn * e), -1, 1)); om = om if ev[2] >= 0 else 2 * math.pi - om
    nu = math.acos(np.clip(ev @ r / (e * rn), -1, 1)); nu = nu if r @ v >= 0 else 2 * math.pi - nu
    return a, e, inc, node, om, nu


def lever(mu, a, e):
    return 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(mu / a)   # m of r_p per m/s, tangential apoapsis burn (stagek_lib.lever_m_per_mps)


def rp_of(r, v, mu):
    a, e = elements(r, v, mu)[:2]; return a * (1 - e)


def dang(x, y):
    return (x - y + math.pi) % (2 * math.pi) - math.pi


def rd(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def load_grid(path):
    d = rd(path); A = np.array([[float(r[k]) for k in ("t_s", "rx", "ry", "rz", "vx", "vy", "vz", "mass_kg")] for r in d]); return A


def row_at(A, t):
    """Grid sample nearest to t (the grids are 5 s inside the burn region, identical times in every run)."""
    i = int(np.argmin(abs(A[:, 0] - t))); return A[i]


def load_events(path):
    return [dict(kind=r["kind"], t=float(r["t_s"]), r=np.array([float(r[k]) for k in ("rx", "ry", "rz")]), v=np.array([float(r[k]) for k in ("vx", "vy", "vz")]), m=float(r["mass_kg"])) for r in rd(path)]


def analyze(dir_, name, ref, mu):
    X = os.path.join(dir_, name); R = os.path.join(dir_, ref)
    b = rd(os.path.join(X, "burn.csv"))[0]; dv = float(b["dv_cmd_mps"]); ws, we = float(b["win_start_s"]), float(b["win_stop_s"]); tc = 0.5 * (ws + we)
    GX, GR = load_grid(os.path.join(X, "grid.csv")), load_grid(os.path.join(R, "grid.csv"))
    EX, ER = load_events(os.path.join(X, "events.csv")), load_events(os.path.join(R, "events.csv"))
    sr = row_at(GR, tc); r_c, v_c = sr[1:4], sr[4:7]; a, e, _, _, om0, nu_c = elements(r_c, v_c, mu)
    vh = v_c / np.linalg.norm(v_c)
    # exact vis-viva impulse of the commanded dv at the burn-center state (tangential, along the inertial velocity)
    ex = elements(r_c, v_c + dv * vh, mu); drp_exact = ex[0] * (1 - ex[1]) - a * (1 - e); dom_exact = dang(ex[4], om0)
    drp_form = lever(mu, a, e) * dv
    # delivered dv by mass and by the velocity difference to the no-burn flight after the burn (first common sample at or after the stop)
    mB = row_at(GX, ws - 10)[7]; mA = row_at(GX, we + 10)[7]; dv_mass = ISP * G0 * math.log(mB / mA)
    ia = int(np.argmax(GX[:, 0] >= we + 10)); tb = GX[ia, 0]; sx = GX[ia]; sn = GR[int(np.argmin(abs(GR[:, 0] - tb)))]
    assert abs(sn[0] - tb) < 1e-6
    dvec = sx[4:7] - sn[4:7]; vn = sn[4:7] / np.linalg.norm(sn[4:7]); hh = np.cross(sn[1:4], sn[4:7]); hh /= np.linalg.norm(hh); rad = np.cross(vn, hh)
    tang, nrm, rdl = dvec @ vn, dvec @ hh, dvec @ rad
    ang = math.degrees(math.acos(np.clip(dvec @ vn / np.linalg.norm(dvec), -1, 1)))
    # osculating change after the burn (same time, burn minus no burn)
    eX = elements(sx[1:4], sx[4:7], mu); eN = elements(sn[1:4], sn[4:7], mu)
    drp_osc = eX[0] * (1 - eX[1]) - eN[0] * (1 - eN[1]); dom_osc = dang(eX[4], eN[4]); de_osc = eX[1] - eN[1]; da_osc = eX[0] - eN[0]
    # next periapsis event (third event: periapsis, apoapsis, periapsis)
    pX = [x for x in EX if x["kind"] == "periapsis"][-1]; pN = [x for x in ER if x["kind"] == "periapsis"][-1]
    drp_ev = np.linalg.norm(pX["r"]) - np.linalg.norm(pN["r"]); dom_ev = dang(elements(pX["r"], pX["v"], mu)[4], elements(pN["r"], pN["v"], mu)[4])
    aN = [x for x in ER if x["kind"] == "apoapsis"][0]; dt_apo = tc - aN["t"]; nu_off = math.degrees(dang(nu_c, math.pi))
    return dict(run=name, dv_cmd=dv, dv_mass=dv_mass, dv_vec_norm=float(np.linalg.norm(dvec)), dv_tang=tang, dv_normal=nrm, dv_radial=rdl, thrust_angle_to_v_deg=ang,
        win_start=ws, win_stop=we, center_minus_apo_s=dt_apo, nu_center_minus_pi_deg=nu_off, a_m=a, e=e, rp_before_m=a * (1 - e),
        drp_exact_m=drp_exact, drp_formula_m=drp_form, drp_osc_m=drp_osc, drp_event_m=drp_ev,
        ratio_osc_exact=drp_osc / drp_exact, ratio_event_exact=drp_ev / drp_exact, ratio_event_formula=drp_ev / drp_form, ratio_osc_formula=drp_osc / drp_form,
        dom_exact_rad=dom_exact, dom_osc_rad=dom_osc, dom_event_rad=dom_ev, da_osc_m=da_osc, de_osc=de_osc)


def selftest():
    mu = 4.2828370371e13; a = 1.6007744e7; rp = 3.4825e6; e = (2 * a - 2 * rp) / (2 * a); ra = a * (1 + e)
    r = np.array([ra, 0, 0]); vt = math.sqrt(mu * (2 / ra - 1 / a)); u = np.array([0, math.cos(1.6), math.sin(1.6)]); v = vt * u   # inclined (node defined)
    for dv in (0.1, 0.25):
        el = elements(r, v + dv * u, mu); d = el[0] * (1 - el[1]) - a * (1 - e)
        assert abs(d / (lever(mu, a, e) * dv) - 1) < 2e-3, (dv, d, lever(mu, a, e) * dv)
        assert abs(dang(el[4], elements(r, v, mu)[4])) < 1e-9
    print("selftest ok: exact impulse = formula lever within 2e-3, tangential burn at apoapsis leaves omega unchanged")


def main(argv):
    if "--selftest" in argv:
        return selftest()
    dir_ = argv[argv.index("--dir") + 1]; mu = float(rd(os.path.join(dir_, "meta.csv"))[0]["mu_harm"])
    cases = [("nom", "none"), ("dv014", "none"), ("tp60", "none"), ("tm60", "none"), ("nom_sun", "none_sun")]
    rows = [analyze(dir_, n, r, mu) for n, r in cases]
    with open(os.path.join(dir_, "analysis.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print("\n".join(f"{k:>24} {v:.6g}" if isinstance(v, float) else f"{k:>24} {v}" for k, v in r.items())); print()
    n, d = rows[0], rows[1]
    print(f"linearity: event drp(0.14)/drp(0.25) = {d['drp_event_m'] / n['drp_event_m']:.4f} (0.14/0.25 = {0.14 / 0.25:.4f}); exact-prediction ratio {d['drp_exact_m'] / n['drp_exact_m']:.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
