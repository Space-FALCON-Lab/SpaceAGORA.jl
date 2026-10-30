"""Self-tests S0-S5 of the apoapsis Pi closed forms against exact Kepler/Gauss quadrature (no truth data).
Criteria: tasks/2026-10-15_periapsis-region-border-targeting/apoapsis_pi_preregistration.md.
Usage: python3 -I apoapsis_pi_selftest.py [--results-root DIR]  -> <root>/apoapsis_pi/selftest.csv
"""
import argparse, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root, REPO
import apoapsis_pi_lib as L

MU_MARS = 4.282837285418775e13   # SOURCED: Mars() struct, planets.jl
H_MARS = 6.308278108e3           # SOURCED: Mars() struct H
MU_EARTH = 3.98600436233e14      # SOURCED: Earth() struct
Q_LO, Q_HI = 0.14e4, 0.32e4      # SOURCED: study_config.toml [band], W/m2


def earth_h_layer(h):
    """Scale height of the corridor US76 fit layer containing h (examples/periapsis_corridor_config.toml)."""
    with open(pathlib.Path(REPO) / "examples" / "periapsis_corridor_config.toml", "rb") as f:
        atm = tomllib.load(f)["atmosphere"]
    th, tr, nodes = np.array(atm["table_altitudes_m"]), np.array(atm["table_density_kg_m3"]), np.array(atm["node_altitudes_m"])
    nr = np.array([tr[np.where(th == n)[0][0]] for n in nodes]); Hs = np.diff(nodes) / np.log(nr[:-1] / nr[1:])
    return Hs[int(np.clip(np.searchsorted(nodes, h, side="right") - 1, 0, len(Hs) - 1))]


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    out = pathlib.Path(results_root(a.results_root)) / "apoapsis_pi"; out.mkdir(parents=True, exist_ok=True)
    rows = []
    def rec(test, case, value, ref, ok, note=""):
        rows.append(dict(test=test, case=case, value=value, reference=ref, passed=bool(ok), note=note))

    H_E = earth_h_layer(140e3)
    planets = {"Mars": (MU_MARS, H_MARS, 3505e3), "Earth": (MU_EARTH, H_E, 6500e3)}
    # S0: Pi_E written through the existing groups equals the direct form (algebra)
    for pl, (mu, H, rp) in planets.items():
        for e in (0.1, 0.5, 0.8):
            a_ = rp / (1 - e); rho = 1e-8; beta = 20.0; Vp = np.sqrt(mu * (1 + e) / rp); Q = 0.5 * rho * Vp**3
            direct = L.pi_E(a_, e, rp, Q, H, beta, mu)
            pidrag = 0.5 * rho * rp * (1 + e) / beta
            via = 2 * pidrag * np.sqrt(2 * np.pi * (H / rp) * (1 + e) / e) / (1 - e)
            rec("S0", f"{pl} e={e}", direct, via, abs(direct / via - 1) < 1e-12)
    # S1: Laplace closed form vs exact Kepler-path quadrature, criterion |ratio-1| <= 0.01 where c = a e / H >= 25
    for pl, (mu, H, rp) in planets.items():
        for e in (0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9):
            a_ = rp / (1 - e); rho = 1e-8; Vp = np.sqrt(mu * (1 + e) / rp)
            r = L.delta_a_quad(a_, e, H, rho, 20.0, mu) / L.delta_a_pass(a_, e, rp, 0.5 * rho * Vp**3, H, 20.0, mu)
            c = a_ * e / H
            rec("S1", f"{pl} e={e} c={c:.1f}", r, 1.0, (abs(r - 1) <= 0.01) if c >= 25 else True, "" if c >= 25 else "c<25: outside validity, reported only")
    # S2: passes and time to go, discrete exact-quadrature map at constant Q and r_p vs closed form (Odyssey-like range, Mars)
    mu, H, rp = planets["Mars"]; beta = 20.0
    for Qw in (0.04e4, 0.13e4, 0.24e4):
        ra, ra_t = 28533e3, 3900e3; n = 0.0; t = 0.0
        e0 = (ra - rp) / (ra + rp); et = (ra_t - rp) / (ra_t + rp)
        while True:
            a_, e = L.elements_from_x(rp, ra / rp); Vp = np.sqrt(mu * (1 + e) / rp)
            dra = 2 * L.delta_a_quad(a_, e, H, 2 * Qw / Vp**3, beta, mu); T = 2 * np.pi * np.sqrt(a_**3 / mu)
            if ra + dra <= ra_t:
                fr = (ra - ra_t) / -dra; n += fr; t += fr * T; break
            ra += dra; n += 1; t += T
        kq = L.k_Q(rp, Qw, H, beta, mu)
        Nc, tc = L.passes_to_go(e0, et, kq), L.time_to_go(rp, e0, et, kq, mu)
        rec("S2", f"passes Q={Qw/1e4}", n, Nc, abs(n / Nc - 1) <= 0.02)
        rec("S2", f"time_s Q={Qw/1e4}", t, tc, abs(t / tc - 1) <= 0.02)
    # S3: drag-driven periapsis drop per pass vs first-order Gauss quadrature; criterion 10 % where c >= 25
    for pl, (mu, H, rp) in planets.items():
        for e in (0.05, 0.1, 0.2, 0.4, 0.8):
            a_ = rp / (1 - e); rho = 1e-8; Vp = np.sqrt(mu * (1 + e) / rp)
            da, de, drp = L.gauss_pass_quad(a_, e, H, rho, 20.0, mu)
            lam = L.lam_rp(e, -da / a_); c = a_ * e / H
            r = (-drp / H) / lam
            rec("S3", f"{pl} e={e} c={c:.1f}", r, 1.0, (abs(r - 1) <= 0.10) if c >= 25 else True)
    # S4: apoapsis-burn periapsis lever vs vis-viva finite difference; spec section 5.1 form reported
    mu = MU_MARS
    for e in (0.05, 0.3, 0.78):
        a_ = 3505e3 / (1 - e); ra = a_ * (1 + e); va = np.sqrt(mu * (2 / ra - 1 / a_)); dv = 1e-3
        a2 = 1 / (2 / ra - (va + dv) ** 2 / mu); fd = ((2 * a2 - ra) - (2 * a_ - ra)) / (dv / np.sqrt(mu / a_))
        rec("S4", f"lever e={e}", fd, L.lever_rp(a_, e), abs(fd / L.lever_rp(a_, e) - 1) <= 1e-4)
        rec("S4-spec", f"spec 2a(1-e) e={e}", 2 * a_ * (1 - e), fd, True, "reported: ratio spec/FD")
    # S5: corridor width in scale heights reproduces the band widths quoted in the study (Mars 5.215 km, Earth 12.09 km)
    w = L.corridor_width_H(Q_LO, Q_HI)
    rec("S5", "Mars band km", w * H_MARS / 1e3, 5.215, abs(w * H_MARS / 1e3 - 5.215) <= 0.01)
    rec("S5", "Earth band km", w * H_E / 1e3, 12.09, abs(w * H_E / 1e3 - 12.09) <= 0.01, f"H(130-150 km layer) = {H_E:.1f} m")
    df = pd.DataFrame(rows); df.to_csv(out / "selftest.csv", index=False)
    pd.set_option("display.width", 200); print(df.to_string(index=False))
    crit = df[~df.test.str.endswith("spec")]
    for t, g in crit.groupby("test"):
        print(f"{t}: {'PASS' if g.passed.all() else 'FAIL'} ({int(g.passed.sum())}/{len(g)})")


if __name__ == "__main__":
    main()
