"""Reference values for the estimator-v4 selftests (a), (b): per-pass tesseral kicks of the Python reference (tesseral_lib.py) at fixed states.

Usage: python3 -I v4_reference_kicks.py [--results-root DIR]
Writes <root>/v4/reference_kicks.csv (one row per state and field subset: Gauss quadrature of 512 nodes `gauss_kick`, the inbound-half r_p kick `term_kicks(half=True)`,
and for the subset field (degree <= 6, order <= 4) the Lagrange-series kick `series_kick`) and <root>/v4/reference_unitY.csv (|Y_nm| of the unit-kick quadrature, 2 states).
The states are fixed here (ASSUMED, Odyssey-like Mars orbits: a, e from the Stage C pass 1; inclinations as the Stage C grid)."""
import argparse, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import tesseral_lib as tl

SPIN = 7.08823596e-5          # SOURCED: planets.jl Mars struct (as tesseral_param.py)
# (a_m, e, inc_deg, raan_deg, omega_deg, theta_p_rad)
STATES = [(16065.2e3, 0.78186, 25.2, 10.0, 30.0, 0.4), (16065.2e3, 0.78186, 25.2, 200.0, 250.0, 2.9), (15500.0e3, 0.7550, 63.4, 100.0, 80.0, 1.7),
          (14800.0e3, 0.7280, 93.522, 300.0, 150.0, 5.1), (16065.2e3, 0.78186, 116.6, 45.0, 310.0, 3.3), (13000.0e3, 0.6400, 45.0, 120.0, 200.0, 0.9),
          (16065.2e3, 0.78186, 80.0, 250.0, 20.0, 4.4), (12000.0e3, 0.5500, 30.0, 330.0, 100.0, 6.0), (9000.0e3, 0.3700, 63.4, 15.0, 60.0, 2.2),
          (16065.2e3, 0.78186, 93.522, 77.0, 288.0, 1.1)]


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    out = pathlib.Path(results_root(a.results_root)) / "v4"; out.mkdir(parents=True, exist_ok=True)
    cs = tl.coefficients("Mars", 20); mu, R = cs["mu"], cs["R"]
    Cn, Sn = tl.select(cs["C"], cs["S"], lambda n, m: m >= 1)
    Cs, Ss = tl.select(cs["C"], cs["S"], lambda n, m: 1 <= m <= min(n, 4) and n <= 6)
    coef = cs["Cbar"] - 1j * cs["Sbar"]; coef[:, 0] = 0
    rows = []; Yrows = []
    for k, (A, e, inc_d, O_d, w_d, thp) in enumerate(STATES):
        inc, O, w = map(math.radians, (inc_d, O_d, w_d))
        for fld, (C, S) in (("full20", (Cn, Sn)), ("sub6x4", (Cs, Ss))):
            g = tl.gauss_kick(A, e, inc, O, w, thp, SPIN, cs, C, S, nf=512)
            Y = tl.term_kicks(A, e, inc, O, w, thp, SPIN, mu, R, 20, nf=256, half=True)
            cm = coef.copy()
            if fld == "sub6x4":
                for n in range(21):
                    for m in range(21):
                        if not (n <= 6 and 1 <= m <= min(n, 4)): cm[n, m] = 0
            half_rp = float(np.real((cm * Y).sum()))
            row = dict(state=k, field=fld, a=A, e=e, inc_deg=inc_d, raan_deg=O_d, omega_deg=w_d, theta_p=thp, **{f"quad_{kk}": v for kk, v in g.items()}, quad_rp_half=half_rp)
            if fld == "sub6x4":
                N = math.sqrt(mu / A ** 3); nu = SPIN / N; tot = dict.fromkeys(tl.KICK_NAMES, 0.0)
                for n in range(2, 7):
                    for m in range(1, min(n, 4) + 1):
                        T = tl.SeriesTerm(n, m, e, inc, NM=8192)
                        sk = tl.series_kick(T, A, O, w, thp, nu, cs)
                        for kk in tot: tot[kk] += sk[kk]
                row.update({f"series_{kk}": v for kk, v in tot.items()})
            rows.append(row)
        if k in (0, 3):
            Y = tl.term_kicks(A, e, inc, 0.0, w, 0.0, SPIN, mu, R, 20, nf=96)
            for n in range(2, 21):
                for m in range(1, n + 1): Yrows.append(dict(state=k, n=n, m=m, absY=float(abs(Y[n, m]))))
    pd.DataFrame(rows).to_csv(out / "reference_kicks.csv", index=False); pd.DataFrame(Yrows).to_csv(out / "reference_unitY.csv", index=False)
    print(pd.DataFrame(rows)[["state", "field", "quad_rp", "quad_rp_half"]].to_string())


if __name__ == "__main__":
    main()
