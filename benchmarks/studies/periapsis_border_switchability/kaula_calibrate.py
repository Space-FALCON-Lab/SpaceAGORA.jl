"""Kaula calibration on the repo gravity files: degree amplitudes, fitted K per planet with scatter, roughness ratios.
Usage: python3 -I kaula_calibrate.py [--results-root DIR] [--gravity-dir DIR]
Writes <root>/kaula/{degree_amplitudes.csv, kaula_fit.csv, kaula_zonals.csv, fig_kaula_fit.png}.
Fit windows fixed in the script before the cell evaluation: W_ALL = [4, nmax], W_HI = [10, nmax]; the pre-registered K is the W_HI geometric mean."""
import argparse, pathlib, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
import kaula_lib as kl

def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--gravity-dir", default=None); a = ap.parse_args()
    out = pathlib.Path(results_root(a.results_root)) / "kaula"; out.mkdir(parents=True, exist_ok=True)
    amps, fits, zon = [], [], []
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for pl in ("Earth", "Mars", "Venus"):
        df, R, gm, meta = kl.read_gravity(pl, a.gravity_dir); nmax = int(df.degree.max())
        n, A = kl.degree_amplitude(df)
        # last degree with a nonzero amplitude (Venus file header says truncated; the rows say otherwise)
        nz = n[np.isfinite(A) & (A > 0)]; nmax_eff = int(nz.max())
        for d, v in zip(n, A): amps.append(dict(planet=pl, n=int(d), A_n=v, K_n=v * d * d))
        J = kl.zonal_J(df)
        for d in range(2, min(nmax_eff, 12) + 1): zon.append(dict(planet=pl, n=d, J_n=J[d], Cbar_n0=-J[d] / np.sqrt(2 * d + 1), A_n=A[d - 2], ratio_Cn0_over_An=abs(J[d] / np.sqrt(2 * d + 1)) / A[d - 2]))
        for lab, lo, hi in (("W_ALL", 4, nmax_eff), ("W_HI", 10, nmax_eff), ("W_20_50", 20, min(50, nmax_eff)), ("W_2_nmax", 2, nmax_eff)):
            K, sc, nu, sl = kl.fit_K(n, A, lo, hi)
            fits.append(dict(planet=pl, window=lab, nlo=lo, nhi=hi, n_used=nu, K=K, scatter_dex=sc, scatter_factor=10 ** sc, free_slope=sl))
        ax.loglog(n, A, ".", ms=4, label=f"{pl} (file to n={nmax_eff})")
        Kh = [f for f in fits if f["planet"] == pl and f["window"] == "W_HI"][0]["K"]
        ax.loglog(n, Kh / n ** 2, "-", lw=0.8, color=ax.lines[-1].get_color())
    ax.set_xlabel("degree n"); ax.set_ylabel("rms coefficient A_n (normalized)"); ax.legend(); ax.grid(alpha=.3, which="both"); ax.set_title("Degree amplitudes and K/n^2 (W_HI fit)")
    fig.tight_layout(); fig.savefig(out / "fig_kaula_fit.png", dpi=130)
    pd.DataFrame(amps).to_csv(out / "degree_amplitudes.csv", index=False); F = pd.DataFrame(fits); F.to_csv(out / "kaula_fit.csv", index=False)
    Z = pd.DataFrame(zon); Z.to_csv(out / "kaula_zonals.csv", index=False)
    pd.set_option("display.width", 200); print(F.to_string()); print(Z[Z.n <= 6].to_string())
main()
