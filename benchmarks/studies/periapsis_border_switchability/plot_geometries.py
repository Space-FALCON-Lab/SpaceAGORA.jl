"""Figure 1: periapsis altitude and peak Q vs time for the chosen geometries, band shaded. Also prints band altitudes vs speed.
Usage: python3 plot_geometries.py [OUT_DIR]"""
import sys, tomllib, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Figure 1: geometries")
_ap.add_argument("out_dir", nargs="?", help="directory to read/write (default: the results root)")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir or results_root(_args.results_root)
study = tomllib.load(open(f"{HERE}/study_config.toml", "rb"))
atm = tomllib.load(open(f"{HERE}/../../../examples/periapsis_corridor_config.toml", "rb"))["atmosphere"]
tab_h, tab_rho = np.array(atm["table_altitudes_m"]), np.array(atm["table_density_kg_m3"])
nodes = np.array(atm["node_altitudes_m"]); nrho = np.array([tab_rho[list(tab_h).index(h)] for h in nodes])
H = np.diff(nodes) / np.log(nrho[:-1] / nrho[1:])
def h_of_rho(rho):  # inverse of the piecewise fit; above the top node the last layer is extended
    k = np.clip(np.searchsorted(-np.log(nrho), -np.log(rho)) - 1, 0, len(H) - 1)
    return nodes[k] - H[k] * np.log(rho / nrho[k])
def band_alt(v, q):  return h_of_rho(2 * q * study["band"]["w_cm2_to_w_m2"] / v ** 3)
qlo, qhi, qh = study["band"]["q_low_w_cm2"], study["band"]["q_high_w_cm2"], study["band"]["q_hard_w_cm2"]
geo = pd.read_csv(f"{OUT}/geometries.csv")
sel = geo[geo.role != "worst"]
fig, ax = plt.subplots(2, len(sel), figsize=(4.2 * len(sel), 7), squeeze=False, constrained_layout=True)
for j, r in enumerate(sel.itertuples()):
    p3 = pd.read_csv(f"{OUT}/grid/{r.tag}/passes_3b.csv"); pc = pd.read_csv(f"{OUT}/grid/{r.tag}/passes_ctl.csv")
    for c, (p, lab, ls) in enumerate([(p3, "Sun+Moon", "-"), (pc, "no third bodies", "--")]):
        ax[0, j].plot(p.time_s / 86400, p.alt_m / 1e3, ls, label=lab, color=f"C{c}")
        ax[1, j].semilogy(p.time_s / 86400, p.q_w_cm2, ls, color=f"C{c}")
    v = p3.vrel_mps.values
    ax[0, j].fill_between(p3.time_s / 86400, band_alt(v, qhi) / 1e3, band_alt(v, qlo) / 1e3, color="green", alpha=0.18, label="band 0.32-0.14 W/cm2")
    ax[0, j].axhline(band_alt(v[0], qh) / 1e3, color="red", lw=0.8, ls=":", label="hard 0.54 W/cm2")
    ax[1, j].axhspan(qlo, qhi, color="green", alpha=0.18); ax[1, j].axhline(qh, color="red", lw=0.8, ls=":")
    ax[0, j].set_title(f"{r.role}: {r.tag}", fontsize=10); ax[1, j].set_xlabel("time (days)")
    ax[1, j].set_ylim(0.03, 3)
ax[0, 0].set_ylabel("geodetic periapsis altitude (km)"); ax[1, 0].set_ylabel("peak Q = 0.5 rho V_rel^3 (W/cm2)")
ax[0, 0].legend(fontsize=7)
fig.savefig(f"{OUT}/fig_geometries_alt_q.png", dpi=130)
# band altitudes at the initial speeds (for the README)
g = pd.read_csv(f"{OUT}/grid_summary.csv")
for v in (g.initial_vrel_mps.min(), g.initial_vrel_mps.median(), g.initial_vrel_mps.max()):
    print(f"V_rel={v:.0f} m/s: h(0.54)={band_alt(v,qh)/1e3:.2f} km  h(0.32)={band_alt(v,qhi)/1e3:.2f}  h(0.14)={band_alt(v,qlo)/1e3:.2f}  width={(band_alt(v,qlo)-band_alt(v,qhi))/1e3:.2f} km")
