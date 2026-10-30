"""Part B analysis: hand-estimate comparison, correlations, heat maps. Usage: python3 analyze_grid.py [OUT_DIR]"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Part B analysis")
_ap.add_argument("out_dir", nargs="?", help="directory to read/write (default: the results root)")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir or results_root(_args.results_root)
g = pd.read_csv(f"{OUT}/grid_summary.csv")
eff = g.third_body_effect_geodetic_m / 1e3
print(f"cells={len(g)}  failed retcodes={(g.retcode_3b!='Success').sum()}+{(g.retcode_ctl!='Success').sum()}")
print(f"third-body effect on geodetic r_p over ~60 d (km): max|.|={eff.abs().max():.1f}  median|.|={eff.abs().median():.1f}  min|.|={eff.abs().min():.1f}")
print(f"no-third-body change (km): max|.|={g.drp_geodetic_ctl_m.abs().max()/1e3:.2f} median|.|={g.drp_geodetic_ctl_m.abs().median()/1e3:.2f}")
ratio = (g.third_body_effect_geodetic_m.abs() / g.drp_geodetic_ctl_m.abs().clip(lower=1.0))
print(f"ratio |3B effect| / |control change|: median={ratio.median():.1f} min={ratio.min():.1f} max={ratio.max():.1f}")
print(f"3B 60 d change itself (km): max|.|={g.drp_geodetic_3b_m.abs().max()/1e3:.1f}")
# correlation with sin(2 omega') : cell-mean of the combined weighted indicator and the double-averaged hand prediction
g["s2_comb"] = g.mean_sin2wp_moon + g.mean_sin2wp_sun
for col, lab in [("s2_comb", "mean sin2w'(Moon)+mean sin2w'(Sun)"), ("mean_sin2wp_moon", "mean sin2w'(Moon)"), ("mean_sin2wp_sun", "mean sin2w'(Sun)")]:
    print(f"corr(3B-ctl r_p effect, {lab}) = {np.corrcoef(g[col], g.third_body_effect_geodetic_m)[0,1]:+.3f}  (expected negative: r_p change = -a de/dt ~ -sin2w')")
print(f"corr(hand double-averaged prediction, measured radius effect) = {np.corrcoef(g.kozai_predicted_radius_change_m, g.third_body_effect_radius_m)[0,1]:+.3f}")
fit = np.polyfit(g.kozai_predicted_radius_change_m / 1e3, g.third_body_effect_radius_m / 1e3, 1)
print(f"measured_radius_effect = {fit[0]:.3f} * predicted + {fit[1]:.2f} km")
print(f"sign agreement hand prediction vs measured: {(np.sign(g.kozai_predicted_radius_change_m)==np.sign(g.third_body_effect_radius_m)).mean()*100:.0f}%")
print("hand bounds: 170-230 km/60 d best geometry; measured max |3B effect| radius: %.1f km, geodetic %.1f km" % (g.third_body_effect_radius_m.abs().max()/1e3, eff.abs().max()))
print(f"in-band days 3B: median={g.in_band_days_3b.median():.1f} max={g.in_band_days_3b.max():.1f}; ctl median={g.in_band_days_ctl.median():.1f}")
print(f"first-exit day (3B): median={g.first_exit_days_3b.median():.1f} min={g.first_exit_days_3b.min():.1f}; edge counts: {g.first_exit_edge_3b.value_counts().to_dict()}")
print(f"cells exceeding Q_hard: {int(g.exceeds_hard_3b.sum())}")
print(g[["tag","third_body_effect_geodetic_m","drp_geodetic_ctl_m","kozai_predicted_radius_change_m","third_body_effect_radius_m","mean_sin2wp_moon","mean_sin2wp_sun","in_band_days_3b","first_exit_days_3b","first_exit_edge_3b","min_alt_3b_m","max_alt_3b_m","apo_end_3b_m","apo_end_ctl_m"]].round(1).to_string())
g.to_csv(f"{OUT}/grid_summary_with_derived.csv", index=False)

# heat map: delta r_p (60 d) vs (RAAN, omega), per epoch; rows = with / without third bodies
raans, omegas = sorted(g.raan_deg.unique()), sorted(g.omega_deg.unique())
epochs = sorted(g.epoch_idx.unique())
fig, ax = plt.subplots(2, len(epochs), figsize=(4.6 * len(epochs), 8), squeeze=False, constrained_layout=True)
vmax = np.nanmax(np.abs(g.drp_geodetic_3b_m)) / 1e3
for j, e in enumerate(epochs):
    for i, (col, lab) in enumerate([("drp_geodetic_3b_m", "with Sun+Moon"), ("drp_geodetic_ctl_m", "no third bodies")]):
        z = np.full((len(raans), len(omegas)), np.nan)
        for r in g[g.epoch_idx == e].itertuples():
            z[raans.index(r.raan_deg), omegas.index(r.omega_deg)] = getattr(r, col) / 1e3
        im = ax[i, j].imshow(z, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax[i, j].set_xticks(range(len(omegas)), [f"{int(o)}" for o in omegas]); ax[i, j].set_yticks(range(len(raans)), [f"{int(r)}" for r in raans])
        ax[i, j].set_xlabel("argument of periapsis (deg)"); ax[i, j].set_ylabel("RAAN (deg)")
        ax[i, j].set_title(f"{g[g.epoch_idx == e].epoch.iloc[0][:10]}  {lab}", fontsize=10)
        for a in range(len(raans)):
            for b in range(len(omegas)):
                if np.isfinite(z[a, b]):
                    cell = g[(g.epoch_idx == e) & (g.raan_deg == raans[a]) & (g.omega_deg == omegas[b])].iloc[0]
                    star = "*" if (i == 0 and cell.retcode_3b != "Success") else ""   # * = run terminated early (orbit decayed), value is at termination
                    ax[i, j].text(b, a, f"{z[a, b]:.1f}{star}", ha="center", va="center", fontsize=8)
fig.colorbar(im, ax=ax, shrink=0.8, label="* = 3B run ended early (periapsis fell, orbit decayed); value at that time.  Measured change of geodetic periapsis altitude over ~60 d (km)")
fig.savefig(f"{OUT}/fig_heatmap_drp60.png", dpi=130)

# scatter: measured vs hand prediction
fig, ax = plt.subplots(figsize=(5.5, 5), constrained_layout=True)
ax.scatter(g.kozai_predicted_radius_change_m / 1e3, g.third_body_effect_radius_m / 1e3, c=g.epoch_idx, cmap="viridis")
lim = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]), abs(ax.get_ylim()[0]), abs(ax.get_ylim()[1]))
ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8); ax.set_xlabel("double-averaged hand prediction, Sun+Moon (km)"); ax.set_ylabel("measured 3B - control, radius change (km)")
ax.axhline(0, c="gray", lw=0.5); ax.axvline(0, c="gray", lw=0.5)
fig.savefig(f"{OUT}/fig_hand_vs_measured.png", dpi=130)
