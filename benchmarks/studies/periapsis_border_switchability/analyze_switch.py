"""Parts C/D analysis: finite-difference table, switch_dv, exit/return times, direct-step comparison, plots.
Usage: python3 analyze_switch.py [OUT_DIR]"""
import sys, tomllib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Parts C/D analysis")
_ap.add_argument("out_dir", nargs="?", help="directory to read/write (default: the results root)")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir or results_root(_args.results_root)
STUDY = tomllib.load(open(__file__.replace("analyze_switch.py", "study_config.toml"), "rb"))
SW = STUDY["switch"]
fd = pd.read_csv(f"{OUT}/switch_fd.csv")
geo = pd.read_csv(f"{OUT}/geometries.csv").set_index("tag")
ref = fd[fd.direction == "reference"].set_index(["tag", "start", "location"])
k = fd[fd.direction != "reference"].copy()
k["per_mps"] = k.delta_rate_km_day / k.dv_mps          # km/day per m/s, one value per probe
k["step_per_mps"] = k.d_first_km / k.dv_mps            # km per m/s of first-pass periapsis step
k["domega_per_mps"] = k.kick_domega_rad / k.dv_mps
k["domega_end_per_mps"] = k.end_domega_rad / k.dv_mps

rows = []
for (tag, start, loc, d), grp in k.groupby(["tag", "start", "location", "direction"]):
    slope = (grp.delta_rate_km_day * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum()   # least squares through origin
    pos = grp[grp.dv_mps > 0].set_index(grp[grp.dv_mps > 0].dv_mps.abs()).per_mps
    neg = grp[grp.dv_mps < 0].set_index(grp[grp.dv_mps < 0].dv_mps.abs()).per_mps
    sizes = sorted(pos.index)
    per_size = [(pos[s] + neg[s]) / 2 for s in sizes]        # odd part: removes the even (second-order) response
    even = [(pos[s] - neg[s]) / 2 for s in sizes]
    rows.append(dict(tag=tag, role=grp.role.iloc[0], start=start, location=loc, direction=d,
        rate_per_mps=slope, rate_small=per_size[0], rate_mid=per_size[1], rate_large=per_size[-1],
        linearity_spread=(max(per_size) - min(per_size)) / max(abs(slope), 1e-12),
        even_over_odd_at_1=abs(even[-1]) / max(abs(per_size[-1]), 1e-12),
        step_km_per_mps=(grp.d_first_km * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        domega_rad_per_mps=(grp.kick_domega_rad * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        domega_end_rad_per_mps=(grp.end_domega_rad * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        da_m_per_mps=(grp.kick_da_m * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        de_per_mps=(grp.kick_de * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        di_rad_per_mps=(grp.kick_di_rad * grp.dv_mps).sum() / (grp.dv_mps ** 2).sum(),
        ref_slope_km_day=ref.loc[(tag, start, loc), "ref_slope_km_day"], ref_q0=ref.loc[(tag, start, loc), "ref_q0"]))
T = pd.DataFrame(rows)
# Step compensation. Every kick also steps the next periapsis (step_km_per_mps). The tangential-at-apoapsis kick moves the periapsis by
# step_tan km per m/s and the drift rate by rate_tan km/day per m/s; cancelling a kick's step with such a tangential kick (superposition,
# linear response assumed) costs |step_i / step_tan| m/s per m/s and adds -step_i * rate_tan / step_tan to the rate. The compensated
# lever is the net rate change per total m/s spent. For the tangential-at-apoapsis kick itself the net is exactly zero.
tan = T[(T.direction == "tangential") & (T.location == "apoapsis")].set_index(["tag", "start"])
T["step_tan"] = [tan.loc[(a, b), "step_km_per_mps"] for a, b in zip(T.tag, T.start)]
T["rate_tan"] = [tan.loc[(a, b), "rate_per_mps"] for a, b in zip(T.tag, T.start)]
T["comp_dv_per_mps"] = (T.step_km_per_mps / T.step_tan).abs()
T["rate_net_per_mps_total"] = (T.rate_per_mps - T.step_km_per_mps * T.rate_tan / T.step_tan) / (1 + T.comp_dv_per_mps)
T.to_csv(f"{OUT}/switch_derivatives.csv", index=False)
print("== periapsis step per m/s of each lever (km per m/s) ==")
print(T.pivot_table(index=["tag", "start"], columns=["location", "direction"], values="step_km_per_mps").map(lambda x: f"{x:+.3f}").to_string())
print("== step-compensated net rate change per total m/s [km/day per m/s] ==")
print(T.pivot_table(index=["tag", "start"], columns=["location", "direction"], values="rate_net_per_mps_total").map(lambda x: f"{x:+.2e}").to_string())
pd.set_option("display.width", 250)
print("== d(drift rate)/d(dv) [km/day per m/s], least squares over +-{0.05,0.2,1.0} m/s; spread = (max-min)/|fit| across sizes of the odd part ==")
print(T.pivot_table(index=["tag", "start"], columns=["location", "direction"], values="rate_per_mps").map(lambda x: f"{x:+.2e}").to_string())
print("== linearity spread (max over all) =", round(T.linearity_spread.max(), 3), " median", round(T.linearity_spread.median(), 3),
      "; even/odd response at 1 m/s: median", round(T.even_over_odd_at_1.median(), 3))
print("== lever table: direct periapsis step per m/s tangential at apoapsis (km per m/s) ==")
print(T[(T.direction == "tangential") & (T.location == "apoapsis")][["tag", "start", "step_km_per_mps", "rate_per_mps", "da_m_per_mps", "de_per_mps"]].to_string())
print("== apse-line rotation lever, immediate osculating (rad per m/s); hand = %.3g ==" % SW["hand_apse_rotation_rad_per_mps"])
print(T[T.direction == "radial"].pivot_table(index=["tag", "start"], columns="location", values="domega_rad_per_mps").to_string())
print("   long-term (omega difference after the window) per m/s:")
print(T[T.direction == "radial"].pivot_table(index=["tag", "start"], columns="location", values="domega_end_rad_per_mps").to_string())
print("== normal kick di (rad per m/s) ==")
print(T[T.direction == "normal"].pivot_table(index=["tag", "start"], columns="location", values="di_rad_per_mps").to_string())

# ---- switchability ---------------------------------------------------------------------------
probe_max = max(SW["probe_dv_mps"])
def time_to(dist, rate):          # time to cover dist (>0) at signed rate; inf if moving away
    return dist / rate if rate > 0 else np.inf
sw_rows = []
for (tag, start), grp in T.groupby(["tag", "start"]):
    best = grp.loc[grp.rate_net_per_mps_total.abs().idxmax()]  # step-compensated best lever
    raw = grp.loc[grp.rate_per_mps.abs().idxmax()]
    r0 = best.ref_slope_km_day
    pas = pd.read_csv(f"{OUT}/kicks/{tag}/{start}/{best.location}/passes_ref.csv")
    h_lo, h_hi = pas.h_qhigh_m.iloc[0] / 1e3, pas.h_qlow_m.iloc[0] / 1e3     # altitude of Q_high (low edge), of Q_low (high edge)
    alt0 = pas.alt_m.iloc[0] / 1e3
    width = h_hi - h_lo
    step_row = grp[(grp.direction == "tangential") & (grp.location == "apoapsis")].iloc[0]
    cost_per_km = 1.0 / abs(step_row.step_km_per_mps)
    for target in SW["target_delta_rate_km_day"]:
        for sign, name in ((+1, "raise"), (-1, "lower")):
            dv_kick = target / abs(best.rate_net_per_mps_total) / (1 + best.comp_dv_per_mps)  # the kick itself
            dv = target / abs(best.rate_net_per_mps_total)   # total spent incl. step compensation
            dv_raw = target / abs(raw.rate_per_mps)
            new_rate = r0 + sign * target
            # distance to the edge being approached, from band centre (actual start state) and from the opposite edge
            d_center = (h_hi - alt0) if sign > 0 else (alt0 - h_lo)
            d_edge = width
            back_rate = r0 - sign * target       # the opposite switch, used for the return
            sw_rows.append(dict(tag=tag, role=best.role, start=start, direction=name, delta_rate_km_day=target,
                best_direction=best.direction, best_location=best.location, rate_per_mps=best.rate_per_mps, rate_net_per_mps_total=best.rate_net_per_mps_total,
                raw_best_direction=raw.direction, raw_best_location=raw.location, raw_switch_dv_mps=dv_raw, kick_dv_mps=dv_kick, comp_dv_mps=dv - dv_kick, step_km_of_kick=best.step_km_per_mps * dv_kick,
                switch_dv_mps=dv, extrapolated=dv_kick > probe_max, extrapolation_factor=dv_kick / probe_max,
                ref_rate_km_day=r0, band_width_km=width, dist_center_to_edge_km=d_center,
                exit_time_center_days=time_to(d_center, sign * new_rate), exit_time_edge_days=time_to(d_edge, sign * new_rate),
                return_dv_mps=dv, return_time_days=time_to(d_center, -sign * back_rate),
                direct_cost_mps_per_km=cost_per_km, direct_step_same_disp_mps=cost_per_km * d_center,
                direct_step_10km_mps=cost_per_km * SW["direct_step_km"]))
S = pd.DataFrame(sw_rows)
S["switch_dv_max_mps"] = SW["switch_dv_max_fraction_of_direct_step"] * S.direct_step_10km_mps
S["switchable"] = S.switch_dv_mps <= S.switch_dv_max_mps
S["ratio_to_direct_10km"] = S.switch_dv_mps / S.direct_step_10km_mps
S.to_csv(f"{OUT}/switchability.csv", index=False)
print("== switchability (best direction/location per geometry; dv linear extrapolation of the FD slope) ==")
print(S[["tag", "start", "direction", "delta_rate_km_day", "best_direction", "best_location", "raw_best_direction", "raw_switch_dv_mps", "kick_dv_mps", "switch_dv_mps", "step_km_of_kick", "extrapolation_factor",
         "ref_rate_km_day", "exit_time_center_days", "exit_time_edge_days", "return_time_days", "direct_step_same_disp_mps",
         "direct_step_10km_mps", "ratio_to_direct_10km", "switchable"]].round(3).to_string())
both = S.groupby(["tag", "start", "delta_rate_km_day"]).switchable.all()
print("== switchable in both directions below switch_dv_max (%.2f x direct 10 km step):" % SW["switch_dv_max_fraction_of_direct_step"], int(both.sum()), "of", len(both))

# ---- plots ---------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
lab = [f"{r.tag}\n{r.start}" for r in S[(S.delta_rate_km_day == SW["target_delta_rate_km_day"][-1]) & (S.direction == "raise")].itertuples()]
x = np.arange(len(lab)); w = 0.2
for j, (rate, direction) in enumerate([(SW["target_delta_rate_km_day"][0], "raise"), (SW["target_delta_rate_km_day"][-1], "raise")]):
    sub = S[(S.delta_rate_km_day == rate) & (S.direction == direction)]
    ax.bar(x + (j - 1.5) * w, sub.switch_dv_mps.values, w, label=f"switch_dv, {rate} km/day (best lever)")
sub = S[(S.delta_rate_km_day == SW["target_delta_rate_km_day"][0]) & (S.direction == "raise")]
ax.bar(x + 0.5 * w, sub.direct_step_same_disp_mps.values, w, color="k", label="direct step, same displacement as the switch buys (half band)")
ax.bar(x + 1.5 * w, sub.direct_step_10km_mps.values, w, color="gray", label=f"direct step {SW['direct_step_km']} km")
ax.set_yscale("log"); ax.set_xticks(x, lab, fontsize=7, rotation=45, ha="right"); ax.set_ylabel("delta-v (m/s)"); ax.legend(fontsize=8)
ax.set_title("Switch delta-v (linear extrapolation of finite differences) vs direct periapsis step")
fig.savefig(f"{OUT}/fig_switchability_bars.png", dpi=130)
