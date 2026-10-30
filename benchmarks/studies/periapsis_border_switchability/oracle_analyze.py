"""Summaries, ranking tables and figures from the CSVs written by oracle_run.py.
Usage:  python3 -I oracle_analyze.py --out-dir DIR      (reads and writes DIR; numpy, pandas, matplotlib only)
"""
import argparse, pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7", grey="#7f7f7f")
POL = [("dp", "DP oracle"), ("lazy_s0", "lazy s=0"), ("lazy_s1", "lazy s=1 km"), ("lazy_s2", "lazy s=2 km"), ("lazy_s4", "lazy s=4 km"),
       ("lazy_s6", "lazy s=6 km"), ("near_passive", "near_passive"), ("reset_opposite", "reset opposite"),
       ("hold_nearest", "hold nearest"), ("hold_center", "hold center"), ("hold_far", "hold far (upstream)")]


def q(s, p): return float(s.quantile(p))


def ranking(D, label):
    rows = []
    for k, name in POL:
        dv = D[k + "_dv"]; b = D[k + "_burns"]; r = dv / D.dp_dv
        rows.append(dict(set=label, policy=k, label=name, n=len(D), dv_median=dv.median(), dv_q1=q(dv, .25), dv_q3=q(dv, .75), dv_min=dv.min(), dv_max=dv.max(),
                         ratio_median=r.median(), ratio_q1=q(r, .25), ratio_q3=q(r, .75), ratio_max=r.max(), burns_median=b.median(), burns_min=b.min(), burns_max=b.max()))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    T = pd.read_csv(out / "oracle_cell_labels.csv"); Bc = pd.read_csv(out / "oracle_burn_budget_curve.csv"); Hh = pd.read_csv(out / "oracle_hold_horizon.csv")
    C = pd.read_csv(out / "oracle_convergence.csv"); R = pd.read_csv(out / "oracle_race.csv"); RP = pd.read_csv(out / "oracle_race_paths.csv")
    comp = T[T.status == "completed"].copy(); grid = comp[comp.inc_deg != 50.0]
    rk = pd.concat([ranking(grid, "completed_border2_grid_172"), ranking(comp, "completed_all_203")]); rk.to_csv(out / "oracle_policy_ranking.csv", index=False)
    pd.set_option("display.width", 250)
    print(rk[["set", "policy", "dv_median", "dv_q1", "dv_q3", "dv_min", "dv_max", "ratio_median", "ratio_max", "burns_median", "burns_min", "burns_max"]].round(3).to_string(index=False))

    # scoping reproduction (res 250 m grid DP, own-pass lever) on the 172 border2 cells
    sc = grid
    print("\nSCOPING REPRO (250 m grid DP): n=%d DP med %.3f min %.3f max %.3f; lazy/DP med %.4f max %.4f; near/DP med %.4f max %.4f; reset/DP med %.4f max %.4f; burns lazy %g s4 %g near %g reset %g" % (
        len(sc), sc.dp_dv_res250.median(), sc.dp_dv_res250.min(), sc.dp_dv_res250.max(), (sc.lazy_s0_dv / sc.dp_dv_res250).median(), (sc.lazy_s0_dv / sc.dp_dv_res250).max(),
        (sc.near_passive_dv / sc.dp_dv_res250).median(), (sc.near_passive_dv / sc.dp_dv_res250).max(), (sc.reset_opposite_dv / sc.dp_dv_res250).median(), (sc.reset_opposite_dv / sc.dp_dv_res250).max(),
        sc.lazy_s0_burns.median(), sc.lazy_s4_burns.median(), sc.near_passive_burns.median(), sc.reset_opposite_burns.median()))
    print("exact-vs-250m DP: median ratio 250/exact %.4f, max %.4f" % ((comp.dp_dv_res250 / comp.dp_dv).median(), (comp.dp_dv_res250 / comp.dp_dv).max()))
    print("lazy_s0 below exact DP anywhere:", int((comp.lazy_s0_dv < comp.dp_dv - 1e-9).sum()), " min lazy/DP", (comp.lazy_s0_dv / comp.dp_dv).min())

    # forecast value vs policy value (completed 203), per-cell differences
    d = pd.DataFrame({"tag": comp.tag})
    d["policy_value_near_minus_dp"] = comp.near_passive_dv - comp.dp_dv
    d["policy_value_reset_minus_dp"] = comp.reset_opposite_dv - comp.dp_dv
    d["forecast_value_dv_lazy_minus_dp"] = comp.lazy_s0_dv - comp.dp_dv
    print("\nPOLICY/FORECAST VALUE (m/s, median [q1,q3] max; completed %d)" % len(comp))
    for c in d.columns[1:]:
        print(f"  {c}: {d[c].median():.3f} [{q(d[c], .25):.3f}, {q(d[c], .75):.3f}] max {d[c].max():.3f}; relative to DP median {100 * (d[c] / comp.dp_dv).median():.1f}%")
    print("  burns: lazy_s0 %g, DP(tie-min) %g, hold_nearest(min count) %g, near_passive %g, reset %g, budget-min %g" % (
        comp.lazy_s0_burns.median(), comp.dp_burns.median(), comp.hold_nearest_burns.median(), comp.near_passive_burns.median(), comp.reset_opposite_burns.median(), comp.budget_min_burns.median()))
    print("  matched-burn oracle saving vs policy (1 - dV_oracle(N)/dV_policy):")
    mb = []
    for p in ["lazy_s0", "lazy_s1", "lazy_s2", "lazy_s4", "near_passive", "reset_opposite"]:
        s_ = 1 - comp[f"budget_dv_at_{p}_burns"] / comp[p + "_dv"]
        mb.append(dict(policy=p, burns_median=comp[p + "_burns"].median(), policy_dv_median=comp[p + "_dv"].median(), oracle_dv_at_same_burns_median=comp[f"budget_dv_at_{p}_burns"].median(),
                       saving_median=s_.median(), saving_q1=q(s_, .25), saving_q3=q(s_, .75), saving_max=s_.max()))
    mb = pd.DataFrame(mb); mb.to_csv(out / "oracle_matched_burn_saving.csv", index=False); print(mb.round(3).to_string(index=False))
    print("  min feasible burn count: median %g range %g-%g; dV at min burns / DP median %.4f max %.4f; hold_nearest burns == min burns in %.0f%% of cells" % (
        comp.budget_min_burns.median(), comp.budget_min_burns.min(), comp.budget_min_burns.max(), (comp.budget_dv_at_min_burns / comp.dp_dv).median(), (comp.budget_dv_at_min_burns / comp.dp_dv).max(),
        100 * (comp.hold_nearest_burns == comp.budget_min_burns).mean()))
    print("  hold gaps (days) median of per-cell medians %.2f, max over cells %.2f" % (comp.hold_nearest_gap_days_median.median(), comp.hold_nearest_gap_days_max.max()))
    hh = Hh[Hh.tag.isin(comp.tag)].merge(comp[["tag", "dp_dv"]], on="tag"); hh["ratio"] = hh.dv / hh.dp_dv
    hs = hh.groupby("forecast_passes").agg(dv_median=("dv", "median"), ratio_median=("ratio", "median"), ratio_max=("ratio", "max"), burns_median=("burns", "median"), burns_max=("burns", "max")).reset_index()
    hs.to_csv(out / "oracle_hold_horizon_summary.csv", index=False); print(hs.round(4).to_string(index=False))

    # lever variants
    print("\nLEVER: DP measured-lever / Kepler: completed med %.4f range %.4f-%.4f; decayed med %.4f range %.4f-%.4f; neighbor-mean a,e: med %.5f range %.4f-%.4f" % (
        (comp.dp_meas_dv / comp.dp_dv).median(), (comp.dp_meas_dv / comp.dp_dv).min(), (comp.dp_meas_dv / comp.dp_dv).max(),
        (T[T.status == "decayed"].dp_meas_dv / T[T.status == "decayed"].dp_dv).median(), (T[T.status == "decayed"].dp_meas_dv / T[T.status == "decayed"].dp_dv).min(), (T[T.status == "decayed"].dp_meas_dv / T[T.status == "decayed"].dp_dv).max(),
        (comp.dp_nbrmean_dv / comp.dp_dv).median(), (comp.dp_nbrmean_dv / comp.dp_dv).min(), (comp.dp_nbrmean_dv / comp.dp_dv).max()))

    # min-burn constraint (grid DP vs the grid DP with no constraint at the same res)
    print("\nMIN BURN dV (completed):")
    mbrows = []
    for dm in (0.01, 0.05, 0.1):
        for p, base in (("dp", "dp_dv_grid_dvmin_res"), ("lazy_s0", "lazy_s0_dv"), ("lazy_s2", "lazy_s2_dv"), ("near_passive", "near_passive_dv"), ("reset_opposite", "reset_opposite_dv")):
            col = f"{p}_dvmin{dm:g}_dv"; bcol = f"{p}_dvmin{dm:g}_burns"; bb = {"dp": "dp_burns"}.get(p, p + "_burns")
            r = comp[col] / comp[base]
            mbrows.append(dict(min_burn_dv_mps=dm, policy=p, dv_median=comp[col].median(), ratio_to_unconstrained_median=r.median(), ratio_max=r.max(), burns_median=comp[bcol].median(), burns_unconstrained_median=comp[bb].median()))
    mbr = pd.DataFrame(mbrows); mbr.to_csv(out / "oracle_min_burn.csv", index=False); print(mbr.round(4).to_string(index=False))
    sm = comp.copy()
    # smallest burn (m/s) of the unconstrained schedules is not stored per burn; report the share of lazy burns below each probe through the dvmin burn-count change
    # decayed
    dec = T[T.status == "decayed"]
    print("\nDECAYED: n=%d (by inc: %s); covered days median %.1f range %.1f-%.1f; cost to termination median %.2f range %.2f-%.2f; dV before first passive 0.54 median %.2f (to day median %.1f)" % (
        len(dec), dec.inc_deg.value_counts().sort_index().to_dict(), dec.covered_days.median(), dec.covered_days.min(), dec.covered_days.max(), dec.cost_to_termination_mps.median(),
        dec.cost_to_termination_mps.min(), dec.cost_to_termination_mps.max(), dec.dp_dv_before_hard.median(), dec.dp_dv_before_hard_days.median()))
    print("passive exceeds 0.54: all %d, completed %d, decayed %d; pass-1 outside shrunk band %d; pass-1 above hard %d; passive stays in shrunk band all passes: %d; DP infeasible: %d" % (
        T.passive_exceeds_hard.sum(), comp.passive_exceeds_hard.sum(), dec.passive_exceeds_hard.sum(), (~T.x0_in_shrunk_band).sum(), (T.passive_hard_first_pass == 1).sum(), T.passive_in_shrunk_band_all.sum(), T.dp_dv.isna().sum()))
    print("lever ratio last/first (completed): median %.3f, q3 %.3f, max %.3f; cells with ratio>1.5: %d" % (
        (comp.lever_kepler_last_mps_per_km / comp.lever_kepler_first_mps_per_km).median(), q(comp.lever_kepler_last_mps_per_km / comp.lever_kepler_first_mps_per_km, .75),
        (comp.lever_kepler_last_mps_per_km / comp.lever_kepler_first_mps_per_km).max(), int(((comp.lever_kepler_last_mps_per_km / comp.lever_kepler_first_mps_per_km) > 1.5).sum())))
    print("\nDP geometry spread: completed all %d: min %.3f (%s) max %.3f (%s) ratio %.1f; grid172 ratio %.1f; median %.3f" % (len(comp), comp.dp_dv.min(), comp.loc[comp.dp_dv.idxmin(), "tag"], comp.dp_dv.max(), comp.loc[comp.dp_dv.idxmax(), "tag"],
        comp.dp_dv.max() / comp.dp_dv.min(), grid.dp_dv.max() / grid.dp_dv.min(), comp.dp_dv.median()))
    print("cheapest 8:", [(r.tag, round(r.dp_dv, 3)) for r in comp.nsmallest(8, "dp_dv").itertuples()])
    print("costliest 5:", [(r.tag, round(r.dp_dv, 3)) for r in comp.nlargest(5, "dp_dv").itertuples()])
    g = comp.groupby("inc_deg").dp_dv.agg(["count", "median", "min", "max"]); print(g.round(3).to_string())
    print(comp.groupby(["inc_deg", "raan_deg"]).dp_dv.agg(["count", "median"]).round(2).unstack(0).to_string())
    print(comp.groupby("epoch").dp_dv.agg(["count", "median"]).round(2).to_string()); print(comp.groupby("omega_deg").dp_dv.agg(["count", "median"]).round(2).to_string())
    print("\nCONVERGENCE (exact = ref):"); cv = C.copy()
    for c in cv.columns[2:-1]: cv[c] = cv[c] / cv["dv_exact"] - 1
    cv["lazy_s0_dv"] = cv["lazy_s0_dv"] / C["dv_exact"] - 1; print((100 * cv.set_index("tag").drop(columns="dv_exact")).round(3).to_string())
    cv.to_csv(out / "oracle_convergence_relative.csv", index=False)

    # race
    Rk = R[(R.lever == "kepler")]
    print("\nRACE (kepler lever): per cell/horizon; constL excess % range over {first,median} x {early,late} ties")
    rr = Rk.groupby(["tag", "horizon_days"]).agg(covered=("covered_days", "first"), lever_ratio=("lever_ratio", "first"), dp=("dp_dv", "first"), dp_burns=("dp_burns", "first"), lazy=("lazy_dv", "first"),
                                                  ex_min=("constL_excess_pct", "min"), ex_max=("constL_excess_pct", "max"), dpday=("dp_mean_burn_day", "first"), cday_min=("constL_mean_burn_day", "min"), cday_max=("constL_mean_burn_day", "max")).reset_index()
    rr["lazy_excess_pct"] = 100 * (rr.lazy / rr.dp - 1); rr.to_csv(out / "oracle_race_summary.csv", index=False); print(rr.round(3).to_string(index=False))
    rm = R[(R.lever == "measured")].groupby(["tag", "horizon_days"]).dp_dv.first().reset_index().merge(rr[["tag", "horizon_days", "dp"]], on=["tag", "horizon_days"]); rm["meas_over_kepler"] = rm.dp_dv / rm.dp
    print("measured-lever race DP / Kepler:", rm.round(4).to_dict("records"))

    figures(out, T, comp, Bc, hs, C, RP, rr, rk)


def figures(out, T, comp, Bc, hs, C, RP, rr, rk):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25, "figure.dpi": 150})
    # --- fig 1 distributions
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw=dict(width_ratios=[1.1, 1.1, 0.8]))
    names = [n for k, n in POL if k != "lazy_s6"]; keys = [k for k, n in POL if k != "lazy_s6"]
    col = [OI["black"] if k == "dp" else OI["blue"] if k.startswith("lazy") else OI["verm"] if k in ("near_passive", "reset_opposite") else OI["green"] for k in keys]
    for j, (a_, yfun, ylab) in enumerate(((ax[0], lambda k: comp[k + "_dv"], "dV per 60 d (m/s)"), (ax[1], lambda k: comp[k + "_dv"] / comp.dp_dv, "dV / exact DP dV"))):
        data = [yfun(k).values for k in keys]
        bp = a_.boxplot(data, vert=False, widths=.6, patch_artist=True, showfliers=False, medianprops=dict(color="white", lw=1.5), whis=(0, 100))
        for patch, c in zip(bp["boxes"], col): patch.set_facecolor(c); patch.set_alpha(.85)
        a_.set_yticks(range(1, len(names) + 1)); a_.set_yticklabels(names if j == 0 else [""] * len(names)); a_.invert_yaxis(); a_.set_xlabel(ylab)
        for i, d_ in enumerate(data): a_.text(1.01 if j else 1.0, i + 1, f" {np.median(d_):.2f}", va="center", ha="left", fontsize=7.5, transform=a_.get_yaxis_transform()) if False else None
        a_.set_title(("(a) " if j == 0 else "(b) ") + ylab + f"; whiskers min-max, n={len(comp)} completed", fontsize=9)
    if True: ax[1].set_xscale("log"); ax[1].axvline(1, color="k", lw=.6)
    bd = [comp[k + "_burns"].values for k in keys]
    bp = ax[2].boxplot(bd, vert=False, widths=.6, patch_artist=True, showfliers=False, medianprops=dict(color="white", lw=1.5), whis=(0, 100))
    for patch, c in zip(bp["boxes"], col): patch.set_facecolor(c); patch.set_alpha(.85)
    ax[2].set_yticks(range(1, len(names) + 1)); ax[2].set_yticklabels([""] * len(names)); ax[2].invert_yaxis(); ax[2].set_xlabel("burns per 60 d"); ax[2].set_title("(c) burn count", fontsize=9)
    fig.suptitle("Policy dV and burns on completed truth histories (perfect forecast, superposition, 1 km margin; exact DP)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "oracle_fig1_policy_distributions.png"); plt.close(fig)

    # --- fig 2 frontier
    A = comp.set_index("tag"); Bw = Bc.pivot(index="tag", columns="burns", values="dv_best_at_most").loc[A.index].div(A.dp_dv.values, axis=0)
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    cs = np.arange(1, 61); med = []; nf = []
    for c in cs:
        v = Bw[c].values; f = np.isfinite(v); nf.append(f.sum()); med.append(np.median(v[f]) if f.any() else np.nan)
    med = np.array(med); nf = np.array(nf); full = nf == len(A)
    ax[0].plot(cs[full], med[full], color=OI["black"], lw=2, label="oracle: min dV with at most N burns (all cells feasible)")
    ax[0].plot(cs[~full], med[~full], color=OI["black"], lw=1.2, ls="--", label="oracle, feasible cells only (n shown)")
    for c in (1, 3, 6, 10): ax[0].annotate(f"n={nf[c - 1]}", (c, med[c - 1]), textcoords="offset points", xytext=(4, 6), fontsize=7)
    pts = [("lazy_s0", "lazy s=0", OI["blue"]), ("lazy_s1", "s=1", OI["blue"]), ("lazy_s2", "s=2", OI["blue"]), ("lazy_s4", "s=4", OI["blue"]), ("near_passive", "near_passive", OI["verm"]),
           ("reset_opposite", "reset opposite", OI["verm"]), ("hold_nearest", "hold nearest", OI["green"]), ("hold_center", "hold center", OI["green"]), ("hold_far", "hold far", OI["green"])]
    xs = [comp[k + "_burns"].median() for k, _, _ in pts]; ys = [(comp[k + "_dv"] / comp.dp_dv).median() for k, _, _ in pts]
    ax[0].plot(xs[:4], ys[:4], "-", color=OI["blue"], lw=1)
    for (k, n, c), x_, y_ in zip(pts, xs, ys):
        ax[0].plot(x_, y_, "o", color=c, ms=7); ax[0].annotate(n, (x_, y_), textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax[0].set_xlabel("burns per 60 d (median over cells)"); ax[0].set_ylabel("dV / exact DP dV (median over cells)"); ax[0].set_xlim(0, 50); ax[0].set_ylim(0.98, 1.75)
    ax[0].set_title(f"(a) dV vs burn-count frontier, medians, n={len(A)} completed", fontsize=9); ax[0].legend(fontsize=7.5, loc="upper right")
    keys2 = ["lazy_s0", "lazy_s1", "lazy_s2", "lazy_s4", "near_passive", "reset_opposite"]; lab2 = ["lazy s=0", "s=1", "s=2", "s=4", "near_passive", "reset opp."]
    sv = [100 * (1 - comp[f"budget_dv_at_{k}_burns"] / comp[k + "_dv"]).values for k in keys2]
    bp = ax[1].boxplot(sv, widths=.55, patch_artist=True, showfliers=True, flierprops=dict(ms=3, alpha=.5), medianprops=dict(color="white", lw=1.5))
    for patch, k in zip(bp["boxes"], keys2): patch.set_facecolor(OI["blue"] if k.startswith("lazy") else OI["verm"]); patch.set_alpha(.85)
    ax[1].set_xticks(range(1, 7)); ax[1].set_xticklabels([f"{n}\n({comp[k + '_burns'].median():.0f} burns)" for n, k in zip(lab2, keys2)], fontsize=8)
    ax[1].set_ylabel("dV saved by the oracle at the same burn count (%)"); ax[1].set_title("(b) what a multi-pass forecast saves at a fixed burn count", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "oracle_fig2_frontier.png"); plt.close(fig)

    # --- fig 3 oracle cost map
    incs = [23.4, 50.0, 63.4, 116.6]; fig, axs = plt.subplots(4, 3, figsize=(10.5, 15))
    vmax = float(T[T.status == "completed"].dp_dv.max()); cmap = plt.get_cmap("viridis_r")
    for i, inc in enumerate(incs):
        for e in (1, 2, 3):
            ax_ = axs[i, e - 1]; d = T[(T.inc_deg == inc) & (T.epoch == e)]
            raans = sorted(d.raan_deg.unique()); oms = sorted(d.omega_deg.unique())
            M = np.full((len(raans), len(oms)), np.nan); dec = np.zeros_like(M, dtype=bool)
            for r in d.itertuples():
                M[raans.index(r.raan_deg), oms.index(r.omega_deg)] = r.dp_dv if r.status == "completed" else np.nan; dec[raans.index(r.raan_deg), oms.index(r.omega_deg)] = r.status != "completed"
            im = ax_.imshow(M, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
            for a_ in range(len(raans)):
                for b_ in range(len(oms)):
                    if dec[a_, b_]: ax_.text(b_, a_, "decayed", ha="center", va="center", fontsize=7, color=OI["grey"]); ax_.add_patch(plt.Rectangle((b_ - .5, a_ - .5), 1, 1, fill=True, color="#e8e8e8", zorder=0))
                    else: ax_.text(b_, a_, f"{M[a_, b_]:.1f}", ha="center", va="center", fontsize=8, color="white" if M[a_, b_] < .55 * vmax else "black")
            ax_.set_xticks(range(len(oms))); ax_.set_xticklabels([f"{int(o)}" for o in oms]); ax_.set_yticks(range(len(raans))); ax_.set_yticklabels([f"{int(r)}" for r in raans])
            ax_.set_title(f"i = {inc:g} deg, epoch {e}", fontsize=9); ax_.grid(False)
            if e == 1: ax_.set_ylabel("RAAN (deg)")
            if i == 3 or inc == 50.0 and False: ax_.set_xlabel("omega (deg)")
            ax_.set_xlabel("omega (deg)")
    fig.suptitle("Exact oracle-DP dV (m/s per 60 d) by geometry; gray = decayed or impacted before 60 d (no label)", fontsize=10)
    fig.subplots_adjust(hspace=.5, top=.94); fig.colorbar(im, ax=axs, shrink=.5, label="oracle dV (m/s)", pad=.02); fig.savefig(out / "oracle_fig3_cost_map.png", bbox_inches="tight"); plt.close(fig)

    # --- fig 4 race timing
    tags = list(rr.tag.unique()); fig, axs = plt.subplots(2, 4, figsize=(16, 7.5)); axs = axs.ravel()
    for ax_, tag in zip(axs, tags):
        d = RP[RP.tag == tag]; t = d.time_days.values
        ax_.plot(t, d.dv_dp.cumsum(), color=OI["black"], lw=2, label="DP, lever(a(t))"); ax_.plot(t, d.dv_const.cumsum(), color=OI["verm"], lw=1.3, ls="--", label="DP, constant lever (first pass; early-tie)")
        ax_.plot(t, d.dv_lazy.cumsum(), color=OI["blue"], lw=1.1, label="lazy s=0")
        ax_.set_title(tag, fontsize=9); ax_.set_xlabel("time (d)"); ax_.set_ylabel("cumulative dV (m/s)")
        a2 = ax_.twinx(); a2.plot(t, d.lever_mps_per_km, color=OI["green"], lw=1, alpha=.8); a2.set_ylabel("lever (m/s per km)", color=OI["green"]); a2.grid(False); a2.spines["right"].set_visible(True)
    axs[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("Race histories (180 d or to decay): cumulative oracle dV under the changing lever vs a constant lever; green = lever", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "oracle_fig4_race_timing.png"); plt.close(fig)

    # --- fig 5 forecast reach and convergence
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    h = hs[hs.forecast_passes > 0]
    ax[0].plot(h.forecast_passes, h.burns_median, "o-", color=OI["green"]); ax[0].axhline(hs[hs.forecast_passes == 0].burns_median.iloc[0], color=OI["green"], ls=":", label="unlimited forecast")
    ax[0].set_xscale("log"); ax[0].set_xlabel("forecast reach (passes, about 0.8 d each)"); ax[0].set_ylabel("burns per 60 d (median)"); ax[0].legend(fontsize=8)
    ax[0].set_title("(a) hold-nearest: burns vs forecast reach (dV ratio to DP stays %.4f-%.4f)" % (h.ratio_median.min(), h.ratio_median.max()), fontsize=8.5)
    resc = [c for c in C.columns if c.startswith("dv_res")]; rs = [float(c[6:]) for c in resc]
    for _, r in C.iterrows(): ax[1].plot(rs, [100 * (r[c] / r["dv_exact"] - 1) for c in resc], "-o", ms=3, lw=1, alpha=.8)
    ax[1].set_xscale("log"); ax[1].set_yscale("symlog", linthresh=.01); ax[1].set_xlabel("offset grid spacing (m)"); ax[1].set_ylabel("grid DP over exact DP (%)"); ax[1].set_title("(b) grid-DP convergence on 12 cells (exact DP = 0)", fontsize=9)
    ax[2].hist((comp.lazy_s0_dv / comp.dp_dv - 1) * 100, bins=40, color=OI["blue"]); ax[2].set_xlabel("lazy s=0 over exact DP (%)"); ax[2].set_ylabel("cells"); ax[2].set_yscale("log"); ax[2].set_title("(c) lazy-to-edge vs oracle, all completed cells", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "oracle_fig5_reach_convergence.png"); plt.close(fig)


if __name__ == "__main__":
    main()
