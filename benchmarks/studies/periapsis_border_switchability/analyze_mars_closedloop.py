"""Stage C3 (Mars): closed-loop results, pre-registered criteria, and the critical-inclination ("quiet set") shortcut test.

Usage: python3 -I analyze_mars_closedloop.py [--results-root DIR]
Reads <root>/mars/{closedloop,closedloop_j3diag}/closedloop_rows.csv and <tag>/passes_cl.csv, closedloop_cells.csv, selection_cell_table.csv, mars_cell_labels.csv.
Writes <root>/mars/{closedloop_summary.csv, closedloop_ratios.csv, closedloop_cells_table.csv, quiet_shortcut.csv, fig_closedloop.png}.
Criteria (mars_generality_plan.md): C3a all runs complete 60 d, no stop/decay/hard breach; C3b every pass of every run in [0.14, 0.32] W/cm2;
C3c median closed-loop dV selected_best / typical <= 0.25 with bootstrap (10,000 draws, seed 20261017) upper limit < 0.5.
Shortcut transfers if (i) median true A1 cost of critical-inclination cells < the all-completed median, AND (ii) B2c top-10% recall >= 0.5,
AND (iii) crit_not_top median closed-loop dV <= 2 x selected_best median.
"""
import argparse, math, pathlib, sys, tomllib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
MC = tomllib.load(open(HERE / "mars_generality_config.toml", "rb"))
QLO, QHI, QHARD = MC["band"]["q_low_w_cm2"], MC["band"]["q_high_w_cm2"], MC["band"]["q_hard_w_cm2"]
CRIT = MC["analysis"]["critical_inclinations_deg"]
STRATA = ["selected_best", "typical", "crit_not_top", "odyssey_ref", "worst_completed"]   # worst_completed: POST HOC addition, frozen predictor only (not in the pre-registered C3a-c)


def boot_ratio(a, b, draws, rng):
    a = np.asarray(a, float); b = np.asarray(b, float)
    r = []
    for _ in range(draws):
        x = rng.choice(a, len(a)); y = rng.choice(b, len(b)); den = np.median(y)
        r.append(np.median(x) / den if den > 0 else np.nan)
    r = np.array(r); r = r[np.isfinite(r)]
    return (float(np.median(a) / np.median(b)) if np.median(b) > 0 else np.nan), (float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))) if len(r) else (np.nan, np.nan)


def load_arm(m, folder):
    R = pd.read_csv(m / folder / "closedloop_rows.csv"); extra = []
    for t in R.tag:
        p = pd.read_csv(m / folder / t / "passes_cl.csv"); tr = pd.read_csv(m / "truth" / t / "passes_3b.csv")
        q = p.q_w_cm2.values
        extra.append(dict(tag=t, n_pass=len(p), n_below_band=int((q < QLO).sum()), n_above_band=int((q > QHI).sum()), q_min=float(q.min()), q_max=float(q.max()),
                          band_edge_shift_km=float((tr.h_qlow_m.iloc[-1] - tr.h_qlow_m.iloc[0]) / 1e3), vrel_first=float(tr.vrel_mps.iloc[0]), vrel_last=float(tr.vrel_mps.iloc[-1]),
                          alt_first_km=float(p.alt_m.iloc[0] / 1e3), alt_last_km=float(p.alt_m.iloc[-1] / 1e3)))
    return R.merge(pd.DataFrame(extra), on="tag")


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    rng = np.random.default_rng(MC["closedloop"]["seed"]); draws = int(MC["closedloop"]["bootstrap_draws"])
    C = pd.read_csv(m / "closedloop_cells.csv").set_index("tag"); T = pd.read_csv(m / "selection_cell_table.csv"); L = pd.read_csv(m / "mars_cell_labels.csv")
    arms = {"frozen": load_arm(m, "closedloop"), "j3diag(UNFROZEN)": load_arm(m, "closedloop_j3diag")}
    worst = load_arm(m, "closedloop_worst") if (m / "closedloop_worst" / "closedloop_rows.csv").exists() else None
    sumrows = []; ratrows = []
    for name, R in arms.items():
        R["open_loop_true_dv"] = R.tag.map(C.true_dv); R["open_loop_pred_dv"] = R.tag.map(C.pred_dv); R["true_class"] = R.tag.map(C.true_class); R["rank_est"] = R.tag.map(C.rank_est)
        R.to_csv(m / f"closedloop_cells_table_{name.split('(')[0]}.csv", index=False)
        c3a = bool(R.completed.all() and (~R.controller_stopped).all() and (~R.exceeds_hard).all() and (R.retcode == "Success").all())
        oob = int((R.n_below_band + R.n_above_band).sum()); runs_oob = int(((R.n_below_band + R.n_above_band) > 0).sum())
        c3b = oob == 0
        for st in STRATA:
            g = R[R.stratum == st]
            if len(g):
                sumrows.append(dict(arm=name, stratum=st, n=len(g), dv_median=g.total_dv_mps.median(), dv_mean=g.total_dv_mps.mean(), dv_min=g.total_dv_mps.min(), dv_max=g.total_dv_mps.max(),
                                    burns_median=g.burns.median(), open_loop_true_dv_median=g.open_loop_true_dv.median(), passes_below_band=int(g.n_below_band.sum()), passes_above_band=int(g.n_above_band.sum()),
                                    runs_with_out_of_band=int(((g.n_below_band + g.n_above_band) > 0).sum()), q_min=g.q_min.min(), q_max=g.q_max.max(), n_true_decayed=int((g.true_class == "decayed").sum())))
        sb = R[R.stratum == "selected_best"].total_dv_mps.values; ty = R[R.stratum == "typical"].total_dv_mps.values
        r, ci = boot_ratio(sb, ty, draws, rng)
        c3c = bool(np.isfinite(r) and r <= 0.25 and ci[1] < 0.5)
        sb2 = R[(R.stratum == "selected_best") & (R.true_class != "decayed")].total_dv_mps.values
        r2, ci2 = boot_ratio(sb2, ty, draws, rng) if len(sb2) else (np.nan, (np.nan, np.nan))
        cn = R[R.stratum == "crit_not_top"].total_dv_mps.values; od = R[R.stratum == "odyssey_ref"].total_dv_mps.values
        rc, cic = boot_ratio(cn, sb, draws, rng) if len(cn) and np.median(sb) > 0 else (np.nan, (np.nan, np.nan))
        ro, cio = boot_ratio(sb, od, draws, rng) if len(od) else (np.nan, (np.nan, np.nan))
        ratrows.append(dict(arm=name, C3a_all_complete_no_stop_no_breach=c3a, C3b_all_passes_in_band=c3b, passes_out_of_band=oob, runs_with_out_of_band=runs_oob, n_runs=len(R),
                            ratio_selected_over_typical=r, ci_lo=ci[0], ci_hi=ci[1], C3c_pass=c3c, ratio_selected_excl_decayed_over_typical=r2, ci2_lo=ci2[0], ci2_hi=ci2[1],
                            ratio_selected_over_odyssey_ref=ro, cio_lo=cio[0], cio_hi=cio[1], ratio_crit_not_top_over_selected=rc, cic_lo=cic[0], cic_hi=cic[1],
                            median_crit_not_top=float(np.median(cn)) if len(cn) else np.nan, median_selected=float(np.median(sb)), median_typical=float(np.median(ty)), median_odyssey_ref=float(np.median(od)) if len(od) else np.nan))
    if worst is not None:   # post hoc stratum, appended after the pre-registered criteria were computed
        cw = pd.read_csv(m / "closedloop_cells_worst.csv").set_index("tag"); worst["open_loop_true_dv"] = worst.tag.map(cw.true_dv); worst["true_class"] = worst.tag.map(cw.true_class)
        fr0 = arms["frozen"]; wv = worst.total_dv_mps.values
        for st in ("selected_best", "typical", "odyssey_ref", "crit_not_top"):
            r_, ci_ = boot_ratio(fr0[fr0.stratum == st].total_dv_mps.values, wv, draws, rng)
            ratrows.append(dict(arm="frozen (POST HOC)", comparison=f"{st} / worst_completed", ratio=r_, ci_lo=ci_[0], ci_hi=ci_[1]))
        sumrows.append(dict(arm="frozen (POST HOC)", stratum="worst_completed", n=len(worst), dv_median=worst.total_dv_mps.median(), dv_mean=worst.total_dv_mps.mean(), dv_min=worst.total_dv_mps.min(), dv_max=worst.total_dv_mps.max(),
                            burns_median=worst.burns.median(), open_loop_true_dv_median=worst.open_loop_true_dv.median(), passes_below_band=int(worst.n_below_band.sum()), passes_above_band=int(worst.n_above_band.sum()),
                            runs_with_out_of_band=int(((worst.n_below_band + worst.n_above_band) > 0).sum()), q_min=worst.q_min.min(), q_max=worst.q_max.max(), n_true_decayed=0))
    S = pd.DataFrame(sumrows); S.to_csv(m / "closedloop_summary.csv", index=False); Rt = pd.DataFrame(ratrows); Rt.to_csv(m / "closedloop_ratios.csv", index=False)
    print(S.round(3).to_string()); print(Rt.T.to_string())
    fr = arms["frozen"]
    print("\nout-of-band runs (frozen):"); print(fr[(fr.n_below_band + fr.n_above_band) > 0][["tag", "stratum", "n_below_band", "n_above_band", "q_min", "q_max", "band_edge_shift_km", "vrel_first", "vrel_last", "true_class"]].round(3).to_string())
    print("band edge shift over the run (km), all flown cells: median %.2f, range %.2f to %.2f; V_rel first/last median %.0f / %.0f m/s" % (fr.band_edge_shift_km.median(), fr.band_edge_shift_km.min(), fr.band_edge_shift_km.max(), fr.vrel_first.median(), fr.vrel_last.median()))
    # --- quiet-set shortcut (open loop, all cells) ---
    comp = T[T.status == "completed"]; med_all = float(comp.true_cost_mps.median()); crit = comp.inc_deg.apply(lambda i: any(abs(i - c) < 0.05 for c in CRIT)).values
    k = max(1, math.ceil(MC["analysis"]["top_fraction"] * len(T)))
    top_true = set(comp.sort_values("true_cost_mps", kind="stable").head(max(1, math.ceil(MC["analysis"]["top_fraction"] * len(T)))).tag)
    crit_all = T.inc_deg.apply(lambda i: any(abs(i - c) < 0.05 for c in CRIT)).values
    isq = {r.tag: bool(c) for r, c in zip(T.itertuples(), crit_all)}
    recall_b2c_det = len([t for t in top_true if isq[t]]) / len(top_true)
    # estimator rank inside the critical-inclination cells (completed)
    from scipy import stats
    cc = comp[crit]; sp_in = stats.spearmanr(np.where(cc.v2_decayed, np.inf, cc.v2_lazy_dv), cc.true_cost_mps)[0]
    byinc = comp.groupby("inc_deg").true_cost_mps.agg(["count", "median", "min", "max", lambda x: int((x == 0).sum())]); byinc.columns = ["n_completed", "median_true_mps", "min", "max", "n_zero"]
    byinc["n_decayed"] = T[T.status == "decayed"].groupby("inc_deg").size(); byinc["n_cells"] = T.groupby("inc_deg").size()
    q = dict(median_true_critical=float(comp[crit].true_cost_mps.median()), median_true_noncritical=float(comp[~crit].true_cost_mps.median()), median_true_all_completed=med_all,
             critical_cells=int(crit_all.sum()), decayed_among_critical=int(((T.status == "decayed").values & crit_all).sum()), decayed_among_noncritical=int(((T.status == "decayed").values & ~crit_all).sum()),
             top10_true_cells=len(top_true), top10_true_that_are_critical=int(sum(isq[t] for t in top_true)), b2c_recall_if_all_critical_picked=recall_b2c_det,
             spearman_est_within_critical=float(sp_in))
    mm = pd.read_csv(m / "selection_metrics.csv"); q["b2c_recall_random_tiebreak_pooled"] = float(mm[(mm.split == "pooled") & (mm.method == "B2c_critical")].recall_top.iloc[0])
    q["est_v2_recall_pooled"] = float(mm[(mm.split == "pooled") & (mm.method == "v2")].recall_top.iloc[0])
    ci_i = q["median_true_critical"] < med_all; ci_ii = q["b2c_recall_random_tiebreak_pooled"] >= 0.5
    cn = arms["frozen"][arms["frozen"].stratum == "crit_not_top"].total_dv_mps.median(); sbm = arms["frozen"][arms["frozen"].stratum == "selected_best"].total_dv_mps.median()
    ci_iii = bool(cn <= 2 * sbm) if np.isfinite(cn) and np.isfinite(sbm) else False
    q.update(shortcut_i_median_cost_lower=bool(ci_i), shortcut_ii_recall_ge_0p5=bool(ci_ii), shortcut_iii_closedloop_crit_not_top_le_2x_selected=ci_iii, shortcut_transfers=bool(ci_i and ci_ii and ci_iii))
    pd.DataFrame([q]).T.to_csv(m / "quiet_shortcut.csv", header=["value"]); byinc.to_csv(m / "quiet_shortcut_by_inclination.csv")
    print("\nquiet shortcut:"); [print(" ", k_, v_) for k_, v_ in q.items()]; print(byinc.round(3))
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2)); xs = np.arange(len(STRATA)); w = 0.35
    for j, (name, R) in enumerate(arms.items()):
        ax[0].bar(xs + j * w, [R[R.stratum == s].total_dv_mps.median() for s in STRATA], w, label=name, color=("#0b6e99", "#2a7f3f")[j])
    ax[0].set_xticks(xs + w / 2); ax[0].set_xticklabels(STRATA, fontsize=8, rotation=15); ax[0].set_ylabel("median closed-loop dV, 60 d (m/s)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2, axis="y")
    ax[1].scatter(fr.open_loop_true_dv, fr.total_dv_mps, c=[("#c4612b" if t == "decayed" else "#0b6e99") for t in fr.true_class], s=25)
    mx = max(fr.open_loop_true_dv.max(), fr.total_dv_mps.max()) * 1.05; ax[1].plot([0, mx], [0, mx], "k", lw=0.8); ax[1].set_xlabel("open-loop A1 label (m/s; decayed cells have none)"); ax[1].set_ylabel("closed-loop dV (m/s)"); ax[1].grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(m / "fig_closedloop.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
