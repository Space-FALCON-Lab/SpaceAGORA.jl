"""Follow-up closed-loop analysis: estimator-selected vs quiet-set cells (Arm B), rank predictiveness (Spearman), Arm C (example predictor, no third bodies) vs Arm B.
python3 -I analyze_followup.py [--results-root R] [--boot 10000] [--seed 20261017]
Inputs under <R>/compare_baselines: closedloop_cells.csv, closedloop_cells_followup.csv, closedloop_cells_all.csv, closedloop_B/, closedloop_B_followup/, closedloop_C_all/, closedloop_A/.
Outputs there: followup_cells_table.csv, followup_estimator_vs_quiet.csv, followup_spearman.csv, followup_armC_vs_armB.csv, followup_armC_failures.csv, fig_followup.png."""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_root import add_results_root_argument, results_root
import compare_baselines as cb
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HARD = 0.54


def rd(base, d):
    return pd.read_csv(os.path.join(base, d, "closedloop_rows.csv"))


def boot_ratio(a, b, n, rng):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ma = np.median(a[rng.integers(0, len(a), (n, len(a)))], axis=1); mb = np.median(b[rng.integers(0, len(b), (n, len(b)))], axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ma / mb
    r = r[np.isfinite(r)]
    return np.quantile(r, .025), np.quantile(r, .975)


def boot_diff(a, b, n, rng):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.median(a[rng.integers(0, len(a), (n, len(a)))], axis=1) - np.median(b[rng.integers(0, len(b), (n, len(b)))], axis=1)
    return np.quantile(d, .025), np.quantile(d, .975)


def spearman(x, y):
    return pd.Series(x).corr(pd.Series(y), method="spearman")


def boot_spearman(x, y, n, rng):
    x, y = np.asarray(x, float), np.asarray(y, float); idx = rng.integers(0, len(x), (n, len(x))); out = []
    for ii in idx:
        v = spearman(x[ii], y[ii])
        if np.isfinite(v):
            out.append(v)
    return np.quantile(out, .025), np.quantile(out, .975)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--boot", type=int, default=10000); ap.add_argument("--seed", type=int, default=20261017)
    a = ap.parse_args(); root = results_root(a.results_root); base = os.path.join(root, "compare_baselines"); rng = np.random.default_rng(a.seed)
    orig = pd.read_csv(os.path.join(base, "closedloop_cells.csv")); new = pd.read_csv(os.path.join(base, "closedloop_cells_followup.csv"))
    L, T, _ = cb.load(root); S, order = cb.build_sets(L, T); rank = {t: i + 1 for i, t in enumerate(order.tag.values)}; Ti = T.set_index("tag")
    quiet = set(S["B-quiet"])
    B = pd.concat([rd(base, "closedloop_B"), rd(base, "closedloop_B_followup")], ignore_index=True)
    allc = pd.concat([orig, new[orig.columns]], ignore_index=True)
    C_ = rd(base, "closedloop_C_all"); A_ = rd(base, "closedloop_A")
    info = allc.set_index("tag").assign(est_rank=lambda d: [rank[t] for t in d.index], est_pred_lazy_dv=lambda d: [float(Ti.loc[t].v1_lazy_dv) for t in d.index],
                                         in_quiet=lambda d: [t in quiet for t in d.index])
    cols = ["tag", "completed", "controller_stopped", "controller_failures", "total_dv_mps", "burns", "first_burn_days", "peak_q_w_cm2", "q_min_w_cm2", "exceeds_hard", "passes", "passes_in_band", "min_alt_km", "t_end_days", "retcode", "predictor_evals", "wall_s"]
    T_ = pd.concat([B[cols].assign(arm="B"), C_[cols].assign(arm="C"), A_[cols].assign(arm="A")]).merge(info[["stratum", "est_rank", "est_pred_lazy_dv", "in_quiet", "open_loop_dp_dv", "open_loop_near_passive_dv", "passive_status"]], left_on="tag", right_index=True)
    T_["breach"] = T_.exceeds_hard.astype(bool) | (T_.peak_q_w_cm2 > HARD); T_["all_in_band"] = T_.passes_in_band == T_.passes
    T_.to_csv(os.path.join(base, "followup_cells_table.csv"), index=False)
    b = T_[T_.arm == "B"]
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 300)
    # 1) estimator-selected vs quiet, Arm B
    def sel(strata, dfx=b):
        return dfx[dfx.stratum.isin(strata) & dfx.completed]
    groups = {"estimator_top_not_quiet (new 6)": ["estimator_top_not_quiet"], "quiet_not_top5 (new 6)": ["quiet_not_top5"], "quiet_in_top (new 6)": ["quiet_in_top"]}
    rows = []
    comps = [("estimator_top_not_quiet", ["estimator_top_not_quiet"], "quiet_not_top5", ["quiet_not_top5"]),
             ("estimator_top_not_quiet", ["estimator_top_not_quiet"], "quiet_in_top", ["quiet_in_top"]),
             ("quiet_in_top", ["quiet_in_top"], "quiet_not_top5", ["quiet_not_top5"]),
             ("estimator_top_not_quiet + selected_best 50deg", None, "all quiet (quiet_not_top5 + quiet_in_top + selected_best 116.6/RAAN0)", None)]
    for n1, s1, n2, s2 in comps:
        if s1 is None:
            g1 = b[((b.stratum == "estimator_top_not_quiet") | ((b.stratum == "selected_best") & ~b.in_quiet)) & b.completed]
            g2 = b[(b.stratum.isin(["quiet_not_top5", "quiet_in_top"]) | ((b.stratum == "selected_best") & b.in_quiet)) & b.completed]
        else:
            g1, g2 = sel(s1), sel(s2)
        for met in ("total_dv_mps", "burns"):
            lo, hi = boot_ratio(g1[met], g2[met], a.boot, rng); dl, dh = boot_diff(g1[met], g2[met], a.boot, rng)
            rows.append(dict(group1=n1, group2=n2, metric=met, n1=len(g1), n2=len(g2), median1=g1[met].median(), median2=g2[met].median(), ratio=g1[met].median() / g2[met].median(), ci95_lo=lo, ci95_hi=hi,
                             diff=g1[met].median() - g2[met].median(), diff_lo=dl, diff_hi=dh))
    E1 = pd.DataFrame(rows); E1.to_csv(os.path.join(base, "followup_estimator_vs_quiet.csv"), index=False); print(E1.round(3).to_string())
    print(b.groupby("stratum").agg(n=("tag", "size"), done=("completed", "sum"), dv_med=("total_dv_mps", "median"), dv_min=("total_dv_mps", "min"), dv_max=("total_dv_mps", "max"), burns_med=("burns", "median"), rank_med=("est_rank", "median"), pred_med=("est_pred_lazy_dv", "median")).round(3).to_string())
    # 2) Spearman: predicted vs closed-loop dV over new + old selected cells (Arm B)
    pool = b[b.stratum.isin(["selected_best", "estimator_top_not_quiet", "quiet_not_top5", "quiet_in_top"])]
    srows = []
    for name, g in (("all new + old selected (24)", pool), ("quiet cells only", pool[pool.in_quiet]), ("non-quiet cells only", pool[~pool.in_quiet])):
        for xn in ("est_pred_lazy_dv", "est_rank"):
            lo, hi = boot_spearman(g[xn], g.total_dv_mps, a.boot, rng) if len(g) > 3 else (np.nan, np.nan)
            srows.append(dict(cells=name, n=len(g), predictor=xn, vs="closed-loop dV Arm B", spearman=spearman(g[xn], g.total_dv_mps), ci95_lo=lo, ci95_hi=hi, pred_range=f"{g[xn].min():.3g}-{g[xn].max():.3g}"))
        g2 = g.copy(); lo, hi = boot_spearman(g2.est_pred_lazy_dv, g2.open_loop_dp_dv, a.boot, rng) if len(g2) > 3 else (np.nan, np.nan)
        srows.append(dict(cells=name, n=len(g2), predictor="est_pred_lazy_dv", vs="open-loop DP dV (Part 1)", spearman=spearman(g2.est_pred_lazy_dv, g2.open_loop_dp_dv), ci95_lo=lo, ci95_hi=hi, pred_range=""))
    pool2 = pool.assign(quiet_flag=pool.in_quiet.astype(float)); srows.append(dict(cells="all new + old selected (24)", n=len(pool2), predictor="in_quiet (0/1)", vs="closed-loop dV Arm B", spearman=spearman(pool2.quiet_flag, pool2.total_dv_mps),
                                                                                  ci95_lo=boot_spearman(pool2.quiet_flag, pool2.total_dv_mps, a.boot, rng)[0], ci95_hi=boot_spearman(pool2.quiet_flag, pool2.total_dv_mps, a.boot, rng)[1], pred_range=""))
    Sp = pd.DataFrame(srows); Sp.to_csv(os.path.join(base, "followup_spearman.csv"), index=False); print(Sp.round(3).to_string())
    # 3) Arm C vs B (and A) on the original 30 (+ all 48 for C)
    m = T_[T_.arm == "C"].merge(T_[T_.arm == "B"][["tag", "total_dv_mps", "burns", "completed", "peak_q_w_cm2"]], on="tag", suffixes=("_C", "_B"))
    m["C_over_B_dv"] = m.total_dv_mps_C / m.total_dv_mps_B; m["burns_C_minus_B"] = m.burns_C - m.burns_B
    m.to_csv(os.path.join(base, "followup_armC_vs_armB.csv"), index=False)
    print(m.groupby("stratum").agg(n=("tag", "size"), C_done=("completed_C", "sum"), B_done=("completed_B", "sum"), C_breach=("breach", "sum"), C_allband=("all_in_band", "sum"), dv_C=("total_dv_mps_C", "median"), dv_B=("total_dv_mps_B", "median"),
                                    ratio_med=("C_over_B_dv", "median"), burns_C=("burns_C", "median"), burns_B=("burns_B", "median"), qmax_C=("peak_q_w_cm2_C", "max"), qmax_B=("peak_q_w_cm2_B", "max")).round(3).to_string())
    o30 = m[m.tag.isin(orig.tag)]
    lo, hi = boot_ratio(o30.total_dv_mps_C[o30.completed_C & o30.completed_B], o30.total_dv_mps_B[o30.completed_C & o30.completed_B], a.boot, rng)
    print("orig 30: C dV total", round(float(o30.total_dv_mps_C.sum()), 2), "B", round(float(o30.total_dv_mps_B.sum()), 2), "median C/B", round(float(o30.C_over_B_dv.median()), 3), "ratio of medians CI", round(float(lo), 3), round(float(hi), 3),
          "C breaches", int(o30.breach.sum()), "C not completed", int((~o30.completed_C).sum()), "C passes out of band", int((~o30.all_in_band).sum()), "pairs", len(o30))
    bad = T_[(T_.arm == "C") & ((~T_.completed) | T_.breach | (~T_.all_in_band))]
    bad[["tag", "stratum", "retcode", "t_end_days", "controller_stopped", "controller_failures", "peak_q_w_cm2", "min_alt_km", "passes", "passes_in_band", "total_dv_mps", "burns"]].to_csv(os.path.join(base, "followup_armC_failures.csv"), index=False)
    print(bad[["tag", "stratum", "retcode", "t_end_days", "controller_failures", "peak_q_w_cm2", "min_alt_km", "passes_in_band", "passes"]].to_string())
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    col = {"selected_best": "C0", "estimator_top_not_quiet": "C1", "quiet_not_top5": "C2", "quiet_in_top": "C3"}
    for s_, g in pool.groupby("stratum"):
        ax[0].scatter(g.est_pred_lazy_dv, g.total_dv_mps, label=s_, color=col[s_])
    ax[0].set_xlabel("estimator-v2 predicted lazy_s0 dV (m/s)"); ax[0].set_ylabel("closed-loop dV, Arm B (m/s)"); ax[0].legend(fontsize=8)
    ax[1].scatter(m.total_dv_mps_B, m.total_dv_mps_C, c=["r" if not c else "k" for c in m.completed_C], s=14); mx = max(m.total_dv_mps_B.max(), m.total_dv_mps_C.max()); ax[1].plot([0, mx], [0, mx], "k", lw=.7)
    ax[1].set_xlabel("Arm B dV (m/s)"); ax[1].set_ylabel("Arm C dV (m/s), red = not completed")
    fig.savefig(os.path.join(base, "fig_followup.png"), dpi=130)


if __name__ == "__main__":
    main()
