"""Part 2: labels, metrics and figures for the blind verification of the estimator-v1 picks.

Usage (numpy, pandas, scipy, matplotlib):
  python3 -I analyze_blind.py [--results-root DIR]

Reads <root>/blind/{blind_cells.csv, blind_run_summary.csv, blind_mc.toml, grid/<tag>/passes_3b.csv} (select_blind.py and
run_blind.jl), <root>/oracle/oracle_cell_labels.csv (A1 sweep labels) and <root>/selection/selection_cell_table.csv (A2 retrospective table).
Truth labels use the A1 pipeline (oracle_lib via analyze_selection_skill.policy_costs): lazy-to-edge s = 0, near_passive, exact DP,
1 km margin per side, Keplerian lever per pass. Writes the tables and PNG figures under <root>/blind.
"""
import argparse, importlib.util, pathlib, sys, tomllib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("ass", HERE / "analyze_selection_skill.py")
ass = importlib.util.module_from_spec(spec); sys.modules["ass"] = ass; spec.loader.exec_module(ass)
lib = ass.lib
CFG = tomllib.load(open(HERE / "blind_config.toml", "rb"))
Q_HARD = ass.Q_HARD
MARGIN = ass.MARGIN
DAY = 86400.0
GROUP_COLOR = {"a_top5": "#0072B2", "b_random_safe": "#E69F00", "c_unsafe": "#999999"}   # Okabe-Ito
GROUP_LABEL = {"a_top5": "(a) predicted cheapest", "b_random_safe": "(b) random predicted-safe", "c_unsafe": "(c) predicted decay/breach"}


def truth_row(root, c, summ):
    d = pd.read_csv(root / "blind" / "grid" / c.tag / "passes_3b.csv")
    x = d.alt_m.values; lo = d.h_qhigh_m.values + MARGIN; hi = d.h_qlow_m.values - MARGIN; L = lib.lever_kepler(d.a_m.values, d.e.values)
    pol = ass.policy_costs(x, d.h_qhigh_m.values, d.h_qlow_m.values, d.a_m.values, d.e.values)
    cost, nb, steps = lib.lazy(x, lo, hi, L, s=0.0)
    u = np.cumsum(steps); pos = x + u
    held = bool(np.all(pos >= lo - 1e-6) and np.all(pos <= hi + 1e-6))
    ok = summ.loc[c.tag, "retcode_3b"] == "Success"
    qmax = float(d.q_w_cm2.max())
    return dict(tag=c.tag, status="completed" if ok else "decayed", t_end_days=float(d.time_s.iloc[-1] / DAY), n_passes=len(d),
                true_dp_mps=pol["dp_dv"], true_lazy_s0_mps=pol["lazy_dv"], true_lazy_burns=pol["lazy_burns"], true_near_mps=pol["near_dv"],
                true_near_burns=pol["near_burns"], true_q_max=qmax, true_exceeds_hard=qmax > Q_HARD,
                true_x0_in_band=bool(lo[0] <= x[0] <= hi[0]), true_alt_min_km=float(x.min() / 1e3), lazy_holds_band=held,
                true_class="decayed" if not ok else ("breach" if qmax > Q_HARD else "safe"))


def build(root):
    B = pd.read_csv(root / "blind" / "blind_cells.csv")
    summ = pd.read_csv(root / "blind" / "blind_run_summary.csv").set_index("tag")
    T = pd.DataFrame([truth_row(root, c, summ) for c in B.itertuples()])
    M = B.merge(T, on="tag")
    M["true_cost_mps"] = np.where(M.status == "completed", M.true_dp_mps, np.nan)   # A1 label: exact DP on completed cells only
    M["pred_unsafe"] = M.pred_class != "safe"
    M["true_unsafe"] = (M.status == "decayed") | M.true_exceeds_hard
    return M


def rank_stats(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return dict(n=int(ok.sum()), spearman=np.nan, kendall=np.nan)
    return dict(n=int(ok.sum()), spearman=float(stats.spearmanr(a[ok], b[ok])[0]), kendall=float(stats.kendalltau(a[ok], b[ok])[0]))


def selection_metrics(M, k=5):
    """Spearman/Kendall (both completed), top-k recall and regret, pooled and per epoch."""
    rows = []
    for name, g in [("all 45", M)] + [(f"epoch {e}", M[M.epoch == e]) for e in (1, 2, 3)]:
        pred = np.where(g.pred_decayed.values, np.inf, g.pred_dv_mps.values)
        true = g.true_cost_mps.values
        r = rank_stats(np.where(np.isfinite(pred), pred, np.nan), true)
        # score for ranking with ties/decay: predicted decayed ranked last; ties by larger predicted band margin
        order_pred = np.lexsort((-g.pred_min_band_margin_m.values, pred))
        top_pred = set(order_pred[:k])
        ok = np.isfinite(true)
        order_true = np.argsort(np.where(ok, true, np.inf), kind="stable")
        top_true = set(order_true[:min(k, ok.sum())])
        recall = len(top_pred & top_true) / max(1, len(top_true))
        pick = order_pred[0]
        best = np.nanmin(true)
        regret = (true[pick] - best) if ok[pick] else np.inf
        # regret over the best predicted-safe pick only
        rows.append(dict(set=name, n=len(g), n_true_complete=int(ok.sum()), n_pred_complete_and_true_complete=r["n"], spearman=r["spearman"],
                         kendall=r["kendall"], k=k, top_k_recall=recall, random_recall=k / len(g), pick_tag=g.tag.values[pick],
                         pick_pred_dv=pred[pick], pick_true_cost=true[pick], best_true_tag=g.tag.values[int(np.nanargmin(true))], best_true_cost=best,
                         regret_mps=regret, pick_true_status=g.status.values[pick]))
    return pd.DataFrame(rows)


def confusion(M):
    labs = ["safe", "breach", "decayed"]
    C = pd.crosstab(M.true_class, M.pred_class).reindex(index=labs, columns=labs, fill_value=0)
    return C


def false_safe_defs(M):
    """(i) literal, (ii) controlled-outcome, plus counts. Predicted safe = pred_class safe."""
    ps = M.pred_class == "safe"
    lit = ps & M.true_unsafe
    ctrl = ps & ~((M.status == "completed") & M.lazy_holds_band)
    return dict(n_pred_safe=int(ps.sum()), literal=int(lit.sum()), controlled=int(ctrl.sum()),
                literal_rate=float(lit.sum() / max(ps.sum(), 1)), controlled_rate=float(ctrl.sum() / max(ps.sum(), 1)),
                literal_tags=";".join(M.tag[lit]), controlled_tags=";".join(M.tag[ctrl]))


def margin_table(pred_q, pred_dec, true_unsafe, true_dec, margins):
    """Predicted safe at margin m: completes and Q_pred (1+m) <= 0.54. False-safe: of those, truth decayed or passive Q > 0.54.
    Miss: truly safe cells (completed and passive Q <= 0.54) declared unsafe at m."""
    rows = []
    true_safe = ~true_unsafe
    for m in margins:
        ps = (~pred_dec) & (pred_q * (1 + m) <= Q_HARD)
        fs = ps & true_unsafe
        miss = (~ps) & true_safe
        rows.append(dict(margin=m, n=len(ps), n_pred_safe=int(ps.sum()), false_safe=int(fs.sum()), false_safe_rate=float(fs.sum() / max(ps.sum(), 1)),
                         false_safe_of_true_unsafe=float(fs.sum() / max(true_unsafe.sum(), 1)),
                         n_true_safe=int(true_safe.sum()), n_missed=int(miss.sum()), miss_rate=float(miss.sum() / max(true_safe.sum(), 1)),
                         false_safe_decayed_only=int((ps & true_dec).sum())))
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(); add_results_root_argument(p); a = p.parse_args()
    root = pathlib.Path(results_root(a.results_root)); out = root / "blind"
    M = build(root); M.to_csv(out / "blind_cell_table.csv", index=False)
    labels = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv")
    comp = labels[labels.status == "completed"].true_cost_label_mps
    sweep = dict(n_completed=int(len(comp)), median=float(comp.median()), q1=float(comp.quantile(.25)), q3=float(comp.quantile(.75)), best=float(comp.min()), worst=float(comp.max()))
    print("sweep (203 completed, A1 DP label):", sweep)
    pd.DataFrame([sweep]).to_csv(out / "sweep_reference.csv", index=False)

    # ---- counts and status
    print(pd.crosstab(M.group, M.status)); print(pd.crosstab(M.group, M.true_class))
    # ---- picks table
    P = M[M.group == "a_top5"][["epoch", "rank_in_group", "tag", "pred_dv_mps", "pred_burns", "pred_q_max", "status", "true_cost_mps", "true_lazy_s0_mps", "true_lazy_burns", "true_q_max", "true_class", "t_end_days"]].copy()
    allc = M[M.status == "completed"].true_cost_mps
    P["true_rank_in_blind_completed"] = [int((allc < v).sum() + 1) if np.isfinite(v) else np.nan for v in P.true_cost_mps]
    P["n_blind_completed"] = len(allc)
    P["true_over_sweep_median"] = P.true_cost_mps / sweep["median"]
    P["true_minus_sweep_best_mps"] = P.true_cost_mps - sweep["best"]
    P["true_pctile_in_sweep"] = [100 * float((comp < v).mean()) if np.isfinite(v) else np.nan for v in P.true_cost_mps]
    P.to_csv(out / "blind_top5_picks.csv", index=False); print(P.round(3).to_string())
    # ---- predicted vs true cost, all groups where both exist
    both = M[(~M.pred_decayed) & (M.status == "completed")]
    cost_err = pd.DataFrame(dict(tag=both.tag, group=both.group, pred=both.pred_dv_mps, true=both.true_lazy_s0_mps, true_dp=both.true_cost_mps))
    cost_err["abs_err"] = cost_err.pred - cost_err.true; cost_err["ratio"] = cost_err.pred / cost_err.true
    cost_err.to_csv(out / "blind_cost_error.csv", index=False)
    print("pred vs true lazy s=0: n", len(cost_err), "median ratio", cost_err.ratio.replace([np.inf], np.nan).median(), "RMS abs err", float(np.sqrt((cost_err.abs_err ** 2).mean())),
          "median abs err", float(cost_err.abs_err.abs().median()), "max abs err", float(cost_err.abs_err.abs().max()))
    by_group = pd.DataFrame([dict(group=gk, n=len(g), rms_abs_err=float(np.sqrt((g.abs_err ** 2).mean())), median_abs_err=float(g.abs_err.abs().median()),
                                  max_abs_err=float(g.abs_err.abs().max()), median_pred=float(g.pred.median()), median_true=float(g.true.median())) for gk, g in cost_err.groupby("group")]).set_index("group")
    print(by_group.round(3)); by_group.to_csv(out / "blind_cost_error_by_group.csv")
    # ---- selection metrics
    SM = selection_metrics(M); SM.to_csv(out / "blind_selection_metrics.csv", index=False); print(SM.round(3).to_string())
    # ---- classification
    C = confusion(M); C.to_csv(out / "blind_confusion.csv"); print(C)
    acc3 = float(np.trace(C.values) / C.values.sum())
    bin_pred_unsafe = M.pred_unsafe.values; bin_true_unsafe = M.true_unsafe.values
    cb = pd.crosstab(pd.Series(np.where(bin_true_unsafe, "unsafe", "safe"), name="truth"), pd.Series(np.where(bin_pred_unsafe, "unsafe", "safe"), name="pred"))
    print("3-class accuracy", acc3); print(cb)
    dec = M[(M.status == "decayed") & M.pred_decayed]
    dt = pd.DataFrame(dict(tag=dec.tag, pred=dec.pred_t_decay_days, true=dec.t_end_days)); dt["err"] = dt.pred - dt.true
    print("decay time error (pred - true, d): n", len(dt), "mean", dt.err.mean(), "RMS", float(np.sqrt((dt.err ** 2).mean())), "range", dt.err.min(), dt.err.max())
    dt.to_csv(out / "blind_decay_time_error.csv", index=False)
    # ---- false-safe definitions
    fs = false_safe_defs(M); print(fs)
    margins = CFG["evaluation"]["q_margins"]
    MT = margin_table(M.pred_q_max.values, M.pred_decayed.values, M.true_unsafe.values, (M.status == "decayed").values, margins)
    fine = np.round(np.arange(0.0, 0.3001, 0.01), 3)
    MT_fine = margin_table(M.pred_q_max.values, M.pred_decayed.values, M.true_unsafe.values, (M.status == "decayed").values, fine)
    # retrospective A2 set (336 cells; selfcal estimator, the scan's input) and pooled
    R = pd.read_csv(root / "selection" / "selection_cell_table.csv")
    Ru = ((R.status == "decayed") | R.true_exceeds_hard).values
    MR = margin_table(R.v1_selfcal_q_max.values, R.v1_selfcal_decayed.values.astype(bool), Ru, (R.status == "decayed").values, margins)
    MR_fine = margin_table(R.v1_selfcal_q_max.values, R.v1_selfcal_decayed.values.astype(bool), Ru, (R.status == "decayed").values, fine)
    pq = np.concatenate([M.pred_q_max.values, R.v1_selfcal_q_max.values]); pd_ = np.concatenate([M.pred_decayed.values, R.v1_selfcal_decayed.values.astype(bool)])
    tu = np.concatenate([M.true_unsafe.values, Ru]); td = np.concatenate([(M.status == "decayed").values, (R.status == "decayed").values])
    MP = margin_table(pq, pd_, tu, td, margins); MP_fine = margin_table(pq, pd_, tu, td, fine)
    for nm, t in (("blind45", MT), ("retro336_selfcal", MR), ("pooled381", MP)):
        t.insert(0, "set", nm)
    pd.concat([MT, MR, MP]).to_csv(out / "false_safe_vs_margin.csv", index=False)
    pd.concat([MT_fine.assign(set="blind45"), MR_fine.assign(set="retro336_selfcal"), MP_fine.assign(set="pooled381")]).to_csv(out / "false_safe_vs_margin_fine.csv", index=False)
    print(pd.concat([MT, MR, MP]).round(4).to_string())
    for nm, t in (("blind45", MT_fine), ("retro336", MR_fine), ("pooled381", MP_fine)):
        z = t[t.false_safe == 0]
        print(nm, "smallest margin on the 0.01 grid with zero false-safe:", None if z.empty else float(z.margin.iloc[0]), "miss rate there:", None if z.empty else float(z.miss_rate.iloc[0]))
    # false-safe with definitions per group, and the same for the retrospective set (literal and controlled; controlled uses A1 lazy_holds on truth)
    for g in ("a_top5", "b_random_safe", "c_unsafe"):
        d = false_safe_defs(M[M.group == g]); print(g, {k: d[k] for k in ("n_pred_safe", "literal", "controlled")})
    rs_ = R.v1_selfcal_class == "safe"
    print("retro336 selfcal: pred safe", int(rs_.sum()), "literal false-safe", int((rs_ & Ru).sum()), "controlled (pred safe & truth decayed)", int((rs_ & (R.status == "decayed")).sum()))
    pd.DataFrame([dict(set="blind45", **{k: v for k, v in fs.items()}), dict(set="retro336_selfcal", n_pred_safe=int(rs_.sum()), literal=int((rs_ & Ru).sum()),
                 controlled=int((rs_ & (R.status == "decayed")).sum()), literal_rate=float((rs_ & Ru).sum() / rs_.sum()), controlled_rate=float((rs_ & (R.status == "decayed")).sum() / rs_.sum()))]
                ).to_csv(out / "false_safe_definitions.csv", index=False)
    # Q prediction error on predicted-complete, truly-completed cells
    q = M[(~M.pred_decayed) & (M.status == "completed")]
    qr = q.pred_q_max / q.true_q_max
    print("Q_pred/Q_true (both complete): n", len(q), "median", qr.median(), "min", qr.min(), "max", qr.max())

    figures(M, P, MT, MR, MP, MT_fine, MR_fine, MP_fine, C, sweep, out)
    M.to_csv(out / "blind_cell_table.csv", index=False)


def figures(M, P, MT, MR, MP, MTf, MRf, MPf, C, sweep, out):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    # ---- fig 1: predicted vs true cost (and peak Q)
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
    both = M[(~M.pred_decayed) & (M.status == "completed")]
    for g, d in both.groupby("group"):
        ax[0].scatter(d.pred_dv_mps + 1e-3, d.true_cost_mps + 1e-3, c=GROUP_COLOR[g], label=GROUP_LABEL[g], s=36, edgecolor="k", linewidth=.4)
    lim = [1e-3, 40]
    ax[0].plot(lim, lim, "k--", lw=.8, label="predicted = true"); ax[0].axhline(sweep["median"], color="#555", lw=.7, ls=":"); ax[0].axhline(sweep["best"], color="#555", lw=.7, ls=":")
    ax[0].text(1.3e-3, sweep["median"] * 1.07, "sweep median", fontsize=8, color="#555"); ax[0].text(1.3e-3, sweep["best"] * 1.1, "sweep best", fontsize=8, color="#555")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlim(lim); ax[0].set_ylim(1e-3, 40)
    ax[0].set_xlabel("predicted dV, lazy s=0 (m/s; +0.001 for the log axis)"); ax[0].set_ylabel("true dV, exact DP on the truth history (m/s; +0.001)")
    ax[0].set_title("Cost: predicted-complete and truly completed cells"); ax[0].legend(fontsize=8, loc="lower right", frameon=False)
    nbad = M[(M.pred_class == "safe") & (M.status == "decayed")]
    if len(nbad):
        ax[0].text(0.98, 0.04, f"{len(nbad)} predicted-safe cell(s) truly decayed (not plotted)", transform=ax[0].transAxes, ha="right", fontsize=8, color="#b00")
    ax[1].scatter(both.pred_q_max, both.true_q_max, c=[GROUP_COLOR[g] for g in both.group], s=36, edgecolor="k", linewidth=.4)
    ax[1].plot([0.15, 1.5], [0.15, 1.5], "k--", lw=.8); ax[1].axvline(Q_HARD, color="#b00", lw=.7); ax[1].axhline(Q_HARD, color="#b00", lw=.7)
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("predicted passive peak Q (W/cm2)"); ax[1].set_ylabel("true passive peak Q (W/cm2)")
    ax[1].set_title("Peak Q (hard limit 0.54 in red)")
    fig.tight_layout(); fig.savefig(out / "fig_blind_pred_vs_true_cost.png", dpi=130); plt.close(fig)
    # ---- fig 2: true dV of the top-5 picks against the sweep distribution
    fig, ax = plt.subplots(figsize=(9, 4.4))
    labs = labels_sorted = P.assign(name=P.tag.str.replace("^I", "", regex=True))
    xs = np.arange(len(P))
    vals = P.true_cost_mps.values
    ax.bar(xs, np.where(np.isfinite(vals), vals, 0), color=[GROUP_COLOR["a_top5"] if np.isfinite(v) else "#b00" for v in vals])
    for i, (v, st) in enumerate(zip(vals, P.status)):
        if not np.isfinite(v):
            ax.text(i, 0.1, "decayed", rotation=90, ha="center", va="bottom", color="#b00", fontsize=8)
    ax.axhline(sweep["median"], color="#555", ls=":", lw=1); ax.text(len(P) - .4, sweep["median"] * 1.02, f"sweep median {sweep['median']:.2f}", ha="right", fontsize=8)
    ax.axhline(sweep["best"], color="#555", ls="--", lw=1); ax.text(len(P) - .4, sweep["best"] + .15, f"sweep best {sweep['best']:.2f}", ha="right", fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels(P.tag.str[1:], rotation=75, ha="right", fontsize=8); ax.set_ylabel("true dV (m/s)")
    ax.set_title("True maintenance cost of the 15 predicted-cheapest picks (blind)")
    fig.tight_layout(); fig.savefig(out / "fig_blind_top5_true_cost.png", dpi=130); plt.close(fig)
    # ---- fig 3: confusion matrix
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    im = ax.imshow(C.values, cmap="Blues")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, int(C.values[i, j]), ha="center", va="center", color="w" if C.values[i, j] > C.values.max() / 2 else "k", fontsize=13)
    ax.set_xticks(range(3)); ax.set_xticklabels(C.columns); ax.set_yticks(range(3)); ax.set_yticklabels(C.index)
    ax.set_xlabel("predicted"); ax.set_ylabel("truth"); ax.set_title("Blind set (45): outcome classes")
    fig.tight_layout(); fig.savefig(out / "fig_blind_confusion.png", dpi=130); plt.close(fig)
    # ---- fig 4: false-safe vs margin
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for t, nm, col, lw, ls in ((MPf, "pooled (381)", "#222222", 4.0, "-"), (MRf, "retrospective (336)", "#D55E00", 2.0, "-"), (MTf, "blind (45)", "#0072B2", 2.0, "--")):
        ax[0].plot(100 * t.margin, t.false_safe, label=nm, color=col, lw=lw, ls=ls)
        ax[1].plot(100 * t.margin, 100 * t.miss_rate, label=nm, color=col, lw=lw, ls=ls)
    for m in CFG["evaluation"]["q_margins"]:
        for a_ in ax:
            a_.axvline(100 * m, color="#ccc", lw=.6, zorder=0)
    ax[0].set_xlabel("margin on predicted peak Q (%)"); ax[0].set_ylabel("false-safe cells (count)"); ax[0].set_title("False-safe, literal definition (i)"); ax[0].legend(frameon=False, fontsize=9)
    ax[1].set_xlabel("margin on predicted peak Q (%)"); ax[1].set_ylabel("truly safe cells rejected (%)"); ax[1].set_title("Miss rate (opportunity lost)")
    fig.tight_layout(); fig.savefig(out / "fig_blind_false_safe_vs_margin.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
