"""Stage E part 2 analysis: closed-loop rows of Arm A (engine predictor + Sun/Moon) and Arm B (estimator v2) against the open-loop Part 1 numbers.
python3 -I analyze_closedloop.py [--results-root R] [--boot 10000] [--seed 20261016]
Inputs under <R>/compare_baselines: closedloop_cells.csv, closedloop_A/closedloop_rows.csv, closedloop_B/closedloop_rows.csv, closedloop_*/<tag>/decisions.csv.
Outputs there: closedloop_cells_table.csv, closedloop_by_stratum.csv, closedloop_selected_vs_baseline.csv, closedloop_vs_openloop.csv, closedloop_failures.csv, fig_closedloop.png.
Cells that did not complete 60 d (decayed or controller stop) have a PARTIAL dV: they are reported in their own columns and excluded from completed-only medians."""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_root import add_results_root_argument, results_root
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HARD = 0.54


def boot_ratio(a, b, n, rng):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    ma = np.median(a[rng.integers(0, len(a), (n, len(a)))], axis=1); mb = np.median(b[rng.integers(0, len(b), (n, len(b)))], axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ma / mb
    r = r[np.isfinite(r)]
    return (np.quantile(r, .025), np.quantile(r, .975)) if len(r) else (np.nan, np.nan)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--boot", type=int, default=10000); ap.add_argument("--seed", type=int, default=20261016)
    a = ap.parse_args(); root = results_root(a.results_root); base = os.path.join(root, "compare_baselines")
    C = pd.read_csv(os.path.join(base, "closedloop_cells.csv"))
    R = pd.concat([pd.read_csv(os.path.join(base, f"closedloop_{x}", "closedloop_rows.csv")).assign(arm=x) for x in "AB"], ignore_index=True)
    R = R.merge(C[["tag", "passive_status", "open_loop_dp_dv", "open_loop_near_passive_dv", "open_loop_near_passive_burns", "inc_deg", "raan_deg", "omega_deg", "epoch"]], on="tag")
    R["breach"] = R.exceeds_hard.astype(bool) | (R.peak_q_w_cm2 > HARD)
    R["status_class"] = np.where(R.completed, "completed_60d", np.where(R.controller_stopped, "controller_stopped", "terminated_other"))
    R["in_band_all_passes"] = R.passes_in_band == R.passes
    R.to_csv(os.path.join(base, "closedloop_cells_table.csv"), index=False)
    rng = np.random.default_rng(a.seed)
    # per stratum x arm
    rows = []
    for (st, arm), g in R.groupby(["stratum", "arm"], sort=False):
        c = g[g.completed]
        rows.append(dict(stratum=st, arm=arm, n=len(g), n_completed_60d=len(c), n_not_completed=len(g) - len(c), n_breach_hard=int(g.breach.sum()),
                         n_all_passes_in_band=int(g.in_band_all_passes.sum()),
                         dv_median_completed=c.total_dv_mps.median(), dv_q1=c.total_dv_mps.quantile(.25), dv_q3=c.total_dv_mps.quantile(.75), dv_min=c.total_dv_mps.min(), dv_max=c.total_dv_mps.max(),
                         burns_median_completed=c.burns.median(), peak_q_median=g.peak_q_w_cm2.median(), peak_q_max=g.peak_q_w_cm2.max(),
                         minutes_in_band_median=g.minutes_in_band.median(), minutes_in_band_frac_median=(g.minutes_in_band / (g.t_end_days * 1440)).median(),
                         dv_partial_median_not_completed=g[~g.completed].total_dv_mps.median(),
                         open_loop_near_passive_dv_median=g.open_loop_near_passive_dv.median(), open_loop_near_passive_burns_median=g.open_loop_near_passive_burns.median()))
    S = pd.DataFrame(rows); S.to_csv(os.path.join(base, "closedloop_by_stratum.csv"), index=False)
    # Selected vs baseline strata, per arm (completed only)
    ratios = []
    for arm in "AB":
        sel = R[(R.arm == arm) & (R.stratum == "selected_best") & R.completed].total_dv_mps.values
        for bs in ("ref50", "median_cost", "worst_completed"):
            b = R[(R.arm == arm) & (R.stratum == bs) & R.completed].total_dv_mps.values
            lo, hi = boot_ratio(sel, b, a.boot, rng)
            ratios.append(dict(arm=arm, selected="selected_best", baseline=bs, n_sel=len(sel), n_base=len(b), median_sel=np.median(sel) if len(sel) else np.nan,
                               median_base=np.median(b) if len(b) else np.nan, ratio=np.median(sel) / np.median(b) if len(b) and len(sel) else np.nan, ci95_lo=lo, ci95_hi=hi))
    # open-loop same-cell ratios
    ol = C[C.passive_status == "completed"].set_index("tag")
    for bs in ("ref50", "median_cost", "worst_completed"):
        s_ = ol.loc[[t for t in C[C.stratum == "selected_best"].tag if t in ol.index]].open_loop_near_passive_dv.values
        b_ = ol.loc[[t for t in C[C.stratum == bs].tag if t in ol.index]].open_loop_near_passive_dv.values
        lo, hi = boot_ratio(s_, b_, a.boot, rng)
        ratios.append(dict(arm="open_loop_near_passive(same cells)", selected="selected_best", baseline=bs, n_sel=len(s_), n_base=len(b_), median_sel=np.median(s_), median_base=np.median(b_),
                           ratio=np.median(s_) / np.median(b_), ci95_lo=lo, ci95_hi=hi))
    pd.DataFrame(ratios).to_csv(os.path.join(base, "closedloop_selected_vs_baseline.csv"), index=False)
    # closed loop vs open loop, per cell
    W = R.pivot(index="tag", columns="arm", values="total_dv_mps").add_prefix("dv_arm").join(R.pivot(index="tag", columns="arm", values="burns").add_prefix("burns_arm"))
    W = W.join(C.set_index("tag")[["stratum", "passive_status", "open_loop_near_passive_dv", "open_loop_near_passive_burns", "open_loop_dp_dv"]])
    W = W.join(R.pivot(index="tag", columns="arm", values="completed").add_prefix("completed_arm"))
    W["armA_over_open_loop"] = W.dv_armA / W.open_loop_near_passive_dv; W["armB_over_open_loop"] = W.dv_armB / W.open_loop_near_passive_dv; W["armB_over_armA"] = W.dv_armB / W.dv_armA
    W.to_csv(os.path.join(base, "closedloop_vs_openloop.csv"))
    # failures: not completed, breach, or not all passes in band
    F = R[(~R.completed) | R.breach | (~R.in_band_all_passes)].copy()
    reasons = []
    for r in F.itertuples():
        d = os.path.join(base, f"closedloop_{r.arm}", r.tag, "decisions.csv")
        last = pd.read_csv(d).iloc[-1] if os.path.isfile(d) and os.path.getsize(d) > 5 else None
        reasons.append("" if last is None else f"last decision t={last.request_time_s/86400:.2f} d reason={last.trigger_reason} feasible={last.feasible}")
    F["last_decision"] = reasons
    F[["arm", "tag", "stratum", "passive_status", "retcode", "t_end_days", "status_class", "controller_failures", "total_dv_mps", "burns", "peak_q_w_cm2", "breach", "passes", "passes_in_band", "min_alt_km", "last_decision"]].to_csv(
        os.path.join(base, "closedloop_failures.csv"), index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 200)
    print(S.round(3).to_string()); print(pd.DataFrame(ratios).round(3).to_string()); print(W.round(3).to_string()); print(F[["arm", "tag", "stratum", "status_class", "controller_failures", "peak_q_w_cm2", "passes_in_band", "passes", "last_decision"]].to_string())
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    strata = list(dict.fromkeys(R.stratum)); w = 0.27
    for k, (lab, col) in enumerate((("open loop near_passive (Part 1)", None), ("closed loop Arm A", "A"), ("closed loop Arm B", "B"))):
        med = [(g.open_loop_near_passive_dv.median() if col is None else g[g.completed].total_dv_mps.median()) for s_ in strata for g in [R[(R.stratum == s_) & (R.arm == "A")] if col is None else R[(R.stratum == s_) & (R.arm == col)]]]
        ax.bar(np.arange(len(strata)) + k * w, med, w, label=lab)
    ax.set_xticks(np.arange(len(strata)) + w); ax.set_xticklabels(strata, rotation=12); ax.set_ylabel("median 60 d dV (m/s); closed loop = completed cells only"); ax.legend()
    fig.savefig(os.path.join(base, "fig_closedloop.png"), dpi=130)


if __name__ == "__main__":
    main()
