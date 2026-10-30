"""Earth analysis on the >= 20x20 truth (python3 -I analyze_earth_g20.py --results-root NEW --mode score|truthdiff ...).

  --mode score --pred-dir DIR      the A2 selection-skill machinery (analyze_selection_skill.py, unchanged) on a prediction directory
                                   (DIR/{pred,pred_selfcal,hand,pred_cells.csv}); the attribution-variant tables are skipped (they need the 3x0 variants).
                                   Writes selection_cell_table.csv, selection_metrics.csv, threshold_table.csv, confusion_*.csv, lookahead_error.csv,
                                   dv_error_table.csv, decay_time_error.csv into DIR.
  --mode truthdiff --old-root OLD  per-cell difference between the earlier (3x0) truth under OLD and the new truth, and the change of the A1 labels.
"""
import argparse, importlib.util, pathlib, sys
import numpy as np, pandas as pd
from scipy import stats
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("analyze_selection_skill", HERE / "analyze_selection_skill.py")
SK = importlib.util.module_from_spec(spec); sys.modules["analyze_selection_skill"] = SK; spec.loader.exec_module(SK)
DAY = 86400.0
_splits_grid = SK.splits_of


def _splits_any(T):
    """The A2 folds; RAAN halves become [0, 180) and [180, 360) when the RAAN values are not on the 45-degree grid (blind cells)."""
    out = _splits_grid(T)
    if not T.raan_deg.isin(sum(SK.CFG["splits"]["raan_halves"], [])).all():
        for k in [k for k in out if k.startswith("raan_")]: del out[k]
        out["raan_0-179"] = (T.raan_deg < 180).values; out["raan_180-359"] = (T.raan_deg >= 180).values
    return out


SK.splits_of = _splits_any


def score(root, d):
    SK.ROOT_LABELS = root / "oracle" / "oracle_cell_labels.csv"
    if not (d / "hand").exists():
        (d / "hand").symlink_to(root / "estimator_v2" / "hand")      # B1 hand model: same inputs for every estimator generation
    T = SK.build_cell_table(root, d); T.to_csv(d / "selection_cell_table.csv", index=False)
    M = SK.evaluate(T); M.to_csv(d / "selection_metrics.csv", index=False)
    TH = SK.threshold_table(M); TH.to_csv(d / "threshold_table.csv", index=False)
    for key in ("v1", "v1_selfcal", "B1_hand"):
        SK.confusion(T, key).to_csv(d / f"confusion_{key}.csv")
    D_all = SK.lookahead_errors(T, root, d); D_comp = D_all[D_all.status == "completed"]
    LA = pd.concat([SK.lookahead_table(D_comp, "completed"), SK.lookahead_table(D_all, "all_cells")]); LA.to_csv(d / "lookahead_error.csv", index=False)
    SK.dv_error_table(T).to_csv(d / "dv_error_table.csv", index=False); SK.decay_time_table(T).to_csv(d / "decay_time_error.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    cols = ["split", "method", "n", "n_complete_true", "spearman", "kendall", "k", "recall_top", "random_recall", "regret_complete_mps", "regret_safe_mps", "false_safe_rate", "false_decay_rate", "n_pred_safe"]
    print(M[cols].round(3).to_string()); print(TH.round(3).to_string()); print(LA.round(3).to_string())
    for key in ("v1", "v1_selfcal"):
        print(key); print(SK.confusion(T, key))
    print(SK.dv_error_table(T).round(3).to_string()); print(SK.decay_time_table(T).round(3).to_string())


def truth_hist(root, tag, inc):
    return pd.read_csv(SK.truth_path(tag, inc, root))


def truthdiff(root, old):
    Ln = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv").set_index("tag"); Lo = pd.read_csv(old / "oracle" / "oracle_cell_labels.csv").set_index("tag")
    rows = []
    for tag, c in Ln.iterrows():
        o = truth_hist(old, tag, c.inc_deg); n = truth_hist(root, tag, c.inc_deg); m = min(len(o), len(n)); d = n.alt_m.values[:m] - o.alt_m.values[:m]
        rows.append(dict(tag=tag, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, epoch=int(c.epoch), n_old=len(o), n_new=len(n), status_old=Lo.loc[tag, "status"], status_new=c.status,
                         pass1_alt_diff_m=n.alt_m.iloc[0] - o.alt_m.iloc[0], end_alt_diff_m=d[-1], rms_alt_diff_m=float(np.sqrt(np.mean(d ** 2))), max_abs_alt_diff_m=float(np.abs(d).max()),
                         t_end_old_d=o.time_s.iloc[-1] / DAY, t_end_new_d=n.time_s.iloc[-1] / DAY,
                         label_old=Lo.loc[tag, "true_cost_label_mps"], label_new=c.true_cost_label_mps, cost_to_term_old=Lo.loc[tag, "cost_to_termination_mps"], cost_to_term_new=c.cost_to_termination_mps,
                         lazy_old=Lo.loc[tag, "lazy_s0_dv"], lazy_new=c.lazy_s0_dv, exceeds_old=bool(Lo.loc[tag, "passive_exceeds_hard"]), exceeds_new=bool(c.passive_exceeds_hard)))
    D = pd.DataFrame(rows); D.to_csv(root / "earth_truth_old_vs_new.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("cells", len(D), "status changes:", int((D.status_old != D.status_new).sum()), "breach-flag changes:", int((D.exceeds_old != D.exceeds_new).sum()))
    print(pd.crosstab(D.status_old, D.status_new))
    for col in ("pass1_alt_diff_m", "end_alt_diff_m", "rms_alt_diff_m", "max_abs_alt_diff_m"):
        print(col, "median %.1f  p90 %.1f  max %.1f m" % (D[col].abs().median(), D[col].abs().quantile(0.9), D[col].abs().max()))
    print("by inclination, RMS altitude difference (m): median / max"); print(D.groupby("inc_deg").rms_alt_diff_m.agg(["median", "max"]).round(1))
    print("by inclination, |end-of-history altitude difference| (m): median / max"); print(D.groupby("inc_deg").end_alt_diff_m.agg(lambda x: x.abs().median()).round(1))
    both = D[(D.status_old == "completed") & (D.status_new == "completed")]
    print("completed in both:", len(both), "label Spearman old-vs-new %.4f" % stats.spearmanr(both.label_old, both.label_new)[0], " label difference (m/s) median %.3f max %.3f, relative median %.3f max %.3f"
          % ((both.label_new - both.label_old).abs().median(), (both.label_new - both.label_old).abs().max(), ((both.label_new - both.label_old).abs() / both.label_old.clip(lower=1e-6)).median(), ((both.label_new - both.label_old).abs() / both.label_old.clip(lower=1e-6)).max()))
    dec = D[(D.status_old == "decayed") & (D.status_new == "decayed")]
    print("decayed in both:", len(dec), "end-time difference (d) median %.3f max %.3f" % ((dec.t_end_new_d - dec.t_end_old_d).abs().median(), (dec.t_end_new_d - dec.t_end_old_d).abs().max()))
    return D


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--mode", required=True, choices=["score", "truthdiff"]); ap.add_argument("--pred-dir"); ap.add_argument("--old-root")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root))
    if a.mode == "score":
        score(root, pathlib.Path(a.pred_dir))
    else:
        truthdiff(root, pathlib.Path(a.old_root))


if __name__ == "__main__":
    main()
