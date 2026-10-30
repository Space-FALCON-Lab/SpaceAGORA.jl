"""A1 labels (exact DP, lazy-to-edge s = 0, near_passive; 1 km margin per side, Keplerian lever per pass) for the blind Earth cells, written in the format of
oracle_cell_labels.csv so that analyze_earth_g20.py runs unchanged. Needs <root>/border2/grid/<tag>/passes_3b.csv with the stored band columns (earth_g20_prepare.py --blind).

Usage: python3 -I earth_blind_labels.py --results-root BLIND_ROOT
"""
import argparse, importlib.util, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("analyze_selection_skill", HERE / "analyze_selection_skill.py")
an = importlib.util.module_from_spec(spec); sys.modules["analyze_selection_skill"] = an; spec.loader.exec_module(an)
lib = an.lib; MARGIN = an.MARGIN; DAY = 86400.0


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root))
    sc = pd.read_csv(root / "border2" / "static_cells.csv"); rs = pd.read_csv(root / "border2" / "static_run_summary.csv").set_index("tag"); rows = []
    for c in sc.itertuples():
        df = pd.read_csv(root / "border2" / "grid" / c.tag / "passes_3b.csv"); x = df.alt_m.values; lo = df.h_qhigh_m.values + MARGIN; hi = df.h_qlow_m.values - MARGIN
        L = lib.lever_kepler(df.a_m.values, df.e.values); q = df.q_w_cm2.values; ok = rs.loc[c.tag, "retcode_3b"] == "Success"
        d = lib.dp_schedule(lo - x, hi - x, L, None, tie=an.ORC["dp"]["tie_lambda_mps"]); dp = np.nan if d is None else d["dv"]
        lz, lzb, _ = lib.lazy(x, lo, hi, L, s=0.0); npv, npb, _ = lib.near_passive(x, lo, hi, L)
        rows.append(dict(tag=c.tag, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, epoch=int(c.epoch_idx), family="blind", status="completed" if ok else "decayed",
                         retcode=rs.loc[c.tag, "retcode_3b"], t_end_days=df.time_s.iloc[-1] / DAY, n_passes=len(df), dp_dv=dp, dp_burns=-1 if d is None else d["burns"], lazy_s0_dv=lz, lazy_s0_burns=lzb,
                         near_passive_dv=npv, near_passive_burns=npb, passive_q_max_w_cm2=float(q.max()), passive_exceeds_hard=bool((q > an.Q_HARD).any()),
                         true_cost_label_mps=dp if ok else np.nan, cost_to_termination_mps=np.nan if ok else dp))
    out = root / "oracle"; out.mkdir(exist_ok=True); T = pd.DataFrame(rows); T.to_csv(out / "oracle_cell_labels.csv", index=False)
    print(len(T), "cells;", T.status.value_counts().to_dict(), "; passive breach:", int(T.passive_exceeds_hard.sum()))


if __name__ == "__main__":
    main()
