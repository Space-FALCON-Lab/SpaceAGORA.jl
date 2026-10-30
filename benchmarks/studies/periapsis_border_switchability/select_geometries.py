"""Pick the Part C geometries from the Part B grid by a rule fixed before any kick result exists.

Rule (all from the passive Sun+Moon 60-day runs, never from kicks; revised once after the first look at the grid, before kicks,
because the first draft's difference-series rule picked edge extrema instead of V shapes):
  crossing: the 3B periapsis-altitude history has an interior turning point (minimum or maximum) between day CROSS_LO and CROSS_HI
            with both legs larger than LEG_MIN_M (a natural V: the drift changes sign mid-campaign). Ranked by 3B in-band days
            (so the chosen V stays near the band); top 2, at most one per (RAAN, omega) pair.
  static  : among cells that are not V-shaped and whose 3B run completed (retcode Success), the longest total 3B in-band time
            , ties -> smallest max |alt - alt0|; top 2.
  worst   : the largest |3B - control| periapsis change of all remaining completed cells (scale only).
Usage: python3 select_geometries.py [OUT_DIR]
"""
import sys, os
import numpy as np, pandas as pd

import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Pick the Part C geometries")
_ap.add_argument("out_dir", nargs="?", help="directory to read/write (default: the results root)")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir or results_root(_args.results_root)
CROSS_LO, CROSS_HI = 10.0, 50.0  # days, ASSUMED window for "mid-campaign"
LEG_MIN_M = 2000.0  # ASSUMED minimum leg size of a V
g = pd.read_csv(f"{OUT}/grid_summary.csv")
rows = []
for r in g.itertuples():
    p3 = pd.read_csv(f"{OUT}/grid/{r.tag}/passes_3b.csv"); pc = pd.read_csv(f"{OUT}/grid/{r.tag}/passes_ctl.csv")
    t3 = p3.time_s.values / 86400; tc = pc.time_s.values / 86400
    d = p3.alt_m.values - np.interp(t3, tc, pc.alt_m.values)
    a3 = p3.alt_m.values
    v = False
    for fn in (np.argmin, np.argmax):
        i = fn(a3)
        if CROSS_LO <= t3[i] <= CROSS_HI and min(abs(a3[0] - a3[i]), abs(a3[-1] - a3[i])) > LEG_MIN_M:
            v = True; t_ext = t3[i]; leg = min(abs(a3[0] - a3[i]), abs(a3[-1] - a3[i]))
    score = leg if v else 0.0
    alt = p3.alt_m.values
    rows.append(dict(tag=r.tag, done=(r.retcode_3b == 'Success'), first_exit=r.first_exit_days_3b, raan_deg=r.raan_deg, omega_deg=r.omega_deg, epoch_idx=r.epoch_idx, v_score_m=score, t_ext_day=t_ext if v else np.nan, d_end_m=d[-1] if 'd' in dir() else 0, in_band_days=r.in_band_days_3b, max_dev_m=np.abs(alt - alt[0]).max()))
s = pd.DataFrame(rows)
s["eff_m"] = g.third_body_effect_geodetic_m.values
cross = s[s.v_score_m > 0].sort_values("in_band_days", ascending=False).head(2).assign(role="crossing")
rest = s[~s.tag.isin(cross.tag) & (s.v_score_m == 0) & s.done]
static = rest.sort_values(["in_band_days", "max_dev_m"], ascending=[False, True]).head(2).assign(role="static")
rest = rest[~rest.tag.isin(static.tag)]
worst = rest.reindex(rest.eff_m.abs().sort_values(ascending=False).index).head(1).assign(role="worst")
sel = pd.concat([cross, static, worst])
sel.to_csv(f"{OUT}/geometries.csv", index=False)
s.to_csv(f"{OUT}/selection_scores.csv", index=False)
print(f"V-shaped cells: {(s.v_score_m > 0).sum()} of {len(s)}")
print(sel.round(1).to_string())
