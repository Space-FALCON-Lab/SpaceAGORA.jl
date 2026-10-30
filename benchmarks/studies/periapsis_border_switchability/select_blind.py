"""Part 2: choose the 45 blind truth cells from the estimator-v1 scan by the rule in blind_config.toml.

Usage: python3 -I select_blind.py [--results-root DIR]
Reads <root>/selection/scan_candidates.csv; writes <root>/blind/blind_cells.csv. numpy and pandas only.
"""
import argparse, pathlib, sys, tomllib
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

p = argparse.ArgumentParser(); add_results_root_argument(p); a = p.parse_args()
root = pathlib.Path(results_root(a.results_root))
cfg = tomllib.load(open(HERE / "blind_config.toml", "rb"))["rule"]
S = pd.read_csv(root / "selection" / "scan_candidates.csv")
key = ["inc_deg", "raan_deg", "omega_deg"]
S = S[S.status.isna()].sort_values(key + ["epoch"]).reset_index(drop=True)    # no truth cell in the 336
rng_s = np.random.default_rng(cfg["seed_random_safe"]); rng_u = np.random.default_rng(cfg["seed_random_unsafe"])
rows = []
for e in (1, 2, 3):
    g = S[S.epoch == e]
    safe = g[g.pred_class == "safe"].sort_values(["lazy_dv", "min_band_margin_m", *key], ascending=[True, False, True, True, True])
    top = safe.head(cfg["n_top"])
    rest = safe.iloc[cfg["n_top"]:].sort_values(key)
    rs = rest.iloc[np.sort(rng_s.choice(len(rest), cfg["n_random_safe"], replace=False))]
    dec = g[g.pred_class == "decayed"].sort_values(key); br = g[g.pred_class == "breach"].sort_values(key)
    ru = pd.concat([dec.iloc[np.sort(rng_u.choice(len(dec), cfg["n_unsafe_decayed"], replace=False))],
                    br.iloc[np.sort(rng_u.choice(len(br), cfg["n_unsafe_breach"], replace=False))]])
    for grp, d in (("a_top5", top), ("b_random_safe", rs), ("c_unsafe", ru)):
        for r, c in enumerate(d.itertuples(), 1):
            rows.append(dict(epoch=e, group=grp, rank_in_group=r, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg,
                             pred_class=c.pred_class, pred_dv_mps=c.lazy_dv, pred_burns=c.lazy_burns, pred_near_dv_mps=c.near_dv,
                             pred_q_max=c.q_max, pred_min_band_margin_m=c.min_band_margin_m, pred_decayed=bool(c.decayed_pred),
                             pred_t_decay_days=c.t_decay_days, pred_n_passes=c.n_pred, seed_m=c.seed_m))
B = pd.DataFrame(rows)
B["tag"] = [f"I{round(10*r.inc_deg):04d}_R{round(r.raan_deg):03d}_W{round(r.omega_deg):03d}_E{r.epoch}" for r in B.itertuples()]
assert B.tag.is_unique and len(B) == 45
out = root / "blind"; out.mkdir(parents=True, exist_ok=True)
B.to_csv(out / "blind_cells.csv", index=False)
print(B[["epoch", "group", "tag", "pred_class", "pred_dv_mps", "pred_q_max"]].to_string())
