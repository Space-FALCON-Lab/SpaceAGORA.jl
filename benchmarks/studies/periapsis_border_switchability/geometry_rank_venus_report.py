"""Summaries of the Venus (and Mars-control) ranking outputs for venus_geometry_rank_results.md. Usage: python3 -I geometry_rank_venus_report.py <venus_rank_data dir>"""
import json, pathlib, sys
import numpy as np, pandas as pd
D = pathlib.Path(sys.argv[1]); pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
CH = ["geometry", "drag", "zonal", "tesseral", "third_body", "speed"]


def classes(d, title):
    print("\n## " + title)
    g = d.groupby("inc_deg").agg(n=("floor_mps", "size"), zero=("floor_mps", lambda x: int((x <= 0).sum())), median=("floor_mps", "median"), p10=("floor_mps", lambda x: x.quantile(.1)), p90=("floor_mps", lambda x: x.quantile(.9)),
                                 best=("floor_mps", "min"), worst=("floor_mps", "max"), zb=("zero_burn_slack", "median"), R=("R", "median"), neutral=("regime", lambda x: int((x == "neutral").sum())))
    print(g.round(3).to_string()); print("zero-floor share %.3f (%d of %d); feasible %.3f; regimes %s" % ((d.floor_mps <= 0).mean(), (d.floor_mps <= 0).sum(), len(d), d.feasible.mean(), d.regime.value_counts().to_dict()))


def top(d, n=10):
    e = d.dropna(subset=["floor_ens_mps"]).sort_values(["floor_ens_mps", "floor_mps"]).head(n)
    print(e[["id", "inc_deg", "omega0_deg", "raan0_deg", "epoch_off_d", "floor_mps", "floor_ens_mps", "floor_ens_sd_mps", "R", "zero_burn_slack", "regime"]].round(4).to_string(index=False))
    e = d.sort_values("floor_mps", ascending=False).head(3); print("worst (det):"); print(e[["id", "inc_deg", "omega0_deg", "raan0_deg", "epoch_off_d", "floor_mps", "R"]].round(3).to_string(index=False))


def channels(d, title):
    print("\n## channels " + title)
    rows = []
    for k in CH:
        if "rng_" + k in d: rows.append((k, d["rng_" + k].median(), d["net_" + k].abs().median(), d["rms_inc_" + k].median(), (d["rng_" + k] ** 2).mean()))
    t = pd.DataFrame(rows, columns=["channel", "median_range_ln", "median_abs_net_ln", "median_rms_per_pass_ln", "mean_range_sq"]); t["share_of_range_sq"] = t.mean_range_sq / t.mean_range_sq.sum()
    print(t.drop(columns="mean_range_sq").round(4).to_string(index=False))


for tag, name in (("rank_cut", "Venus cut-256"), ("rank_fine", "Venus fine (19 incl, full omega/Omega/epoch)"), ("rank_cut_nosun", "Venus cut-256 WITHOUT Sun"), ("mars_control", "Mars control (no mean field) cut-256")):
    p = D / tag / "geometry_rank.csv"
    if not p.exists(): continue
    d = pd.read_csv(p); classes(d, name); channels(d, name)
    a = pd.read_csv(D / tag / "anchor.csv").iloc[0]; print("anchor: det %.3f ens %.3f +- %.3f, R %.2f, regime %s, q_req %.4f W_eff %.3f days %.1f passes %d" % (a.floor_det_mps, a.floor_ens_mps, a.floor_ens_sd_mps, a.R, a.regime, a.q_req, a.W_eff, a.days, a.n_passes))
    if d.floor_ens_mps.notna().any(): print("top by ensemble:"); top(d)
print("\n## sensitivities (deterministic single-realization floors; the rule is applied to the deterministic anchor and best floors)")
rows = []
for p in sorted(D.glob("sens*/geometry_rank.csv")) + [D / "rank_cut" / "geometry_rank.csv", D / "rank_cut_nosun" / "geometry_rank.csv"]:
    d = pd.read_csv(p); a = pd.read_csv(p.parent / "anchor.csv").iloc[0]; m = d.groupby("inc_deg").floor_mps.median(); b = d.loc[d.floor_mps.idxmin()]
    dec = json.loads((p.parent / "geometry_decision.json").read_text())["campaign"]
    rows.append(dict(run=p.parent.name, f_anchor_det=a.floor_det_mps, f_best_det=b.floor_mps, rule_passes=bool(b.floor_mps <= 0.7 * a.floor_det_mps), best_inc=b.inc_deg, zero_share=(d.floor_mps <= 0).mean(),
                     med_crit=float(np.mean([m[63.4], m[116.6]])), med_polar=m[89.876], med_39=m[39.2], q_req_anchor=a.q_req, W_eff=a.W_eff, ra_target_km=dec["ra_target_m"] / 1e3, d_ref=dec["d_ref_days"]))
print(pd.DataFrame(rows).round(3).to_string(index=False))
