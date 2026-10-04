#!/usr/bin/env python3
"""Summarise a dispersed Odyssey campaign results directory.

    python3 analyze_campaign.py <results_dir> <telemetry_dir> <out_dir> [n_passes=40]

The runs end at the apoapsis after pass n_passes+1 (run_campaign.jl --orbits=n+2),
so apoapsis 1..n+1 and periapsis/pass 1..n are kept: decay over n passes is
apo[1] - apo[n+1]. Percentiles are numpy's default (linear interpolation).
Writes per_sample_summary.csv, per_orbit_histories.csv, per_pass_r.csv,
nominal_series.csv, flight_series.csv, spread.csv, r_stats.csv.
"""
import sys, os, glob, tomllib
import numpy as np
import pandas as pd

res, tele, out = sys.argv[1:4]
N = int(sys.argv[4]) if len(sys.argv) > 4 else 40
os.makedirs(out, exist_ok=True)
samples = pd.read_csv(os.path.join(res, "samples.csv"))
with open(os.path.join(res, "campaign.toml"), "rb") as f:
    camp = tomllib.load(f)

def member(d):
    po = pd.read_csv(os.path.join(d, "per_orbit.csv"))
    apo = po.apo_km.to_numpy()[: N + 1]
    pp = po.iloc[:N]
    pr_path = os.path.join(d, "per_pass_r.csv")
    pr = pd.read_csv(pr_path).iloc[:N] if os.path.isfile(pr_path) else None
    s = dict(n_apo_kept=len(apo), n_pass_kept=len(pp),
             apo_first_km=apo[0], apo_final_km=apo[-1], apo_decay_km=apo[0] - apo[-1],
             peri_min_km=pp.peri_km.min(), heat_load_40_passes_Jcm2=pp.pass_heat_load_Jcm2.sum(),
             peak_heat_rate_Wcm2=pp.pass_peak_heat_rate_Wcm2.max())
    if pr is not None:
        s.update(r_density_weighted_mean=pr.r_density_weighted.mean(),
                 r_density_weighted_std=pr.r_density_weighted.std(ddof=1))
    hist = po.iloc[: N + 1].copy()
    hist.loc[hist.index >= N, ["peri_time_s", "peri_km", "pass_heat_load_Jcm2", "pass_peak_heat_rate_Wcm2"]] = np.nan
    return s, hist, pr

rows, hists, prs = [], [], []
for _, r in samples.sort_values("seed").iterrows():
    d = os.path.join(res, f"sample_seed{int(r.seed)}")
    base = dict(seed=int(r.seed), termination=r.get("retcode", "ERROR"), dispatch_success=r.dispatch_success,
                error=r.get("error", ""), solve_s=r.get("solve_s", np.nan),
                dispatch_elapsed_s=r.dispatch_elapsed_s, worker_pid=r.get("pid", np.nan),
                naccept=r.get("naccept", np.nan), nreject=r.get("nreject", np.nan), nf=r.get("nf", np.nan))
    if os.path.isfile(os.path.join(d, "per_orbit.csv")):
        s, h, pr = member(d)
        base.update(s)
        h.insert(0, "seed", int(r.seed)); hists.append(h)
        if pr is not None:
            pr.insert(0, "seed", int(r.seed)); prs.append(pr)
    rows.append(base)
summ = pd.DataFrame(rows)
summ.to_csv(os.path.join(out, "per_sample_summary.csv"), index=False)
pd.concat(hists).to_csv(os.path.join(out, "per_orbit_histories.csv"), index=False)
prall = pd.concat(prs); prall.to_csv(os.path.join(out, "per_pass_r.csv"), index=False)

nom_s, nom_h, _ = member(os.path.join(res, "nominal"))
nom_h.to_csv(os.path.join(out, "nominal_series.csv"), index=False)
nom_sum = pd.read_csv(os.path.join(res, "nominal_summary.csv")).iloc[0]

fa = pd.read_feather(os.path.join(tele, "True_Odyssey_apoapsis_alts_kernel.feather"))
fp = pd.read_feather(os.path.join(tele, "True_Odyssey_periapsis_alts_kernel.feather"))
orb = nom_h.flight_orbit.to_numpy()
fa = fa.set_index("orbit"); fp = fp.set_index("orbit")
flight = pd.DataFrame(dict(flight_orbit=orb, apo_km=[fa.altitude.get(o, np.nan) for o in orb],
                           peri_km=[fp.altitude.get(o, np.nan) if i < N else np.nan for i, o in enumerate(orb)]))
flight.to_csv(os.path.join(out, "flight_series.csv"), index=False)
flight_vals = dict(apo_decay_km=flight.apo_km.iloc[0] - flight.apo_km.iloc[N],
                   peri_min_km=flight.peri_km.iloc[:N].min(), apo_final_km=flight.apo_km.iloc[N])

ok = summ[summ.termination == "Terminated"]
spread = []
for q in ("apo_decay_km", "peri_min_km", "apo_final_km", "heat_load_40_passes_Jcm2", "peak_heat_rate_Wcm2"):
    x = ok[q].to_numpy()
    row = dict(quantity=q, n=len(x), min=x.min(), p5=np.percentile(x, 5), p50=np.percentile(x, 50),
               p95=np.percentile(x, 95), max=x.max(), mean=x.mean(), std=x.std(ddof=1),
               nominal=nom_s[q], nominal_rank_below=int((x < nom_s[q]).sum()))
    if q in flight_vals:
        row.update(flight=flight_vals[q], flight_rank_below=int((x < flight_vals[q]).sum()))
    spread.append(row)
pd.DataFrame(spread).to_csv(os.path.join(out, "spread.csv"), index=False)

rdw = prall.r_density_weighted
pd.DataFrame([dict(passes=len(prall), members=prall.seed.nunique(),
                   r_dw_mean=rdw.mean(), r_dw_median=rdw.median(), r_dw_std=rdw.std(ddof=1),
                   r_dw_p5=np.percentile(rdw, 5), r_dw_p95=np.percentile(rdw, 95),
                   r_std_below130_mean=prall.r_std_below130.mean(),
                   sigma_below130_median=prall.sigma_median_below130.median(),
                   sigma_density_weighted_median=prall.sigma_density_weighted.median(),
                   low_side_sigma_equiv=1 - 1 / (1 + prall.sigma_median_below130.median()))]
             ).to_csv(os.path.join(out, "r_stats.csv"), index=False)

with open(os.path.join(out, "campaign_facts.txt"), "w") as f:
    for k in ("route", "consumers", "local_slots", "wall_s", "sum_member_solve_s", "sum_dispatch_elapsed_s",
              "distinct_member_pids", "n_samples", "n_failed", "started_utc", "finished_utc", "commit"):
        f.write(f"{k} = {camp.get(k)}\n")
    f.write(f"nominal_solve_s = {nom_sum.solve_s}\nnominal_retcode = {nom_sum.retcode}\n")
print(open(os.path.join(out, "campaign_facts.txt")).read())
print(pd.DataFrame(spread).to_string())
print(pd.read_csv(os.path.join(out, "r_stats.csv")).T.to_string())
