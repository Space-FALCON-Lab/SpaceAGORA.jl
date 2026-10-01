#!/usr/bin/env python3
"""Analyse a perturbed-density comparison results directory.

    python3 analyze.py <results_dir> <telemetry_dir> <out_dir>

<results_dir> holds one subdirectory per run tag (run_summary.toml, extrema.csv,
simulation_results.csv, perturbation_log.csv) as written by run_mode.jl.
Writes to <out_dir>:
  runs.csv               one row per run: solver stats, wall time, final apsides,
                         decay, deltas vs nominal, bit-identity flags
  per_orbit.csv          per orbit and run: apo/peri, deltas vs nominal,
                         per-pass decay delta, flight reference (context)
  per_pass_r.csv         per pass and run: realized r stats vs GRAM sigma
  sensitivity.csv        tight / dt01 / rep vs their base, beside nominal's own
Every number is computed from the files; nothing is rounded before writing.
"""
import sys, os, tomllib
import numpy as np
import pandas as pd

res, tele, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
tags = sorted(d for d in os.listdir(res) if os.path.isfile(os.path.join(res, d, "run_summary.toml")))

def load(tag):
    with open(os.path.join(res, tag, "run_summary.toml"), "rb") as f:
        s = tomllib.load(f)
    ex_path = os.path.join(res, tag, "extrema.csv")
    ex = pd.read_csv(ex_path) if os.path.isfile(ex_path) else None
    return s, ex

runs = {t: load(t) for t in tags}
nom_s, nom_ex = runs["nominal"]
nom_apo = nom_ex[nom_ex.event == "apo"].altitude_km.to_numpy()
nom_peri = nom_ex[nom_ex.event == "peri"].altitude_km.to_numpy()

fl_apo = pd.read_feather(os.path.join(tele, "True_Odyssey_apoapsis_alts_kernel.feather"))
fl_peri = pd.read_feather(os.path.join(tele, "True_Odyssey_periapsis_alts_kernel.feather"))
fl_apo = dict(zip(fl_apo.orbit, fl_apo.altitude))
fl_peri = dict(zip(fl_peri.orbit, fl_peri.altitude))

def results_equal(a, b):
    pa = os.path.join(res, a, "simulation_results.csv"); pb = os.path.join(res, b, "simulation_results.csv")
    if not (os.path.isfile(pa) and os.path.isfile(pb)):
        return None
    with open(pa, "rb") as fa, open(pb, "rb") as fb:
        return fa.read() == fb.read()

rows, orbit_rows = [], []
for t, (s, ex) in runs.items():
    r = dict(tag=t, mode=s["mode"], seed=s["gram_seed"], tight=s["tight"],
             dt_max_atmosphere=s["dt_max_atmosphere"], reltol_orbit=s["reltol_orbit"],
             retcode=s["retcode"], naccept=s["naccept"], nreject=s["nreject"], nf=s["nf"],
             wall_s=s["wall_s"], solver_sequence=s["solver_sequence"], commit=s["commit"])
    if ex is not None:
        apo = ex[ex.event == "apo"].altitude_km.to_numpy()
        peri = ex[ex.event == "peri"].altitude_km.to_numpy()
        n = min(len(apo), len(nom_apo)); m = min(len(peri), len(nom_peri))
        r.update(n_apo=len(apo), n_peri=len(peri), apo_first_km=apo[0], apo_final_km=apo[-1],
                 apo_decay_km=apo[0] - apo[-1], peri_min_km=peri.min(), peri_final_km=peri[-1],
                 d_apo_final_km=apo[n - 1] - nom_apo[n - 1],
                 d_decay_km=(apo[0] - apo[n - 1]) - (nom_apo[0] - nom_apo[n - 1]),
                 max_abs_d_apo_km=np.max(np.abs(apo[:n] - nom_apo[:n])),
                 max_abs_d_peri_km=np.max(np.abs(peri[:m] - nom_peri[:m])),
                 rms_d_peri_km=np.sqrt(np.mean((peri[:m] - nom_peri[:m]) ** 2)))
        for k in range(max(len(apo), len(peri))):
            orb = 19.0 + k
            o = dict(tag=t, index=k + 1, flight_orbit=orb)
            if k < len(apo):
                o["apo_km"] = apo[k]
                if k < len(nom_apo):
                    o["d_apo_km"] = apo[k] - nom_apo[k]
                    if k >= 1:
                        o["d_pass_decay_km"] = (apo[k - 1] - apo[k]) - (nom_apo[k - 1] - nom_apo[k])
            if k < len(peri):
                o["peri_km"] = peri[k]
                if k < len(nom_peri):
                    o["d_peri_km"] = peri[k] - nom_peri[k]
            o["flight_apo_km"] = fl_apo.get(orb, np.nan)
            o["flight_peri_km"] = fl_peri.get(orb, np.nan)
            orbit_rows.append(o)
    rows.append(r)
runs_df = pd.DataFrame(rows)
for t in tags:
    if t.endswith("_rep"):
        runs_df.loc[runs_df.tag == t, "bit_identical_to_base"] = results_equal(t, t[:-4])
runs_df.to_csv(os.path.join(out, "runs.csv"), index=False)
pd.DataFrame(orbit_rows).to_csv(os.path.join(out, "per_orbit.csv"), index=False)

# Realized r per pass. A: kind 1 (one walk sample per accepted step below EI).
# B: kind 2 (knots, carry GRAM mean and sigma) and kind 3 (factor applied at
# accepted steps). Density-weighted mean uses GRAM's mean density at the same
# sample, which is what sets the drag impulse.
pr = []
for t in tags:
    lp = os.path.join(res, t, "perturbation_log.csv")
    if not os.path.isfile(lp):
        continue
    lg = pd.read_csv(lp).rename(columns={"pass": "ps"})
    for kind, label in ((1, "A_step_sample"), (2, "B_knot"), (3, "B_applied")):
        g0 = lg[lg.kind == kind]
        for ps, g in g0.groupby("ps"):
            row = dict(tag=t, kind=label, pass_=ps, n=len(g), r_mean=g.r.mean(), r_std=g.r.std(ddof=1),
                       r_min=g.r.min(), r_max=g.r.max())
            if kind in (1, 2):
                w = g.mean_density_kgm3.to_numpy()
                row.update(r_density_weighted=np.sum(w * g.r) / np.sum(w),
                           gram_sigma_median=g.sigma_frac.median(),
                           gram_sigma_density_weighted=np.sum(w * g.sigma_frac) / np.sum(w))
                low = g[g.alt_m < 130e3]
                row.update(n_below130=len(low), r_std_below130=low.r.std(ddof=1) if len(low) > 1 else np.nan,
                           gram_sigma_median_below130=low.sigma_frac.median() if len(low) else np.nan)
            if kind == 3:
                row.update(max_abs_pred_minus_actual_alt_m=g.aux.abs().max())
            pr.append(row)
pr_df = pd.DataFrame(pr).rename(columns={"pass_": "pass"})
pr_df.to_csv(os.path.join(out, "per_pass_r.csv"), index=False)

# Sensitivity: each variant against its base, beside nominal's own response.
sens = []
def final_apo(t):
    ex = runs[t][1]
    return ex[ex.event == "apo"].altitude_km.to_numpy()
def peri_series(t):
    ex = runs[t][1]
    return ex[ex.event == "peri"].altitude_km.to_numpy()
for base in ("nominal", "A_s11", "B_s11"):
    for suf in ("_tight", "_dt01", "_rep"):
        v = base + suf
        if v not in runs or runs[v][1] is None or runs[base][1] is None:
            continue
        a, b = final_apo(v), final_apo(base)
        p, q = peri_series(v), peri_series(base)
        n = min(len(a), len(b)); m = min(len(p), len(q))
        sens.append(dict(variant=v, base=base, d_apo_final_km=a[n - 1] - b[n - 1],
                         max_abs_d_apo_km=np.max(np.abs(a[:n] - b[:n])),
                         max_abs_d_peri_km=np.max(np.abs(p[:m] - q[:m])),
                         bit_identical=results_equal(v, base)))
pd.DataFrame(sens).to_csv(os.path.join(out, "sensitivity.csv"), index=False)
print(runs_df.to_string())
print(pd.DataFrame(sens).to_string())
