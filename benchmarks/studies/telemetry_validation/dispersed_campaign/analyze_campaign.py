#!/usr/bin/env python3
"""Summarise a dispersed Odyssey campaign results directory.

    python3 analyze_campaign.py <results_dir> <telemetry_dir> <out_dir> [n_passes=40]

The analysis horizon is n_passes passes: apoapses 1..n+1 and periapses (passes)
1..n, so decay over the horizon is apo[1] - apo[n+1]. run_campaign.jl
--orbits=n+2 reaches it. Each member is classified before any statistic:

  dispatch_failure   the dispatcher holds no member value
  solver_failure     the member's solve threw or returned an unsuccessful retcode
  early_termination  the member's apsides do not cover the horizon (Terminated
                     covers both the orbit count and impact)
  invalid_output     a member that covers the horizon but whose output is
                     missing, out of event order or not finite there
  complete           the rest

Only complete members enter spread.csv, r_stats.csv and per_orbit_histories.csv.
Excluded members are listed with their status and reason in
per_sample_summary.csv (their histories in excluded_member_histories.csv), the
counts are in member_status_counts.csv and campaign_facts.txt, and every
statistic reports the number of members it used. The nominal member must be
complete and the flight telemetry must cover the horizon, or the analysis stops
with ValueError before writing anything. Percentiles are numpy's default
(linear interpolation).

Writes per_sample_summary.csv, member_status_counts.csv, per_orbit_histories.csv,
excluded_member_histories.csv, per_pass_r.csv, nominal_series.csv,
flight_series.csv, spread.csv, r_stats.csv and campaign_facts.txt.
"""
import os
import sys

import numpy as np
import pandas as pd

STATUSES = ("complete", "early_termination", "invalid_output", "solver_failure", "dispatch_failure")
SOLVER_OK = ("Success", "Terminated")
SPREAD_QUANTITIES = ("apo_decay_km", "peri_min_km", "apo_final_km", "heat_load_horizon_Jcm2",
                     "peak_heat_rate_Wcm2")
PASS_COLUMNS = ["peri_time_s", "peri_km", "pass_heat_load_Jcm2", "pass_peak_heat_rate_Wcm2"]
FACT_KEYS = ("route", "consumers", "local_slots", "wall_s", "sum_member_solve_s", "sum_dispatch_elapsed_s",
             "distinct_member_pids", "n_samples", "horizon_passes", "n_complete", "n_early_terminated",
             "n_solver_failed", "n_dispatch_failed", "n_failed", "started_utc", "finished_utc", "commit")


def load_toml(path):
    try:
        import tomllib
    except ModuleNotFoundError:  # Python < 3.11
        import tomli as tomllib
    with open(path, "rb") as f:
        return tomllib.load(f)


def _is_true(x):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    return str(x).strip().lower() == "true"


def _text(x):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x)


def run_status(row):
    """(status, reason) from a member's summary row, or (None, "") when its run succeeded."""
    if not _is_true(row.get("dispatch_success", False)):
        return "dispatch_failure", _text(row.get("error", "")) or "the dispatcher holds no member value"
    err = _text(row.get("error", "")).strip()
    if err:
        return "solver_failure", err
    retcode = _text(row.get("retcode", "")) or "ERROR"
    if retcode not in SOLVER_OK:
        return "solver_failure", f"retcode {retcode}"
    return None, ""


def horizon_record(po, n):
    """(status, reason, metrics) of one member's per-orbit series over an n-pass horizon."""
    apo = po.apo_km.to_numpy(dtype=float)
    peri = po.peri_km.to_numpy(dtype=float)
    n_apo, n_peri = int(np.isfinite(apo).sum()), int(np.isfinite(peri).sum())
    if n_apo < n + 1 or n_peri < n:
        return ("early_termination",
                f"{n_apo} apoapses and {n_peri} passes; the {n}-pass horizon needs {n + 1} and {n}", None)
    h_apo = po[["apo_time_s", "apo_km"]].to_numpy(dtype=float)[: n + 1]
    h_pass = po[PASS_COLUMNS].to_numpy(dtype=float)[:n]
    if not (np.isfinite(h_apo).all() and np.isfinite(h_pass).all()):
        return "invalid_output", "a non-finite value within the horizon", None
    apo_t, peri_t = h_apo[:, 0], h_pass[:, 0]
    if not (np.all(apo_t[:-1] < peri_t) and np.all(peri_t < apo_t[1:])):
        return "invalid_output", "apsides out of order (pass k must lie between apoapses k and k+1)", None
    a = h_apo[:, 1]
    metrics = dict(apo_first_km=a[0], apo_final_km=a[n], apo_decay_km=a[0] - a[n],
                   peri_min_km=h_pass[:, 1].min(), heat_load_horizon_Jcm2=h_pass[:, 2].sum(),
                   peak_heat_rate_Wcm2=h_pass[:, 3].max())
    return "complete", "", metrics


def horizon_history(po, n):
    """The horizon's rows: apoapses 1..n+1, passes 1..n (the pass after the horizon blanked)."""
    hist = po.iloc[: n + 1].copy()
    hist.loc[hist.index >= n, PASS_COLUMNS] = np.nan
    return hist


def member_record(row, member_dir, n):
    """(status, reason, metrics, per_orbit, per_pass_r) for one member."""
    status, reason = run_status(row)
    if status is not None:
        return status, reason, None, None, None
    path = os.path.join(member_dir, "per_orbit.csv")
    if not os.path.isfile(path):
        return "invalid_output", "per_orbit.csv is missing", None, None, None
    po = pd.read_csv(path)
    status, reason, metrics = horizon_record(po, n)
    pr_path = os.path.join(member_dir, "per_pass_r.csv")
    pr = pd.read_csv(pr_path) if os.path.isfile(pr_path) else None
    return status, reason, metrics, po, pr


def spread_row(q, x, n_members, nominal, flight=None):
    x = np.asarray(x, dtype=float)
    row = dict(quantity=q, n=len(x), n_members=n_members, n_excluded=n_members - len(x),
               min=x.min(), p5=np.percentile(x, 5), p50=np.percentile(x, 50), p95=np.percentile(x, 95),
               max=x.max(), mean=x.mean(), std=x.std(ddof=1) if len(x) > 1 else np.nan,
               nominal=nominal, nominal_rank_below=int((x < nominal).sum()))
    if flight is not None:
        row.update(flight=flight, flight_rank_below=int((x < flight).sum()))
    return row


def flight_series(tele, orbits, n):
    fa = pd.read_feather(os.path.join(tele, "True_Odyssey_apoapsis_alts_kernel.feather")).set_index("orbit")
    fp = pd.read_feather(os.path.join(tele, "True_Odyssey_periapsis_alts_kernel.feather")).set_index("orbit")
    flight = pd.DataFrame(dict(
        flight_orbit=orbits,
        apo_km=[fa.altitude.get(o, np.nan) for o in orbits],
        peri_km=[fp.altitude.get(o, np.nan) if i < n else np.nan for i, o in enumerate(orbits)]))
    missing = [o for o, v in zip(orbits, flight.apo_km) if not np.isfinite(v)]
    missing += [o for o, v in zip(orbits[:n], flight.peri_km[:n]) if not np.isfinite(v)]
    if missing:
        raise ValueError(f"the flight telemetry does not cover the {n}-pass horizon: no apsis for orbit(s) "
                         f"{sorted(set(missing))}")
    vals = dict(apo_decay_km=flight.apo_km.iloc[0] - flight.apo_km.iloc[n],
                peri_min_km=flight.peri_km.iloc[:n].min(), apo_final_km=flight.apo_km.iloc[n])
    return flight, vals


def analyze(res, tele, n):
    """Every output table, computed and validated before anything is written."""
    if n < 1:
        raise ValueError(f"n_passes must be at least 1, got {n}")
    samples = pd.read_csv(os.path.join(res, "samples.csv"))
    camp = load_toml(os.path.join(res, "campaign.toml"))

    rows, hists, excluded, prs = [], [], [], []
    for _, r in samples.sort_values("seed").iterrows():
        seed = int(r.seed)
        status, reason, metrics, po, pr = member_record(r, os.path.join(res, f"sample_seed{seed}"), n)
        base = dict(seed=seed, status=status, status_reason=reason, termination=_text(r.get("retcode", "")),
                    termination_cause=_text(r.get("termination_cause", "")),
                    completed_orbits=r.get("completed_orbits", np.nan), driver_status=_text(r.get("status", "")),
                    dispatch_success=r.get("dispatch_success"), error=_text(r.get("error", "")),
                    solve_s=r.get("solve_s", np.nan), dispatch_elapsed_s=r.get("dispatch_elapsed_s", np.nan),
                    worker_pid=r.get("pid", np.nan), naccept=r.get("naccept", np.nan),
                    nreject=r.get("nreject", np.nan), nf=r.get("nf", np.nan))
        if po is not None:
            base.update(n_apo=int(np.isfinite(po.apo_km.to_numpy(dtype=float)).sum()),
                        n_peri=int(np.isfinite(po.peri_km.to_numpy(dtype=float)).sum()))
        if status == "complete":
            base.update(metrics)
            h = horizon_history(po, n); h.insert(0, "seed", seed); hists.append(h)
            if pr is not None:
                pr = pr.iloc[:n].copy()
                base.update(n_r_passes=len(pr), r_density_weighted_mean=pr.r_density_weighted.mean(),
                            r_density_weighted_std=pr.r_density_weighted.std(ddof=1))
                pr.insert(0, "seed", seed); prs.append(pr)
        elif po is not None:
            h = po.copy(); h.insert(0, "seed", seed); h.insert(1, "status", status); excluded.append(h)
        rows.append(base)
    summ = pd.DataFrame(rows)
    counts = pd.DataFrame(dict(status=list(STATUSES),
                               count=[int((summ.status == s).sum()) for s in STATUSES]))
    ok = summ[summ.status == "complete"]
    if ok.empty:
        raise ValueError(f"no complete member over the {n}-pass horizon; member statuses: "
                         f"{dict(zip(counts.status, counts['count']))}")

    nom_row = pd.read_csv(os.path.join(res, "nominal_summary.csv")).iloc[0].to_dict()
    nom_row["dispatch_success"] = True    # run directly, not dispatched
    nom_status, nom_reason, nom_metrics, nom_po, _ = member_record(nom_row, os.path.join(res, "nominal"), n)
    if nom_status != "complete":
        raise ValueError(f"the nominal run is {nom_status} ({nom_reason}); the {n}-pass analysis needs "
                         f"{n + 1} apoapses and {n} passes (run_campaign.jl --orbits={n + 2})")
    nom_h = horizon_history(nom_po, n)
    flight, flight_vals = flight_series(tele, nom_h.flight_orbit.to_numpy(), n)

    n_members = len(summ)
    spread = pd.DataFrame([spread_row(q, ok[q], n_members, nom_metrics[q], flight_vals.get(q))
                           for q in SPREAD_QUANTITIES])
    spread.insert(1, "horizon_passes", n)

    prall = pd.concat(prs) if prs else pd.DataFrame(
        columns=["seed", "pass", "r_density_weighted", "r_std_below130", "sigma_median_below130",
                 "sigma_density_weighted"])
    rdw = prall.r_density_weighted.to_numpy(dtype=float)
    r_stats = pd.DataFrame([dict(
        horizon_passes=n, passes=len(prall), members=prall.seed.nunique(),
        r_dw_mean=rdw.mean() if len(rdw) else np.nan, r_dw_median=np.median(rdw) if len(rdw) else np.nan,
        r_dw_std=rdw.std(ddof=1) if len(rdw) > 1 else np.nan,
        r_dw_p5=np.percentile(rdw, 5) if len(rdw) else np.nan,
        r_dw_p95=np.percentile(rdw, 95) if len(rdw) else np.nan,
        r_std_below130_mean=prall.r_std_below130.mean(),
        sigma_below130_median=prall.sigma_median_below130.median(),
        sigma_density_weighted_median=prall.sigma_density_weighted.median(),
        low_side_sigma_equiv=1 - 1 / (1 + prall.sigma_median_below130.median()))])

    facts = [f"{k} = {camp[k]}" for k in FACT_KEYS if k in camp]
    facts += [f"analysis_horizon_passes = {n}", f"analysis_n_members = {n_members}"]
    facts += [f"analysis_n_{s} = {c}" for s, c in zip(counts.status, counts["count"])]
    facts += [f"nominal_solve_s = {nom_row.get('solve_s')}", f"nominal_retcode = {nom_row.get('retcode')}"]

    excl = pd.concat(excluded) if excluded else pd.DataFrame(columns=["seed", "status"])
    return dict(per_sample_summary=summ, member_status_counts=counts, per_orbit_histories=pd.concat(hists),
                excluded_member_histories=excl, per_pass_r=prall, nominal_series=nom_h, flight_series=flight,
                spread=spread, r_stats=r_stats, facts=facts)


def main(argv):
    res, tele, out = argv[:3]
    n = int(argv[3]) if len(argv) > 3 else 40
    tables = analyze(res, tele, n)
    os.makedirs(out, exist_ok=True)
    for name, frame in tables.items():
        if name != "facts":
            frame.to_csv(os.path.join(out, f"{name}.csv"), index=False)
    with open(os.path.join(out, "campaign_facts.txt"), "w") as f:
        f.write("\n".join(tables["facts"]) + "\n")
    print("\n".join(tables["facts"]))
    print(tables["spread"].to_string())
    print(tables["r_stats"].T.to_string())


if __name__ == "__main__":
    main(sys.argv[1:])
