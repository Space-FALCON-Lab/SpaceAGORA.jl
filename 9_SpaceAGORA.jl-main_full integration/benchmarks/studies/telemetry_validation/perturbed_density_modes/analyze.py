#!/usr/bin/env python3
"""Analyse a perturbed-density comparison results directory.

    python3 analyze.py <results_dir> <telemetry_dir> <out_dir>

<results_dir> holds one subdirectory per run tag (run_summary.toml, extrema.csv,
simulation_results.csv, perturbation_log.csv) as written by run_mode.jl.

A run is used only when it is verified: its summary says the attempt completed
(attempt_record.jl), every output the summary names matches its recorded
SHA-256, its runtime orbit counter reached the request with a consistent
termination cause, and its saved apsides cover the run. Every run is listed in
runs.csv with whether it was verified and why not; an unverified run enters no
comparison, so a failed rerun beside an earlier run's files cannot report
agreement. Two runs are compared only when both are verified and requested the
same number of orbits. The nominal run must be verified.

Writes to <out_dir>:
  runs.csv               one row per run: verification, solver stats, wall time,
                         final apsides, decay, deltas vs nominal, bit-identity flags
  per_orbit.csv          per orbit and verified run: apo/peri, deltas vs nominal,
                         per-pass decay delta, flight reference (context)
  per_pass_r.csv         per pass and verified run: realized r stats vs GRAM sigma
  sensitivity.csv        tight / dt01 / rep vs their base, beside nominal's own;
                         comparable=False rows carry the reason and no numbers
Every number is computed from the files; nothing is rounded before writing.
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd

VARIANT_SUFFIXES = ("_tight", "_dt01", "_rep")
REQUIRED_OUTPUTS = ("extrema.csv", "simulation_results.csv")


def load_toml(path):
    try:
        import tomllib
    except ModuleNotFoundError:  # Python < 3.11
        import tomli as tomllib
    with open(path, "rb") as f:
        return tomllib.load(f)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_run(run_dir, summary):
    """(verified, reason, extrema) for one run tag's directory and summary."""
    status = summary.get("status")
    if status is None:
        return False, "no attempt record (the summary predates attempt_record.jl)", None
    if status != "complete":
        why = summary.get("status_reason", "")
        return False, f"attempt {status}" + (f": {why}" if why else ""), None
    retcode = summary.get("retcode")
    error = summary.get("error")
    if not isinstance(error, str) or error.strip() or retcode not in ("Success", "Terminated"):
        return False, "missing or unsuccessful solver outcome", None
    orbits = summary.get("orbits_requested")
    completed = summary.get("completed_orbits")
    if not isinstance(orbits, int) or isinstance(orbits, bool) or orbits < 1:
        return False, "missing or invalid requested orbit count", None
    if not isinstance(completed, int) or isinstance(completed, bool) or completed < 0:
        return False, "missing or invalid completed orbit count", None
    expected_cause = ("end_of_time_span" if retcode == "Success" else
                      "orbit_count" if completed >= orbits else "terminated_before_orbit_count")
    if summary.get("termination_cause") != expected_cause:
        return False, "missing or inconsistent termination cause", None
    if completed < orbits:
        return False, f"only {completed} completed orbit events for {orbits} requested", None
    hashes = summary.get("output_sha256") or {}
    for name in REQUIRED_OUTPUTS:
        if name not in hashes:
            return False, f"{name} has no recorded identity", None
    for name, digest in sorted(hashes.items()):
        path = os.path.join(run_dir, name)
        if not os.path.isfile(path):
            return False, f"{name} is missing", None
        if sha256(path) != digest:
            return False, f"{name} does not match its recorded identity", None
    ex = pd.read_csv(os.path.join(run_dir, "extrema.csv"))
    n_apo = int((ex.event == "apo").sum())
    n_peri = int((ex.event == "peri").sum())
    if n_apo < orbits - 1 or n_peri < orbits - 1:
        return False, f"apsides cover {n_apo} apoapses and {n_peri} periapses of {orbits} orbits", None
    if not np.isfinite(ex.altitude_km.to_numpy(dtype=float)).all():
        return False, "non-finite apsis altitude", None
    return True, "", ex


def load_runs(res):
    tags = sorted(d for d in os.listdir(res) if os.path.isfile(os.path.join(res, d, "run_summary.toml")))
    runs = {}
    for t in tags:
        s = load_toml(os.path.join(res, t, "run_summary.toml"))
        ok, why, ex = verify_run(os.path.join(res, t), s)
        runs[t] = dict(summary=s, verified=ok, reason=why, extrema=ex)
    return runs


def comparable(runs, a, b):
    """(comparable, reason) for comparing run a with run b."""
    for t in (a, b):
        if t not in runs:
            return False, f"{t} has no run"
        if not runs[t]["verified"]:
            return False, f"{t} is not verified ({runs[t]['reason']})"
    if int(runs[a]["summary"]["orbits_requested"]) != int(runs[b]["summary"]["orbits_requested"]):
        return False, "different requested orbit counts"
    return True, ""


def apsides(run, event):
    ex = run["extrema"]
    return ex[ex.event == event].altitude_km.to_numpy()


def bit_identical(runs, a, b):
    # Both runs are verified: the recorded identities are those of the files.
    return (runs[a]["summary"]["output_sha256"]["simulation_results.csv"] ==
            runs[b]["summary"]["output_sha256"]["simulation_results.csv"])


def run_rows(runs, fl_apo, fl_peri):
    nom = runs["nominal"]
    nom_apo, nom_peri = apsides(nom, "apo"), apsides(nom, "peri")
    rows, orbit_rows = [], []
    for t, run in runs.items():
        s = run["summary"]
        r = dict(tag=t, verified=run["verified"], unverified_reason=run["reason"],
                 status=s.get("status", ""), attempt_id=s.get("attempt_id", ""),
                 termination_cause=s.get("termination_cause", ""), orbits_requested=s.get("orbits_requested"),
                 mode=s.get("mode"), seed=s.get("gram_seed"), tight=s.get("tight"),
                 dt_max_atmosphere=s.get("dt_max_atmosphere"), reltol_orbit=s.get("reltol_orbit"),
                 retcode=s.get("retcode"), naccept=s.get("naccept"), nreject=s.get("nreject"), nf=s.get("nf"),
                 wall_s=s.get("wall_s"), solver_sequence=s.get("solver_sequence"), commit=s.get("commit"))
        if run["verified"]:
            apo, peri = apsides(run, "apo"), apsides(run, "peri")
            r.update(n_apo=len(apo), n_peri=len(peri), apo_first_km=apo[0], apo_final_km=apo[-1],
                     apo_decay_km=apo[0] - apo[-1], peri_min_km=peri.min(), peri_final_km=peri[-1])
            ok, why = comparable(runs, t, "nominal")
            r.update(comparable_to_nominal=ok, not_comparable_reason=why)
            if ok:
                n = min(len(apo), len(nom_apo)); m = min(len(peri), len(nom_peri))
                r.update(d_apo_final_km=apo[n - 1] - nom_apo[n - 1],
                         d_decay_km=(apo[0] - apo[n - 1]) - (nom_apo[0] - nom_apo[n - 1]),
                         max_abs_d_apo_km=np.max(np.abs(apo[:n] - nom_apo[:n])),
                         max_abs_d_peri_km=np.max(np.abs(peri[:m] - nom_peri[:m])),
                         rms_d_peri_km=np.sqrt(np.mean((peri[:m] - nom_peri[:m]) ** 2)))
            for k in range(max(len(apo), len(peri))):
                orb = 19.0 + k
                o = dict(tag=t, index=k + 1, flight_orbit=orb)
                if k < len(apo):
                    o["apo_km"] = apo[k]
                    if ok and k < len(nom_apo):
                        o["d_apo_km"] = apo[k] - nom_apo[k]
                        if k >= 1:
                            o["d_pass_decay_km"] = (apo[k - 1] - apo[k]) - (nom_apo[k - 1] - nom_apo[k])
                if k < len(peri):
                    o["peri_km"] = peri[k]
                    if ok and k < len(nom_peri):
                        o["d_peri_km"] = peri[k] - nom_peri[k]
                o["flight_apo_km"] = fl_apo.get(orb, np.nan)
                o["flight_peri_km"] = fl_peri.get(orb, np.nan)
                orbit_rows.append(o)
        if t.endswith("_rep"):
            ok, _ = comparable(runs, t, t[:-4])
            r["bit_identical_to_base"] = str(bit_identical(runs, t, t[:-4])) if ok else "not_comparable"
        rows.append(r)
    return pd.DataFrame(rows), pd.DataFrame(orbit_rows)


def per_pass_rows(res, runs):
    # Realized r per pass. A: kind 1 (one walk sample per accepted step below EI).
    # B: kind 2 (knots, carry GRAM mean and sigma) and kind 3 (factor applied at
    # accepted steps). Density-weighted mean uses GRAM's mean density at the same
    # sample, which is what sets the drag impulse. Only a log that is part of a
    # verified attempt (its identity recorded and matched) is read.
    pr = []
    for t, run in runs.items():
        if not run["verified"] or "perturbation_log.csv" not in run["summary"]["output_sha256"]:
            continue
        lg = pd.read_csv(os.path.join(res, t, "perturbation_log.csv")).rename(columns={"pass": "ps"})
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
    return pd.DataFrame(pr).rename(columns={"pass_": "pass"})


def sensitivity_rows(runs):
    # Each variant against its base, beside nominal's own response.
    sens = []
    for base in sorted({t for t in runs if not t.endswith(VARIANT_SUFFIXES)}):
        for suf in VARIANT_SUFFIXES:
            v = base + suf
            if v not in runs:
                continue
            ok, why = comparable(runs, v, base)
            row = dict(variant=v, base=base, comparable=ok, not_comparable_reason=why,
                       d_apo_final_km=np.nan, max_abs_d_apo_km=np.nan, max_abs_d_peri_km=np.nan,
                       bit_identical=None)
            if ok:
                a, b = apsides(runs[v], "apo"), apsides(runs[base], "apo")
                p, q = apsides(runs[v], "peri"), apsides(runs[base], "peri")
                n = min(len(a), len(b)); m = min(len(p), len(q))
                row.update(d_apo_final_km=a[n - 1] - b[n - 1],
                           max_abs_d_apo_km=np.max(np.abs(a[:n] - b[:n])),
                           max_abs_d_peri_km=np.max(np.abs(p[:m] - q[:m])),
                           bit_identical=bit_identical(runs, v, base))
            sens.append(row)
    return pd.DataFrame(sens, columns=["variant", "base", "comparable", "not_comparable_reason", "d_apo_final_km",
                                       "max_abs_d_apo_km", "max_abs_d_peri_km", "bit_identical"])


def main(argv):
    res, tele, out = argv[:3]
    runs = load_runs(res)
    if "nominal" not in runs:
        raise ValueError(f"no nominal run in {res}; every comparison is against it")
    if not runs["nominal"]["verified"]:
        raise ValueError(f"the nominal run is not verified ({runs['nominal']['reason']}); "
                         "every comparison is against it")
    fl_apo = pd.read_feather(os.path.join(tele, "True_Odyssey_apoapsis_alts_kernel.feather"))
    fl_peri = pd.read_feather(os.path.join(tele, "True_Odyssey_periapsis_alts_kernel.feather"))
    fl_apo = dict(zip(fl_apo.orbit, fl_apo.altitude))
    fl_peri = dict(zip(fl_peri.orbit, fl_peri.altitude))

    runs_df, orbit_df = run_rows(runs, fl_apo, fl_peri)
    pr_df = per_pass_rows(res, runs)
    sens_df = sensitivity_rows(runs)

    os.makedirs(out, exist_ok=True)
    runs_df.to_csv(os.path.join(out, "runs.csv"), index=False)
    orbit_df.to_csv(os.path.join(out, "per_orbit.csv"), index=False)
    pr_df.to_csv(os.path.join(out, "per_pass_r.csv"), index=False)
    sens_df.to_csv(os.path.join(out, "sensitivity.csv"), index=False)
    print(runs_df.to_string())
    print(sens_df.to_string())


if __name__ == "__main__":
    main(sys.argv[1:])
