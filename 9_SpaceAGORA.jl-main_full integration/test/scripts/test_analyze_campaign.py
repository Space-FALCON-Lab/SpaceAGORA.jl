"""Tests for the dispersed-campaign analyzer
(benchmarks/studies/telemetry_validation/dispersed_campaign/analyze_campaign.py)
on fabricated campaign directories. The numbers are invented to exercise the
member classification and the fixed-horizon statistics; none is a measurement.
Needs numpy, pandas, pyarrow (feather) and tomllib or tomli."""
import contextlib
import importlib.util
import io
import os
import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[2]
TOOL = REPO / "benchmarks" / "studies" / "telemetry_validation" / "dispersed_campaign" / "analyze_campaign.py"
spec = importlib.util.spec_from_file_location("analyze_campaign", TOOL)
ac = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ac)

N = 40
OFFSET = 19


def per_orbit(n_apo, n_peri, apo0=30000.0, decay=1.0, heat=None, peri=None):
    """A per_orbit.csv frame as odyssey_sample.jl writes it: rows to the larger
    count, the shorter series NaN-padded, pass k bounded by apoapses k and k+1."""
    n_rows = max(n_apo, n_peri)
    k = np.arange(n_rows, dtype=float)

    def pad(v, m):
        return np.concatenate([np.asarray(v, dtype=float)[:m], np.full(n_rows - m, np.nan)])
    m = min(n_apo, n_peri)
    return pd.DataFrame(dict(
        index=np.arange(1, n_rows + 1), flight_orbit=OFFSET + np.arange(n_rows),
        apo_time_s=pad(k * 1000.0, n_apo), apo_km=pad(apo0 - decay * k, n_apo),
        peri_time_s=pad(k * 1000.0 + 500.0, n_peri),
        peri_km=pad(np.full(n_rows, 100.0) if peri is None else peri, n_peri),
        pass_heat_load_Jcm2=pad(np.ones(n_rows) if heat is None else heat, m),
        pass_peak_heat_rate_Wcm2=pad(np.ones(n_rows) if heat is None else heat, m)))


def write_campaign(root, members, nominal=None, flight_orbits=80):
    """members: seed -> (summary row fields, per_orbit frame or None)."""
    root = pathlib.Path(root)
    rows = []
    for seed, (fields, po) in members.items():
        row = dict(seed=seed, dispatch_success=True, dispatch_elapsed_s=1.0, error="", retcode="Terminated",
                   termination_cause="orbit_count", completed_orbits=N + 2, solve_s=1.0, pid=1,
                   naccept=1, nreject=0, nf=1)
        row.update(fields)
        rows.append(row)
        if po is not None:
            d = root / f"sample_seed{seed}"
            d.mkdir(parents=True)
            po.to_csv(d / "per_orbit.csv", index=False)
            n = len(po)
            pd.DataFrame({"pass": np.arange(1, n + 1), "knots": 6, "r_density_weighted": np.ones(n),
                          "sigma_density_weighted": np.full(n, 0.35), "r_std": np.full(n, 0.2),
                          "r_std_below130": np.full(n, 0.2), "sigma_median_below130": np.full(n, 0.35)}
                         ).to_csv(d / "per_pass_r.csv", index=False)
    pd.DataFrame(rows).to_csv(root / "samples.csv", index=False)
    (root / "campaign.toml").write_text('route = "threads"\nn_samples = %d\nhorizon_passes = %d\n' % (len(rows), N))
    (root / "nominal").mkdir()
    (per_orbit(N + 2, N + 1) if nominal is None else nominal).to_csv(root / "nominal" / "per_orbit.csv", index=False)
    pd.DataFrame([dict(seed=1001, retcode="Terminated", error="", solve_s=2.0)]).to_csv(
        root / "nominal_summary.csv", index=False)
    tele = root / "telemetry"
    tele.mkdir()
    orbits = OFFSET + np.arange(flight_orbits)
    pd.DataFrame(dict(orbit=orbits, altitude=30000.0 - 1.5 * np.arange(flight_orbits))).to_feather(
        tele / "True_Odyssey_apoapsis_alts_kernel.feather")
    pd.DataFrame(dict(orbit=orbits, altitude=np.full(flight_orbits, 99.0))).to_feather(
        tele / "True_Odyssey_periapsis_alts_kernel.feather")
    return str(root), str(tele)


class MemberClassification(unittest.TestCase):
    def test_run_status(self):
        self.assertEqual(ac.run_status(dict(dispatch_success=True, error=np.nan, retcode="Terminated")), (None, ""))
        self.assertEqual(ac.run_status(dict(dispatch_success="true", error="", retcode="Success")), (None, ""))
        self.assertEqual(ac.run_status(dict(dispatch_success=False, error="lost worker"))[0], "dispatch_failure")
        self.assertEqual(ac.run_status(dict(dispatch_success=True, error="DomainError", retcode="ERROR"))[0],
                         "solver_failure")
        self.assertEqual(ac.run_status(dict(dispatch_success=True, error=np.nan, retcode="Unstable")),
                         ("solver_failure", "retcode Unstable"))

    def test_horizon_record(self):
        status, _, m = ac.horizon_record(per_orbit(N + 2, N + 1), N)
        self.assertEqual(status, "complete")
        self.assertEqual(m["apo_decay_km"], 40.0)
        self.assertEqual(ac.horizon_record(per_orbit(N + 1, N), N)[0], "complete")
        # The old writer's shape for --orbits=41: apoapses truncated to the passes.
        self.assertEqual(ac.horizon_record(per_orbit(N, N), N)[0], "early_termination")
        self.assertEqual(ac.horizon_record(per_orbit(2, 2), N)[0], "early_termination")
        po = per_orbit(N + 2, N + 1)
        po.loc[5, "pass_heat_load_Jcm2"] = np.nan
        self.assertEqual(ac.horizon_record(po, N)[0], "invalid_output")
        po = per_orbit(N + 2, N + 1)
        po.loc[3, "peri_time_s"] = po.loc[4, "apo_time_s"] + 1.0
        self.assertEqual(ac.horizon_record(po, N)[0], "invalid_output")

    def test_horizon_metrics_use_only_the_horizon(self):
        # Pass N+1 has the largest heat, the deepest periapsis: outside the horizon.
        heat = np.arange(1.0, N + 3)
        peri = np.full(N + 2, 100.0)
        peri[N] = 50.0
        status, _, m = ac.horizon_record(per_orbit(N + 2, N + 1, heat=heat, peri=peri), N)
        self.assertEqual(status, "complete")
        self.assertEqual(m["heat_load_horizon_Jcm2"], sum(range(1, N + 1)))
        self.assertEqual(m["peak_heat_rate_Wcm2"], float(N))
        self.assertEqual(m["peri_min_km"], 100.0)
        self.assertEqual(m["apo_final_km"], 30000.0 - N)


class CampaignAnalysis(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_main(self, res, tele):
        out = os.path.join(self.tmp.name, "out")
        with contextlib.redirect_stdout(io.StringIO()):
            ac.main([res, tele, out])
        return out

    def test_incomplete_member_is_excluded_and_counted(self):
        # The review's counterexample: a two-pass Terminated member beside a
        # complete one must not enter the 40-pass spread.
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {
            101: (dict(), per_orbit(N + 2, N + 1)),
            102: (dict(termination_cause="terminated_before_orbit_count", completed_orbits=2), per_orbit(2, 2)),
        })
        out = self.run_main(res, tele)
        spread = pd.read_csv(os.path.join(out, "spread.csv")).set_index("quantity")
        self.assertEqual(spread.loc["apo_decay_km", "n"], 1)
        self.assertEqual(spread.loc["apo_decay_km", "n_members"], 2)
        self.assertEqual(spread.loc["apo_decay_km", "n_excluded"], 1)
        self.assertEqual(spread.loc["apo_decay_km", "mean"], 40.0)
        self.assertEqual(spread.loc["heat_load_horizon_Jcm2", "mean"], 40.0)
        self.assertTrue((spread.horizon_passes == N).all())
        summ = pd.read_csv(os.path.join(out, "per_sample_summary.csv")).set_index("seed")
        self.assertEqual(summ.loc[102, "status"], "early_termination")
        self.assertEqual(summ.loc[102, "termination_cause"], "terminated_before_orbit_count")
        self.assertTrue(np.isnan(summ.loc[102, "apo_decay_km"]))
        hist = pd.read_csv(os.path.join(out, "per_orbit_histories.csv"))
        self.assertEqual(sorted(hist.seed.unique()), [101])
        self.assertEqual(len(hist), N + 1)
        excluded = pd.read_csv(os.path.join(out, "excluded_member_histories.csv"))
        self.assertEqual(sorted(excluded.seed.unique()), [102])
        r = pd.read_csv(os.path.join(out, "r_stats.csv")).iloc[0]
        self.assertEqual((r.members, r.passes), (1, N))
        counts = pd.read_csv(os.path.join(out, "member_status_counts.csv")).set_index("status")["count"]
        self.assertEqual(counts.to_dict(), dict(complete=1, early_termination=1, invalid_output=0,
                                                solver_failure=0, dispatch_failure=0))
        facts = open(os.path.join(out, "campaign_facts.txt")).read()
        self.assertIn("analysis_n_complete = 1", facts)
        self.assertIn("analysis_n_early_termination = 1", facts)
        self.assertIn("analysis_n_members = 2", facts)

    def test_failures_are_reported_separately(self):
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {
            101: (dict(), per_orbit(N + 2, N + 1)),
            102: (dict(retcode="ERROR", error="DomainError(-1.0)", termination_cause="error"), None),
            103: (dict(retcode="Unstable"), per_orbit(N + 2, N + 1)),
            104: (dict(dispatch_success=False, error="ProcessExitedException(3)", retcode=np.nan), None),
            105: (dict(), None),
        })
        out = self.run_main(res, tele)
        summ = pd.read_csv(os.path.join(out, "per_sample_summary.csv")).set_index("seed")
        self.assertEqual(summ.status.to_dict(), {101: "complete", 102: "solver_failure", 103: "solver_failure",
                                                 104: "dispatch_failure", 105: "invalid_output"})
        self.assertEqual(summ.loc[105, "status_reason"], "per_orbit.csv is missing")
        spread = pd.read_csv(os.path.join(out, "spread.csv"))
        self.assertTrue((spread.n == 1).all() and (spread.n_members == 5).all())
        facts = open(os.path.join(out, "campaign_facts.txt")).read()
        for line in ("analysis_n_solver_failure = 2", "analysis_n_dispatch_failure = 1",
                     "analysis_n_invalid_output = 1", "analysis_n_complete = 1", "route = threads"):
            self.assertIn(line, facts)

    def test_short_nominal_stops_before_writing(self):
        # The old default (--orbits=41 with apoapses truncated) gave the nominal
        # 40 apoapses: the 40-pass analysis needs 41.
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {101: (dict(), per_orbit(N + 2, N + 1))},
                                   nominal=per_orbit(N, N))
        out = os.path.join(self.tmp.name, "out")
        with self.assertRaisesRegex(ValueError, r"nominal run is early_termination.*--orbits=42"):
            ac.main([res, tele, out])
        self.assertFalse(os.path.exists(out))

    def test_short_flight_series_stops_before_writing(self):
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {101: (dict(), per_orbit(N + 2, N + 1))},
                                   flight_orbits=N)
        out = os.path.join(self.tmp.name, "out")
        with self.assertRaisesRegex(ValueError, "flight telemetry does not cover the 40-pass horizon"):
            ac.main([res, tele, out])
        self.assertFalse(os.path.exists(out))

    def test_no_complete_member_stops_before_writing(self):
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {101: (dict(), per_orbit(2, 2))})
        out = os.path.join(self.tmp.name, "out")
        with self.assertRaisesRegex(ValueError, "no complete member"):
            ac.main([res, tele, out])
        self.assertFalse(os.path.exists(out))

    def test_flight_comparison_over_the_horizon(self):
        res, tele = write_campaign(os.path.join(self.tmp.name, "c"), {101: (dict(), per_orbit(N + 2, N + 1))})
        out = self.run_main(res, tele)
        spread = pd.read_csv(os.path.join(out, "spread.csv")).set_index("quantity")
        self.assertEqual(spread.loc["apo_decay_km", "flight"], 1.5 * N)
        self.assertEqual(spread.loc["apo_final_km", "flight"], 30000.0 - 1.5 * N)
        self.assertEqual(spread.loc["peri_min_km", "flight"], 99.0)
        flight = pd.read_csv(os.path.join(out, "flight_series.csv"))
        self.assertEqual(len(flight), N + 1)
        self.assertTrue(np.isnan(flight.peri_km.iloc[N]))


if __name__ == "__main__":
    unittest.main()
