"""Tests for the perturbed-density comparison analyzer
(benchmarks/studies/telemetry_validation/perturbed_density_modes/analyze.py) on
fabricated run directories written the way run_mode.jl and attempt_record.jl
write them. The numbers are invented; none is a measurement. Needs numpy,
pandas, pyarrow (feather) and tomllib or tomli."""
import contextlib
import hashlib
import importlib.util
import io
import os
import pathlib
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[2]
TOOL = REPO / "benchmarks" / "studies" / "telemetry_validation" / "perturbed_density_modes" / "analyze.py"
spec = importlib.util.spec_from_file_location("perturbed_modes_analyze", TOOL)
an = importlib.util.module_from_spec(spec)
spec.loader.exec_module(an)

ORBITS = 5


def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    return '"%s"' % v


def write_summary(d, fields, hashes=None):
    lines = ["%s = %s" % (k, toml_value(v)) for k, v in fields.items()]
    if hashes is not None:
        lines.append("")
        lines.append("[output_sha256]")
        lines += ['"%s" = "%s"' % (k, v) for k, v in sorted(hashes.items())]
    (d / "run_summary.toml").write_text("\n".join(lines) + "\n")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_run(res, tag, apo0=30000.0, n_apo=ORBITS - 1, n_peri=ORBITS - 1, status="complete", orbits=ORBITS,
              trajectory="time,x\n0,1\n", record=True):
    """A run tag directory; with record=True the summary names its outputs, as
    finish_attempt! writes it. status=None writes a summary from before the
    attempt record existed."""
    d = pathlib.Path(res) / tag
    d.mkdir(parents=True, exist_ok=True)
    (d / "simulation_results.csv").write_text(trajectory)
    pd.DataFrame(dict(
        event=["peri"] * n_peri + ["apo"] * n_apo,
        index=list(range(1, n_peri + 1)) + list(range(1, n_apo + 1)),
        flight_orbit=[19.0 + k for k in range(n_peri)] + [19.0 + k for k in range(n_apo)],
        time_s=[1000.0 * k + 500.0 for k in range(n_peri)] + [1000.0 * k for k in range(n_apo)],
        altitude_km=[100.0] * n_peri + [apo0 - k for k in range(n_apo)])).to_csv(d / "extrema.csv", index=False)
    fields = dict(tag=tag, mode="pass", gram_seed=11, tight=False, dt_max_atmosphere=0.2, reltol_orbit=1e-9,
                  retcode="Terminated", error="", naccept=10, nreject=0, nf=60, wall_s=1.0,
                  solver_sequence="Tsit5", commit="test", orbits_requested=orbits,
                  completed_orbits=orbits, termination_cause="orbit_count")
    if status is not None:
        fields.update(status=status, status_reason="", attempt_id="attempt-" + tag)
    hashes = {n: sha256(d / n) for n in ("extrema.csv", "simulation_results.csv")} if record else None
    write_summary(d, fields, hashes)
    return d


def write_failed_rerun(res, tag):
    """A rerun that failed before writing a trajectory, leaving the earlier
    attempt's files in place (what the guard-less producer did)."""
    d = pathlib.Path(res) / tag
    write_summary(d, dict(tag=tag, mode="pass", gram_seed=11, tight=False, dt_max_atmosphere=0.2,
                          reltol_orbit=1e-9, retcode="ERROR", error="solver failure", naccept=-1, nreject=-1,
                          nf=-1, wall_s=1.0, solver_sequence="", commit="test", orbits_requested=ORBITS,
                          status="failed", status_reason="the solve threw", attempt_id="attempt-2"), {})


def write_telemetry(tele):
    tele = pathlib.Path(tele)
    tele.mkdir(parents=True)
    orbits = 19.0 + np.arange(10)
    pd.DataFrame(dict(orbit=orbits, altitude=30000.0 - np.arange(10))).to_feather(
        tele / "True_Odyssey_apoapsis_alts_kernel.feather")
    pd.DataFrame(dict(orbit=orbits, altitude=np.full(10, 100.0))).to_feather(
        tele / "True_Odyssey_periapsis_alts_kernel.feather")
    return str(tele)


class PerturbedModesAnalysis(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = pathlib.Path(tmp.name)
        self.res = str(self.root / "res")
        self.tele = write_telemetry(self.root / "tele")
        write_run(self.res, "nominal")
        write_run(self.res, "B_s11", apo0=30001.0)

    def analyze(self):
        out = str(self.root / "out")
        with contextlib.redirect_stdout(io.StringIO()):
            an.main([self.res, self.tele, out])
        runs = pd.read_csv(os.path.join(out, "runs.csv")).set_index("tag")
        sens = pd.read_csv(os.path.join(out, "sensitivity.csv")).set_index("variant")
        return runs, sens

    def test_identical_rerun_agrees(self):
        write_run(self.res, "B_s11_rep", apo0=30001.0)
        runs, sens = self.analyze()
        self.assertTrue(sens.loc["B_s11_rep", "comparable"])
        self.assertTrue(sens.loc["B_s11_rep", "bit_identical"])
        self.assertEqual(sens.loc["B_s11_rep", "max_abs_d_apo_km"], 0.0)
        self.assertEqual(str(runs.loc["B_s11_rep", "bit_identical_to_base"]), "True")
        self.assertEqual(runs.loc["B_s11", "d_apo_final_km"], 1.0)

    def test_failed_rerun_cannot_report_agreement(self):
        # The review's regression: a successful tag rerun unsuccessfully, with the
        # earlier attempt's identical files still on disk.
        write_run(self.res, "B_s11_rep", apo0=30001.0)
        write_failed_rerun(self.res, "B_s11_rep")
        runs, sens = self.analyze()
        row = sens.loc["B_s11_rep"]
        self.assertFalse(row.comparable)
        self.assertIn("attempt failed", row.not_comparable_reason)
        self.assertTrue(pd.isna(row.bit_identical))
        self.assertTrue(np.isnan(row.max_abs_d_apo_km) and np.isnan(row.d_apo_final_km))
        self.assertFalse(runs.loc["B_s11_rep", "verified"])
        self.assertEqual(runs.loc["B_s11_rep", "bit_identical_to_base"], "not_comparable")
        per_orbit = pd.read_csv(self.root / "out" / "per_orbit.csv")
        self.assertNotIn("B_s11_rep", set(per_orbit.tag))

    def test_summary_without_attempt_record_is_not_compared(self):
        write_run(self.res, "B_s11_rep", apo0=30001.0, status=None, record=False)
        runs, sens = self.analyze()
        self.assertFalse(sens.loc["B_s11_rep", "comparable"])
        self.assertIn("predates", runs.loc["B_s11_rep", "unverified_reason"])

    def test_changed_output_is_not_compared(self):
        d = write_run(self.res, "B_s11_rep", apo0=30001.0)
        (d / "simulation_results.csv").write_text("time,x\n0,2\n")
        runs, sens = self.analyze()
        self.assertFalse(sens.loc["B_s11_rep", "comparable"])
        self.assertIn("does not match its recorded identity", runs.loc["B_s11_rep", "unverified_reason"])

    def test_short_coverage_is_not_compared(self):
        write_run(self.res, "B_s11_rep", apo0=30001.0, n_apo=2, n_peri=2)
        runs, sens = self.analyze()
        self.assertFalse(sens.loc["B_s11_rep", "comparable"])
        self.assertIn("2 apoapses and 2 periapses of 5 orbits", runs.loc["B_s11_rep", "unverified_reason"])

    def test_different_orbit_counts_are_not_compared(self):
        write_run(self.res, "B_s11_tight", apo0=30001.0, orbits=ORBITS + 1, n_apo=ORBITS, n_peri=ORBITS)
        _, sens = self.analyze()
        self.assertFalse(sens.loc["B_s11_tight", "comparable"])
        self.assertEqual(sens.loc["B_s11_tight", "not_comparable_reason"], "different requested orbit counts")

    def test_near_end_early_termination_cannot_report_agreement(self):
        # Four apoapses and periapses can precede the fifth requested apoapsis.
        # Hashes and "complete" alone must not override the runtime counter.
        d = write_run(self.res, "B_s11_rep", apo0=30001.0)
        fields = an.load_toml(d / "run_summary.toml")
        hashes = fields.pop("output_sha256")
        fields.update(completed_orbits=ORBITS - 1, termination_cause="terminated_before_orbit_count")
        write_summary(d, fields, hashes)
        runs, sens = self.analyze()
        self.assertFalse(runs.loc["B_s11_rep", "verified"])
        self.assertIn("4 completed orbit events for 5 requested", runs.loc["B_s11_rep", "unverified_reason"])
        self.assertFalse(sens.loc["B_s11_rep", "comparable"])
        self.assertTrue(pd.isna(sens.loc["B_s11_rep", "bit_identical"]))
        self.assertTrue(np.isnan(sens.loc["B_s11_rep", "d_apo_final_km"]))

    def test_missing_or_malformed_completion_metadata_is_rejected(self):
        d = write_run(self.res, "B_s11_rep", apo0=30001.0)
        original = an.load_toml(d / "run_summary.toml")
        invalid = [
            ("completed_orbits", None), ("completed_orbits", "5"),
            ("completed_orbits", 5.0), ("completed_orbits", True), ("completed_orbits", -1),
            ("orbits_requested", None), ("orbits_requested", 0), ("orbits_requested", "5"),
            ("orbits_requested", True), ("termination_cause", None),
            ("termination_cause", "terminated_unknown"), ("termination_cause", "end_of_time_span"),
            ("termination_cause", "terminated_before_orbit_count"), ("retcode", "ERROR"),
            ("retcode", "Success"), ("error", "solver failure"), ("error", None),
        ]
        for key, value in invalid:
            with self.subTest(key=key, value=value):
                fields = original.copy()
                if value is None:
                    fields.pop(key)
                else:
                    fields[key] = value
                self.assertFalse(an.verify_run(str(d), fields)[0])
        complete_at_time_limit = dict(original, retcode="Success", termination_cause="end_of_time_span")
        self.assertTrue(an.verify_run(str(d), complete_at_time_limit)[0])
        short_at_time_limit = dict(complete_at_time_limit, completed_orbits=ORBITS - 1)
        self.assertFalse(an.verify_run(str(d), short_at_time_limit)[0])

    def test_unverified_nominal_stops_the_analysis(self):
        write_failed_rerun(self.res, "nominal")
        with self.assertRaisesRegex(ValueError, "nominal run is not verified"):
            an.main([self.res, self.tele, str(self.root / "out")])
        self.assertFalse((self.root / "out").exists())


if __name__ == "__main__":
    unittest.main()
