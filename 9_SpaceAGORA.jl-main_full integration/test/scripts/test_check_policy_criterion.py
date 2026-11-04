"""Tests for scripts/check_policy_criterion.py's full-budget (P5f) rule, on
fabricated raw CSVs. The numbers below are invented to exercise the rule; none
is a measurement."""
import importlib.util
import io
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[2]
TOOL = REPO / "scripts" / "check_policy_criterion.py"
spec = importlib.util.spec_from_file_location("check_policy_criterion", TOOL)
cpc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cpc)

BUDGET = 12
SPLITS = [(w, BUDGET // w) for w in range(1, BUDGET + 1) if BUDGET % w == 0]
STATIC = ["outer_threads", "outer_process", "outer_inner_static"]
CASE = "mcgrid_16sat_8mc"


def rows(phase, case, mode, w, t, times, success=True, allocation=""):
    return [
        {"phase_id": phase, "case": case, "mode": mode, "thread_count": t,
         "process_workers": w, "mc_samples": 8, "repeat": i + 1,
         "wall_time_s": x, "success": success, "adaptive_allocation": allocation}
        for i, x in enumerate(times)
    ]


def p5f_frame(static_time, adaptive_times, adaptive_success=True, extra=()):
    """Every static route at every split of BUDGET, serial at every split, and
    the adaptive mode once at the full budget (BUDGET workers x BUDGET threads)."""
    out = []
    for w, t in SPLITS:
        out += rows("P5f", CASE, "serial", w, t, [10.0, 10.0, 10.0])
        for m in STATIC:
            x = static_time(m, w, t)
            out += rows("P5f", CASE, m, w, t, [x, x * 1.01, x * 0.99])
    out += rows("P5f", CASE, "predictive", BUDGET, BUDGET, adaptive_times,
                success=adaptive_success, allocation="threads:w8+l0:b1")
    out += list(extra)
    return pd.DataFrame(out)


def spread(m, w, t):
    # The fastest static allocation is outer_threads at 2x6, 2.0 s; everything
    # else is at least 25% slower, so the bias correction takes the minimum.
    if (m, w, t) == ("outer_threads", 2, 6):
        return 2.0
    return 2.5 + 0.1 * w + (0.2 if m == "outer_process" else 0.0)


class FullBudgetRule(unittest.TestCase):
    def evaluate(self, df, adaptive="predictive"):
        with redirect_stdout(io.StringIO()):
            return cpc.evaluate_run(df, 0.10, 0.03, adaptive)

    def test_baseline_is_best_static_over_every_split(self):
        recs = self.evaluate(p5f_frame(spread, [2.1, 2.15, 2.05]))
        self.assertEqual(len(recs), 1)
        r = recs[0]
        self.assertEqual(r["phase"], "P5f")
        self.assertAlmostEqual(r["baseline"], 2.0)
        self.assertEqual(r["baseline_kind"], "min: outer_threads 2x6")
        self.assertAlmostEqual(r["ratio"], 2.1 / 2.0)
        self.assertEqual(r["verdict"], "PASS")
        # The axis names the full budget and the allocation the row recorded.
        self.assertEqual(r["axis"], "full 12x12 [threads:w8+l0:b1]")

    def test_slower_than_tolerance_fails(self):
        r = self.evaluate(p5f_frame(spread, [2.4, 2.4, 2.4]))[0]
        self.assertEqual(r["verdict"], "FAIL")
        self.assertTrue(r["reason"].startswith("ratio"))

    def test_tied_statics_take_the_mean(self):
        # Every (route, split) median inside the 3% band: the baseline is their
        # mean, as for any other point.
        tied = lambda m, w, t: 2.0 + 0.001 * w
        r = self.evaluate(p5f_frame(tied, [2.0, 2.0, 2.0]))[0]
        medians = [2.0 + 0.001 * w for w, _ in SPLITS for _ in STATIC]
        self.assertAlmostEqual(r["baseline"], sum(medians) / len(medians))
        self.assertTrue(r["baseline_kind"].startswith("mean of 18"))

    def test_failed_campaign_fails(self):
        r = self.evaluate(p5f_frame(spread, [2.0, 2.0, 2.0], adaptive_success=False))[0]
        self.assertEqual(r["verdict"], "FAIL")
        self.assertTrue(r["reason"].startswith("failed campaigns"))

    def test_statics_at_another_budget_are_ignored(self):
        # A faster static route at a different total budget is not a split of
        # this one and must not become the baseline.
        extra = rows("P5f", CASE, "outer_threads", 2, 12, [0.5, 0.5, 0.5])
        r = self.evaluate(p5f_frame(spread, [2.1, 2.1, 2.1], extra=extra))[0]
        self.assertAlmostEqual(r["baseline"], 2.0)

    def test_p5f_static_rows_are_not_scored_per_split(self):
        # P5f's per-split points hold no adaptive row, and P5 in the same run is
        # still scored point by point at its own splits.
        p5 = []
        for m in STATIC:
            p5 += rows("P5", CASE, m, 2, 6, [3.0, 3.0, 3.0])
        p5 += rows("P5", CASE, "predictive", 2, 6, [3.1, 3.1, 3.1])
        recs = self.evaluate(p5f_frame(spread, [2.1, 2.1, 2.1], extra=p5))
        self.assertEqual(sorted(r["phase"] for r in recs), ["P5", "P5f"])
        p5rec = next(r for r in recs if r["phase"] == "P5")
        self.assertEqual(p5rec["axis"], "2x6")

    def test_policy_v2_has_no_p5f_point(self):
        self.assertEqual(self.evaluate(p5f_frame(spread, [2.1]), adaptive="policy_v2"), [])

    def test_cli_exit_codes(self):
        with tempfile.TemporaryDirectory() as d:
            for times, code in (([2.1, 2.1, 2.1], 0), ([2.5, 2.5, 2.5], 1)):
                path = os.path.join(d, "paper_benchmarks_raw_20260101_000000.csv")
                p5f_frame(spread, times).to_csv(path, index=False)
                res = subprocess.run([sys.executable, str(TOOL), d, "--adaptive", "predictive"],
                                     capture_output=True, text=True)
                self.assertEqual(res.returncode, code, res.stdout + res.stderr)
                self.assertIn("full 12x12 [threads:w8+l0:b1]", res.stdout)


if __name__ == "__main__":
    unittest.main()
