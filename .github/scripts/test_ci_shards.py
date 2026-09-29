#!/usr/bin/env python3
"""Self-test of the CI shard plan, its matrix shape and its verifier.

Run by the `plan` job before the matrices are used. It checks that:
  - every test item and example is planned into exactly one shard per mode;
  - each matrix is `{"include": [...]}` and nothing else, so GitHub makes one
    job per entry (an `include` entry next to another matrix key is merged
    into the existing combinations instead, which once silently dropped three
    suites);
  - julia-ci.yml feeds each shard job exactly that planner output, and each
    aggregator waits for its shard job even when it fails;
  - ci_shard_verify.py rejects a missing shard, a missing item and an item run
    twice, and accepts a complete set.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import ci_shard_plan  # noqa: E402


def plan():
    out = subprocess.run([sys.executable, str(HERE / "ci_shard_plan.py")], capture_output=True,
                         text=True, check=True, env={k: v for k, v in os.environ.items() if k != "GITHUB_OUTPUT"})
    return {k: json.loads(v) for k, v in (line.split("=", 1) for line in out.stdout.splitlines())}


class PlanTests(unittest.TestCase):
    def test_every_item_once(self):
        p = plan()
        suites = ci_shard_plan.suite_files()
        tests = ([f"suite:{s[:2]}" for s in suites if s != "09_probe_drivers.jl"]
                 + [f"probe:{x}" for x in ci_shard_plan.probe_files()]
                 + [f"unit:{x}" for x in ci_shard_plan.unit_files()])
        for mode, extra in (("plain", ["grid", "depwarn"]), ("coverage", ["grid"])):
            planned = [i for e in p[f"{mode}_matrix"]["include"] for i in e["items"].split()]
            self.assertEqual(sorted(planned), sorted(tests + extra), mode)
        planned = [i for e in p["examples_matrix"]["include"] for i in e["items"].split()]
        self.assertEqual(sorted(planned), ci_shard_plan.example_files())

    def test_matrix_shape(self):
        for name, matrix in plan().items():
            self.assertEqual(list(matrix), ["include"], name)
            shards = [e["shard"] for e in matrix["include"]]
            self.assertEqual(len(shards), len(set(shards)), name)
            self.assertTrue(all(e["items"].split() for e in matrix["include"]), name)

    def test_workflow_wiring(self):
        wf = (ROOT / ".github/workflows/julia-ci.yml").read_text()
        for job, output, aggregator in (("test-shards", "plain_matrix", "tests"),
                                        ("coverage-shards", "coverage_matrix", "coverage-quality-gate"),
                                        ("example-shards", "examples_matrix", "example-suite-smoke")):
            body = re.search(rf"^  {job}:\n(.*?)(?=^  \S)", wf, re.S | re.M)
            self.assertIsNotNone(body, job)
            self.assertIn(f"matrix: ${{{{ fromJSON(needs.plan.outputs.{output}) }}}}", body.group(1), job)
            agg = re.search(rf"^  {aggregator}:\n(.*?)(?=^  \S)", wf, re.S | re.M)
            self.assertIsNotNone(agg, aggregator)
            self.assertIn("if: always()", agg.group(1), aggregator)
            self.assertRegex(agg.group(1), rf"needs: \[plan, {job}\]", aggregator)
            self.assertIn(f"SHARD_MATRIX: ${{{{ needs.plan.outputs.{output} }}}}", agg.group(1), aggregator)


class VerifyTests(unittest.TestCase):
    MATRIX = {"include": [{"shard": "t-1", "items": "suite:01 unit:a.jl grid"},
                          {"shard": "t-2", "items": "probe:p.jl"}]}
    UNIVERSE = ["probe:p.jl", "suite:01", "unit:a.jl"]

    def write(self, root, shard, name, attempted, universe=None, requested=None):
        d = root / shard
        d.mkdir(parents=True, exist_ok=True)
        (d / "job.ok").touch()
        uni = universe or self.UNIVERSE
        req = requested or next(e["items"] for e in self.MATRIX["include"] if e["shard"] == shard).split()
        lines = [f"universe = {json.dumps(uni)}", f"requested = {json.dumps(req)}",
                 f"attempted = {json.dumps(attempted)}", "failed = []", "[timings]"]
        (d / f"{name}.toml").write_text("\n".join(lines) + "\n")

    def verify(self, root):
        env = dict(os.environ, SHARD_MATRIX=json.dumps(self.MATRIX))
        return subprocess.run([sys.executable, str(HERE / "ci_shard_verify.py"), "tests", str(root)],
                              env=env, capture_output=True, text=True)

    def complete(self, root):
        self.write(root, "t-1", "integration", ["suite:01"])
        self.write(root, "t-1", "unit", ["unit:a.jl"], universe=["unit:a.jl"])
        (root / "t-1" / "grid.done").touch()
        self.write(root, "t-2", "integration", ["probe:p.jl"])

    def test_complete_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.complete(Path(tmp))
            r = self.verify(Path(tmp))
            self.assertEqual(r.returncode, 0, r.stdout)

    def test_missing_shard_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.complete(Path(tmp))
            for f in (Path(tmp) / "t-2").iterdir():
                f.unlink()
            self.assertNotEqual(self.verify(Path(tmp)).returncode, 0)

    def test_missing_item_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.complete(Path(tmp))
            self.write(Path(tmp), "t-2", "integration", [])
            self.assertNotEqual(self.verify(Path(tmp)).returncode, 0)

    def test_duplicate_item_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.complete(Path(tmp))
            self.write(Path(tmp), "t-2", "integration", ["probe:p.jl", "suite:01"])
            self.assertNotEqual(self.verify(Path(tmp)).returncode, 0)

    def test_workflow_item_not_done_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.complete(Path(tmp))
            (Path(tmp) / "t-1" / "grid.done").unlink()
            self.assertNotEqual(self.verify(Path(tmp)).returncode, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
