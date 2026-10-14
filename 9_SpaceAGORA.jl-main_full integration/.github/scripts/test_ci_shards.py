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
import textwrap
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


def job_env_runner_references(workflow):
    """Find runner expressions in this repository's job-level env blocks.

    GitHub permits runner in step contexts, but not jobs.<job_id>.env.
    Keep this bounded check stdlib-only, like the other workflow wiring tests;
    it supplements YAML parsing rather than claiming full Actions validation.
    """
    invalid = []
    for job in re.finditer(r"^  ([\w-]+):\n(.*?)(?=^  [\w-]+:|\Z)", workflow, re.M | re.S):
        for block in re.finditer(r"^    env:([^\n]*(?:\n(?: {6,}[^\n]*|[ \t]*))*)", job.group(2), re.M):
            value = "\n".join(line for line in block.group(1).splitlines()
                              if not line.lstrip().startswith("#"))
            for expression in re.findall(r"\$\{\{(.*?)\}\}", value, re.S):
                if re.search(r"\brunner\s*(?:\.|\[)", expression):
                    invalid.append(job.group(1))
    return invalid


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
            expected_needs = ["plan", job]
            if aggregator in ("tests", "coverage-quality-gate"):
                expected_needs.append("optional-hypr-installation")
                self.assertIn("needs.optional-hypr-installation.result", agg.group(1))
                if aggregator == "tests":
                    expected_needs.append("optional-mujoco-installation")
                    self.assertIn("needs.optional-mujoco-installation.result", agg.group(1))
                expected_needs.append("optional-visualizer-installation")
                self.assertIn("needs.optional-visualizer-installation.result", agg.group(1))
                expected_needs.append("optional-edg-installation")
                self.assertIn("needs.optional-edg-installation.result", agg.group(1))
            needs = re.search(r"needs: \[(.*?)\]", agg.group(1))
            self.assertIsNotNone(needs, aggregator)
            self.assertEqual([x.strip() for x in needs.group(1).split(",")], expected_needs)
            self.assertIn(f"SHARD_MATRIX: ${{{{ needs.plan.outputs.{output} }}}}", agg.group(1), aggregator)

    def test_edg_preparation_fails_before_export_on_installation_failure(self):
        action = (ROOT / ".github/actions/setup/action.yml").read_text()
        step = action.split("    - name: Prepare mandatory EDG provider for runtime tests\n", 1)[1].split("    - name:", 1)[0]
        script = textwrap.dedent(step.split("      run: |\n", 1)[1]).strip()
        with tempfile.TemporaryDirectory(prefix="edg workflow ") as temporary:
            root = Path(temporary)
            executable = root / "julia"
            executable.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$STUB_JULIA_ARGS"\nexit "$STUB_JULIA_STATUS"\n')
            executable.chmod(0o755)
            exported, arguments = root / "GitHub environment", root / "Julia arguments"
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ["PATH"],
                       RUNNER_TEMP=str(root), GITHUB_ENV=str(exported), STUB_JULIA_ARGS=str(arguments))
            for status in (0, 23):
                exported.write_text("PREVIOUS=value\n")
                result = subprocess.run(["bash", "-c", script], env=dict(env, STUB_JULIA_STATUS=str(status)),
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, status, result.stderr)
                self.assertEqual(arguments.read_text().splitlines(), ["--startup-file=no", "--project=@stdlib",
                                 "scripts/setup_edg.jl", str(root / "spaceagora-edg-env")])
                expected = "PREVIOUS=value\n"
                if status == 0:
                    expected += f"SPACEAGORA_EDG_TEST_ENV={root / 'spaceagora-edg-env'}\n"
                    expected += f"SPACEAGORA_EDG_ENV={root / 'spaceagora-edg-env'}\n"
                    expected += f"JULIA_LOAD_PATH=@:{root / 'spaceagora-edg-env'}:@stdlib\n"
                self.assertEqual(exported.read_text(), expected)

    def test_job_environment_rejects_runner_context(self):
        valid = """jobs:
  sample:
    env:
      ROOT: ${{ github.workspace }}
    steps:
      - run: echo '${{ runner.temp }}'
        env:
          STEP_ROOT: ${{ runner.temp }}
        with:
          path: ${{ runner.temp }}
"""
        self.assertEqual(job_env_runner_references(valid), [])
        for value in ("${{ runner.temp }}", "${{ runner['temp'] }}", "|\n        ${{ runner.temp }}"):
            invalid = valid.replace("ROOT: ${{ github.workspace }}", "ROOT: " + value)
            self.assertEqual(job_env_runner_references(invalid), ["sample"])
        inline = valid.replace("env:\n      ROOT: ${{ github.workspace }}", "env: {ROOT: '${{ runner.temp }}'}")
        self.assertEqual(job_env_runner_references(inline), ["sample"])
        for workflow in sorted((ROOT / ".github/workflows").glob("*.yml")):
            self.assertEqual(job_env_runner_references(workflow.read_text()), [], workflow.name)

    def test_visualizer_preparation_exports_only_a_successful_environment(self):
        wf = (ROOT / ".github/workflows/julia-ci.yml").read_text()
        job = re.search(r"^  optional-visualizer-installation:\n(.*?)(?=^  \S|\Z)", wf, re.S | re.M).group(1)
        preparation = re.search(r"^      - name: Prepare independent core and pinned renderer projects\n"
                                r"        run: \|\n((?:          [^\n]*\n)+)", job, re.M)
        self.assertIsNotNone(preparation)
        script = textwrap.dedent(preparation.group(1))
        with tempfile.TemporaryDirectory(prefix="visualizer workflow ") as temporary:
            root = Path(temporary)
            executable = root / "julia"
            executable.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$STUB_JULIA_ARGS"\nexit "$STUB_JULIA_STATUS"\n')
            executable.chmod(0o755)
            exported, arguments = root / "GitHub environment", root / "Julia arguments"
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ["PATH"],
                       RUNNER_TEMP=str(root), GITHUB_ENV=str(exported), STUB_JULIA_ARGS=str(arguments))
            for status in (0, 23):
                exported.write_text("PREVIOUS=value\n")
                result = subprocess.run(["bash", "-c", script], env=dict(env, STUB_JULIA_STATUS=str(status)),
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, status, result.stderr)
                self.assertEqual(arguments.read_text().splitlines(), ["--startup-file=no", "--project=.",
                                 "test/optional_visualizer/setup.jl", str(root / "visualizer-installations")])
                expected = "PREVIOUS=value\n"
                if status == 0:
                    expected += f"SPACEAGORA_VISUALIZER_TEST_ENV={root / 'visualizer-installations/visualizer'}\n"
                self.assertEqual(exported.read_text(), expected)

    def test_external_mujoco_companion_is_mandatory(self):
        wf = (ROOT / ".github/workflows/julia-ci.yml").read_text()
        match = re.search(r"^  optional-mujoco-installation:\n(.*?)(?=^  \S|\Z)", wf, re.S | re.M)
        self.assertIsNotNone(match)
        job = match.group(1)
        # An external checkout must stay pinned; testing its default core pin
        # would fail to qualify the SpaceAGORA revision under review.
        checkout = re.search(r"repository: Space-FALCON-Lab/SpaceAGORAMuJoCo\.jl\n"
                             r"(?P<settings>(?:          [^\n]*\n)+)", job)
        self.assertIsNotNone(checkout)
        self.assertIn("ref: ${{ steps.mujoco-source.outputs.revision }}", checkout.group("settings"))
        self.assertIn("git rev-parse HEAD:packages/robotics/SpaceAGORAMuJoCo.jl", job)
        self.assertIn('[[ "$revision" =~ ^[0-9a-f]{40}$ ]] || exit 1', job)
        self.assertIn("scripts/optional_packages.jl check mujoco", job)
        self.assertIn('scripts/optional_packages.jl setup mujoco "$RUNNER_TEMP/mujoco-env"', job)
        self.assertIn("path: packages/robotics/SpaceAGORAMuJoCo.jl", checkout.group("settings"))
        self.assertIn("persist-credentials: false", checkout.group("settings"))
        self.assertIn("SPACEAGORA_CORE_PATH: ${{ github.workspace }}", job)
        self.assertIn('SPACEAGORA_CORE_REV: ""', job)
        self.assertNotIn("packages/SpaceAGORAMuJoCo", job)
        self.assertIn("hashFiles('packages/robotics/SpaceAGORAMuJoCo.jl/Artifacts.toml')", job)
        for required in ('SPACEAGORA_MUJOCO_REQUIRED: "1"', 'JULIA_LOAD_PATH: "@:@stdlib"',
                         'JULIA_NUM_THREADS: "4"', '@assert Threads.nthreads() == 4'):
            self.assertIn(required, job)
        # Every isolation/load-order/native probe is a separate fresh process.
        commands = re.findall(r"^        run: (julia [^\n]+)$", job, re.M)
        for script in ("test/installation_plan_tests.jl", "test/source_acquisition_tests.jl",
                       "test/core_only.jl", "test/incompatible_core.jl",
                       "test/load_order.jl core-first", "test/load_order.jl companion-first",
                       "test/runtests.jl"):
            self.assertEqual(sum(f"packages/robotics/SpaceAGORAMuJoCo.jl/{script}" in command
                                 for command in commands), 1, script)
        self.assertNotRegex(job, r"(?m)^\s*(if|continue-on-error):")
        self.assertNotIn("|| true", job)
        self.assertIn("git -C packages/robotics/SpaceAGORAMuJoCo.jl diff --exit-code", job)


    def test_domain_packages_use_selected_setup(self):
        wf = (ROOT / ".github/workflows/julia-ci.yml").read_text()
        job = re.search(r"^  optional-hypr-installation:\n(.*?)(?=^  \S|\Z)", wf, re.S | re.M).group(1)
        self.assertIn("test/optional_packages_tests.jl", job)
        setup = (ROOT / "test/optional_hypr/setup.jl").read_text()
        self.assertIn('OptionalPackages.setup("hypr",', setup)
        self.assertIn('with_hypr=false', setup)
        self.assertNotIn("submodules: recursive", job)

    def test_external_visualizer_is_mandatory_and_core_isolated(self):
        wf = (ROOT / ".github/workflows/julia-ci.yml").read_text()
        match = re.search(r"^  optional-visualizer-installation:\n(.*?)(?=^  \S|\Z)", wf, re.S | re.M)
        self.assertIsNotNone(match)
        job = match.group(1)
        pin = (ROOT / "packages/VisualizerSource.toml").read_text()
        self.assertRegex(pin, r'(?m)^rev = "[0-9a-f]{40}"$')
        self.assertIn('url = "https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl.git"', pin)
        self.assertIn('JULIA_LOAD_PATH: "@:@stdlib"', job)
        self.assertIn('SPACEAGORA_VISUALIZER_PATH: ""', job)
        for probe in ("installation_plan_tests.jl", "setup.jl", "core.jl",
                      "optional.jl core-first", "optional.jl visualizer-first"):
            self.assertEqual(job.count("test/optional_visualizer/" + probe), 1, probe)
        self.assertIn("Pkg.test(\"SpaceAGORAVisualizer\")", job)
        self.assertIn("test/unit/visualization/terrain_export_tests.jl", job)
        self.assertIn('"$terrain_fixture_dir/fixtures.json" "$SPACEAGORA_VISUALIZER_DIR"', job)
        self.assertIn("viewer_channels_models_tests.mjs", job)
        self.assertIn("test/scripts", job)
        self.assertNotIn("continue-on-error:", job)
        self.assertNotIn("|| true", job)
        self.assertNotRegex(job, r"(?m)^    if:")
        self.assertEqual(re.findall(r"(?m)^        if: (.*)$", job), ["always()", "always()"])
        coverage = re.search(r"^  coverage-quality-gate:\n(.*?)(?=^  \S)", wf, re.S | re.M).group(1)
        self.assertIn("name: visualizer-installation-coverage", coverage)
        self.assertIn("tar -xzmf visualizer-installation-coverage/coverage.tgz", coverage)
        setup = (ROOT / ".github/actions/setup/action.yml").read_text()
        self.assertNotIn("setup_visualizer", setup)
        self.assertNotIn("SpaceAGORAVisualizer", setup)



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
