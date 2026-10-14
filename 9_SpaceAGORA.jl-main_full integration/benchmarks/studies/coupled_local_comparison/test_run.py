"""Bounded runner tests: mock locked points, never launch Julia or a workload."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "source"
        harness = self.source / "benchmarks/studies/coupled_local_comparison"
        harness.mkdir(parents=True)
        self.runner = harness / "run.sh"
        shutil.copy2(Path(__file__).with_name("run.sh"), self.runner)
        self.runner.chmod(0o755)
        for tree in (self.source, self.root / "pre", self.root / "post"):
            files = tree / "benchmarks/studies/parallelization_performance"
            files.mkdir(parents=True)
            for name in ("cases.jl", "cli.jl", "execution.jl"):
                (files / name).write_text("# e6_actuated_saved mock harness\n")
        self.bin = self.root / "bin"
        self.bin.mkdir()
        # GNU realpath -m is not available on every developer host.
        realpath = self.bin / "realpath"
        realpath.write_text("#!/bin/sh\n[ \"$1\" != -m ] || shift\nprintf '%s\\n' \"$1\"\n")
        realpath.chmod(0o755)
        date = self.bin / "date"
        date.write_text("#!/bin/sh\nprintf '%s\\n' '2026-10-12T00:00:00Z'\n")
        date.chmod(0o755)
        flock = self.bin / "flock"
        flock.write_text("""#!/bin/sh
shift
mode=$6
printf '%s\n' "$mode" >> "$MOCK_CALLS"
if [ "${MOCK_WRITE_ARTIFACTS:-0}" = 1 ]; then
    dir=$4
    mkdir -p "$dir/states" "$dir/worker_rows"
    printf '%s\n' "$MOCK_ATTEMPT" > "$dir/states/state_c1.bin"
    if [ "$MOCK_ATTEMPT" = failed ]; then
        printf 'stale high counter\n' > "$dir/states/state_c99.bin"
    fi
    printf '%s\n' "$MOCK_ATTEMPT" > "$dir/parallelization_performance_raw_$MOCK_ATTEMPT.csv"
    printf '%s\n' "$MOCK_ATTEMPT" > "$dir/worker_rows/$MOCK_ATTEMPT.csv"
    printf 'attempt=%s\n' "$MOCK_ATTEMPT"
fi
if [ "$mode" = serial ]; then exit "$MOCK_SERIAL_RC"; fi
exit 0
""")
        flock.chmod(0o755)
        self.out = self.root / "out"
        self.calls = self.root / "calls"
        self.env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}",
                        PRE_TREE=str(self.root / "pre"), POST_TREE=str(self.root / "post"),
                        LOCKS=str(self.root / "mock.lock"), SWEEP_N="2", SWEEP_THREADS="2",
                        ROUTES="rhs_flat", QUIET_WAIT_S="0", QUIET_RETRY_S="0",
                        MOCK_CALLS=str(self.calls), MOCK_SERIAL_RC="0")

    def run_stage(self, rc):
        self.env["MOCK_SERIAL_RC"] = str(rc)
        return subprocess.run(["bash", str(self.runner), "B", str(self.out)],
                              env=self.env, text=True, capture_output=True, timeout=10)

    def test_success_and_resume(self):
        first = self.run_stage(0)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        self.assertIn("[done]", first.stdout)
        self.assertEqual(self.calls.read_text().splitlines(), ["serial", "rhs_flat"])
        second = self.run_stage(23)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertEqual(second.stdout.count("[skip]"), 2)
        self.assertEqual(self.calls.read_text().splitlines(), ["serial", "rhs_flat"])

    def test_failure_continues_but_fails_stage(self):
        result = self.run_stage(23)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("exit=23", result.stdout)
        self.assertNotIn("[done]", result.stdout)
        self.assertEqual(self.calls.read_text().splitlines(), ["serial", "rhs_flat"])
        self.assertFalse((self.out / "post/stack2_e6_actuated_saved_serial_t1/done").exists())
        self.assertTrue((self.out / "post/stack2_e6_actuated_saved_rhs_flat_t2/done").exists())
        resumed = self.run_stage(0)
        self.assertEqual(resumed.returncode, 0, resumed.stdout + resumed.stderr)
        self.assertEqual(resumed.stdout.count("[skip]"), 1)

    def test_retry_isolates_partial_artifacts_and_preserves_completed_point(self):
        self.env.update(MOCK_WRITE_ARTIFACTS="1", MOCK_ATTEMPT="failed")
        first = self.run_stage(23)
        self.assertNotEqual(first.returncode, 0, first.stdout + first.stderr)
        partial = self.out / "post/stack2_e6_actuated_saved_serial_t1"
        completed = self.out / "post/stack2_e6_actuated_saved_rhs_flat_t2"

        def snapshot(path):
            return {str(p.relative_to(path)): (p.read_bytes(), p.stat().st_mtime_ns)
                    for p in path.rglob("*") if p.is_file()}

        failed_files = snapshot(partial)
        completed_files = snapshot(completed)
        self.env["MOCK_ATTEMPT"] = "success"
        second = self.run_stage(0)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertEqual(second.stdout.count("[skip]"), 1)
        self.assertTrue((partial / "done").is_file())
        self.assertEqual([p.name for p in partial.glob("parallelization_performance_raw_*")],
                         ["parallelization_performance_raw_success.csv"])
        self.assertEqual([p.name for p in (partial / "states").iterdir()], ["state_c1.bin"])
        self.assertEqual([p.name for p in (partial / "worker_rows").iterdir()], ["success.csv"])
        archives = list((self.out / ".incomplete").glob("post_stack2_e6_actuated_saved_serial_t1.*/point"))
        self.assertEqual(len(archives), 1)
        self.assertEqual(snapshot(archives[0]), failed_files)
        self.assertEqual(snapshot(completed), completed_files)
        self.assertEqual(self.calls.read_text().splitlines(), ["serial", "rhs_flat", "serial"])

    def test_symlink_point_is_refused_without_touching_target(self):
        target = self.root / "external-point"
        target.mkdir()
        (target / "diagnostic").write_text("keep")
        point = self.out / "post/stack2_e6_actuated_saved_serial_t1"
        point.parent.mkdir(parents=True)
        point.symlink_to(target, target_is_directory=True)
        result = self.run_stage(0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing ambiguous point path", result.stderr)
        self.assertTrue(point.is_symlink())
        self.assertEqual(list(target.iterdir()), [target / "diagnostic"])
        self.assertEqual((target / "diagnostic").read_text(), "keep")
        self.assertEqual(self.calls.read_text().splitlines(), ["rhs_flat"])

    def test_non_directory_point_is_refused(self):
        point = self.out / "post/stack2_e6_actuated_saved_serial_t1"
        point.parent.mkdir(parents=True)
        point.write_text("keep")
        result = self.run_stage(0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing ambiguous point path", result.stderr)
        self.assertEqual(point.read_text(), "keep")
        self.assertEqual(self.calls.read_text().splitlines(), ["rhs_flat"])

    def test_symlink_archive_is_refused_without_moving_partial_output(self):
        point = self.out / "post/stack2_e6_actuated_saved_serial_t1"
        point.mkdir(parents=True)
        (point / "diagnostic").write_text("keep")
        target = self.root / "external-archive"
        target.mkdir()
        (self.out / ".incomplete").symlink_to(target, target_is_directory=True)
        result = self.run_stage(0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing ambiguous archive path", result.stderr)
        self.assertEqual((point / "diagnostic").read_text(), "keep")
        self.assertEqual(list(target.iterdir()), [])
        self.assertEqual(self.calls.read_text().splitlines(), ["rhs_flat"])

    def test_quiet_timeout_continues_but_fails_stage(self):
        result = self.run_stage(75)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("machine not quiet", result.stdout)
        self.assertNotIn("[done]", result.stdout)
        self.assertEqual(self.calls.read_text().splitlines(), ["serial", "rhs_flat"])
        self.assertFalse((self.out / "post/stack2_e6_actuated_saved_serial_t1/done").exists())
        self.assertTrue((self.out / "post/stack2_e6_actuated_saved_rhs_flat_t2/done").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
