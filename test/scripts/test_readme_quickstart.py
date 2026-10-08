"""Regression boundaries for the README smoke runner, without simulations."""

import csv
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import zlib


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("readme_smoke", ROOT / "test/smoke/ci_readme_quickstart.py")
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def png_chunk(kind, payload):
    return (struct.pack(">I", len(payload)) + kind + payload
            + struct.pack(">I", zlib.crc32(kind + payload)))


def tiny_png(width=1):
    return (b"\x89PNG\r\n\x1a\n"
            + png_chunk(b"IHDR", struct.pack(">IIBBBBB", width, 1, 8, 2, 0, 0, 0))
            + png_chunk(b"IDAT", zlib.compress(b"\x00\x00\x00\x00"))
            + png_chunk(b"IEND", b""))


def write_csv(output, times=(0.0, 60.0, 120.0), *, columns=smoke.REQUIRED_COLUMNS, bad=None):
    with (output / "simulation_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for time in times:
            row = {name: time if name == "time" else 1.0 for name in columns}
            if bad:
                row.update(bad)
            writer.writerow(row)


def fixture(output):
    (output / "plots").mkdir(parents=True)
    write_csv(output)
    for name in smoke.PLOTS:
        (output / "plots" / name).write_bytes(tiny_png())


class ArtifactTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="readme-artifacts-")
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name) / "output"
        fixture(self.output)

    def test_valid_artifacts_accept_different_sampling(self):
        write_csv(self.output, times=(0.0, 7.0, 51.0, 120.0))
        result = smoke.verify_outputs(self.output)
        self.assertEqual(result["rows"], 4)
        self.assertEqual(set(result["plots"]), set(smoke.PLOTS))
        self.assertEqual(set(result["plots"].values()), {(1, 1)})

    def test_missing_csv(self):
        (self.output / "simulation_results.csv").unlink()
        with self.assertRaises(FileNotFoundError):
            smoke.verify_outputs(self.output)

    def test_empty_or_header_only_csv(self):
        for content in ("", ",".join(smoke.REQUIRED_COLUMNS) + "\n"):
            with self.subTest(content=content):
                (self.output / "simulation_results.csv").write_text(content)
                with self.assertRaises(ValueError):
                    smoke.verify_outputs(self.output)

    def test_each_chart_input_column_is_required(self):
        for name in smoke.REQUIRED_COLUMNS:
            with self.subTest(column=name):
                write_csv(self.output, columns=tuple(c for c in smoke.REQUIRED_COLUMNS if c != name))
                with self.assertRaisesRegex(ValueError, "columns"):
                    smoke.verify_outputs(self.output)

    def test_duplicate_columns(self):
        write_csv(self.output, columns=smoke.REQUIRED_COLUMNS + ("time",))
        with self.assertRaisesRegex(ValueError, "columns"):
            smoke.verify_outputs(self.output)

    def test_invalid_chart_values(self):
        for name in smoke.REQUIRED_COLUMNS:
            for value in ("NaN", "Inf", "-Inf", "", "not-a-number"):
                with self.subTest(column=name, value=value):
                    write_csv(self.output, bad={name: value})
                    with self.assertRaises(ValueError):
                        smoke.verify_outputs(self.output)

    def test_ragged_csv(self):
        path = self.output / "simulation_results.csv"
        for row in ("0,1\n", ",".join(["1"] * (len(smoke.REQUIRED_COLUMNS) + 1)) + "\n"):
            with self.subTest(row=row):
                path.write_text(",".join(smoke.REQUIRED_COLUMNS) + "\n" + row)
                with self.assertRaises(ValueError):
                    smoke.verify_outputs(self.output)

    def test_incomplete_or_unordered_time(self):
        for times in ((0.0,), (0.0, 119.0), (1.0, 120.0),
                      (0.0, 0.0, 120.0), (0.0, 61.0, 60.0, 120.0)):
            with self.subTest(times=times):
                write_csv(self.output, times=times)
                with self.assertRaises(ValueError):
                    smoke.verify_outputs(self.output)

    def test_each_plot_is_required(self):
        for name in smoke.PLOTS:
            with self.subTest(plot=name):
                path = self.output / "plots" / name
                path.unlink()
                with self.assertRaises(FileNotFoundError):
                    smoke.verify_outputs(self.output)
                path.write_bytes(tiny_png())

    def test_invalid_png(self):
        png = tiny_png()
        bad_crc = bytearray(png)
        bad_crc[29] ^= 1
        no_pixels = png[:33] + png_chunk(b"IEND", b"")
        cases = (b"", b"not a PNG", png[:16], png[:-12], png + b"trailing",
                 bytes(bad_crc), tiny_png(width=0), no_pixels)
        for data in cases:
            with self.subTest(data=data):
                (self.output / "plots" / smoke.PLOTS[0]).write_bytes(data)
                with self.assertRaises(ValueError):
                    smoke.verify_outputs(self.output)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="readme runner ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.repo = self.root / "repo with spaces"
        self.repo.mkdir()
        self.source = self.root / "fixture"
        fixture(self.source)
        # Old, valid repository outputs must never satisfy a new run's checks.
        fixture(self.repo / "output")
        self.record = self.root / "launch.json"
        self.child = self.root / "fake julia"

    def child_program(self, body):
        self.child.write_text(
            f"#!{sys.executable}\nimport os, sys, json, shutil\nfrom pathlib import Path\n"
            f"Path({str(self.record)!r}).write_text(json.dumps({{'argv':sys.argv[1:], 'cwd':os.getcwd(), "
            "'override':os.environ.get('SPACEAGORA_CLI_OUTPUT_DIR'), "
            "'smoke':os.environ.get('SPACEAGORA_EXAMPLE_SMOKE'), "
            "'results':os.environ.get('SPACEAGORA_EXAMPLE_SMOKE_RESULTS'), "
            "'duration':os.environ.get('SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME'), "
            "'parallel':os.environ.get('SPACEAGORA_EXAMPLE_PARALLEL')}))\n"
            + body + "\n")
        self.child.chmod(0o755)

    def test_documented_cli_and_fresh_outputs(self):
        self.child_program(f"shutil.copytree({str(self.source)!r}, 'output')")
        inherited = {"SPACEAGORA_CLI_OUTPUT_DIR": str(self.repo / "output"),
                     "SPACEAGORA_EXAMPLE_SMOKE": "1", "SPACEAGORA_EXAMPLE_SMOKE_RESULTS": "1",
                     "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME": "1", "SPACEAGORA_EXAMPLE_PARALLEL": "1"}
        with patch.dict(os.environ, inherited):
            result = smoke.run_smoke(str(self.child), repo_root=self.repo)
        launch = json.loads(self.record.read_text())
        self.assertEqual(result["rows"], 3)
        self.assertEqual(launch["argv"], ["--startup-file=no", "--compiled-modules=existing",
                         f"--project={self.repo}", str(self.repo / "src/cli/main.jl"),
                         "run", "--example=AGORA_Basic_Quickstart.jl", "--smoke"])
        self.assertNotEqual(launch["cwd"], str(self.repo))
        self.assertFalse(Path(launch["cwd"]).exists())
        self.assertIsNone(launch["override"])
        self.assertIsNone(launch["smoke"])
        self.assertIsNone(launch["results"])
        self.assertEqual(launch["duration"], "120.0")
        self.assertEqual(launch["parallel"], "0")
        self.assertTrue((self.repo / "output/simulation_results.csv").is_file())
        documented = "julia --project=. src/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --smoke"
        self.assertIn(documented, (ROOT / "README.md").read_text())

    def test_success_without_outputs_rejects_stale_repository_files(self):
        self.child_program("sys.exit(0)")
        with self.assertRaises(FileNotFoundError):
            smoke.run_smoke(str(self.child), repo_root=self.repo)
        self.assertFalse(Path(json.loads(self.record.read_text())["cwd"]).exists())
        self.assertTrue((self.repo / "output/simulation_results.csv").is_file())

    def test_timeout_stops_simulation_descendant_before_cleanup(self):
        lock_path = self.root / "descendant.lock"
        ready_path = self.root / "descendant.ready"
        child_code = ("import fcntl, time\nfrom pathlib import Path\n"
                      f"lock = open({str(lock_path)!r}, 'w')\n"
                      "fcntl.flock(lock, fcntl.LOCK_EX)\n"
                      f"Path({str(ready_path)!r}).write_text('ready')\n"
                      "time.sleep(60)\n")
        self.child_program("import subprocess, time\n"
                           f"subprocess.Popen([sys.executable, '-c', {child_code!r}])\n"
                           "time.sleep(60)")
        with self.assertRaises(subprocess.TimeoutExpired):
            smoke.run_smoke(str(self.child), repo_root=self.repo, timeout=2)
        self.assertTrue(ready_path.is_file(), "descendant must hold the lock before timeout")
        # An unreaped killed child can briefly be a zombie; its released file
        # lock tests that it has actually stopped, without depending on init.
        with lock_path.open("a") as lock:
            deadline = time.monotonic() + 2
            while True:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    if time.monotonic() >= deadline:
                        self.fail("simulation descendant survived the runner timeout")
                    time.sleep(0.01)
        self.assertFalse(Path(json.loads(self.record.read_text())["cwd"]).exists())

    def test_child_failure_is_not_hidden_by_valid_outputs(self):
        self.child_program(f"shutil.copytree({str(self.source)!r}, 'output')\nsys.exit(7)")
        with self.assertRaises(subprocess.CalledProcessError) as failure:
            smoke.run_smoke(str(self.child), repo_root=self.repo)
        self.assertEqual(failure.exception.returncode, 7)
        self.assertFalse(Path(json.loads(self.record.read_text())["cwd"]).exists())


if __name__ == "__main__":
    unittest.main()
