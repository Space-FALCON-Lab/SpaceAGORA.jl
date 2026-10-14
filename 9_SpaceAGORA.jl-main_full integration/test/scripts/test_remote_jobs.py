"""Offline transport-failure regressions for the remote-job shell functions."""

import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/remote/spaceagora-remote"


class RemoteJobTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="remote-job-tests-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.functions = self.root / "functions.sh"
        # Exercise the real functions without main/config loading. Configuration
        # discovery is mocked, so Bash 3 on macOS needs no associative array.
        source = SCRIPT.read_text()
        self.assertTrue(source.endswith('main "$@"\n'))
        self.functions.write_text(source.removesuffix('main "$@"\n').replace(
            "declare -A REMOTES", ": # configuration is mocked"))
        self.prelude = f"""
source {shlex.quote(str(self.functions))}
REMOTE_NAME=mock
REMOTE_BASE={shlex.quote(str(self.root / 'remote'))}
SSH_ALIAS=mock
RUNS_DIR={shlex.quote(str(self.root / 'runs'))}
resolve_remote() {{ return 0; }}
ensure_remote_dirs() {{ return 0; }}
new_job_id() {{ printf job; }}
spawn_watcher() {{ echo UNEXPECTED_WATCHER; }}
"""

    def run_shell(self, body, *, prelude=None):
        return subprocess.run(
            ["bash", "-c", (self.prelude if prelude is None else prelude) + body],
            text=True, capture_output=True, timeout=10,
            cwd=self.root,
            env={**os.environ, "TMPDIR": str(self.root)},
        )

    def test_shell_syntax(self):
        result = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_run_does_not_launch_after_either_upload_fails(self):
        for target in ("run.sh", "job.meta"):
            for failure_code in (1, 255):
                with self.subTest(target=target, code=failure_code):
                    result = self.run_shell(f"""
rsh() {{
  case "$1" in
    *"cat >"*"{target}"*) cat >/dev/null; return {failure_code} ;;
    *tmux*) echo TMUX_LAUNCHED ;;
    *) cat >/dev/null ;;
  esac
}}
cmd_run --remote mock --no-watch -- true
""")
                    self.assertEqual(result.returncode, failure_code, result.stderr)
                    self.assertNotIn("TMUX_LAUNCHED", result.stdout)
                    self.assertNotIn("job id:", result.stdout)
                    self.assertNotIn("UNEXPECTED_WATCHER", result.stdout)

    def test_run_preserves_literal_cwd_arguments_and_exit_code(self):
        for cwd in (None, "directory with spaces", "apostrophe's $literal; directory", "-leading"):
            with self.subTest(cwd=cwd):
                if cwd is not None:
                    (self.root / cwd).mkdir()
                args = ["--remote", "mock", "--no-watch"]
                if cwd is not None:
                    args += ["--cwd", cwd]
                command_args = ["two words", "apostrophe's", "$(printf unexpected)", "a; exit 9"]
                args += ["--", "python3", "-c",
                         "import json,os,sys; print(json.dumps([os.getcwd(),sys.argv[1:]])); sys.exit(3)",
                         *command_args]
                result = self.run_shell("""
rsh() {
  case "$1" in
    tmux*) bash "$REMOTE_BASE/jobs/job/run.sh" ;;
    *) eval "$1" ;;
  esac
}
cmd_run """ + shlex.join(args))
                self.assertEqual(result.returncode, 0, result.stderr)
                job = self.root / "remote/jobs/job"
                expected_cwd = Path.home().resolve() if cwd is None else (self.root / cwd).resolve()
                self.assertEqual(json.loads((job / "job.log").read_text()), [str(expected_cwd), command_args])
                self.assertIn("EXIT_CODE=3\n", (job / "job.meta").read_text())

    def test_pull_preserves_previous_file_on_interrupted_transfer(self):
        local = self.root / "runs/mock/job"
        local.mkdir(parents=True)
        for target in ("job.log", "job.meta"):
            with self.subTest(target=target):
                (local / target).write_text("previous successful pull\n")
                result = self.run_shell(f"""
rsh() {{
  case "$1" in
    "[[ -d "*) return 1 ;;
    "[[ -f "*) return 0 ;;
    *"cat "*"{target}"*) printf 'partial transfer'; return 255 ;;
    *"cat "*) printf 'complete transfer\\n' ;;
    *) return 99 ;;
  esac
}}
do_pull job
""")
                self.assertEqual(result.returncode, 255, result.stderr)
                self.assertEqual((local / target).read_text(), "previous successful pull\n")
                self.assertNotIn("pulled to", result.stdout)
                self.assertEqual(list(local.glob(f".{target}.*")), [])

    def test_pull_propagates_directory_probe_and_rsync_failures(self):
        for probe_code, transfer_code in ((255, 0), (0, 12)):
            with self.subTest(probe=probe_code, transfer=transfer_code):
                result = self.run_shell(f"""
rsh() {{ return {probe_code}; }}
rsync() {{ return {transfer_code}; }}
do_pull job
""")
                self.assertEqual(result.returncode, probe_code or transfer_code, result.stderr)
                self.assertNotIn("pulled to", result.stdout)

    def test_died_job_without_log_still_pulls_metadata(self):
        result = self.run_shell("""
rsh() {
  case "$1" in
    "[[ -d "*|"[[ -f "*) return 1 ;;
    *"cat "*"job.meta"*) printf 'JOB_ID=job\\n' ;;
    *) return 99 ;;
  esac
}
do_pull job
""")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.root / "runs/mock/job/job.meta").read_text(), "JOB_ID=job\n")

    def test_status_distinguishes_lifecycle_states_and_outage(self):
        jobs = self.root / "remote/jobs"
        for name in ("running", "finished", "died"):
            (jobs / name).mkdir(parents=True)
            (jobs / name / "job.meta").write_text("JOB_ID=" + name + "\n")
        (jobs / "finished/job.meta").write_text("EXIT_CODE=3\n")
        for name, expected in (("running", "running"), ("finished", "finished (EXIT_CODE=3)"),
                               ("died", "died (no EXIT_CODE and no tmux session)"),
                               ("missing", "unknown job")):
            with self.subTest(state=name):
                result = self.run_shell(f"""
tmux() {{ [[ "$3" == spaceagora-running ]]; }}
rsh() {{ eval "$1"; }}
cmd_status --remote mock {name}
""")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines()[0], "state: " + expected)
        result = self.run_shell("rsh() { return 255; }; cmd_status --remote mock running")
        self.assertEqual(result.returncode, 255)
        self.assertIn("state: unreachable", result.stdout)

    def test_watcher_retries_outage_and_failed_pull_before_completing(self):
        cli_dir = self.root / "cli"
        cli_dir.mkdir()
        state_file = self.root / "calls.json"
        state_file.write_text(json.dumps({"status": 0, "pull": 0}))
        cli = cli_dir / "spaceagora-remote"
        cli.write_text(f"""#!/usr/bin/env python3
import json, sys
from pathlib import Path
p = Path({str(state_file)!r})
state = json.loads(p.read_text())
cmd = sys.argv[1]
state[cmd] += 1
p.write_text(json.dumps(state))
if cmd == 'status':
    n = state[cmd]
    print('state: unreachable' if n == 1 else 'state: running' if n == 2 else 'state: finished (EXIT_CODE=0)')
    sys.exit(255 if n == 1 else 0)
sys.exit(255 if state['pull'] == 1 else 0)
""")
        cli.chmod(0o755)
        # Use the production watcher generator, but suppress its background
        # launch and run the generated file synchronously with a hard timeout.
        prelude = self.prelude.replace('spawn_watcher() { echo UNEXPECTED_WATCHER; }', '')
        result = self.run_shell(f"""
SCRIPT_DIR={shlex.quote(str(cli_dir))}
setsid() {{ return 0; }}
spawn_watcher job 0.01
bash "$RUNS_DIR/mock/job/.watch.sh"
""", prelude=prelude)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(state_file.read_text()), {"status": 4, "pull": 2})
        self.assertEqual(result.stdout.count("watch complete:"), 1)
        self.assertIn("pull failed; will retry", result.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
