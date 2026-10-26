#!/usr/bin/env python3
"""Run the extracted, pinned Mars RAAN smoke using this job's existing native assets.

No licensed data or native source is fetched, copied, cached, or uploaded here.
The public package wrapper and its low-level Julia wrapper come from the public
GRAM pin; only Build, SPICE and planet data are linked from existing CI inputs.
"""
from __future__ import annotations

import argparse
import base64
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import tomllib
import unittest
from unittest.mock import patch

STUDY_URL = "https://github.com/Space-FALCON-Lab/SpaceAGORAMissionReconstructions.git"
STUDY_REV = "ff67a512face495afa6dec10171806209a62033d"
CORE_URL = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git"
CORE_REV = "bfa32f6a44c5f242714b8cf182a451e863119a05"
GRAM_URL = "https://github.com/Space-FALCON-Lab/GRAMSuite.jl.git"
GRAM_REV = "b6675517b2f7080044d7bd40ba5502dc463710f2"
READ_TOKEN_ENV = "SPACEAGORA_RESEARCH_READ_TOKEN"
NAN_MARKERS = ("init_NaN", "First function call produced NaNs")


def require_pin(revision: str) -> None:
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("An immutable full commit SHA is required")


def native_assets(workspace: Path) -> Path:
    root = workspace / "data/GRAMSuite.jl/GRAM Suite 2.0"
    for relative in ("Build/lib/libGRAM.so", "SPICE", "Mars/data", "Earth/data"):
        path = root / relative
        valid = path.is_file() if relative.endswith(".so") else path.is_dir()
        if not valid:
            raise ValueError(f"Existing native CI input is missing: {path}")
    return root.resolve()


def study_fetch_environment(url: str, environment: dict[str, str], token: str) -> dict[str, str]:
    if url != STUDY_URL:
        raise ValueError("Research credentials are restricted to the pinned Mission repository")
    result = {key: value for key, value in environment.items()
              if key != READ_TOKEN_ENV and not key.startswith("GIT_TRACE")
              and key != "GIT_CURL_VERBOSE"}
    count = int(result.get("GIT_CONFIG_COUNT", "0"))
    if count < 0:
        raise ValueError("Invalid inherited Git configuration count")
    result["GIT_CONFIG_COUNT"] = str(count + 1)
    result[f"GIT_CONFIG_KEY_{count}"] = f"http.{STUDY_URL}.extraheader"
    encoded = base64.b64encode(f"x-access-token:{token}".encode()).decode()
    result[f"GIT_CONFIG_VALUE_{count}"] = f"AUTHORIZATION: basic {encoded}"
    return result


def checkout(url: str, revision: str, destination: Path, environment: dict[str, str],
             read_token: str | None = None) -> None:
    require_pin(revision)
    fetch_environment = (study_fetch_environment(url, environment, read_token)
                         if read_token else environment)
    if destination.exists():
        raise ValueError(f"Checkout destination must be new: {destination}")
    git = ["git", "-c", "credential.helper=", "-c", "filter.lfs.smudge=cat",
           "-c", "filter.lfs.required=false"]
    subprocess.run(git + ["init", str(destination)], check=True, env=environment)
    subprocess.run(git + ["-C", str(destination), "remote", "add", "origin", url], check=True, env=environment)
    subprocess.run(git + ["-C", str(destination), "fetch", "--depth=1", "--no-tags", "origin", revision], check=True, env=fetch_environment)
    subprocess.run(git + ["-C", str(destination), "checkout", "--detach", revision], check=True, env=environment)
    actual = subprocess.run(git + ["-C", str(destination), "rev-parse", "HEAD"],
                            check=True, env=environment, text=True, capture_output=True).stdout.strip()
    if actual != revision:
        raise ValueError(f"Checkout identity differs: expected {revision}, found {actual}")


def verify_study_pins(study: Path) -> None:
    for name, expected_url, expected_revision in (
        ("SpaceAGORASource.toml", CORE_URL, CORE_REV),
        ("GRAMSuiteSource.toml", GRAM_URL, GRAM_REV),
    ):
        with (study / name).open("rb") as stream:
            pin = tomllib.load(stream)
        if (pin.get("url"), pin.get("rev")) != (expected_url, expected_revision):
            raise ValueError(f"Published study pin differs from the core CI contract: {name}")


def make_overlay(public_gram: Path, assets: Path, destination: Path) -> Path:
    julia_wrapper = public_gram / "GRAM Suite 2.0/Julia"
    if not (julia_wrapper / "GRAM.jl").is_file():
        raise ValueError("Pinned public GRAM Julia wrapper is missing")
    destination.mkdir()  # refuse an existing overlay rather than replace anything
    (destination / "Julia").symlink_to(julia_wrapper, target_is_directory=True)
    for name in ("Build", "SPICE"):
        (destination / name).symlink_to(assets / name, target_is_directory=True)
    for name in ("Earth", "Mars"):
        (destination / name).mkdir()
        (destination / name / "data").symlink_to(assets / name / "data", target_is_directory=True)
    return destination


def link_core_asset_root(core: Path, overlay: Path) -> None:
    # Core's pinned SPICE constants still resolve this default path. Only the
    # newly checked-out core's uninitialized gitlink directory may be populated.
    parent = core / "data/GRAMSuite.jl"
    if parent.is_symlink() or (parent.exists() and
                               (not parent.is_dir() or any(parent.iterdir()))):
        raise ValueError(f"Core GRAM asset directory must be new or empty: {parent}")
    if not overlay.is_dir():
        raise ValueError(f"Native overlay is missing: {overlay}")
    parent.mkdir(parents=True, exist_ok=True)
    (parent / "GRAM Suite 2.0").symlink_to(overlay, target_is_directory=True)


def smoke_environment(base: dict[str, str], core: Path, gram: Path, overlay: Path,
                      environment: Path) -> dict[str, str]:
    result = dict(base)
    result.update({
        "GIT_LFS_SKIP_SMUDGE": "1", "GIT_TERMINAL_PROMPT": "0",
        "JULIA_PKG_PRECOMPILE_AUTO": "0", "JULIA_LOAD_PATH": "@:@stdlib",
        "SPACEAGORA_CORE_ROOT": str(core), "SPACEAGORA_GRAM_SOURCE": str(gram),
        "SPACEAGORA_DATA_ROOT": str(core / "data"),
        "SPACEAGORA_GRAM_DATA_ROOT": str(overlay),
        "SPACEAGORA_MISSION_ENV": str(environment),
        "GRAM_ROOT": str(overlay), "GRAM_LIB": str(overlay / "Build/lib/libGRAM.so"),
        # Match the former automatic example smoke, including output production.
        "SPACEAGORA_EXAMPLE_SMOKE": "1", "SPACEAGORA_EXAMPLE_SMOKE_RESULTS": "1",
        "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME": "120.0",
        "SPACEAGORA_WARN_DEPRECATED_CONFIG": "0", "SPACEAGORA_WARN_NORMALIZE": "0",
    })
    return result


def run_smoke(julia: str, study: Path, environment: Path, variables: dict[str, str], log: Path) -> None:
    # check_mars_raan.jl deliberately disables results; use the same driver directly
    # to retain the old automatic suite's results=true branch and NaN rejection.
    command = [julia, "--startup-file=no", "--compiled-modules=existing", "--depwarn=error",
               f"--project={environment}", str(study / "scripts/run.jl"),
               "examples/AGORA_Mars_RAAN_Scenario.jl"]
    result = subprocess.run(command, cwd=study, env=variables, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.write_text(result.stdout)
    print(result.stdout, end="", flush=True)
    if result.returncode != 0 or any(marker in result.stdout for marker in NAN_MARKERS):
        raise RuntimeError("Pinned Mars RAAN native smoke failed; see its process output")


def run(workspace: Path, scratch: Path, julia: str) -> None:
    base = dict(os.environ, GIT_LFS_SKIP_SMUDGE="1", GIT_TERMINAL_PROMPT="0",
                JULIA_PKG_PRECOMPILE_AUTO="0")
    read_token = base.pop(READ_TOKEN_ENV, "")
    if not read_token:
        message = ("SKIP: pinned Mars RAAN native coverage not exercised; "
                   "private study read credential is unavailable.")
        print(message, flush=True)
        if summary := base.get("GITHUB_STEP_SUMMARY"):
            with Path(summary).open("a") as stream:
                stream.write(message + "\n")
        return
    workspace = workspace.resolve()
    # With credentials supplied, missing inputs and all later failures are errors.
    assets = native_assets(workspace)
    if scratch.exists():
        raise ValueError(f"Choose a new scratch directory: {scratch}")
    scratch.mkdir(parents=True)
    scratch = scratch.resolve()
    study, core, gram = (scratch / name for name in ("study", "core", "public-gram"))
    for url, revision, destination in ((STUDY_URL, STUDY_REV, study), (CORE_URL, CORE_REV, core),
                                       (GRAM_URL, GRAM_REV, gram)):
        checkout(url, revision, destination, base, read_token if url == STUDY_URL else None)
    verify_study_pins(study)
    gravity = core / "data/planets/gravity/Mars50c.csv"
    if not gravity.is_file() or gravity.read_bytes().startswith(b"version https://git-lfs.github.com/spec/v1"):
        raise ValueError("Pinned public Mars gravity input is missing or an unhydrated LFS pointer")
    overlay = make_overlay(gram, assets, scratch / "native-overlay")
    link_core_asset_root(core, overlay)
    environment = scratch / "environment"
    variables = smoke_environment(base, core, gram, overlay, environment)
    subprocess.run([julia, "--startup-file=no", "--compiled-modules=existing", "--project=@stdlib",
                    str(study / "scripts/setup.jl"), "--with-gram"],
                   cwd=study, env=variables, check=True)
    run_smoke(julia, study, environment, variables, scratch / "mars-raan-smoke.log")
    print(f"Pinned Mars RAAN smoke passed: study={STUDY_REV} core={CORE_REV} GRAMSuite={GRAM_REV}")


class RouteControls(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.assets = self.root / "data/GRAMSuite.jl/GRAM Suite 2.0"
        for path in ("Build/lib", "SPICE", "Mars/data", "Earth/data"):
            (self.assets / path).mkdir(parents=True, exist_ok=True)
        (self.assets / "Build/lib/libGRAM.so").write_text("mock native input; never loaded")
        self.gram = self.root / "public-gram"
        wrapper = self.gram / "GRAM Suite 2.0/Julia"
        wrapper.mkdir(parents=True)
        (wrapper / "GRAM.jl").write_text("# exact public wrapper fixture")

    def test_missing_assets_refuse_before_fetch(self):
        (self.assets / "Build/lib/libGRAM.so").unlink()
        with patch.dict(os.environ, {READ_TOKEN_ENV: "dummy"}), \
                patch.object(sys.modules[__name__], "checkout") as fetch:
            with self.assertRaisesRegex(ValueError, "Existing native CI input"):
                run(self.root, self.root / "new-scratch", "julia")
            fetch.assert_not_called()
        self.assertFalse((self.root / "new-scratch").exists())

    def test_absent_or_empty_credentials_skip_before_work(self):
        summary = self.root / "summary"
        for token in (None, ""):
            variables = {"GITHUB_STEP_SUMMARY": str(summary)}
            if token is not None:
                variables[READ_TOKEN_ENV] = token
            with patch.dict(os.environ, variables, clear=True), \
                    patch.object(sys.modules[__name__], "native_assets") as assets, \
                    patch.object(sys.modules[__name__], "checkout") as fetch:
                run(self.root / "missing", self.root / "new-scratch", "julia")
                assets.assert_not_called()
                fetch.assert_not_called()
            self.assertFalse((self.root / "new-scratch").exists())
        self.assertEqual(summary.read_text().count("native coverage not exercised"), 2)

    def test_auth_is_scoped_to_one_fetch_and_never_arguments(self):
        token = "dummy-read-token"
        encoded = base64.b64encode(f"x-access-token:{token}".encode()).decode()
        base = {"GIT_CONFIG_COUNT": "1", "GIT_CONFIG_KEY_0": "core.autocrlf",
                "GIT_CONFIG_VALUE_0": "false", "GIT_TRACE": "1",
                "GIT_TRACE_CURL": "1", "GIT_CURL_VERBOSE": "1"}
        original = dict(base)
        for url, revision in ((STUDY_URL, STUDY_REV), (CORE_URL, CORE_REV), (GRAM_URL, GRAM_REV)):
            fake = subprocess.CompletedProcess([], 0, stdout=revision + "\n")
            with patch("subprocess.run", return_value=fake) as call:
                checkout(url, revision, self.root / "checkout", base,
                         token if url == STUDY_URL else None)
            authenticated = 0
            for invocation in call.call_args_list:
                command, = invocation.args
                variables = invocation.kwargs["env"]
                self.assertNotIn(token, str(command))
                self.assertNotIn(encoded, str(command))
                self.assertNotIn(READ_TOKEN_ENV, variables)
                if url == STUDY_URL and "fetch" in command:
                    authenticated += 1
                    self.assertEqual(variables["GIT_CONFIG_COUNT"], "2")
                    self.assertEqual(variables["GIT_CONFIG_KEY_1"], f"http.{STUDY_URL}.extraheader")
                    self.assertEqual(variables["GIT_CONFIG_VALUE_1"], f"AUTHORIZATION: basic {encoded}")
                    self.assertEqual(variables["GIT_CONFIG_VALUE_0"], "false")
                    self.assertFalse(any(key.startswith("GIT_TRACE") or key == "GIT_CURL_VERBOSE"
                                         for key in variables))
                else:
                    self.assertNotIn(encoded, str(variables))
                    self.assertEqual(variables, original)
            self.assertEqual(authenticated, int(url == STUDY_URL))
        self.assertEqual(base, original)
        with self.assertRaisesRegex(ValueError, "restricted"):
            study_fetch_environment(CORE_URL, {}, token)

    def test_git_url_matching_does_not_authenticate_other_repositories(self):
        variables = study_fetch_environment(STUDY_URL, dict(os.environ), "dummy-read-token")
        for url in (STUDY_URL, CORE_URL, GRAM_URL,
                    "https://github.com/Space-FALCON-Lab/SpaceAGORAMissionReconstructionsSibling.git"):
            result = subprocess.run(["git", "config", "--get-urlmatch", "http.extraheader", url],
                                    env=variables, text=True, capture_output=True)
            if url == STUDY_URL:
                self.assertEqual(result.returncode, 0)
                self.assertIn("AUTHORIZATION: basic ", result.stdout)
            else:
                self.assertNotIn("AUTHORIZATION: basic ", result.stdout)

    def test_run_removes_token_from_setup_and_runtime_environments(self):
        token = "dummy-read-token"
        seen = []
        def fake_checkout(url, revision, destination, environment, read_token=None):
            seen.append((url, read_token, dict(environment)))
            destination.mkdir()
            if url == CORE_URL:
                gravity = destination / "data/planets/gravity/Mars50c.csv"
                gravity.parent.mkdir(parents=True)
                gravity.write_text("gravity fixture")
        with patch.dict(os.environ, {READ_TOKEN_ENV: token}), \
                patch.object(sys.modules[__name__], "checkout", side_effect=fake_checkout), \
                patch.object(sys.modules[__name__], "verify_study_pins"), \
                patch.object(sys.modules[__name__], "make_overlay", return_value=self.root / "overlay"), \
                patch.object(sys.modules[__name__], "link_core_asset_root"), \
                patch("subprocess.run") as setup, \
                patch.object(sys.modules[__name__], "run_smoke") as smoke:
            run(self.root, self.root / "scratch", "julia")
        self.assertEqual([(url, credential) for url, credential, _ in seen],
                         [(STUDY_URL, token), (CORE_URL, None), (GRAM_URL, None)])
        environments = [env for _, _, env in seen] + [setup.call_args.kwargs["env"], smoke.call_args.args[3]]
        for environment in environments:
            self.assertNotIn(READ_TOKEN_ENV, environment)
            self.assertNotIn(token, str(environment))

    def test_supplied_credential_failure_is_not_a_skip(self):
        summary = self.root / "summary"
        failure = subprocess.CalledProcessError(128, ["git", "fetch"])
        with patch.dict(os.environ, {READ_TOKEN_ENV: "dummy", "GITHUB_STEP_SUMMARY": str(summary)}), \
                patch.object(sys.modules[__name__], "checkout", side_effect=failure):
            with self.assertRaises(subprocess.CalledProcessError):
                run(self.root, self.root / "failed-scratch", "julia")
        self.assertFalse(summary.exists())

    def test_each_required_asset_is_mandatory(self):
        for relative in ("SPICE", "Mars/data", "Earth/data"):
            path = self.assets / relative
            path.rmdir()
            with self.assertRaisesRegex(ValueError, "Existing native CI input"):
                native_assets(self.root)
            path.mkdir()

    def test_overlay_uses_public_julia_without_copying_inputs(self):
        overlay = make_overlay(self.gram, native_assets(self.root), self.root / "overlay")
        self.assertEqual((overlay / "Julia").resolve(), (self.gram / "GRAM Suite 2.0/Julia").resolve())
        for relative in ("Build", "SPICE", "Mars/data", "Earth/data"):
            self.assertTrue((overlay / relative).is_symlink())
            self.assertEqual((overlay / relative).resolve(), (self.assets / relative).resolve())
        with self.assertRaises(FileExistsError):
            make_overlay(self.gram, self.assets, overlay)

    def test_core_default_spice_resolves_existing_assets(self):
        overlay = make_overlay(self.gram, native_assets(self.root), self.root / "overlay")
        core = self.root / "core"
        # A gitlink is commonly represented by an empty directory after checkout.
        (core / "data/GRAMSuite.jl").mkdir(parents=True)
        link_core_asset_root(core, overlay)
        default_root = core / "data/GRAMSuite.jl/GRAM Suite 2.0"
        self.assertTrue(default_root.is_symlink())
        self.assertEqual((default_root / "SPICE").resolve(), (self.assets / "SPICE").resolve())
        self.assertEqual((default_root / "Julia").resolve(), (self.gram / "GRAM Suite 2.0/Julia").resolve())
        self.assertFalse((self.assets / "Julia").exists())
        with self.assertRaisesRegex(ValueError, "must be new or empty"):
            link_core_asset_root(core, overlay)
        fresh_core = self.root / "fresh-core"
        link_core_asset_root(fresh_core, overlay)
        self.assertEqual((fresh_core / "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE").resolve(),
                         (self.assets / "SPICE").resolve())

    def test_core_asset_parent_refuses_existing_inputs(self):
        overlay = make_overlay(self.gram, native_assets(self.root), self.root / "overlay")
        core = self.root / "core"
        parent = core / "data/GRAMSuite.jl"
        parent.mkdir(parents=True)
        sentinel = parent / "existing"
        sentinel.write_text("preserve")
        with self.assertRaisesRegex(ValueError, "must be new or empty"):
            link_core_asset_root(core, overlay)
        self.assertEqual(sentinel.read_text(), "preserve")
        sentinel.unlink()
        parent.rmdir()
        parent.symlink_to(self.assets, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "must be new or empty"):
            link_core_asset_root(core, overlay)
        self.assertFalse((self.assets / "GRAM Suite 2.0").exists())

    def test_missing_public_wrapper_refuses(self):
        (self.gram / "GRAM Suite 2.0/Julia/GRAM.jl").unlink()
        with self.assertRaisesRegex(ValueError, "public GRAM Julia wrapper"):
            make_overlay(self.gram, self.assets, self.root / "overlay")
        self.assertFalse((self.root / "overlay").exists())

    def test_wrong_pin_is_not_a_fallback(self):
        study = self.root / "study"
        study.mkdir()
        for name, url, rev in (("SpaceAGORASource.toml", CORE_URL, CORE_REV),
                               ("GRAMSuiteSource.toml", GRAM_URL, GRAM_REV)):
            (study / name).write_text(f'url = "{url}"\nrev = "{rev}"\n')
        verify_study_pins(study)
        for name in ("SpaceAGORASource.toml", "GRAMSuiteSource.toml"):
            original = (study / name).read_text()
            (study / name).write_text(original.replace('rev = "', 'rev = "0'))
            with self.assertRaisesRegex(ValueError, "Published study pin differs"):
                verify_study_pins(study)
            (study / name).write_text(original)
        for value in ("main", "v0.5.1", CORE_REV[:7], CORE_REV + "x"):
            with self.assertRaises(ValueError):
                require_pin(value)

    def test_checkout_head_mismatch_refuses(self):
        fake = subprocess.CompletedProcess([], 0, stdout="0" * 40 + "\n")
        with patch("subprocess.run", return_value=fake):
            with self.assertRaisesRegex(ValueError, "Checkout identity differs"):
                checkout(CORE_URL, CORE_REV, self.root / "new-checkout", {})

    def test_smoke_keeps_results_and_refuses_exit_or_nan(self):
        env = smoke_environment({"JULIA_NUM_THREADS": "4"}, self.root / "core", self.gram,
                                self.root / "overlay", self.root / "env")
        self.assertEqual(env["SPACEAGORA_EXAMPLE_SMOKE_RESULTS"], "1")
        self.assertEqual(env["SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME"], "120.0")
        self.assertEqual(env["JULIA_LOAD_PATH"], "@:@stdlib")
        self.assertEqual(env["JULIA_NUM_THREADS"], "4")
        self.assertEqual(env["GRAM_ROOT"], env["SPACEAGORA_GRAM_DATA_ROOT"])
        for returncode, output, failed in ((0, "native mock passed\n", False), (1, "stopped\n", True),
                                           (0, "init_NaN\n", True), (0, "First function call produced NaNs\n", True)):
            with patch("subprocess.run", return_value=subprocess.CompletedProcess([], returncode, stdout=output)) as call:
                if failed:
                    with self.assertRaises(RuntimeError):
                        run_smoke("julia", self.root, self.root / "env", env, self.root / "log")
                else:
                    run_smoke("julia", self.root, self.root / "env", env, self.root / "log")
                self.assertIn("--compiled-modules=existing", call.call_args.args[0])
                self.assertEqual(call.call_args.kwargs["env"], env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--julia", default="julia")
    options = parser.parse_args()
    if options.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RouteControls))
        return 0 if result.wasSuccessful() else 1
    if options.workspace is None or options.scratch is None:
        parser.error("--workspace and --scratch are required")
    run(options.workspace, options.scratch, options.julia)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
