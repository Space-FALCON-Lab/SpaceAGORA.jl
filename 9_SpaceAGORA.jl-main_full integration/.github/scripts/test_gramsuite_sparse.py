#!/usr/bin/env python3
"""Sparse GRAMSuite retrieval against a local fixture (no network, no GRAM data).

Builds a throwaway "dev" wrapper repository with a GRAMSuite-like layout, points
the production URL at it, and runs .github/scripts/update_gramsuite_submodule.sh
with each sparse profile in .github/gramsuite-sparse/. Checks that the pinned
commit is checked out, that the files the profile selects are there and the
excluded ones are not (and were never downloaded), and that a profile naming a
path the revision lacks fails instead of quietly checking out less.
"""
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / ".github/scripts/update_gramsuite_submodule.sh"
PROFILES = ROOT / ".github/gramsuite-sparse"

# One representative file per pattern the profiles include or exclude.
LAYOUT = {
    "Project.toml": True, "src/GRAMSuite.jl": True,
    "GRAM Suite 2.0/README.md": "native",
    "GRAM Suite 2.0/Julia/GRAM.jl": True,
    "GRAM Suite 2.0/SPICE/pck/pck00011.tpc": "source+",
    "GRAM Suite 2.0/SPICE/spk/planets/de430.bsp": "source+",
    "GRAM Suite 2.0/SPICE/spk/missions/m01_ab_v2.bsp": "native",
    "GRAM Suite 2.0/Build/makefile": "native",
    "GRAM Suite 2.0/common/source/a.cpp": "native",
    "GRAM Suite 2.0/GRAM/examples/x.txt": "native",
    "GRAM Suite 2.0/simulation/GRAM/build_gram.sh": "native",
    "GRAM Suite 2.0/simulation/GRAM/static_grids/grid.bin": False,
    "GRAM Suite 2.0/Earth/source/e.cpp": "native",
    "GRAM Suite 2.0/Earth/data/modeldata/topo.txt": "native",
    "GRAM Suite 2.0/Earth/data/MERRA2data/MERRA2info.txt": "native",
    "GRAM Suite 2.0/Earth/data/MERRA2data/All Mean/MERRA2All_01.bin": "native",
    "GRAM Suite 2.0/Earth/data/MERRA2data/00Z/MERRA2_3hr_00Z_01.bin": False,
    "GRAM Suite 2.0/Earth/earth_surrogate.jls": False,
    "GRAM Suite 2.0/Mars/data/MOLA_data.bin": "native",
    "GRAM Suite 2.0/Venus/source/v.cpp": "native",
    "GRAM Suite 2.0/Jupiter/source/j.cpp": "native",
    "GRAM Suite 2.0/Titan/source/t.cpp": "native",
    "GRAM Suite 2.0/Uranus/source/u.cpp": "native",
    "GRAM Suite 2.0/Uranus/uranus_surrogate.jls": False,
    "GRAM Suite 2.0/Neptune/source/n.cpp": "native",
    "GRAM Suite 2.0/Documentation/manual.pdf": False,
}


def expected(profile, rule):
    if rule is True or rule is False:
        return rule
    if rule == "source+":
        return profile in ("source", "native")
    return profile == rule


class SparseRetrieval(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.env = {k: v for k, v in os.environ.items()
                   if not k.startswith("GIT_") and k not in ("GH_TOKEN", "GITHUB_TOKEN")}
        cls.env.update(GIT_CONFIG_GLOBAL=str(cls.tmp / "gitconfig"), GIT_CONFIG_NOSYSTEM="1",
                       GIT_TERMINAL_PROMPT="0", GIT_ALLOW_PROTOCOL="file", LC_ALL="C")
        remote = cls.tmp / "dev-remote"
        cls.git("init", "--quiet", "--initial-branch=main", str(remote))
        for rel in LAYOUT:
            path = remote / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"fixture {rel}\n")
        cls.git("add", ".", cwd=remote)
        cls.git("-c", "user.name=f", "-c", "user.email=f@example.invalid", "commit", "--quiet", "-m", "pin", cwd=remote)
        cls.pin = cls.git("rev-parse", "HEAD", cwd=remote).strip()
        cls.git("config", "uploadpack.allowFilter", "true", cwd=remote)
        cls.git("config", "uploadpack.allowAnySHA1InWant", "true", cwd=remote)
        cls.git("config", "--global", f"url.{remote.as_uri()}.insteadOf",
                "https://github.com/Space-FALCON-Lab/dev-GRAMSuite.jl.git")
        cls.git("config", "--global", "protocol.file.allow", "always")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def git(cls, *args, cwd=None):
        return subprocess.run(["git", *args], cwd=cwd, env=cls.env, check=True,
                              capture_output=True, text=True).stdout

    def superproject(self, name, profile_dir=None):
        repo = self.tmp / name
        self.git("init", "--quiet", "--initial-branch=main", str(repo))
        (repo / ".github/scripts").mkdir(parents=True)
        shutil.copy(SCRIPT, repo / ".github/scripts")
        shutil.copytree(profile_dir or PROFILES, repo / ".github/gramsuite-sparse")
        (repo / ".github/gramsuite-revisions").write_text(f"regular {'0' * 40}\ndev {self.pin}\n")
        (repo / ".gitmodules").write_text(
            '[submodule "GRAMSuite.jl"]\n\tpath = data/GRAMSuite.jl\n'
            '\turl = https://github.com/Space-FALCON-Lab/GRAMSuite.jl.git\n')
        self.git("add", ".", cwd=repo)
        self.git("update-index", "--add", "--cacheinfo", f"160000,{self.pin},data/GRAMSuite.jl", cwd=repo)
        self.git("-c", "user.name=f", "-c", "user.email=f@example.invalid", "commit", "--quiet", "-m", "s", cwd=repo)
        return repo

    def run_script(self, repo, profile):
        env = dict(self.env, GRAMSUITE_SPARSE_PROFILE=profile)
        return subprocess.run(["bash", ".github/scripts/update_gramsuite_submodule.sh", "dev"],
                              cwd=repo, env=env, capture_output=True, text=True)

    def test_profiles(self):
        for profile in ("wrapper", "source", "native"):
            with self.subTest(profile=profile):
                repo = self.superproject(f"super-{profile}")
                r = self.run_script(repo, profile)
                self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                sub = repo / "data/GRAMSuite.jl"
                self.assertEqual(self.git("rev-parse", "HEAD", cwd=sub).strip(), self.pin)
                for rel, rule in LAYOUT.items():
                    self.assertEqual((sub / rel).exists(), expected(profile, rule), f"{profile}: {rel}")
                # Excluded blobs were never fetched, not just left unchecked-out.
                missing = self.git("rev-list", "--objects", "--missing=print", "HEAD", cwd=sub)
                self.assertIn("?", missing, profile)

    def test_missing_required_path_fails(self):
        profiles = self.tmp / "bad-profiles"
        profiles.mkdir()
        (profiles / "bad.txt").write_text("/src/\n/GRAM Suite 2.0/NoSuchDir/\n")
        repo = self.superproject("super-bad", profiles)
        r = self.run_script(repo, "bad")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("required path missing", r.stderr)

    def test_unknown_profile_fails(self):
        repo = self.superproject("super-unknown")
        self.assertNotEqual(self.run_script(repo, "nosuch").returncode, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
