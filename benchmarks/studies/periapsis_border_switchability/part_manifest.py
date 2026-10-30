"""Manifest for the blind and superposition parts: config, kernel hashes, git summary, versions, route/workers per run.
Usage: python3 -I part_manifest.py [--results-root DIR]   (writes <root>/blind/manifest.toml and <root>/superposition/manifest.toml)"""
import argparse, hashlib, pathlib, subprocess, sys
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
REPO = HERE.parents[2]; SPICE = REPO / "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE"
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() if pathlib.Path(p).is_file() else "missing"
def run(*a):
    try: return subprocess.run(a, capture_output=True, text=True, cwd=REPO).stdout.strip()
    except Exception as e: return f"unavailable ({e})"
p = argparse.ArgumentParser(); add_results_root_argument(p); root = pathlib.Path(results_root(p.parse_args().results_root))
parts = {"blind": (["blind_config.toml", "run_blind.jl", "select_blind.py", "analyze_blind.py", "run_static.jl", "common.jl"], ["blind"]),
         "superposition": (["superposition_config.toml", "run_superposition.jl", "analyze_superposition.py", "run_static.jl", "common.jl"], ["null", "replay_it0", "replay_it1", "replay_it2", "finite"])}
for name, (files, mcs) in parts.items():
    out = root / name; L = ["[run]", f'julia_version = "{run("julia", "--version")}"', f'git_head = "{run("git", "rev-parse", "HEAD")}"',
         f'git_branch = "{run("git", "rev-parse", "--abbrev-ref", "HEAD")}"', f'manifest_toml_sha256 = "{sha(REPO / "Manifest.toml")}"',
         f'corridor_config_sha256 = "{sha(REPO / "examples/periapsis_corridor_config.toml")}"', f'study_config_sha256 = "{sha(HERE / "study_config.toml")}"',
         f'gravity_file_sha256 = "{sha(REPO / "data/planets/gravity/EarthGGM05C.csv")}"', ""]
    for rel in ("spk/planets/de430.bsp", "lsk/naif0012.tls", "pck/pck00011.tpc", "pck/earth_latest_high_prec.bpc", "tf/earth_assoc_itrf93.tf"):
        L += ["[[spice_kernel]]", f'path = "{rel}"', f'sha256 = "{sha(SPICE / rel)}"', ""]
    L += ["[file_sha256]"] + [f'"{f}" = "{sha(HERE / f)}"' for f in files] + [""]
    for m in mcs:
        f = out / f"{m}_mc.toml"
        if f.exists(): L += [f"[mc_{m}]", f.read_text()]
    L += ["[config]", (HERE / files[0]).read_text()]
    (out / "manifest.toml").write_text("\n".join(L)); print("wrote", out / "manifest.toml")
