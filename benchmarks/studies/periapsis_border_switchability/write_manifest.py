"""Write <OUT>/manifest.toml: config, kernel hashes, git summary, versions, Julia wall time. Usage: python3 write_manifest.py [OUT_DIR]"""
import sys, os, hashlib, subprocess, glob, re, tomllib, json
import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Write <OUT>/manifest.toml")
_ap.add_argument("out_dir", nargs="?", help="directory to read/write (default: the results root)")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir or results_root(_args.results_root)
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(f"{HERE}/../../..")
SPICE = f"{REPO}/data/GRAMSuite.jl/GRAM Suite 2.0/SPICE"
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def run(*a):
    try: return subprocess.run(a, capture_output=True, text=True, cwd=REPO).stdout.strip()
    except Exception as e: return f"unavailable ({e})"
study = tomllib.load(open(f"{HERE}/study_config.toml", "rb"))
files = {"spk/planets/de430.bsp": "SPICE planetary ephemeris (Sun, Moon relative to Earth)", "lsk/naif0012.tls": "leap seconds",
         "pck/pck00011.tpc": "planetary constants", "pck/earth_latest_high_prec.bpc": "Earth orientation (ITRF93)", "tf/earth_assoc_itrf93.tf": "ITRF93 frame kernel"}
wall = 0.0; walls = {}
for f in glob.glob(f"{OUT}/logs/*.time") + glob.glob(f"{OUT}/logs/*.time.txt"):
    m = re.search(r"real\s+(\d+)m([\d.]+)s", open(f).read())
    if m: walls[os.path.basename(f)] = int(m[1]) * 60 + float(m[2]); wall += walls[os.path.basename(f)]
vals = lambda d: json.dumps(d, indent=0)
lines = ["[study]", 'name = "periapsis_border_switchability"', f'julia_version = "{run("julia", "--version")}"',
         f'git_head = "{run("git", "rev-parse", "HEAD")}"', f'git_branch = "{run("git", "rev-parse", "--abbrev-ref", "HEAD")}"',
         f'git_status_summary = {json.dumps(run("git", "status", "--short"))}',
         f'corridor_config_sha256 = "{sha(f"{REPO}/examples/periapsis_corridor_config.toml")}"',
         f'study_config_sha256 = "{sha(f"{HERE}/study_config.toml")}"',
         f'gravity_file_sha256 = "{sha(f"{REPO}/data/planets/gravity/EarthGGM05C.csv")}"',
         f'manifest_toml_sha256 = "{sha(f"{REPO}/Manifest.toml")}"',
         f'julia_wall_seconds_by_phase = {json.dumps(walls)}', f"julia_wall_seconds_total_recorded = {wall:.1f}",
         'ephemerides_note = "SpiceEphemeridesModel for every run (engine refuses N-body under SimpleEphemeridesModel); Earth frame ITRF93"', ""]
for rel, what in files.items():
    p = f"{SPICE}/{rel}"
    lines += [f'[[spice_kernel]]', f'path = "{rel}"', f'purpose = "{what}"', f'sha256 = "{sha(p) if os.path.isfile(p) else "missing"}"', f'bytes = {os.path.getsize(p) if os.path.isfile(p) else 0}', ""]
extra = ""
if os.path.exists(f"{HERE}/static_config.toml") and os.path.exists(f"{OUT}/grid_mc.toml"):  # static-region sweep: config, route/workers per phase, kernel hash of the sweep driver
    extra = "\n[static_region_sweep_config]\n" + open(f"{HERE}/static_config.toml").read()
    for ph in ("grid", "control", "race"):
        if os.path.exists(f"{OUT}/{ph}_mc.toml"): extra += f"\n[static_region_mc_{ph}]\n" + open(f"{OUT}/{ph}_mc.toml").read()
    extra += f'\n[static_region_files]\nrun_static_jl_sha256 = "{sha(f"{HERE}/run_static.jl")}"\ncommon_jl_sha256 = "{sha(f"{HERE}/common.jl")}"\nanalyze_static_py_sha256 = "{sha(f"{HERE}/analyze_static.py")}"\n'
if os.path.exists(f"{OUT}/attribution_mc.toml"):  # attribution study: config, route/workers, file hashes
    extra += "\n[attribution_config]\n" + open(f"{HERE}/attribution_config.toml").read() + "\n[attribution_mc]\n" + open(f"{OUT}/attribution_mc.toml").read()
    extra += "\n[attribution_files]\n" + "".join(f'{n}_sha256 = "{sha(f"{HERE}/{n}")}"\n'.replace(".", "_", 1) for n in ("run_attribution.jl", "estimator_check.jl", "analyze_attribution.py", "analyze_estimator.py", "attribution_config.toml"))
open(f"{OUT}/manifest.toml", "w").write("\n".join(lines) + "\n" + "\n[effective_study_config]\n" + open(f"{HERE}/study_config.toml").read() + extra)
print(open(f"{OUT}/manifest.toml").read()[:1800])
