"""WP-A0: copy the volatile /tmp tables that the analytic validation reads into tasks/.../analytic_data/snapshot/, with a sha256 manifest.

Usage: ANALYTIC_TMP_ROOT=<dir with pc20, pcstagek, pcstagef> python3 -I analytic_snapshot.py [--dest DIR]   (idempotent; refuses to overwrite a file whose hash differs)
Sources are the Stage K pair tables, the apoapsis-Pi dev pairs and the F20 open-loop pass tables (Mars, Earth) and the Mars mean-field truth ablation tables (3b, zonal, deg3). The Stage O dev tables are already local.
"""
import argparse, hashlib, pathlib, shutil

import os, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from analytic_pairs import TASKS
TMP = pathlib.Path(os.environ.get("ANALYTIC_TMP_ROOT", "/tmp"))   # directory holding the volatile results roots pc20, pcstagek, pcstagef of earlier study runs (set ANALYTIC_TMP_ROOT)
DEFAULT = TASKS / "analytic_data" / "snapshot"


def sources():
    out = []
    for f in ("pairs_full_mars.csv", "pairs_full_mars_tess.csv", "pairs_full_earth.csv", "pairs_full_venus.csv", "k_results.json", "k_tess_results.json"):
        out.append((TMP / "pcstagek" / f, f"stagek/{f}"))
    for f in ("dev_pairs.csv", "odyssey_phases.csv", "criteria.csv", "diagnostics.csv"):
        out.append((TMP / "pc20/apoapsis_pi" / f, f"apoapsis_pi/{f}"))
    for tag_dir in sorted((TMP / "pc20/mars/truth").iterdir()):
        for f in ("passes_3b.csv", "passes_zonal.csv", "passes_deg3.csv"):
            out.append((tag_dir / f, f"mars_truth/{tag_dir.name}/{f}"))
    for pl in ("mars", "earth"):
        root = TMP / "pcstagef/prod/stagef20" / pl
        out.append((root / "cells_base.csv", f"stagef20/{pl}/cells_base.csv"))
        for d in sorted((root / "runs").glob("open_mean_*")):
            out.append((d / "passes.csv", f"stagef20/{pl}/runs/{d.name}/passes.csv"))
    return out


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dest", default=str(DEFAULT)); dest = pathlib.Path(ap.parse_args().dest)
    rows = []
    for src, rel in sources():
        if not src.exists(): rows.append((None, rel, "MISSING", 0)); continue
        dst = dest / rel; dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() and sha(dst) != sha(src): raise SystemExit(f"hash differs, refusing to overwrite {dst}")
        if not dst.exists(): shutil.copy2(src, dst)
        n = sum(1 for _ in open(dst, "rb")) - (1 if dst.suffix == ".csv" else 0) if dst.suffix == ".csv" else 0
        rows.append((sha(dst), rel, "ok", n))
    (dest / "manifest.sha256").write_text("".join(f"{h}  {rel}\n" for h, rel, st, n in rows if h))
    (dest / "row_counts.txt").write_text("".join(f"{rel}\t{st}\t{n}\n" for h, rel, st, n in rows))
    print(f"{sum(1 for r in rows if r[0])} files, {sum(1 for r in rows if not r[0])} missing -> {dest}")
    for h, rel, st, n in rows:
        if not h: print("MISSING", rel)


if __name__ == "__main__":
    main()
