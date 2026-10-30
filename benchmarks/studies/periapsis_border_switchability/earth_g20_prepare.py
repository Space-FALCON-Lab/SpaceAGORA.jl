"""Lay the >= 20x20 Earth truth (run_earth_g20.jl, <root>/earth_g20/grid/<tag>) out in the directory structure the earlier analysis scripts read, so they run
unchanged on the new results root: <root>/border/grid/R..W..E.. (50 deg), <root>/border2/grid/I... (the other inclinations), <root>/border2/ref50,
static_cells.csv, static_run_summary.csv. The race and direct-step files exist only for the earlier 3x0 truth; oracle_run.py reads them to build a measured-lever
table that the cell labels do not use, so they are linked from --old-root (the labels, A1, are computed from the new truth only).

Usage: python3 -I earth_g20_prepare.py --results-root NEW --old-root OLD
"""
import argparse, os, pathlib, sys
import pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def link(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    os.symlink(src, dst)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--old-root", required=True); ap.add_argument("--blind", action="store_true", help="layout for the blind Earth cells (earth_cells_blind.csv); no 50 deg cells, no race/direct-step links")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)).resolve(); old = pathlib.Path(a.old_root).resolve()
    S = pd.read_csv(root / "earth_g20" / "grid_run_summary.csv"); cells = pd.read_csv(root / ("earth_cells_blind.csv" if a.blind else "earth_cells_336.csv"))
    assert len(S) == len(cells) and set(S.tag) == set(cells.tag), f"expected {len(cells)} completed cells, found {len(S)}"
    fam = pd.read_csv(old / "border2" / "static_cells.csv").set_index("tag").family
    for t in S.tag:
        src = root / "earth_g20" / "grid" / t
        if t.startswith("I0500_"):
            link(src, root / "border" / "grid" / t[6:]); link(src, root / "border2" / "ref50" / t[6:])
        else:
            link(src, root / "border2" / "grid" / t)
    sc = S[["tag", "inc_deg", "raan_deg", "omega_deg", "epoch_idx"]].copy(); sc["family"] = sc.tag.map(fam) if not a.blind else "blind"
    sc["ended_early"] = (S.retcode_3b != "Success").values
    (root / "border2").mkdir(exist_ok=True); sc.to_csv(root / "border2" / "static_cells.csv", index=False); S.to_csv(root / "border2" / "static_run_summary.csv", index=False)
    if a.blind:
        print(f"linked {len(S)} blind cells; status counts: {S.retcode_3b.value_counts().to_dict()}"); return
    for f in ("race_cells.csv", "race_run_summary.csv"):
        link(old / "border2" / f, root / "border2" / f)
    link(old / "border2" / "race", root / "border2" / "race"); link(old / "border" / "direct_step.csv", root / "border" / "direct_step.csv")
    print(f"linked {len(S)} cells; status counts: {S.retcode_3b.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
