"""Stage F at 20x20 cell draw (copy of stagef_cells.py; same seed and grid, so the cells equal the Stage F cells) (written and run BEFORE any production run). For each planet, epoch and inclination, draw n_per_group distinct (RAAN, omega)
combinations with numpy default_rng(cell_seed); the first n_dev_per_group of each group are development cells, the rest held out.
Usage: python3 -I stagef_cells.py --results-root DIR -> DIR/stagef/<planet>/cells_base.csv"""
import argparse, pathlib, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
C = tomllib.load(open(HERE / "stagef20_config.toml", "rb"))
S = C["stagef"]

def tag(i, r, w, e): return f"I{round(10*i):04d}_R{round(r):03d}_W{round(w):03d}_E{e}"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); a = ap.parse_args()
    rng = np.random.default_rng(S["cell_seed"])
    for pl in S["planets"]:
        g = C["planets"][pl]["grid"]; combos = [(r, w) for r in g["raan_deg"] for w in g["omega_deg"]]; rows = []
        for e in range(1, S["n_epochs"] + 1):
            for inc in g["inclination_deg"]:
                pick = rng.choice(len(combos), size=S["n_per_group"], replace=False)
                for k, j in enumerate(pick):
                    r, w = combos[j]
                    rows.append(dict(tag=tag(inc, r, w, e), inc=inc, raan=r, omega=w, epoch_idx=e, split="dev" if k < S["n_dev_per_group"] else "heldout"))
        out = pathlib.Path(a.results_root) / S["results_subdir"] / pl; out.mkdir(parents=True, exist_ok=True)
        D = pd.DataFrame(rows); D.to_csv(out / "cells_base.csv", index=False); print(pl, len(D), D.split.value_counts().to_dict())
main()
