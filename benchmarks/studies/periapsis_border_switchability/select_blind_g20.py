"""Draw the blind cells of the >= 20x20 estimator-v3 test (rule and seeds in blind_g20_config.toml; run before any blind truth run).

Usage: python3 -I select_blind_g20.py --results-root BLIND_ROOT [--mars-config mars_generality_config_g20.toml]
Writes  <root>/mars/cells.csv            (tag, inc, raan, omega, epoch_idx)  and  <root>/mars_blind_config.toml  (the g20 Mars config with the new epochs)
        <root>/earth_cells_blind.csv     (tag, inc_deg, raan_deg, omega_deg, epoch_idx)  and  <root>/earth_epochs.txt (comma-separated ISO epochs for run_earth_g20.jl --epochs)
Disjointness from the development sets is asserted (different epochs, RAAN/omega off the development grids).
"""
import argparse, pathlib, re, sys, tomllib
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def draw_epochs(rng, base, n, lo, hi, excl):
    offs = set()
    while len(offs) < n:
        d = int(rng.integers(lo, hi + 1))
        if all(abs(d - x) > 3 for x in excl):
            offs.add(d)
    offs = sorted(offs); b = datetime.fromisoformat(base)
    return offs, [(b + timedelta(days=d)).strftime("%Y-%m-%dT%H:%M:%S") for d in offs]


def tag(inc, raan, om, e):
    return f"I{round(10 * inc):04d}_R{round(raan):03d}_W{round(om):03d}_E{e}"


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--mars-config", default=str(HERE / "mars_generality_config_g20.toml"))
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); root.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.load(open(HERE / "blind_g20_config.toml", "rb"))
    # ---- Mars
    m = cfg["mars"]; rng = np.random.default_rng(m["seed"])
    offs, epochs = draw_epochs(rng, m["base_epoch_utc"], m["n_epochs"], *m["epoch_offset_range_days"], excl=[0, 45])
    rows = []
    for inc, n in zip(m["inclination_deg"], m["cells_per_inclination"]):
        seen = set()
        while len(seen) < n:
            r, w = int(rng.integers(0, 360)), int(rng.integers(0, 360)); e = int(rng.integers(1, m["n_epochs"] + 1))
            if r % 90 == 0 or w % 90 == 0 or (r, w, e) in seen: continue
            seen.add((r, w, e)); rows.append(dict(tag=tag(inc, r, w, e), inc=inc, raan=float(r), omega=float(w), epoch_idx=e))
    C = pd.DataFrame(rows); assert C.tag.is_unique and len(C) == sum(m["cells_per_inclination"])
    (root / "mars").mkdir(exist_ok=True); C.to_csv(root / "mars" / "cells.csv", index=False)
    txt = open(a.mars_config).read(); new = "epochs_utc = [" + ", ".join(f'"{x}"' for x in epochs) + "]"
    txt2 = re.sub(r"^epochs_utc = \[.*\]$", new, txt, flags=re.M); assert txt2 != txt
    (root / "mars_blind_config.toml").write_text(txt2)
    dev = tomllib.load(open(a.mars_config, "rb"))["grid"]
    assert not set(epochs) & set(dev["epochs_utc"]) and not (set(C.raan) & set(dev["raan_deg"])) and not (set(C.omega) & set(dev["omega_deg"]))
    print(f"Mars: {len(C)} cells; epochs {epochs} (offsets {offs} d)"); print(C.groupby("inc").size().to_string())
    # ---- Earth
    e = cfg["earth"]; rng = np.random.default_rng(e["seed"])
    offs, epochs = draw_epochs(rng, e["base_epoch_utc"], e["n_epochs"], *e["epoch_offset_range_days"], excl=[0, 45, 91])
    rows = []
    for inc, n in zip(e["inclination_deg"], e["cells_per_inclination"]):
        seen = set()
        while len(seen) < n:
            r, w, ep = int(rng.integers(0, 360)), int(rng.integers(0, 360)), int(rng.integers(1, e["n_epochs"] + 1))
            if r % 45 == 0 or w % 90 == 0 or (r, w, ep) in seen: continue
            seen.add((r, w, ep)); rows.append(dict(tag=tag(inc, r, w, ep), inc_deg=inc, raan_deg=float(r), omega_deg=float(w), epoch_idx=ep))
    E = pd.DataFrame(rows); assert E.tag.is_unique and len(E) == sum(e["cells_per_inclination"])
    E.to_csv(root / "earth_cells_blind.csv", index=False); (root / "earth_epochs.txt").write_text(",".join(epochs))
    print(f"Earth: {len(E)} cells; epochs {epochs} (offsets {offs} d)"); print(E.groupby("inc_deg").size().to_string())


if __name__ == "__main__":
    main()
