"""Draw the second blind Mars set of the estimator-v4 test (rule and seed in blind_v4_config.toml; run before any truth run on it).

Usage: python3 -I select_blind_v4.py --results-root V4BLIND_ROOT --exclude-root V3BLIND_ROOT [--mars-config mars_generality_config_g20.toml]
Writes  <root>/mars/cells.csv  (tag, inc, raan, omega, epoch_idx)  and  <root>/mars_blind_config.toml  (the g20 Mars config with the new epochs).
Reads the v3 blind set (<exclude-root>/mars/cells.csv and <exclude-root>/mars_blind_config.toml) only to keep the new cells disjoint from it.
Asserted: epochs differ from the development and v3-blind epochs by more than 3 d, RAAN and omega values are disjoint from the v3-blind values and are not multiples of 90.
"""
import argparse, pathlib, re, sys, tomllib
from datetime import datetime
import numpy as np
import pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
from select_blind_g20 import draw_epochs, tag


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--exclude-root", required=True)
    ap.add_argument("--mars-config", default=str(HERE / "mars_generality_config_g20.toml"))
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); root.mkdir(parents=True, exist_ok=True)
    m = tomllib.load(open(HERE / "blind_v4_config.toml", "rb"))["mars"]; rng = np.random.default_rng(m["seed"])
    ex = pathlib.Path(a.exclude_root); old = pd.read_csv(ex / "mars" / "cells.csv"); old_ep = tomllib.load(open(ex / "mars_blind_config.toml", "rb"))["grid"]["epochs_utc"]
    dev_ep = tomllib.load(open(a.mars_config, "rb"))["grid"]["epochs_utc"]
    base = datetime.fromisoformat(m["base_epoch_utc"])
    excl = [0, 45] + [round((datetime.fromisoformat(x) - base).total_seconds() / 86400) for x in old_ep]
    offs, epochs = draw_epochs(rng, m["base_epoch_utc"], m["n_epochs"], *m["epoch_offset_range_days"], excl=excl)
    rows = []; used_r, used_w = set(old.raan.astype(int)), set(old.omega.astype(int))
    for inc, n in zip(m["inclination_deg"], m["cells_per_inclination"]):
        seen = set()
        while len(seen) < n:
            r, w, e = int(rng.integers(0, 360)), int(rng.integers(0, 360)), int(rng.integers(1, m["n_epochs"] + 1))
            if r % 90 == 0 or w % 90 == 0 or r in used_r or w in used_w or (r, w, e) in seen: continue
            seen.add((r, w, e)); rows.append(dict(tag=tag(inc, r, w, e), inc=inc, raan=float(r), omega=float(w), epoch_idx=e))
    C = pd.DataFrame(rows); assert C.tag.is_unique and len(C) == sum(m["cells_per_inclination"])
    assert not (set(C.raan.astype(int)) & set(old.raan.astype(int))) and not (set(C.omega.astype(int)) & set(old.omega.astype(int)))
    assert not (set(epochs) & (set(old_ep) | set(dev_ep)))
    (root / "mars").mkdir(exist_ok=True); C.to_csv(root / "mars" / "cells.csv", index=False)
    txt = open(a.mars_config).read(); new = "epochs_utc = [" + ", ".join(f'"{x}"' for x in epochs) + "]"
    txt2 = re.sub(r"^epochs_utc = \[.*\]$", new, txt, flags=re.M); assert txt2 != txt
    (root / "mars_blind_config.toml").write_text(txt2)
    print(f"Mars: {len(C)} cells; epochs {epochs} (offsets {offs} d; excluded offsets {excl})"); print(C.groupby("inc").size().to_string())


if __name__ == "__main__":
    main()
