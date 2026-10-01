"""Collect scripts/xval_fullarc.jl output into the paper's CSV deliverables.

Usage:
    python3 scripts/xval_fullarc_export.py <xval outdir> <paper repo root>

Reads <outdir>/gmat_committed and <outdir>/stk_committed (the primary runs) and
any other <target>_<variant> run directories (sensitivity runs), and writes:

    <paper>/data/cross_validation_rmse_fullarc_source.csv
        one row per case x reference (primary runs only)
    <paper>/data/cross_validation_rmse_fullarc_variants.csv
        every run directory, including the sensitivity variants
    <paper>/figures/validation_results/cross_validation_error_series.csv
        per-epoch |dr| for the primary runs, decimated to at most ~1000 points
        per case by keeping every k-th sample (k = ceil(n / 1000)) plus the last
"""

import math
import os
import sys

import pandas as pd
import pyarrow.feather as feather

PRIMARY = {"gmat_committed": "GMAT", "stk_committed": "STK"}
BODY_ORDER = ["earth", "mars", "venus", "moon"]
GRAV_ORDER = ["j0", "j2", "j50"]


def split_scenario(name):
    body, grav, tb = name.split("_")
    return body, grav, ("on" if tb == "tbtrue" else "off")


def main(outdir, paper):
    variant_rows = []
    primary_rows = []
    series_frames = []
    for run in sorted(os.listdir(outdir)):
        res_path = os.path.join(outdir, run, "results.csv")
        if not os.path.isfile(res_path):
            continue
        res = pd.read_csv(res_path)
        for _, r in res.iterrows():
            body, grav, tb = split_scenario(r["scenario"])
            row = {
                "body": body,
                "gravity": grav,
                "third_body": tb,
                "reference": r["target"].upper(),
                "variant": r["variant"],
                "rms_m": repr(float(r["rms_m"])),
                "max_m": repr(float(r["max_m"])),
                "n_points": int(r["n_points"]),
                "arc_end_s": repr(float(r["arc_end_s"])),
                "rms_first10k_m": repr(float(r["rms_first10k_m"])),
                "gravity_degree": int(r["gravity_degree"]),
                "gravity_order": int(r["gravity_order"]),
                "gravity_file": r["gravity_file"],
                "gm_override_m3s2": r["gm_override_m3s2"],
                "third_bodies": "" if pd.isna(r["nbody_bodies"]) else r["nbody_bodies"],
                "reference_file": r["reference_path"],
            }
            variant_rows.append(row)
            if run in PRIMARY:
                primary_rows.append(row)
                s = feather.read_table(os.path.join(outdir, run, r["scenario"], "series.arrow")).to_pandas()
                k = max(1, math.ceil(len(s) / 1000))
                idx = list(range(0, len(s), k))
                if idx[-1] != len(s) - 1:
                    idx.append(len(s) - 1)
                d = s.iloc[idx]
                series_frames.append(pd.DataFrame({
                    "body": body, "gravity": grav, "third_body": tb,
                    "reference": PRIMARY[run],
                    "t_s": [repr(float(x)) for x in d["t_s"]],
                    "err_m": [repr(float(x)) for x in d["err_m"]],
                }))

    def order(df):
        df = df.copy()
        df["_b"] = df["body"].map(BODY_ORDER.index)
        df["_g"] = df["gravity"].map(GRAV_ORDER.index)
        sort_cols = ["_b", "_g", "third_body", "reference"] + (["variant"] if "variant" in df else [])
        return df.sort_values(sort_cols, kind="stable").drop(columns=["_b", "_g"])

    prim = order(pd.DataFrame(primary_rows))
    if len(prim) != 48:
        raise SystemExit(f"expected 48 primary rows, got {len(prim)}")
    prim[["body", "gravity", "third_body", "reference", "rms_m", "max_m", "n_points",
          "arc_end_s", "reference_file"]].to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_source.csv"), index=False)
    order(pd.DataFrame(variant_rows)).to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_variants.csv"), index=False)
    series = pd.concat(series_frames, ignore_index=True)
    series.to_csv(os.path.join(paper, "figures", "validation_results",
                               "cross_validation_error_series.csv"), index=False)
    print(f"wrote {len(prim)} primary rows, {len(variant_rows)} variant rows, {len(series)} series rows")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
