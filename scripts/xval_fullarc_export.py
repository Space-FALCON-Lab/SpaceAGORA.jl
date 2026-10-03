"""Collect scripts/xval_fullarc.jl output into the paper's CSV deliverables.

Usage:
    python3 scripts/xval_fullarc_export.py <xval outdir> <paper repo root>

Validates run_info.json and retained output hashes before accepting
<outdir>/gmat_committed and <outdir>/stk_committed as primary runs, and reads
any other <target>_<variant> run directories (sensitivity runs), and writes:

    <paper>/data/cross_validation_rmse_fullarc_source.csv
        one row per case x reference (primary runs only)
    <paper>/data/cross_validation_rmse_fullarc_variants.csv
        every run directory, including the sensitivity variants
    <paper>/figures/validation_results/cross_validation_error_series.csv
        per-epoch |dr| for the primary runs, decimated to at most ~1000 points
        per case by keeping every k-th sample (k = ceil(n / 1000)) plus the last
"""

import csv
import hashlib
import json
import math
from pathlib import Path
import re
import os
import sys


PRIMARY = {"gmat_committed": "GMAT", "stk_committed": "STK"}
BODY_ORDER = ["earth", "mars", "venus", "moon"]
GRAV_ORDER = ["j0", "j2", "j50"]


def split_scenario(name):
    body, grav, tb = name.split("_")
    return body, grav, ("on" if tb == "tbtrue" else "off")


# These affect selection or logging, not the model or reference inputs.
PRIMARY_CONTROLS = {"XVAL_SCENARIOS", "SPACEAGORA_WARN_NORMALIZE",
                    "SPACEAGORA_WARN_DEPRECATED_CONFIG"}
SCENARIOS = {f"{b}_{g}_{tb}" for b in BODY_ORDER for g in GRAV_ORDER
             for tb in ("tbfalse", "tbtrue")}


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_run(runpath, *, primary=False):
    """Validate recorded identity and retained outputs before writing any export.

    A hash proves file identity, not the reference tool's generation settings.
    Those settings remain explicitly unverified in every exported row.
    """
    runpath = Path(runpath)
    with (runpath / "run_info.json").open() as stream:
        info = json.load(stream)
    require(info.get("schema_version") == 1 and info.get("status") == "complete",
            f"{runpath.name}: missing or incomplete supported run record")
    target, variant = info["target"], info["variant"]
    require(runpath.name == f"{target}_{variant}", "Run directory and recorded target/variant differ")
    require(info["source"] == info["source_end"], "Source changed during run")
    require(info["reference_provenance"] == "unverified_generation_settings",
            "Unknown reference-provenance status; obtain a separately reviewed provenance record")
    require(info["results_sha256"] == digest(runpath / "results.csv"), "Results checksum mismatch")
    with (runpath / "results.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    names = [row["scenario"] for row in rows]
    require(bool(names) and len(names) == len(set(names)), "Empty or duplicate scenarios")
    require(set(names) <= SCENARIOS and set(names) == set(info["scenarios"]) == set(info["cases"]),
            "Scenario coverage disagrees with run record")
    for row in rows:
        require(row["target"] == target and row["variant"] == variant, "Mislabelled result row")
        case = info["cases"][row["scenario"]]
        for suffix, key in (("manifest.toml", "manifest_sha256"), ("series.arrow", "series_sha256")):
            require(case[key] == digest(runpath / row["scenario"] / suffix), f"{suffix} checksum mismatch")
        require(bool(case["inputs"]) and bool(case["kernels"]), "Missing input/kernel identities")
        for item in case["inputs"] + case["kernels"]:
            require(bool(item["path"]) and re.fullmatch(r"[a-f0-9]{64}", item["sha256"]),
                    "Invalid input identity")
        require(case["inputs"][0]["path"] == row["reference_path"], "Reference identity mismatch")
        require(case["planetary_kernel"] in case["inputs"] and
                case["planetary_kernel"] in case["kernels"], "Planetary-kernel identity mismatch")
        require(bool(case["solver_environment"]) and bool(case["model"]), "Missing effective configuration")
        require(case["solver_retcode"] == row["solver_retcode"], "Solver status mismatch")
        if primary:
            require(row["solver_retcode"] in ("Success", "ReturnCode.Success"), "Unsuccessful primary solve")
    if primary:
        require(target in ("gmat", "stk") and variant == "committed", "Sensitivity run cannot be primary")
        require(info["source"]["dirty"] is False, "Dirty source cannot be primary")
        for key in ("commit", "tree"):
            require(re.fullmatch(r"[a-f0-9]{40}", info["source"][key]), "Missing source identity")
        forbidden = set(info["controls"]) - PRIMARY_CONTROLS
        require(not forbidden, f"Scientific/environment overrides cannot be primary: {sorted(forbidden)}")
        require(set(names) == SCENARIOS, "Primary run must contain the exact 24-case matrix")
    return info


def validate_campaign(outdir):
    outdir = Path(outdir)
    for name in PRIMARY:
        require((outdir / name / "results.csv").is_file(), f"Missing primary run {name}")
    records = {p.name: validate_run(p, primary=p.name in PRIMARY)
               for p in sorted(outdir.iterdir()) if (p / "results.csv").is_file()}
    sources = [records[name]["source"] for name in PRIMARY]
    require(sources[0] == sources[1], "Primary runs use different source revisions")
    return records


def main(outdir, paper):
    records = validate_campaign(outdir)  # Refuse before creating any deliverable.
    import pandas as pd
    import pyarrow.feather as feather
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
                "variant": run,  # run directory: <target>_<variant>
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
                "source_commit": records[run]["source"]["commit"],
                "source_dirty": records[run]["source"]["dirty"],
                "run_record_sha256": digest(Path(outdir) / run / "run_info.json"),
                "reference_sha256": records[run]["cases"][r["scenario"]]["inputs"][0]["sha256"],
                "reference_provenance": records[run]["reference_provenance"],
            }
            variant_rows.append(row)
            if run in PRIMARY:
                primary_rows.append(row)
                s = feather.read_table(os.path.join(outdir, run, r["scenario"], "series.arrow")).to_pandas()
                require(len(s) == int(r["n_points"]) and len(s) > 0, "Series sample count mismatch")
                require(float(s["t_s"].iloc[-1]) == float(r["arc_end_s"]), "Series endpoint mismatch")
                k = max(1, math.ceil(len(s) / 1000))
                idx = list(range(0, len(s), k))
                if idx[-1] != len(s) - 1:
                    idx.append(len(s) - 1)
                d = s.iloc[idx]
                series_frames.append(pd.DataFrame({
                    "body": body, "gravity": grav, "third_body": tb,
                    "reference": PRIMARY[run], "variant": run,
                    "source_commit": row["source_commit"],
                    "run_record_sha256": row["run_record_sha256"],
                    "reference_provenance": row["reference_provenance"],
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
    prim.to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_source.csv"), index=False)
    order(pd.DataFrame(variant_rows)).to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_variants.csv"), index=False)
    series = pd.concat(series_frames, ignore_index=True)
    series.to_csv(os.path.join(paper, "figures", "validation_results",
                               "cross_validation_error_series.csv"), index=False)
    print(f"wrote {len(prim)} primary rows, {len(variant_rows)} variant rows, {len(series)} series rows")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
