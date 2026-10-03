"""Collect scripts/xval_fullarc.jl output into the paper's CSV deliverables.

Requires Python 3.11+ (tomllib), pandas, and pyarrow. Validation tests need only Python.

Usage:
    python3 scripts/xval_fullarc_export.py <xval outdir> <paper repo root>

Reads <outdir>/gmat_committed and <outdir>/stk_committed (the primary runs) and
any other <target>_<variant> run directories (sensitivity runs), and writes:

    <paper>/data/cross_validation_rmse_fullarc_source.csv
        one row per case x reference (primary runs only), including model/input identities
        and the digest of the completed run_info.toml. Legacy outputs must be rerun.
    <paper>/data/cross_validation_rmse_fullarc_variants.csv
        every run directory, including the sensitivity variants
    <paper>/figures/validation_results/cross_validation_error_series.csv
        per-epoch |dr| for the primary runs, decimated to at most ~1000 points
        per case by keeping every k-th sample (k = ceil(n / 1000)) plus the last
"""

import csv
import hashlib
import math
import os
from pathlib import Path
import re
import sys
import tomllib

PRIMARY = {"gmat_committed": "GMAT", "stk_committed": "STK"}
BODY_ORDER = ["earth", "mars", "venus", "moon"]
GRAV_ORDER = ["j0", "j2", "j50"]


def split_scenario(name):
    body, grav, tb = name.split("_")
    return body, grav, ("on" if tb == "tbtrue" else "off")


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_primary_run(run_dir, target):
    """Fail closed on legacy, sensitivity, incomplete, or mixed run artifacts.

    This checks execution identity, not the scientific provenance of references.
    It deliberately uses only the Python standard library so it can be tested
    without the optional pandas/pyarrow paper-export dependencies.
    """
    run_dir = Path(run_dir)

    def require(condition, message):
        if not condition:
            raise ValueError(f"{run_dir}: {message}; rerun an override-free committed campaign")

    info_path = run_dir / "run_info.toml"
    require(info_path.is_file(), "missing versioned primary provenance")
    with info_path.open("rb") as stream:
        info = tomllib.load(stream)
    require(type(info.get("provenance_schema")) is int and info["provenance_schema"] == 1,
            "unsupported or legacy provenance schema")
    require(info.get("status") == "complete", "campaign is not complete")
    require(info.get("target") == target and info.get("variant") == "committed",
            "run identity does not match the primary directory")
    require(info.get("primary_eligible") is True and info.get("input_overrides") == {},
            "model/reference overrides are not eligible for primary export")
    require(info.get("src_test_scripts_dirty") is False, "primary source tree was not clean")
    require(re.fullmatch(r"[0-9a-f]{40}", str(info.get("commit", ""))) is not None,
            "missing source commit")
    expected = {f"{body}_{gravity}_{tb}" for body in BODY_ORDER for gravity in GRAV_ORDER
                for tb in ("tbfalse", "tbtrue")}
    scenarios = info.get("scenarios", [])
    require(isinstance(scenarios, list) and len(scenarios) == len(expected)
            and all(isinstance(item, str) for item in scenarios) and set(scenarios) == expected,
            "primary campaign must contain each of the 24 scenarios exactly once")
    inputs = info.get("model_inputs", {})
    require(isinstance(inputs, dict) and bool(inputs.get("planetary_kernel"))
            and re.fullmatch(r"[0-9a-f]{64}", str(inputs.get("planetary_kernel_sha256", "")))
            and isinstance(inputs.get("solver_environment"), dict) and bool(inputs["solver_environment"]),
            "missing effective kernel or solver inputs")
    kernels = info.get("spice_kernels", [])
    require(isinstance(kernels, list) and bool(kernels)
            and all(isinstance(kernel, dict) and bool(kernel.get("path"))
                    and bool(kernel.get("kind")) and re.fullmatch(r"[0-9a-f]{64}", str(kernel.get("sha256", "")))
                    for kernel in kernels)
            and any(kernel["path"] == inputs["planetary_kernel"]
                    and kernel["sha256"] == inputs["planetary_kernel_sha256"] for kernel in kernels),
            "missing loaded-kernel identity or planetary-kernel mismatch")
    required = {"results.csv"} | {f"{case}/{name}" for case in expected
                                     for name in ("manifest.toml", "series.arrow")}
    artifacts = info.get("artifacts_sha256", {})
    require(isinstance(artifacts, dict) and set(artifacts) == required,
            "incomplete artifact identity")
    for relative in sorted(required):
        path = run_dir / relative
        require(path.is_file() and file_digest(path) == artifacts[relative],
                f"artifact differs from the completed run: {relative}")
    with (run_dir / "results.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == len(expected) and {r.get("scenario") for r in rows} == expected,
            "results do not contain the complete unique primary matrix")
    require(all(r.get("target") == target and r.get("variant") == "committed" for r in rows),
            "results contain a different target or variant")
    require(all(re.fullmatch(r"[0-9a-f]{64}", r.get(key, "")) for r in rows
                for key in ("reference_sha256", "gravity_file_sha256")),
            "missing reference or gravity-file identity")
    return {**info, "run_info_sha256": file_digest(info_path)}


def main(outdir, paper):
    # Validate before writing anything, including when an interrupted or legacy
    # run happens to leave exactly 48 rows behind.
    primary_info = {run: validate_primary_run(Path(outdir) / run, target.lower())
                    for run, target in PRIMARY.items()}
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
                "reference_sha256": r.get("reference_sha256", ""),
                "gravity_file_sha256": r.get("gravity_file_sha256", ""),
                "run_commit": primary_info[run]["commit"] if run in primary_info else "",
                "run_info_sha256": primary_info[run]["run_info_sha256"] if run in primary_info else "",
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
          "arc_end_s", "reference_file", "variant", "gravity_degree", "gravity_order",
          "gravity_file", "gm_override_m3s2", "third_bodies", "reference_sha256",
          "gravity_file_sha256", "run_commit", "run_info_sha256"]].to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_source.csv"), index=False)
    order(pd.DataFrame(variant_rows)).to_csv(
        os.path.join(paper, "data", "cross_validation_rmse_fullarc_variants.csv"), index=False)
    series = pd.concat(series_frames, ignore_index=True)
    series.to_csv(os.path.join(paper, "figures", "validation_results",
                               "cross_validation_error_series.csv"), index=False)
    print(f"wrote {len(prim)} primary rows, {len(variant_rows)} variant rows, {len(series)} series rows")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
