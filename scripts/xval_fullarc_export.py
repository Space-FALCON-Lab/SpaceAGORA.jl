"""Collect scripts/xval_fullarc.jl output into the paper's CSV deliverables.

Requires Python 3.11+ (stdlib tomllib), or tomli on older Python.
The CSV/Arrow export additionally requires pandas and pyarrow.

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

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError as exc:
        raise ImportError("Manifest validation requires Python 3.11+ or the tomli package") from exc


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


def declared_input(declared_path, inputs, label):
    """Bind a declared file to one recorded identity without requiring old inputs.

    Schema 1 records absolute input paths but no original repository root. For a
    relative declaration, compare whole trailing path components, independently
    of the exporter's current directory. Ambiguous matches and parent traversal
    are rejected. This checks recorded consistency, not unavailable input bytes.
    """
    require(isinstance(declared_path, str) and declared_path and "\0" not in declared_path,
            f"Invalid {label} path")
    path = Path(declared_path)
    if path.is_absolute():
        normalized = os.path.normpath(declared_path)
        matches = [item for item in inputs if os.path.normpath(item["path"]) == normalized]
    else:
        require(path.parts and ".." not in path.parts, f"Invalid relative {label} path")
        matches = [item for item in inputs if Path(item["path"]).parts[-len(path.parts):] == path.parts]
    require(len(matches) == 1, f"Missing or ambiguous {label} input identity")
    return matches[0]


def validate_configuration(case, row, manifest_path):
    """Reconcile the producer's effective configuration and retained evidence."""
    model = case["model"]
    model_keys = {"name", "planet", "gravity_model", "gravity_harmonics_degree",
                  "gravity_harmonics_order", "gravity_harmonics_file", "nbody_bodies",
                  "orbit_altitude_mode", "kind", "comparison_mode", "events", "telemetry",
                  "telemetry_columns", "max_points_quick", "max_points_full", "min_eval_points",
                  "units", "tolerances_quick", "tolerances_full", "initial_time", "spacecraft",
                  "atmosphere_truth", "calibration", "drag_enabled", "EI_km"}
    require(isinstance(model, dict) and model_keys <= model.keys(),
            "Missing substantive model configuration")
    for key in ("name", "planet", "gravity_model", "orbit_altitude_mode", "kind",
                "comparison_mode", "telemetry"):
        require(isinstance(model[key], str) and bool(model[key]), f"Invalid model {key}")
    for key in ("telemetry_columns", "units", "tolerances_quick", "tolerances_full",
                "initial_time", "spacecraft", "atmosphere_truth", "calibration"):
        require(isinstance(model[key], dict) and bool(model[key]), f"Invalid model {key}")
    for key in ("max_points_quick", "max_points_full", "min_eval_points"):
        require(type(model[key]) is int and model[key] > 0, f"Invalid model {key}")
    require(type(model["drag_enabled"]) is bool, "Invalid model drag_enabled")
    require(type(model["EI_km"]) in (int, float) and math.isfinite(model["EI_km"]) and model["EI_km"] >= 0,
            "Invalid model EI_km")
    events = model["events"]
    require(isinstance(events, list) and events and all(isinstance(event, str) and event for event in events),
            "Invalid model events")
    require(model["name"] == row["scenario"] and model["planet"] == split_scenario(row["scenario"])[0],
            "Model scenario/planet mismatch")
    degree, order = model["gravity_harmonics_degree"], model["gravity_harmonics_order"]
    require(type(degree) is int and type(order) is int and 0 <= order <= degree,
            "Invalid model gravity degree/order")
    field = model["gravity_harmonics_file"]
    require(isinstance(field, str), "Invalid model gravity file")
    bodies = model["nbody_bodies"]
    require(isinstance(bodies, list) and all(isinstance(body, str) and body for body in bodies),
            "Invalid model third bodies")
    with manifest_path.open("rb") as stream:
        manifest = tomllib.load(stream)
    scenarios = manifest.get("scenarios")
    require(isinstance(scenarios, list) and len(scenarios) == 1 and scenarios[0] == model,
            "Manifest and recorded model differ")
    require(int(row["gravity_degree"]) == degree and int(row["gravity_order"]) == order,
            "Result and model gravity degree/order differ")
    require(row["gravity_file"] == field, "Result and model gravity file differ")
    require(row["nbody_bodies"] == "+".join(bodies), "Result and model third bodies differ")
    gm_key = "gravity_harmonics_gm_override_m3s2"
    row_gm = float(row["gm_override_m3s2"]) if row["gm_override_m3s2"] else math.nan
    if gm_key in model:
        gm = model[gm_key]
        require(type(gm) in (int, float) and math.isfinite(gm) and gm > 0,
                "Invalid model GM override")
        require(row_gm == gm, "Result and model GM override differ")
    else:
        require(math.isnan(row_gm), "Result has an unrecorded GM override")
    # Telemetry verification needs the harmonics file for nonzero degree OR an
    # explicit GM override. Only a point mass without that override is exempt.
    if field:
        declared_input(field, case["inputs"], "gravity file")
    else:
        require(degree == 0 and gm_key not in model, "Model requires a declared gravity file")

    environment = case["solver_environment"]
    numeric_keys = {"SPACEAGORA_TELEMETRY_DT_MAX_ORBIT", "SPACEAGORA_TELEMETRY_RELTOL_ORBIT",
                    "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT", "SPACEAGORA_TELEMETRY_RELTOL_ATM",
                    "SPACEAGORA_TELEMETRY_ABSTOL_ATM"}
    mode_key = "SPACEAGORA_TELEMETRY_SOLVER_MODE"
    kernel_key = "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH"
    require(isinstance(environment, dict) and numeric_keys | {mode_key, kernel_key} <= environment.keys(),
            "Missing substantive solver configuration")
    for key in numeric_keys | {mode_key, kernel_key}:
        require(isinstance(environment[key], str) and bool(environment[key]), f"Invalid solver {key}")
    for key in numeric_keys:
        value = float(environment[key])
        require(math.isfinite(value) and value > 0, f"Invalid solver {key}")
    require(declared_input(environment[kernel_key], case["inputs"], "planetary kernel") ==
            case["planetary_kernel"], "Solver planetary-kernel selector mismatch")


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
            require(isinstance(item["path"], str) and Path(item["path"]).is_absolute() and
                    "\0" not in item["path"] and isinstance(item["sha256"], str) and
                    re.fullmatch(r"[a-f0-9]{64}", item["sha256"]),
                    "Invalid input identity")
        require(case["inputs"][0]["path"] == row["reference_path"], "Reference identity mismatch")
        require(case["planetary_kernel"] in case["inputs"] and
                case["planetary_kernel"] in case["kernels"], "Planetary-kernel identity mismatch")
        validate_configuration(case, row, runpath / row["scenario"] / "manifest.toml")
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
