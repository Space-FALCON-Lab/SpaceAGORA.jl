#!/usr/bin/env python3
"""Archive a benchmark run into the paper repository's ``data/raw`` tree.

Benchmark runs land in gitignored ``output/`` directories inside simulator
worktrees, which means a reboot, a ``git clean`` or a deleted worktree takes
them with it.  Every figure in the manuscript has to be regenerable from files
the paper repository actually carries, so each run gets copied into

    <paper repo>/data/raw/<run_id>/

with a ``manifest.toml`` describing it, a row in ``index.csv`` and an entry in
``PROVENANCE.md``.

``run_id`` is ``<machine>_<harness>_<store>_<YYYYMMDD_HHMMSS>``:

* ``machine``  -- ``trx50``, ``workstation`` or ``macbook``
* ``harness``  -- ``ppb`` (benchmarks/studies/paper_parallelization_benchmarks)
                  or ``ps`` (benchmarks/studies/paper_scenarios)
* ``store``    -- calibration-store state: ``cold``, ``converged``, ``mixed``
                  or ``na`` for harnesses that do not use the store
* stamp        -- from the run directory name, else from the run's timestamps

Every manifest field is derived from the run's own CSV columns wherever the
harness records it.  A field the data does not contain is written as
``"unknown"`` rather than guessed.  The two host facts the CSVs never record
(CPU model, physical core count, installed memory) come from the table in
``HOST_FACTS`` below, and each entry names the file it was read from; a host
that is not in the table gets ``"unknown"``.

Usage::

    python3 scripts/archive_paper_run.py <run_dir> --archive <paper>/data/raw \\
        --store converged [--machine trx50] [--notes "..."] [--dry-run]

    python3 scripts/archive_paper_run.py --verify --archive <paper>/data/raw

``--store`` is required for ``ppb`` runs: the calibration store's state is not
recorded in the CSV and only the person who launched the run knows it.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import os
import re
import shutil
import sys

import pandas as pd

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - older interpreters
    tomllib = None


# --------------------------------------------------------------------------
# Host facts the harness CSVs do not record.
#
# `cpu_threads` is deliberately absent here: every harness records it per row,
# so it is derived from the data instead.  Each entry names where the value was
# read from, and that source string is copied into the manifest.
# --------------------------------------------------------------------------
HOST_FACTS = {
    "space-falcon-1": {
        "machine": "workstation",
        "cpu_model": "AMD Ryzen 9 9900X 12-Core Processor",
        "physical_cores": 12,
        "memory_gb": 60,
        "source": (
            "lscpu on the host; "
            "docs/architecture/paper_routing_figures_20260915.md "
            '("space-falcon-1, Ryzen 9 9900X, 12 physical cores / 24 threads, 60 GB")'
        ),
    },
    "space-falcon-lab-TRX50-AERO-D": {
        "machine": "trx50",
        "cpu_model": "AMD Ryzen Threadripper PRO 9985WX",
        "physical_cores": 64,
        "memory_gb": 250,
        "source": (
            "docs/adaptive_policy_r6_nested_campaign_record_20260906.md "
            '("TRX50 (64-core Threadripper PRO 9985WX, 250 GB)"); '
            "physical core count also in docs/architecture/paper_routing_figures_20260915.md"
        ),
    },
}

MACHINE_NAMES = {"trx50", "workstation", "macbook"}
STORE_STATES = {"cold", "converged", "mixed", "na"}

PPB_RAW_GLOB = "paper_benchmarks_raw_*.csv"
PPB_AGG_GLOB = "paper_benchmarks_aggregated_*.csv"
PS_CSV_RE = re.compile(r"^s[1-9]_.*\.csv$")

STAMP_RE = re.compile(r"(\d{8}_\d{6})")

# The paper repository is shared with a co-author and with Overleaf, and the
# local tooling layout of this checkout is not part of a measurement's
# provenance.  Hidden per-tool worktree directories are therefore elided from
# the paths recorded there; everything else in the path is kept verbatim.
AGENT_WORKTREE_RE = re.compile(r"/\.[A-Za-z0-9_.-]+/worktrees/")
AGENT_WORKTREE_PLACEHOLDER = "/<agent-worktrees>/"


def recorded_path(path: str) -> str:
    """The absolute source path as it is written into the archive."""
    return AGENT_WORKTREE_RE.sub(AGENT_WORKTREE_PLACEHOLDER, os.path.abspath(path))


INDEX_COLUMNS = [
    "run_id",
    "machine",
    "harness",
    "store",
    "launched_utc",
    "simulator_commit",
    "phases",
    "n_rows",
    "path",
]

MANIFEST_ORDER = [
    "run_id",
    "machine",
    "hostname",
    "cpu_model",
    "physical_cores",
    "cpu_threads",
    "memory_gb",
    "harness",
    "simulator_commit",
    "simulator_branch",
    "phases",
    "cases",
    "modes",
    "repeats",
    "thread_counts",
    "process_workers",
    "mc_samples",
    "calibration_store",
    "launched_utc",
    "finished_utc",
    "source_path",
    "command",
    "notes",
    # Derived extras, beyond the required set.
    "simulator_commits",
    "n_rows",
    "raw_csv",
    "aggregate_csv",
    "timestamps_source",
    "host_facts_source",
    "archived_utc",
]

UNKNOWN = "unknown"


# --------------------------------------------------------------------------
# Small TOML writer.  tomli_w is not installed here, and the manifests only
# need strings, integers and flat lists of those.
# --------------------------------------------------------------------------
def _toml_escape(text: str) -> str:
    out = []
    for ch in text:
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif ord(ch) < 0x20:
            out.append("\\u%04x" % ord(ch))
        else:
            out.append(ch)
    return "".join(out)


def _toml_value(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    return '"' + _toml_escape(str(value)) + '"'


def write_toml(path: str, data: dict) -> None:
    lines = [
        "# Generated by scripts/archive_paper_run.py. Every field is derived from",
        "# the run's own files except cpu_model, physical_cores and memory_gb,",
        "# whose source is named in host_facts_source.",
        "",
    ]
    for key in MANIFEST_ORDER:
        if key in data:
            lines.append(f"{key} = {_toml_value(data[key])}")
    for key in sorted(k for k in data if k not in MANIFEST_ORDER):
        lines.append(f"{key} = {_toml_value(data[key])}")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def read_toml(path: str) -> dict:
    if tomllib is None:
        raise SystemExit("tomllib is unavailable; Python 3.11 or newer is required to verify")
    with open(path, "rb") as handle:
        return tomllib.load(handle)


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------
def find_files(run_dir: str, pattern: str) -> list[str]:
    return sorted(glob.glob(os.path.join(run_dir, "**", pattern), recursive=True))


def detect_harness(run_dir: str) -> str:
    if find_files(run_dir, PPB_RAW_GLOB):
        return "ppb"
    for name in sorted(os.listdir(run_dir)):
        if PS_CSV_RE.match(name):
            return "ps"
    # A run directory whose raw CSV was renamed by hand (the `_cold_store`
    # siblings are the known case) still belongs to the ppb harness.
    if find_files(run_dir, "*_raw*.csv"):
        return "ppb"
    raise SystemExit(f"cannot tell which harness produced {run_dir}")


def ppb_csvs(run_dir: str, raw_override: str | None) -> tuple[str, str | None]:
    if raw_override:
        raw = os.path.abspath(raw_override)
        if not os.path.isfile(raw):
            raise SystemExit(f"--raw {raw} does not exist")
    else:
        hits = find_files(run_dir, PPB_RAW_GLOB)
        if not hits:
            hits = [p for p in find_files(run_dir, "*_raw*.csv") if os.path.dirname(p) == run_dir]
        if not hits:
            raise SystemExit(f"no raw benchmark CSV found under {run_dir}")
        if len(hits) > 1:
            listed = "\n  ".join(hits)
            raise SystemExit(f"several raw CSVs under {run_dir}; pass --raw:\n  {listed}")
        raw = hits[0]

    agg_hits = find_files(run_dir, PPB_AGG_GLOB)
    if not agg_hits:
        base = os.path.basename(raw)
        guess = os.path.join(os.path.dirname(raw), base.replace("_raw", "_aggregated"))
        agg_hits = [guess] if os.path.isfile(guess) else []
    agg = agg_hits[0] if len(agg_hits) == 1 else (agg_hits[0] if agg_hits else None)
    return raw, agg


def ps_csvs(run_dir: str) -> list[str]:
    hits = [
        os.path.join(run_dir, name)
        for name in sorted(os.listdir(run_dir))
        if PS_CSV_RE.match(name)
    ]
    if not hits:
        raise SystemExit(f"no scenario CSVs (s<N>_*.csv) found in {run_dir}")
    return hits


# --------------------------------------------------------------------------
# Derivation helpers
# --------------------------------------------------------------------------
def distinct(frame: pd.DataFrame, column: str):
    if column not in frame.columns:
        return None
    values = frame[column].dropna().unique().tolist()
    cleaned = []
    for value in values:
        if isinstance(value, float) and value.is_integer():
            cleaned.append(int(value))
        elif hasattr(value, "item"):
            cleaned.append(value.item())
        else:
            cleaned.append(value)
    try:
        return sorted(cleaned)
    except TypeError:
        return sorted(str(v) for v in cleaned)


def distinct_or_unknown(frame: pd.DataFrame, column: str):
    values = distinct(frame, column)
    return UNKNOWN if values is None else values


def first_column(frame: pd.DataFrame, candidates) -> str | None:
    for name in candidates:
        if name in frame.columns:
            return name
    return None


def utc_stamp(epoch: float) -> str:
    return dt.datetime.fromtimestamp(epoch, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")


def newest_and_oldest(paths) -> tuple[float, float]:
    times = [os.path.getmtime(p) for p in paths]
    return min(times), max(times)


def stamp_from(run_dir: str, frame: pd.DataFrame, paths) -> str:
    match = STAMP_RE.search(os.path.basename(run_dir))
    if match:
        return match.group(1)
    if "timestamp_utc" in frame.columns and not frame["timestamp_utc"].dropna().empty:
        first = str(sorted(frame["timestamp_utc"].dropna().astype(str))[0])
        cleaned = re.sub(r"[^0-9]", "", first)
        if len(cleaned) >= 14:
            return cleaned[:8] + "_" + cleaned[8:14]
    oldest, _ = newest_and_oldest(paths)
    return dt.datetime.fromtimestamp(oldest, dt.timezone.utc).strftime("%Y%m%d_%H%M%S")


def commit_fields(frame: pd.DataFrame) -> tuple[str, list[str]]:
    if "git_commit" not in frame.columns:
        return UNKNOWN, []
    commits = sorted(str(c) for c in frame["git_commit"].dropna().unique())
    # The harness writes the literal string "unknown" when it cannot read the
    # hash, which is not a commit.
    commits = [c for c in commits if c and c.lower() != UNKNOWN]
    if not commits:
        return UNKNOWN, []
    if len(commits) == 1:
        return commits[0], commits
    return "multiple", commits


def repeat_counts_ppb(frame: pd.DataFrame) -> list[int]:
    keys = [
        c
        for c in ("phase_id", "case", "mode", "thread_count", "process_workers", "mc_samples")
        if c in frame.columns
    ]
    if not keys:
        return distinct(frame, "repeat") or []
    counts = frame.groupby(keys, dropna=False).size().unique().tolist()
    return sorted(int(c) for c in counts)


def repeat_counts_ps(frame: pd.DataFrame) -> list[int]:
    if "times_s" not in frame.columns:
        return []
    counts = {
        len(str(value).split("|"))
        for value in frame["times_s"].dropna()
        if str(value).strip()
    }
    return sorted(counts)


def command_from_logs(run_dir: str) -> str:
    """Recover the launch command if a log or report records one verbatim."""
    patterns = ("*.log", "*_report_*.md")
    for pattern in patterns:
        for path in find_files(run_dir, pattern):
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as handle:
                    head = [next(handle, "") for _ in range(80)]
            except OSError:
                continue
            for line in head:
                stripped = line.strip().lstrip("$ ").strip()
                if stripped.startswith("julia ") and "--project" in stripped:
                    return stripped
                if stripped.startswith("bin/ppb") or stripped.startswith("./bin/ppb"):
                    return stripped
    return UNKNOWN


# --------------------------------------------------------------------------
# Manifest construction
# --------------------------------------------------------------------------
def build_manifest_ppb(run_dir: str, raw: str, agg: str | None, store: str,
                       machine_override: str | None, notes: str) -> tuple[dict, pd.DataFrame]:
    frame = pd.read_csv(raw, low_memory=False)
    hostname = UNKNOWN
    host_col = first_column(frame, ("machine", "hostname"))
    if host_col is not None:
        hosts = sorted(str(h) for h in frame[host_col].dropna().unique())
        if len(hosts) == 1:
            hostname = hosts[0]
        elif hosts:
            raise SystemExit(f"{raw} mixes hosts {hosts}; it is not one run")

    facts = HOST_FACTS.get(hostname, {})
    machine = machine_override or facts.get("machine", UNKNOWN)

    cpu_threads = distinct(frame, "cpu_threads")
    cpu_threads_value = cpu_threads[0] if cpu_threads and len(cpu_threads) == 1 else (
        cpu_threads if cpu_threads else UNKNOWN
    )

    commit, commits = commit_fields(frame)

    if "timestamp_utc" in frame.columns and not frame["timestamp_utc"].dropna().empty:
        stamps = sorted(str(t) for t in frame["timestamp_utc"].dropna())
        launched, finished = stamps[0], stamps[-1]
        ts_source = "raw CSV timestamp_utc column"
    else:
        oldest, newest = newest_and_oldest([raw])
        launched, finished = utc_stamp(oldest), utc_stamp(newest)
        ts_source = "file modification time (the CSV has no timestamp column)"

    manifest = {
        "run_id": "",  # filled by the caller
        "machine": machine,
        "hostname": hostname,
        "cpu_model": facts.get("cpu_model", UNKNOWN),
        "physical_cores": facts.get("physical_cores", UNKNOWN),
        "cpu_threads": cpu_threads_value,
        "memory_gb": facts.get("memory_gb", UNKNOWN),
        "harness": "ppb",
        "simulator_commit": commit,
        "simulator_commits": commits,
        "simulator_branch": UNKNOWN,
        "phases": distinct_or_unknown(frame, "phase_id"),
        "cases": distinct_or_unknown(frame, "case"),
        "modes": distinct_or_unknown(frame, "mode"),
        "repeats": repeat_counts_ppb(frame),
        "thread_counts": distinct_or_unknown(frame, "thread_count"),
        "process_workers": distinct_or_unknown(frame, "process_workers"),
        "mc_samples": distinct_or_unknown(frame, "mc_samples"),
        "calibration_store": store,
        "launched_utc": launched,
        "finished_utc": finished,
        "source_path": recorded_path(run_dir),
        "command": command_from_logs(run_dir),
        "notes": notes or "",
        "n_rows": int(len(frame)),
        "raw_csv": os.path.basename(raw),
        "aggregate_csv": os.path.basename(agg) if agg else UNKNOWN,
        "timestamps_source": ts_source,
        "host_facts_source": facts.get("source", UNKNOWN),
    }
    return manifest, frame


def build_manifest_ps(run_dir: str, csv_paths: list[str], store: str,
                      machine_override: str | None, notes: str) -> tuple[dict, pd.DataFrame]:
    frames = [pd.read_csv(p, low_memory=False) for p in csv_paths]
    frame = pd.concat(frames, ignore_index=True, sort=False)

    hostname = UNKNOWN
    hosts = sorted(str(h) for h in frame.get("hostname", pd.Series(dtype=str)).dropna().unique())
    if len(hosts) == 1:
        hostname = hosts[0]
    elif hosts:
        raise SystemExit(f"{run_dir} mixes hosts {hosts}; it is not one run")

    facts = HOST_FACTS.get(hostname, {})
    machine = machine_override or facts.get("machine", UNKNOWN)

    cpu_threads = distinct(frame, "cpu_threads")
    cpu_threads_value = cpu_threads[0] if cpu_threads and len(cpu_threads) == 1 else (
        cpu_threads if cpu_threads else UNKNOWN
    )

    modes: list[str] = []
    for column in ("mode", "backend", "profile"):
        values = distinct(frame, column)
        if values:
            modes.extend(str(v) for v in values)
    modes = sorted(set(modes))

    workers_column = first_column(frame, ("proc_workers", "workers", "outer_workers"))
    cases_column = first_column(frame, ("case", "workload"))

    oldest, newest = newest_and_oldest(csv_paths)

    manifest = {
        "run_id": "",
        "machine": machine,
        "hostname": hostname,
        "cpu_model": facts.get("cpu_model", UNKNOWN),
        "physical_cores": facts.get("physical_cores", UNKNOWN),
        "cpu_threads": cpu_threads_value,
        "memory_gb": facts.get("memory_gb", UNKNOWN),
        "harness": "ps",
        "simulator_commit": UNKNOWN,
        "simulator_commits": [],
        "simulator_branch": UNKNOWN,
        "phases": distinct_or_unknown(frame, "scenario"),
        "cases": distinct(frame, cases_column) if cases_column else [],
        "modes": modes,
        "repeats": repeat_counts_ps(frame),
        "thread_counts": distinct_or_unknown(frame, "julia_threads"),
        "process_workers": distinct(frame, workers_column) if workers_column else [],
        "mc_samples": distinct_or_unknown(frame, "samples"),
        "calibration_store": store,
        # The scenario harness writes no timestamp column, and the file
        # modification times record when results were copied between machines
        # rather than when the run started or finished, so neither is known.
        "launched_utc": UNKNOWN,
        "finished_utc": UNKNOWN,
        "source_path": recorded_path(run_dir),
        "command": command_from_logs(run_dir),
        "notes": notes or "",
        "n_rows": int(len(frame)),
        "raw_csv": ";".join(os.path.basename(p) for p in csv_paths),
        "aggregate_csv": UNKNOWN,
        "source_mtime_earliest_utc": utc_stamp(oldest),
        "source_mtime_latest_utc": utc_stamp(newest),
        "timestamps_source": (
            "none: the scenario CSVs carry no timestamp column, and their "
            "modification times record when the results were copied between "
            "machines; the run_id stamp is the earliest of those mtimes"
        ),
        "host_facts_source": facts.get("source", UNKNOWN),
    }
    return manifest, frame


# --------------------------------------------------------------------------
# Copying
# --------------------------------------------------------------------------
def copy_run(run_dir: str, dest: str, excludes: list[str]) -> tuple[int, int]:
    """Copy the run directory, routing *.log files into ``logs/``."""
    files = 0
    total = 0
    for root, dirnames, filenames in os.walk(run_dir):
        dirnames[:] = [d for d in sorted(dirnames) if d not in excludes]
        rel_dir = os.path.relpath(root, run_dir)
        for name in sorted(filenames):
            source = os.path.join(root, name)
            if name.endswith(".log"):
                rel = os.path.join("logs", "" if rel_dir == "." else rel_dir, name)
            else:
                rel = name if rel_dir == "." else os.path.join(rel_dir, name)
            target = os.path.join(dest, os.path.normpath(rel))
            os.makedirs(os.path.dirname(target), exist_ok=True)
            shutil.copy2(source, target)
            files += 1
            total += os.path.getsize(source)
    return files, total


# --------------------------------------------------------------------------
# index.csv and PROVENANCE.md
# --------------------------------------------------------------------------
def append_index(archive: str, manifest: dict) -> None:
    path = os.path.join(archive, "index.csv")
    rows = []
    if os.path.isfile(path):
        with open(path, newline="", encoding="utf-8") as handle:
            rows = [row for row in csv.DictReader(handle) if row.get("run_id") != manifest["run_id"]]
    phases = manifest["phases"]
    row = {
        "run_id": manifest["run_id"],
        "machine": manifest["machine"],
        "harness": manifest["harness"],
        "store": manifest["calibration_store"],
        "launched_utc": manifest["launched_utc"],
        "simulator_commit": manifest["simulator_commit"],
        "phases": ";".join(str(p) for p in phases) if isinstance(phases, list) else str(phases),
        "n_rows": manifest["n_rows"],
        "path": manifest["run_id"],
    }
    rows.append(row)
    rows.sort(key=lambda r: (r["machine"], r["harness"], str(r["launched_utc"]), r["run_id"]))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=INDEX_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


PROVENANCE_MARKER = "<!-- archived runs: entries below are appended by scripts/archive_paper_run.py -->"


def append_provenance(archive: str, manifest: dict) -> None:
    path = os.path.join(archive, "PROVENANCE.md")
    text = ""
    if os.path.isfile(path):
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
    if PROVENANCE_MARKER not in text:
        text = text.rstrip("\n") + "\n\n## Archived benchmark runs\n\n" + PROVENANCE_MARKER + "\n"
    heading = f"### `{manifest['run_id']}`"
    if heading in text:
        return
    phases = manifest["phases"]
    phase_text = ", ".join(str(p) for p in phases) if isinstance(phases, list) else str(phases)
    if UNKNOWN in (str(manifest["launched_utc"]), str(manifest["finished_utc"])):
        when = f"- Raw rows: {manifest['n_rows']}; measurement time not recorded ({manifest['timestamps_source']})"
    else:
        when = (
            f"- Raw rows: {manifest['n_rows']}; measured {manifest['launched_utc']} to "
            f"{manifest['finished_utc']} ({manifest['timestamps_source']})"
        )
    commit_line = f"- Simulator commit: `{manifest['simulator_commit']}`"
    if manifest["simulator_commit"] == "multiple":
        commit_line += " — " + ", ".join(f"`{c}`" for c in manifest.get("simulator_commits", []))
    entry = [
        "",
        heading,
        "",
        f"- Source: `{manifest['source_path']}`",
        f"- Host: `{manifest['hostname']}` ({manifest['cpu_model']}, "
        f"{manifest['physical_cores']} physical cores, {manifest['cpu_threads']} hardware threads, "
        f"{manifest['memory_gb']} GB)",
        f"- Harness: `{manifest['harness']}`; calibration store `{manifest['calibration_store']}`",
        commit_line,
        f"- Phases: {phase_text}",
        when,
    ]
    if manifest.get("notes"):
        entry.append(f"- {manifest['notes']}")
    entry.append("")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text.rstrip("\n") + "\n" + "\n".join(entry))


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def archive(args) -> int:
    run_dir = os.path.abspath(args.run_dir)
    if not os.path.isdir(run_dir):
        raise SystemExit(f"{run_dir} is not a directory")
    archive_dir = os.path.abspath(args.archive)

    harness = args.harness or detect_harness(run_dir)
    if harness == "ppb" and not args.store:
        raise SystemExit(
            "--store is required for ppb runs: the calibration store's state is not in the CSV"
        )
    store = args.store or ("na" if harness == "ps" else None)
    if store not in STORE_STATES:
        raise SystemExit(f"--store must be one of {sorted(STORE_STATES)}")

    if harness == "ppb":
        raw, agg = ppb_csvs(run_dir, args.raw)
        manifest, frame = build_manifest_ppb(run_dir, raw, agg, store, args.machine, args.notes)
        data_paths = [raw]
    else:
        csv_paths = ps_csvs(run_dir)
        manifest, frame = build_manifest_ps(run_dir, csv_paths, store, args.machine, args.notes)
        raw, agg = csv_paths[0], None
        data_paths = csv_paths

    if manifest["machine"] not in MACHINE_NAMES:
        raise SystemExit(
            f"machine for host {manifest['hostname']!r} is unknown; pass --machine "
            f"(one of {sorted(MACHINE_NAMES)})"
        )

    stamp = args.stamp or stamp_from(run_dir, frame, data_paths)
    run_id = f"{manifest['machine']}_{harness}_{store}_{stamp}"
    manifest["run_id"] = run_id
    manifest["archived_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")

    dest = os.path.join(archive_dir, run_id)
    if os.path.exists(dest) and not args.force:
        raise SystemExit(f"{dest} already exists; pass --force to replace it")

    if args.dry_run:
        print(f"[dry run] would archive {run_dir} as {run_id}")
        for key in MANIFEST_ORDER:
            if key in manifest:
                print(f"  {key} = {manifest[key]}")
        return 0

    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.makedirs(dest, exist_ok=True)
    excludes = list(args.exclude or [])
    if harness == "ps":
        excludes.append("plots")
    files, total = copy_run(run_dir, dest, excludes)

    # WS9b's figure scripts glob `paper_benchmarks_raw_*.csv`.  A run whose raw
    # CSV was renamed by hand gets a copy under the canonical name so exactly
    # one file matches that glob; the original is kept beside it.
    canonical = f"paper_benchmarks_raw_{stamp}.csv"
    if harness == "ppb" and not glob.glob(os.path.join(dest, "**", PPB_RAW_GLOB), recursive=True):
        shutil.copy2(raw, os.path.join(dest, canonical))
        manifest["raw_csv"] = canonical
        manifest["notes"] = (
            (manifest["notes"] + " " if manifest["notes"] else "")
            + f"The raw CSV was copied to {canonical} from {os.path.basename(raw)}; "
            "the original name is kept in the run directory."
        ).strip()
        files += 1
        total += os.path.getsize(raw)
        if agg and not glob.glob(os.path.join(dest, "**", PPB_AGG_GLOB), recursive=True):
            canonical_agg = f"paper_benchmarks_aggregated_{stamp}.csv"
            shutil.copy2(agg, os.path.join(dest, canonical_agg))
            manifest["aggregate_csv"] = canonical_agg
            files += 1
            total += os.path.getsize(agg)

    write_toml(os.path.join(dest, "manifest.toml"), manifest)
    append_index(archive_dir, manifest)
    append_provenance(archive_dir, manifest)

    print(f"{run_id}\t{files} files\t{total / 1048576:.1f} MB\t{manifest['n_rows']} raw rows")
    return 0


def verify(args) -> int:
    archive_dir = os.path.abspath(args.archive)
    index_path = os.path.join(archive_dir, "index.csv")
    if not os.path.isfile(index_path):
        raise SystemExit(f"no index.csv in {archive_dir}")
    with open(index_path, newline="", encoding="utf-8") as handle:
        index = {row["run_id"]: row for row in csv.DictReader(handle)}

    run_dirs = sorted(
        name
        for name in os.listdir(archive_dir)
        if os.path.isdir(os.path.join(archive_dir, name))
    )
    problems = []
    for name in run_dirs:
        manifest_path = os.path.join(archive_dir, name, "manifest.toml")
        if not os.path.isfile(manifest_path):
            problems.append(f"{name}: no manifest.toml")
            continue
        manifest = read_toml(manifest_path)
        if manifest.get("run_id") != name:
            problems.append(f"{name}: manifest run_id is {manifest.get('run_id')!r}")
        if name not in index:
            problems.append(f"{name}: missing from index.csv")

        run_dir = os.path.join(archive_dir, name)
        if manifest.get("harness") == "ppb":
            hits = sorted(glob.glob(os.path.join(run_dir, PPB_RAW_GLOB)))
            if len(hits) != 1:
                problems.append(f"{name}: {len(hits)} files match {PPB_RAW_GLOB} at the run root")
                continue
            frame = pd.read_csv(hits[0], low_memory=False)
        else:
            hits = [
                os.path.join(run_dir, f)
                for f in sorted(os.listdir(run_dir))
                if PS_CSV_RE.match(f)
            ]
            if not hits:
                problems.append(f"{name}: no scenario CSVs")
                continue
            frame = pd.concat(
                [pd.read_csv(p, low_memory=False) for p in hits], ignore_index=True, sort=False
            )

        if int(manifest.get("n_rows", -1)) != len(frame):
            problems.append(
                f"{name}: manifest n_rows={manifest.get('n_rows')} but the CSV has {len(frame)}"
            )
        commit, _ = commit_fields(frame)
        if str(manifest.get("simulator_commit")) != commit:
            problems.append(
                f"{name}: manifest simulator_commit={manifest.get('simulator_commit')!r} "
                f"but the CSV says {commit!r}"
            )
        row = index.get(name)
        if row and int(row["n_rows"]) != len(frame):
            problems.append(f"{name}: index.csv n_rows={row['n_rows']} but the CSV has {len(frame)}")
        if row and row["simulator_commit"] != commit:
            problems.append(
                f"{name}: index.csv simulator_commit={row['simulator_commit']!r} "
                f"but the CSV says {commit!r}"
            )
        print(f"ok  {name}  {len(frame)} rows  {commit}")

    for run_id in index:
        if run_id not in run_dirs:
            problems.append(f"{run_id}: in index.csv but no directory")

    if problems:
        print("\nFAILED", file=sys.stderr)
        for problem in problems:
            print("  " + problem, file=sys.stderr)
        return 1
    print(f"\n{len(run_dirs)} runs verified")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", nargs="?", help="benchmark run directory to archive")
    parser.add_argument("--archive", required=True, help="the paper repository's data/raw directory")
    parser.add_argument("--machine", choices=sorted(MACHINE_NAMES),
                        help="override the machine name derived from the hostname")
    parser.add_argument("--store", choices=sorted(STORE_STATES),
                        help="calibration store state; required for ppb runs")
    parser.add_argument("--harness", choices=("ppb", "ps"), help="override harness detection")
    parser.add_argument("--raw", help="the raw CSV, when a run directory holds more than one")
    parser.add_argument("--stamp", help="override the YYYYMMDD_HHMMSS stamp in the run_id")
    parser.add_argument("--notes", default="", help="free text copied into the manifest and PROVENANCE")
    parser.add_argument("--exclude", action="append", default=[],
                        help="directory name to skip while copying (repeatable)")
    parser.add_argument("--force", action="store_true", help="replace an existing run_id")
    parser.add_argument("--dry-run", action="store_true", help="derive and print, copy nothing")
    parser.add_argument("--verify", action="store_true",
                        help="re-read every archived run and check it against its manifest")
    args = parser.parse_args()

    if args.verify:
        return verify(args)
    if not args.run_dir:
        parser.error("a run directory is required unless --verify is given")
    return archive(args)


if __name__ == "__main__":
    raise SystemExit(main())
