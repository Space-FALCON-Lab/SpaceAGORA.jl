#!/usr/bin/env python3
"""Check that a set of CI shards ran everything, once.

Usage: ci_shard_verify.py tests|examples REPORTS_DIR  (matrix JSON in $SHARD_MATRIX)

REPORTS_DIR holds one directory per shard (named by its `shard` key), each with
the TOML reports that .github/scripts/ci_shard_hooks.jl wrote, the `*.done`
markers of the workflow-step items (grid, depwarn) and a `job.ok` marker the
shard writes as its last step. Fails when a planned shard left no reports or did
not finish, when shards disagree about the universe of items, or when any item
of that universe was not attempted exactly once. The universe comes from the
Julia side's parse of the sources, not from the planner, so an item the planner
missed is caught here.

Prints each item's measured minutes, for re-balancing .github/ci-shards.toml
(ci_shard_weights.py), and a warning annotation for each test item that took
far longer than its planned weight and each shard whose start-up before its
first item took long (start-up runs in every shard).
"""
import json
import os
import re
import sys
import tomllib
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ci_shard_plan  # noqa: E402

PREAMBLE_WARN_MIN = 6.0


def load(path):
    return tomllib.loads(path.read_text()) if path.is_file() else None


def slow_timings(matrix, timings):
    """(row, minutes, planned) for each test item whose measured minutes exceed
    both twice and five minutes more than the weight the planner used for it,
    and for each shard preamble over PREAMBLE_WARN_MIN (planned None)."""
    modes = {e["shard"].rsplit("-", 1)[0] for e in matrix}
    cfg = tomllib.loads(ci_shard_plan.CONFIG.read_text())["tests"]
    mode = modes.pop() if len(modes) == 1 else None
    weights = cfg[mode].get("weights", {}) if mode in cfg else None
    defaults = cfg[mode]["default_weight"] if mode in cfg else None
    slow = []
    for row, secs in sorted(timings.items()):
        if re.fullmatch(r"[\w-]+-\d+:preamble", row):
            if secs / 60 > PREAMBLE_WARN_MIN:
                slow.append((row, secs / 60, None))
        elif weights is not None and ci_shard_plan.kind(row) in ("suite", "probe", "unit"):
            planned = float(weights.get(row, defaults.get(ci_shard_plan.kind(row), defaults.get("other", 1.0))))
            if secs / 60 > max(2 * planned, planned + 5):
                slow.append((row, secs / 60, planned))
    return slow


def main():
    kind, reports = sys.argv[1], Path(sys.argv[2])
    matrix = json.loads(os.environ["SHARD_MATRIX"])["include"]
    problems = []
    universes = {}
    attempted = Counter()
    timings = {}
    planned_workflow_items = Counter()
    done_workflow_items = Counter()

    for entry in matrix:
        shard = entry["shard"]
        sdir = reports / shard
        if not (sdir / "job.ok").is_file():
            problems.append(f"{shard}: did not finish (no job.ok marker)")
        names = ["examples"] if kind == "examples" else ["integration", "unit"]
        items = entry["items"].split()
        wants_unit = any(i.startswith("unit:") for i in items)
        for name in names:
            rep = load(sdir / f"{name}.toml")
            if rep is None:
                if name == "unit" and not wants_unit:
                    continue
                problems.append(f"{shard}: no {name}.toml report")
                continue
            universes.setdefault(name, {})[shard] = tuple(rep["universe"])
            requested = sorted(i for i in rep["requested"]
                               if (name != "unit" or i.startswith("unit:")))
            expected = sorted(i for i in items if name != "unit" or i.startswith("unit:"))
            if requested != expected:
                problems.append(f"{shard}: {name} report was asked for {requested}, planned {expected}")
            attempted.update(rep["attempted"])
            for failed in rep.get("failed", []):
                problems.append(f"{shard}: {failed} failed")
            for item, secs in rep.get("timings", {}).items():
                timings[f"{shard}:{item}" if item in ("preamble", "suite:09") else item] = secs
        for item in items:
            if ":" not in item and kind == "tests":
                planned_workflow_items[item] += 1
                if (sdir / f"{item}.done").is_file():
                    done_workflow_items[item] += 1

    universe = set()
    for name, per_shard in universes.items():
        distinct = set(per_shard.values())
        if len(distinct) > 1:
            problems.append(f"shards disagree about the {name} universe")
        for u in distinct:
            universe.update(u)
    if not universe:
        problems.append("no shard reported a universe")

    for item in sorted(universe):
        n = attempted[item]
        if n != 1:
            problems.append(f"{item}: attempted {n} times, expected once")
    for item in sorted(set(attempted) - universe):
        problems.append(f"{item}: attempted but not in the universe")
    for item, n in planned_workflow_items.items():
        if n != 1 or done_workflow_items[item] != 1:
            problems.append(f"workflow item {item}: planned {n} times, completed {done_workflow_items[item]} times")

    print(f"{kind}: {len(matrix)} shards, {len(universe)} items in the universe")
    for item, secs in sorted(timings.items(), key=lambda kv: -kv[1]):
        print(f"  {secs / 60:6.2f} min  {item}")
    if kind == "tests":
        slow = slow_timings(matrix, timings)
        starts = sorted((m, row.split(":")[0]) for row, m, planned in slow if planned is None)
        if starts:
            print(f"::warning title=Slow CI shard start-up::{len(starts)} of {len(matrix)} shards took over "
                  f"{PREAMBLE_WARN_MIN:g} min before their first item (longest {starts[-1][0]:.1f} min, "
                  f"{starts[-1][1]}); tests run there repeat in every shard")
        for row, minutes, planned in slow:
            if planned is not None:
                print(f"::warning title=Stale CI shard weight::{row} took {minutes:.1f} min, planned "
                      f"{planned:.1f} min; refresh .github/ci-shards.toml with .github/scripts/ci_shard_weights.py")
    if problems:
        print("\nFAILED:")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    print("every item ran exactly once")


if __name__ == "__main__":
    main()
