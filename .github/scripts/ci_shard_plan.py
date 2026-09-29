#!/usr/bin/env python3
"""Plan the Julia CI shards.

Finds every test item (suites, suite 09's probes, the unit files) and every
example the smoke runs, and splits each set across a fixed number of shards,
balanced by the measured minutes in .github/ci-shards.toml. Items without a
measurement get their kind's default, so a new probe, unit file or example is
assigned to some shard without anyone editing the plan.

Writes the job matrices to $GITHUB_OUTPUT (or prints them) as
{"include": [...]}, one entry per shard and no other matrix key, so every
entry is a job of its own.

The Julia side (.github/scripts/ci_shard_hooks.jl) rediscovers the items from
the parsed source and refuses any item it does not know, and
ci_shard_verify.py checks the shards' reports against that universe, so a
disagreement between this planner's text scan and Julia's parse fails CI
instead of dropping a test.
"""
import json
import os
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / ".github" / "ci-shards.toml"


def suite_files():
    text = (ROOT / "test/integration/runtests.jl").read_text()
    block = re.search(r"const _ALL_SUITES = \[(.*?)\]", text, re.S)
    if not block:
        sys.exit("ci_shard_plan: _ALL_SUITES not found in test/integration/runtests.jl")
    suites = re.findall(r'"(\d\d_[A-Za-z0-9_]+\.jl)"', block.group(1))
    on_disk = sorted(p.name for p in (ROOT / "test/suites").glob("[0-9][0-9]_*.jl"))
    if sorted(suites) != on_disk:
        sys.exit(f"ci_shard_plan: _ALL_SUITES {suites} does not match test/suites {on_disk}")
    return suites


def probe_files():
    text = (ROOT / "test/suites/09_probe_drivers.jl").read_text()
    block = re.findall(r"probe_files\s*=\s*\[(.*?)\]", text, re.S)
    if len(block) != 1:
        sys.exit("ci_shard_plan: expected one probe_files list in 09_probe_drivers.jl")
    lines = [l.split("#", 1)[0] for l in block[0].splitlines()]
    return re.findall(r'"([^"]+)"', "\n".join(lines))


def unit_files():
    labels = []
    for line in (ROOT / "test/unit/runtests.jl").read_text().splitlines():
        m = re.match(r"^include\((.*)\)\s*(#.*)?$", line.strip())
        if m:
            labels.append("/".join(re.findall(r'"([^"]*)"', m.group(1))))
    return labels


def example_files():
    # Mirrors list_examples() in test/smoke/ci_examples_suite_smoke.jl.
    skip = {"common.jl", "aerobraking_mission_plot_utils.jl", "odyssey_surrogate.jl"}
    return sorted(p.name for p in (ROOT / "examples").glob("*.jl") if p.name not in skip)


def kind(item):
    return item.split(":", 1)[0] if ":" in item else item


def partition(items, shards, weights, defaults, unit_overhead=0.0):
    """Longest-processing-time assignment. The unit files a shard gets run in
    one extra Julia process, so the first one also pays that process's
    start-up."""
    def w(item):
        return float(weights.get(item, defaults.get(kind(item), defaults.get("other", 1.0))))

    loads = [0.0] * shards
    has_unit = [False] * shards
    groups = [[] for _ in range(shards)]
    for item in sorted(items, key=lambda i: (-w(i), i)):
        def cost(s):
            extra = unit_overhead if kind(item) == "unit" and not has_unit[s] else 0.0
            return loads[s] + w(item) + extra
        best = min(range(shards), key=lambda s: (cost(s), s))
        loads[best] = cost(best)
        has_unit[best] = has_unit[best] or kind(item) == "unit"
        groups[best].append(item)
    order = {item: n for n, item in enumerate(items)}
    return [(sorted(g, key=order.get), round(l, 1)) for g, l in zip(groups, loads)]


def main():
    cfg = tomllib.loads(CONFIG.read_text())
    suites = suite_files()
    tests = ([f"suite:{s[:2]}" for s in suites if s != "09_probe_drivers.jl"]
             + [f"probe:{p}" for p in probe_files()]
             + [f"unit:{u}" for u in unit_files()])
    if len(set(tests)) != len(tests):
        sys.exit("ci_shard_plan: duplicate test items")

    outputs = {}
    for mode in ("plain", "coverage"):
        mcfg = cfg["tests"][mode]
        items = tests + list(mcfg.get("workflow_items", []))
        shards = partition(items, int(mcfg["shards"]), mcfg.get("weights", {}),
                           mcfg["default_weight"], float(mcfg.get("unit_process_overhead", 0.0)))
        entries = []
        for n, (group, load) in enumerate(shards, start=1):
            if not group:
                sys.exit(f"ci_shard_plan: {mode} shard {n} is empty; lower the shard count")
            entries.append({
                "shard": f"{mode}-{n}",
                "name": f"{mode} {n}/{len(shards)}",
                "items": " ".join(group),
                "grid": "grid" in group,
                "depwarn": "depwarn" in group,
                "estimate_min": load,
            })
            print(f"{mode} {n}/{len(shards)} ~{load} min: {' '.join(group)}", file=sys.stderr)
        outputs[f"{mode}_matrix"] = {"include": entries}

    ecfg = cfg["examples"]
    shards = partition(example_files(), int(ecfg["shards"]), ecfg.get("weights", {}), ecfg["default_weight"])
    entries = []
    for n, (group, load) in enumerate(shards, start=1):
        entries.append({"shard": f"examples-{n}", "name": f"{n}/{len(shards)}",
                        "items": " ".join(group), "estimate_min": load})
        print(f"examples {n}/{len(shards)} ~{load} min: {' '.join(group)}", file=sys.stderr)
    outputs["examples_matrix"] = {"include": entries}

    out = os.environ.get("GITHUB_OUTPUT")
    lines = [f"{k}={json.dumps(v, separators=(',', ':'))}" for k, v in outputs.items()]
    if out:
        with open(out, "a") as fh:
            fh.write("\n".join(lines) + "\n")
    else:
        print("\n".join(lines))


if __name__ == "__main__":
    main()
