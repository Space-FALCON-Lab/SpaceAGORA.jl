#!/usr/bin/env python3
"""Refresh the measured weights in .github/ci-shards.toml.

Usage: ci_shard_weights.py plain|coverage LOG [LOG ...]

Each LOG is the output of ci_shard_verify.py for that mode: the log of the
`tests` job (plain) or of the `coverage-quality-gate` job (coverage), for
instance from
  gh api repos/Space-FALCON-Lab/SpaceAGORA.jl/actions/jobs/<job id>/logs > LOG
The weight of every item measured in at least one LOG becomes its mean minutes
over the LOGs that measured it. Items the LOGs do not measure keep their
current weight, and weights of items the planner no longer finds are dropped,
so the table always matches the current items. The table is rewritten in
place; everything outside it is left as it is. The example shards record no
per-example times, so [examples.weights] is still kept by hand.
"""
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ci_shard_plan  # noqa: E402

# "2026-10-15T14:48:41.6786877Z    37.17 min  unit:visualization/articulation_tests.jl"
LINE = re.compile(r"^(?:\S+Z )?\s*(\d+(?:\.\d+)?) min  (\S+)\s*$")
SECTIONS = {"plain": "tests.plain.weights", "coverage": "tests.coverage.weights"}
KIND_ORDER = {"suite": 0, "probe": 1, "unit": 2}


def measured(logs):
    """Mean minutes per item. Per-shard rows (`coverage-3:preamble`,
    `plain-1:suite:09`) describe a shard, not an item, and are skipped."""
    samples = {}
    for log in logs:
        for line in Path(log).read_text(errors="replace").splitlines():
            m = LINE.match(line)
            if m and not re.match(r"^(plain|coverage)-\d+:", m.group(2)):
                samples.setdefault(m.group(2), []).append(float(m.group(1)))
    return {item: statistics.mean(v) for item, v in samples.items()}


def universe(mode):
    cfg = ci_shard_plan.tomllib.loads(ci_shard_plan.CONFIG.read_text())
    suites = ci_shard_plan.suite_files()
    return set([f"suite:{s[:2]}" for s in suites if s != "09_probe_drivers.jl"]
               + [f"probe:{p}" for p in ci_shard_plan.probe_files()]
               + [f"unit:{u}" for u in ci_shard_plan.unit_files()]
               + list(cfg["tests"][mode].get("workflow_items", [])))


def toml_key(item):
    return item if re.fullmatch(r"[A-Za-z0-9_-]+", item) else '"' + item + '"'


def refresh(text, mode, measurements, items):
    """Return `text` with the mode's weight table replaced, and a summary."""
    header = f"[{SECTIONS[mode]}]"
    lines = text.splitlines(keepends=True)
    try:
        start = next(i for i, l in enumerate(lines) if l.strip() == header)
    except StopIteration:
        sys.exit(f"ci_shard_weights: {header} not found in {ci_shard_plan.CONFIG}")
    end = next((i for i in range(start + 1, len(lines)) if lines[i].lstrip().startswith("[")), len(lines))
    while end > start + 1 and not lines[end - 1].strip():
        end -= 1
    current = ci_shard_plan.tomllib.loads("".join(lines[start:end]))
    for part in SECTIONS[mode].split("."):
        current = current[part]

    weights = {i: w for i, w in current.items() if i in items}
    weights.update({i: w for i, w in measurements.items() if i in items})
    order = sorted(weights, key=lambda i: (KIND_ORDER.get(ci_shard_plan.kind(i), 3), i))
    body = [header + "\n"] + [f"{toml_key(i)} = {weights[i]:.2f}\n" for i in order]
    summary = {
        "updated": sorted(i for i in measurements if i in items),
        "kept": sorted(i for i in weights if i not in measurements),
        "dropped": sorted(i for i in current if i not in items),
        "unknown": sorted(i for i in measurements if i not in items),
        "default": sorted(items - set(weights)),
    }
    return "".join(lines[:start] + body + lines[end:]), summary


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in SECTIONS:
        sys.exit(__doc__)
    mode, logs = sys.argv[1], sys.argv[2:]
    measurements = measured(logs)
    if not measurements:
        sys.exit("ci_shard_weights: no `<minutes> min  <item>` lines in the logs")
    text, summary = refresh(ci_shard_plan.CONFIG.read_text(), mode, measurements, universe(mode))
    ci_shard_plan.CONFIG.write_text(text)
    for key, items in summary.items():
        print(f"{key}: {len(items)}" + (f"  {' '.join(items)}" if key != "updated" and items else ""))


if __name__ == "__main__":
    main()
