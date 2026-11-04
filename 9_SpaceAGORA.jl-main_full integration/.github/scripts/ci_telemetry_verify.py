#!/usr/bin/env python3
"""Check that the telemetry quick-regression shards graded every scenario once.

Usage: ci_telemetry_verify.py SUMMARIES_DIR

Reads every *_summary.csv under SUMMARIES_DIR (one per shard) and compares the
scenarios they report with the scenarios test/telemetry_benchmark_manifest.toml
defines, which is what an unsharded `quick` run grades. A scenario added to the
manifest but to no shard in julia-ci.yml fails here instead of going ungraded.
"""
import csv
import sys
import tomllib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    manifest = tomllib.loads((ROOT / "test/telemetry_benchmark_manifest.toml").read_text())
    expected = {s["name"].lower() for s in manifest["scenarios"]}
    seen = Counter()
    for path in sorted(Path(sys.argv[1]).rglob("*_summary.csv")):
        with path.open(newline="") as fh:
            names = {row["scenario"].lower() for row in csv.DictReader(fh)}
        print(f"{path}: {', '.join(sorted(names))}")
        seen.update(names)
    problems = [f"{s}: graded {seen[s]} times, expected once" for s in sorted(expected) if seen[s] != 1]
    problems += [f"{s}: graded but not in the manifest" for s in sorted(set(seen) - expected)]
    if problems:
        print("FAILED:\n  " + "\n  ".join(problems))
        sys.exit(1)
    print(f"all {len(expected)} telemetry scenarios graded exactly once")


if __name__ == "__main__":
    main()
