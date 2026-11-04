#!/usr/bin/env python3
"""Compare two P-series runs point by point: did a code change move the numbers?

Joins two aggregated CSVs on (phase, case, mode, threads, workers, mc_samples)
and reports the ratio new/old at every point both measured. The question this
answers is not "is the new run fast" but "is it the same measurement", so the
output leads with how many points moved further than run-to-run noise.

A mode one run measured and the other did not (a run taken before `predictive`
joined the ladder, say) is not an error: it simply has no counterpart to be
compared against, and is reported as such rather than dropped silently.

The noise band defaults to 8%, which is the worst static-route drift observed
between two runs of identical code on this harness (median 1.6%); pass
--band to set your own.

Usage:
    python3 scripts/compare_paper_routing_runs.py OLD NEW [--band 0.08]
        [--phases P1,P5] [--label-old NAME --label-new NAME]
"""
from __future__ import annotations

import argparse
import statistics
import sys

from make_paper_routing_tables import load

KEY = ["phase_id", "case", "mode", "thread_count", "process_workers", "mc_samples"]


def index(df):
    out = {}
    for _, r in df.iterrows():
        t = r.get("wall_time_median_s")
        if t is None or t != t or float(t) <= 0:
            continue
        out[tuple(str(r[k]) for k in KEY)] = float(t)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--band", type=float, default=0.08,
                    help="fractional move treated as run-to-run noise (default 0.08)")
    ap.add_argument("--phases", default=None, help="comma-separated subset, e.g. P1,P5")
    ap.add_argument("--label-old", default="old")
    ap.add_argument("--label-new", default="new")
    args = ap.parse_args()

    a, b = index(load(args.old)), index(load(args.new))
    want = set(args.phases.split(",")) if args.phases else None

    shared = [k for k in a if k in b and (want is None or k[0] in want)]
    only_a = [k for k in a if k not in b and (want is None or k[0] in want)]
    only_b = [k for k in b if k not in a and (want is None or k[0] in want)]
    if not shared:
        raise SystemExit("no points measured by both runs")

    rows = sorted(((b[k] / a[k], k) for k in shared), key=lambda x: -abs(x[0] - 1))
    moved = [(r, k) for r, k in rows if abs(r - 1) > args.band]

    print(f"{args.label_old}: {args.old}")
    print(f"{args.label_new}: {args.new}")
    print(f"\n{len(shared)} points measured by both"
          + (f"; {len(only_a)} only in {args.label_old}, {len(only_b)} only in {args.label_new}"
             if only_a or only_b else ""))
    # Name the modes that exist on one side only. Without this, adding a mode to
    # the ladder shows up as nothing more than a larger "only in new" count.
    modes_a = {k[2] for k in a}
    modes_b = {k[2] for k in b}
    for side, extra in ((args.label_old, modes_a - modes_b), (args.label_new, modes_b - modes_a)):
        if extra:
            print(f"  modes only in {side}: {', '.join(sorted(extra))} "
                  "(not compared: no counterpart run)")
    ratios = [r for r, _ in rows]
    print(f"ratio {args.label_new}/{args.label_old}: median {statistics.median(ratios):.3f}, "
          f"range {min(ratios):.3f}-{max(ratios):.3f}")
    print(f"outside the +/-{args.band:.0%} noise band: {len(moved)} of {len(shared)}")

    if moved:
        print(f"\n{'phase':<6}{'case':<34}{'mode':<20}{'thr':>4}{'wk':>4}"
              f"{args.label_old:>10}{args.label_new:>10}{'ratio':>8}")
        for r, k in moved:
            ph, case, mode, thr, wk, _mc = k
            print(f"{ph:<6}{case[:33]:<34}{mode:<20}{thr:>4}{wk:>4}"
                  f"{a[k]:>10.3f}{b[k]:>10.3f}{r:>8.3f}")
    else:
        print("\nNo point moved further than the noise band: the two runs are "
              "the same measurement.")

    print("\nper phase:")
    for ph in sorted({k[0] for k in shared}):
        rs = [b[k] / a[k] for k in shared if k[0] == ph]
        out = sum(1 for x in rs if abs(x - 1) > args.band)
        print(f"  {ph}: n={len(rs):>3}  median {statistics.median(rs):.3f}  "
              f"range {min(rs):.3f}-{max(rs):.3f}  outside band: {out}")

    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    raise SystemExit(main())
