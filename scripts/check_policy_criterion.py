#!/usr/bin/env python3
"""Check the R6 (policy_v2) adaptive-routing acceptance criterion against raw
per-repeat benchmark data.

Criterion: at every measured point, the median wall time of R6 over ALL
repeats must be no more than `--tolerance` (default 10%) slower than the
bias-corrected best pinned static route at that point, AND no R6 campaign at
any point may have failed.

This reads the RAW per-repeat CSV (`paper_benchmarks_raw_<stamp>.csv`), not
the aggregated one that `scripts/make_paper_routing_tables.py` reads: the
harness's aggregation step drops rows with `success=false` before computing
medians, which hides a campaign that ran 3500x slower than baseline because it
failed outright rather than merely running slow. Working from raw repeats is
the only way to see that failure in the median and in the failed-campaign
count.

Usage:
    python3 scripts/check_policy_criterion.py RUN_DIR [RUN_DIR ...] \
        [--label NAME ...] [--tolerance 0.10] [--equiv-band 0.03] [--markdown]

Each RUN_DIR is searched recursively for `paper_benchmarks_raw_*.csv`; every
match found under it is loaded. --label assigns a display name to each
RUN_DIR by position (repeatable); the default label is the directory
basename. Exit code is 1 if any point in any run fails the criterion.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import pandas as pd

STATIC_PARALLEL = ["outer_threads", "outer_process", "inner_only", "outer_inner_static"]
ADAPTIVE = "policy_v2"

POINT_KEYS = ["phase_id", "case", "thread_count", "process_workers", "mc_samples"]

# Phases whose axis is a single already-present column; P1's axis (spacecraft
# count) is parsed out of the case name instead, and P5's ("workers x
# threads") is built from two columns. A phase not in this table still gets a
# usable axis label from the generic fallback below -- it is not excluded.
PHASE_AXIS_COL = {
    "P2": ("thread_count", "threads"),
    "P3": ("process_workers", "process_workers"),
    "P4": ("process_workers", "process_workers"),
}


def _n_sat_from_case(case: str):
    m = re.search(r"(\d+)sat", str(case))
    return int(m.group(1)) if m else None


def axis_for(row) -> tuple[str, object]:
    """Return (label, sort_key) for one point's row of key columns."""
    phase = row["phase_id"]
    if phase == "P1":
        n_sat = _n_sat_from_case(row["case"])
        if n_sat is not None:
            return str(n_sat), n_sat
    if phase == "P5":
        w, t = int(row["process_workers"]), int(row["thread_count"])
        return f"{w}x{t}", w
    if phase in PHASE_AXIS_COL:
        col, _ = PHASE_AXIS_COL[phase]
        val = row[col]
        return str(int(val)), int(val)
    # Generic fallback for any phase outside P1-P5: show the full config so
    # nothing is silently dropped, ordered by (threads, workers, mc_samples).
    tc, pw, mc = int(row["thread_count"]), int(row["process_workers"]), int(row["mc_samples"])
    return f"tc={tc},pw={pw},mc={mc}", (tc, pw, mc)


def normalize_success(series: pd.Series) -> pd.Series:
    """Coerce a `success` column (bool, or string 'True'/'False' etc.) to bool."""

    def conv(v):
        if isinstance(v, bool):
            return v
        if pd.isna(v):
            return True
        if isinstance(v, (int, float)):
            return bool(v)
        return str(v).strip().lower() in ("true", "1", "yes")

    return series.map(conv)


def load_raw(run_dir: str) -> pd.DataFrame:
    hits = sorted(glob.glob(os.path.join(run_dir, "**", "paper_benchmarks_raw_*.csv"), recursive=True))
    if not hits:
        raise SystemExit(f"no paper_benchmarks_raw_*.csv found under {run_dir}")
    frames = []
    for path in hits:
        df = pd.read_csv(path)
        df["_source"] = path
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def mode_stats(grp: pd.DataFrame) -> dict:
    """Per-mode {median, n, n_failed} for one point's rows."""
    out = {}
    for mode, mgrp in grp.groupby("mode"):
        n = len(mgrp)
        if "success" in mgrp.columns:
            ok = normalize_success(mgrp["success"])
            n_failed = int((~ok).sum())
        else:
            n_failed = 0
        out[mode] = {
            "median": float(mgrp.wall_time_s.median()),
            "n": n,
            "n_failed": n_failed,
        }
    return out


def evaluate_point(key: tuple, grp: pd.DataFrame, tolerance: float, equiv_band: float):
    stats = mode_stats(grp)
    statics = {m: s for m, s in stats.items() if m in STATIC_PARALLEL}
    if not statics:
        return None  # no pinned static route measured at this point -- can't baseline it

    medians = {m: s["median"] for m, s in statics.items()}
    lo, hi = min(medians.values()), max(medians.values())
    if lo > 0 and (hi / lo - 1.0) <= equiv_band:
        baseline, kind = sum(medians.values()) / len(medians), "mean"
    else:
        baseline, kind = lo, "min"

    if ADAPTIVE not in stats:
        return None  # R6 did not run at this point

    r6 = stats[ADAPTIVE]
    ratio = r6["median"] / baseline if baseline else float("inf")

    row0 = grp.iloc[0]
    axis, order = axis_for(row0)

    if r6["n_failed"] > 0:
        verdict, reason = "FAIL", f"failed campaigns: {r6['n_failed']}/{r6['n']}"
    elif ratio > 1.0 + tolerance:
        verdict, reason = "FAIL", f"ratio {ratio:.3f}"
    else:
        verdict, reason = "PASS", ""

    return {
        "phase": row0["phase_id"],
        "case": row0["case"],
        "axis": axis,
        "_order": order,
        "r6_median": r6["median"],
        "baseline": baseline,
        "baseline_kind": kind,
        "ratio": ratio,
        "n_failed": r6["n_failed"],
        "n": r6["n"],
        "verdict": verdict,
        "reason": reason,
    }


def evaluate_run(df: pd.DataFrame, tolerance: float, equiv_band: float) -> list[dict]:
    results = []
    skipped_no_static = 0
    skipped_no_r6 = 0
    for key, grp in df.groupby(POINT_KEYS, dropna=False):
        rec = evaluate_point(key, grp, tolerance, equiv_band)
        if rec is None:
            statics_present = any(m in grp["mode"].unique() for m in STATIC_PARALLEL)
            if not statics_present:
                skipped_no_static += 1
            elif ADAPTIVE not in grp["mode"].unique():
                skipped_no_r6 += 1
            continue
        results.append(rec)
    results.sort(key=lambda r: (r["phase"], r["case"], r["_order"]))
    if skipped_no_static:
        print(f"  ({skipped_no_static} point(s) skipped: no pinned static route present)")
    if skipped_no_r6:
        print(f"  ({skipped_no_r6} point(s) skipped: no policy_v2 (R6) data present)")
    return results


def fmt(x, digits=3):
    return f"{x:.{digits}f}"


def print_table(label: str, rows: list[dict], markdown: bool) -> None:
    headers = [
        "phase", "case", "axis", "R6 median (s)", "baseline (s)", "baseline kind",
        "ratio", "n_failed/n", "verdict", "reason",
    ]
    cells = [headers]
    for r in rows:
        cells.append([
            r["phase"],
            r["case"],
            str(r["axis"]),
            fmt(r["r6_median"]),
            fmt(r["baseline"]),
            r["baseline_kind"],
            fmt(r["ratio"]),
            f"{r['n_failed']}/{r['n']}",
            r["verdict"],
            r["reason"] or "-",
        ])

    print(f"\n== {label} ==")
    if markdown:
        print("| " + " | ".join(headers) + " |")
        print("|" + "|".join(["---"] * len(headers)) + "|")
        for row in cells[1:]:
            print("| " + " | ".join(row) + " |")
    else:
        widths = [max(len(row[i]) for row in cells) for i in range(len(headers))]
        for row in cells:
            print("  ".join(c.ljust(w) for c, w in zip(row, widths)))

    n_fail_campaigns = sum(1 for r in rows if r["verdict"] == "FAIL" and r["reason"].startswith("failed campaigns"))
    n_fail_ratio = sum(1 for r in rows if r["verdict"] == "FAIL" and r["reason"].startswith("ratio"))
    n_failing = n_fail_campaigns + n_fail_ratio
    print(
        f"\n{label}: {len(rows)} points, {n_failing} failing "
        f"({n_fail_campaigns} by failed campaigns, {n_fail_ratio} by ratio)"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dirs", nargs="+", help="Run directories containing paper_benchmarks_raw_*.csv")
    parser.add_argument("--label", action="append", default=[], help="Display label for each run dir, by position (repeatable)")
    parser.add_argument("--tolerance", type=float, default=0.10, help="Max allowed R6/baseline overhead (default 0.10 = 10%%)")
    parser.add_argument("--equiv-band", type=float, default=0.03, help="Static routes within this relative band are treated as tied (default 0.03)")
    parser.add_argument("--markdown", action="store_true", help="Emit GitHub-flavored markdown tables instead of plain text")
    args = parser.parse_args()

    labels = list(args.label) + [os.path.basename(os.path.normpath(d)) for d in args.run_dirs[len(args.label):]]

    any_failing = False
    for run_dir, label in zip(args.run_dirs, labels):
        df = load_raw(run_dir)
        rows = evaluate_run(df, args.tolerance, args.equiv_band)
        print_table(label, rows, args.markdown)
        if any(r["verdict"] == "FAIL" for r in rows):
            any_failing = True

    return 1 if any_failing else 0


if __name__ == "__main__":
    sys.exit(main())
