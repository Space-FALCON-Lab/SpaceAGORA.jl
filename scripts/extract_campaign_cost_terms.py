#!/usr/bin/env python3
"""Extract the predictive planner's campaign cost terms from an archived
targeted run, and rank the measured plans under the planner's model.

Two constants come out of the per-mode CSVs and dispatch traces of a run made
with `SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1`:

`round_tail`, in units of one sample time. For every timed campaign that ran
the pure process route (`process@wW+l0`), the campaign wall in sample units
minus the round count the list schedule predicts, per unit of the final
round's expected excess:

    round_tail = (wall / mean_sample - R) / (H_k - 1),   R = ceil(n / W),
                                                          k = n - W * (R - 1)

with `H_k` the k-th harmonic number (`H_k - 1` is the expected excess of the
slowest of k exponential-tailed samples over their mean). The reported value is
the median over every such row of the run. The same rows also give a per-round
figure, `(wall - sum_of_samples / W) / R / mean_sample`, which is printed so
the two readings can be compared: see docs/architecture/predictive_routing_r7.md
("The final round").

`pool_startup_s`, in seconds. For the FIRST campaign of each mode's process
(the only one that meets a cold pool), each pool worker's first take plus its
first occupancy minus its first work, i.e. how long after the dispatch started
that worker's first result arrived beyond the work it did. Median per campaign,
then the median across the run's cold campaigns.

With `--rank`, every plan the run measured at each point is priced by the
planner's makespan model (a Python mirror of `predictive_plan_makespan` in
src/simulation/campaigns/predictive_planner.jl) under four variants -- the
model without the new terms, with the round tail, with the round tail and a
scaled heap term, and with the round tail and the local slots' slowdown as the
run's own traces measured it -- and ranked against the measured medians.

With `--write-constants PATH`, the `[campaign]` table of the machine constants
file at PATH is replaced (or appended) with the two values and their source.
Only do that with a run measured on the machine PATH belongs to.

Usage:
    python3 scripts/extract_campaign_cost_terms.py RUN_DIR [--rank]
        [--usl-alpha A --usl-beta B] [--write-constants PATH]
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import statistics as st
import sys

TRACE_ARRAY = r"{key}=\[([^\]]*)\]"


def _array(line: str, key: str) -> list[float] | None:
    m = re.search(TRACE_ARRAY.format(key=key), line)
    if not m or not m.group(1).strip():
        return None
    return [float(x) for x in m.group(1).split(",")]


def _rows(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _log_lines(run_dir: str, label: str) -> list[str]:
    path = os.path.join(run_dir, "logs", f"{label}.log")
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        return fh.read().splitlines()


def _plans_in_log(lines: list[str]) -> list[tuple[str, int, int]]:
    """(route, workers, local_slots) per campaign, in campaign order."""
    plans = []
    for line in lines:
        m = re.search(r"\[predictive\] chosen (\w+)@w(\d+)\+l(\d+)", line)
        if m:
            plans.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return plans


def _bandit_plans(lines: list[str]) -> list[tuple[str, int, int]]:
    plans = []
    pending = None
    for line in lines:
        m = re.search(r"\[dispatch-trace\] plan=\S+ route=(\w+) threads=(\d+)", line)
        if m:
            pending = [m.group(1), int(m.group(2)), 0]
            if m.group(1) != "process":
                plans.append(tuple(pending))
                pending = None
            continue
        m = re.search(r"ensure_process_workers=\S+ workers=(\d+) local_slots=(\d+)", line)
        if m and pending is not None:
            pending[1], pending[2] = int(m.group(1)), int(m.group(2))
            plans.append(tuple(pending))
            pending = None
    return plans


def round_excess(k: int) -> float:
    """H_k - 1."""
    return sum(1.0 / i for i in range(2, k + 1))


def tail_rows(run_dir: str):
    """Yield (label, round_tail, per_round) for every pure-process timed row."""
    for path in sorted(glob.glob(os.path.join(run_dir, "*.csv"))):
        label = os.path.basename(path)[:-4]
        rows = _rows(path)
        if not rows:
            continue
        mode = rows[0]["mode"]
        if mode == "predictive":
            plans = _plans_in_log(_log_lines(run_dir, label))
            if not plans or any(p[0] != "process" or p[2] != 0 for p in plans):
                continue
        elif mode != "outer_process":
            continue
        for r in rows:
            if r["success"] != "true" or r["outer_backend_actual"] != "process":
                continue
            n = int(r["mc_samples"])
            w = int(r["process_workers"])
            wall = float(r["wall_time_s"])
            mean = float(r["mean_sample_wall_time_s"])
            total = float(r["sample_wall_time_sum_s"])
            rounds = math.ceil(n / w)
            k_last = n - w * (rounds - 1)
            if k_last < 2:
                continue
            tail = (wall / mean - rounds) / round_excess(k_last)
            per_round = (wall - total / w) / rounds / mean
            yield label, tail, per_round


def startup_campaigns(run_dir: str):
    """Yield (label, median_ms) for the first traced campaign of each log."""
    for path in sorted(glob.glob(os.path.join(run_dir, "logs", "*.log"))):
        label = os.path.basename(path)[:-4]
        for line in open(path):
            if "worker first_take=" not in line:
                continue
            take = _array(line, "first_take")
            occ = _array(line, "first_occupancy")
            work = _array(line, "first_work")
            if take and occ and work:
                yield label, st.median(t + o - w for t, o, w in zip(take, occ, work))
            break


# ── The planner's makespan model (mirror of predictive_planner.jl) ────────────

def usl_heap(alpha: float, beta: float, k: int) -> float:
    if k <= 1:
        return 1.0
    denom = 1.0 + alpha * (k - 1) + beta * k * (k - 1)
    speedup = max(1.0, k / denom) if denom > 0 else 1.0
    return k / speedup


def plan_makespan(n: int, workers: int, local_slots: int, heap: float,
                  tail: float = 0.0, startup: float = 0.0) -> float:
    slow = [1.0] * workers + [heap] * local_slots
    free = [startup] * workers + [0.0] * local_slots
    taken = [0] * len(slow)
    for _ in range(n):
        j = min(range(len(free)), key=lambda i: (free[i], i))
        free[j] += slow[j]
        taken[j] += 1
    finish = 0.0
    for lo, hi in ((0, workers), (workers, workers + local_slots)):
        if lo >= hi:
            continue
        top = max(taken[lo:hi])
        if top == 0:
            continue
        k = taken[lo:hi].count(top)
        finish = max(finish, max(free[lo:hi]) + tail * round_excess(k))
    return finish


def observed_heap(run_dir: str) -> dict[tuple[int, int], float]:
    """(workers, local_slots) -> median over campaigns of the local slots' mean
    first work over the pool workers' mean first work, from mixed traces."""
    ratios: dict[tuple[int, int], list[float]] = {}
    for path in sorted(glob.glob(os.path.join(run_dir, "logs", "*.log"))):
        worker = None
        for line in open(path):
            if "worker first_take=" in line:
                worker = _array(line, "first_work")
            elif "local first_take=" in line and worker:
                local = _array(line, "first_work")
                if local:
                    key = (len(worker), len(local))
                    ratios.setdefault(key, []).append(st.mean(local) / st.mean(worker))
                worker = None
    return {k: st.median(v) for k, v in ratios.items()}


def rank(run_dir: str, alpha: float, beta: float, tail: float, heap_scale: float):
    measured_heap = observed_heap(run_dir)
    points: dict[str, dict[tuple, list[float]]] = {}
    for path in sorted(glob.glob(os.path.join(run_dir, "*.csv"))):
        label = os.path.basename(path)[:-4]
        rows = _rows(path)
        if not rows:
            continue
        mode = rows[0]["mode"]
        point = label[: -len(mode) - 1] if label.endswith("_" + mode) else label
        n = int(rows[0]["mc_samples"])
        w = int(rows[0]["process_workers"])
        lines = _log_lines(run_dir, label)
        if mode == "outer_process":
            plans = [("process", w, 0)] * len(rows)
        elif mode == "predictive":
            plans = _plans_in_log(lines)[-len(rows):]
        else:
            plans = _bandit_plans(lines)[-len(rows):]
        if len(plans) != len(rows):
            print(f"  skip {label}: {len(plans)} plans for {len(rows)} rows", file=sys.stderr)
            continue
        for plan, r in zip(plans, rows):
            if plan[0] != "process" or r["success"] != "true":
                continue
            key = (n, plan[1], plan[2])
            points.setdefault(point, {}).setdefault(key, []).append(float(r["wall_time_s"]))
    for point, plans in sorted(points.items()):
        print(f"\n{point}")
        print("  plan          measured_s  rank |  model v1  rank |  +tail  rank | +tail,heap x{:.3f}  rank |"
              " +tail,heap observed  rank".format(heap_scale))
        table = []
        for (n, w, l), walls in plans.items():
            h = usl_heap(alpha, beta, l) if l > 0 else 1.0
            h_c = max(1.0, heap_scale * h) if l > 0 else 1.0
            h_o = measured_heap.get((w, l), float("nan")) if l > 0 else 1.0
            table.append(((w, l), st.median(walls),
                          plan_makespan(n, w, l, h),
                          plan_makespan(n, w, l, h, tail),
                          plan_makespan(n, w, l, h_c, tail),
                          plan_makespan(n, w, l, h_o, tail) if h_o == h_o else float("nan"),
                          h_o))

        def ranks(col):
            order = sorted(range(len(table)), key=lambda i: round(table[i][col], 9))
            out = [0] * len(table)
            for pos, i in enumerate(order):
                out[i] = pos + 1
            return out

        r1, r2, r3, r4, r5 = ranks(1), ranks(2), ranks(3), ranks(4), ranks(5)
        for i, ((w, l), m, a, b, c, d, h_o) in enumerate(table):
            print(f"  process@w{w}+l{l:<3} {m:9.3f}  {r1[i]:4d} | {a:8.3f}  {r2[i]:4d} | {b:6.3f}  {r3[i]:4d} |"
                  f" {c:18.3f}  {r4[i]:4d} | {d:9.3f} (x{h_o:.2f})  {r5[i]:4d}")


def write_constants(path: str, tail: float, startup_s: float, source: str) -> None:
    with open(path) as fh:
        text = fh.read()
    table = ("[campaign]\n"
             f"round_tail = {tail:.6g}\n"
             f"pool_startup_s = {startup_s:.6g}\n"
             f"source = \"{source}\"\n")
    pattern = re.compile(r"^\[campaign\]\n(?:(?!\[).*\n?)*", re.MULTILINE)
    text = pattern.sub(table, text) if pattern.search(text) else text.rstrip("\n") + "\n\n" + table
    with open(path + ".tmp", "w") as fh:
        fh.write(text)
    os.replace(path + ".tmp", path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--rank", action="store_true")
    ap.add_argument("--usl-alpha", type=float, default=0.0)
    ap.add_argument("--usl-beta", type=float, default=0.0)
    ap.add_argument("--heap-scale", type=float, default=1.0,
                    help="heap correction to apply in the last ranking column")
    ap.add_argument("--write-constants", metavar="PATH")
    args = ap.parse_args()

    tails = list(tail_rows(args.run_dir))
    if not tails:
        print("no pure-process rows found", file=sys.stderr)
        return 1
    print("pure-process rows (round tail, per-round excess), sample units:")
    for label in sorted({d[0] for d in tails}):
        vals = [d for d in tails if d[0] == label]
        print(f"  {label:24s} n={len(vals):2d}  round_tail median {st.median(v[1] for v in vals):.4f}"
              f"  per-round median {st.median(v[2] for v in vals):.3f}")
    tail = st.median(d[1] for d in tails)
    print(f"round_tail = {tail:.4f}  (median of {len(tails)} rows)")

    starts = list(startup_campaigns(args.run_dir))
    for label, ms in starts:
        print(f"  cold campaign {label:24s} take+occupancy-work median {ms:.1f} ms")
    startup_s = st.median(ms for _, ms in starts) / 1000.0 if starts else 0.0
    print(f"pool_startup_s = {startup_s:.4f}  (median of {len(starts)} cold campaigns)")

    if args.rank:
        rank(args.run_dir, args.usl_alpha, args.usl_beta, tail, args.heap_scale)

    if args.write_constants:
        source = (f"{os.path.basename(os.path.normpath(args.run_dir))}: round_tail = median over "
                  f"{len(tails)} pure-process rows; pool_startup_s = median over {len(starts)} cold campaigns")
        write_constants(args.write_constants, tail, startup_s, source)
        print(f"wrote [campaign] to {args.write_constants}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
