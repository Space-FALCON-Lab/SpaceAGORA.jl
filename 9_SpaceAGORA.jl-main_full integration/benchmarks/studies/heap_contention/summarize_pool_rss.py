#!/usr/bin/env python3
"""Summarize peak RSS per process and peak total RSS from pool_rss_probe.sh's
`ps`-sampling logs, and report the wall-time ratio between two arms.

Usage:
    summarize_pool_rss.py <nohint_rss_log> <hint_rss_log> <nohint_run_log> <hint_run_log>

Each rss log line is: "<unix_ts> <pid> <ppid> <rss_kb> <args...>" (one line
per julia-named process per 2 s sample). Each run log is the julia script's
own stdout, which prints "MARKER wall_s=<float>" once the timed batch ends.
"""
import re
import sys
from collections import defaultdict


def parse_rss_log(path):
    """Parse the raw `ps` sample log and keep only THIS run's own processes:
    the coordinator (identified by "pool_rss_probe.jl" appearing in its own
    argv, which is unique to it) plus its direct Distributed pool workers
    (ppid == coordinator pid and "--worker" in argv). The plain `grep julia`
    the sampling loop uses also catches unrelated julia-named processes on a
    shared machine -- another agent's job, the VS Code language server, a
    precompile subprocess -- which is not this run and must not count toward
    its peak RSS.
    """
    rows = []  # (ts, pid, ppid, rss_kb, args)
    with open(path) as f:
        for line in f:
            fields = line.split(None, 4)
            if len(fields) < 5:
                continue
            try:
                ts = float(fields[0])
                pid = int(fields[1])
                ppid = int(fields[2])
                rss_kb = int(fields[3])
            except ValueError:
                continue
            rows.append((ts, pid, ppid, rss_kb, fields[4]))

    coordinator_pids = {pid for (_, pid, _, _, args) in rows if "pool_rss_probe.jl" in args}
    if not coordinator_pids:
        return {}, {}, None
    coordinator_pid = min(coordinator_pids)  # stable across ticks; only one per log in practice

    def is_this_run(pid, ppid, args):
        if pid == coordinator_pid:
            return True
        return ppid == coordinator_pid and "--worker" in args

    per_pid = defaultdict(list)
    per_ts = defaultdict(int)
    for ts, pid, ppid, rss_kb, args in rows:
        if not is_this_run(pid, ppid, args):
            continue
        per_pid[pid].append((ts, rss_kb))
        per_ts[ts] += rss_kb
    return per_pid, per_ts, coordinator_pid


def summarize_arm(label, rss_log):
    per_pid, per_ts, coordinator_pid = parse_rss_log(rss_log)
    if not per_pid:
        print(f"{label}: no samples parsed from {rss_log}")
        return
    peak_per_pid = {pid: max(rss for _, rss in samples) for pid, samples in per_pid.items()}
    top5 = sorted(peak_per_pid.items(), key=lambda kv: -kv[1])[:5]
    peak_total_kb = max(per_ts.values()) if per_ts else 0
    n_samples = len(per_ts)
    print(f"{label}: coordinator pid={coordinator_pid}, {len(per_pid)} pids belonging to this run "
          f"(coordinator + pool workers, other julia-named processes on the shared machine excluded), "
          f"{n_samples} sampling ticks")
    print(f"{label}: peak TOTAL RSS across this run's own processes at one tick = {peak_total_kb/1024/1024:.2f} GiB")
    print(f"{label}: per-process peaks (GiB), highest first:")
    for pid, rss_kb in top5:
        tag = "coordinator" if pid == coordinator_pid else "worker"
        print(f"{label}:   pid={pid} ({tag}) peak_rss={rss_kb/1024/1024:.2f} GiB")


def wall_time(run_log):
    with open(run_log) as f:
        text = f.read()
    m = re.search(r"MARKER wall_s=([0-9.]+)", text)
    return float(m.group(1)) if m else None


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    nohint_rss, hint_rss, nohint_run, hint_run = sys.argv[1:5]

    summarize_arm("no-hint", nohint_rss)
    print()
    summarize_arm("hint", hint_rss)
    print()

    w_nohint = wall_time(nohint_run)
    w_hint = wall_time(hint_run)
    print(f"wall_s no-hint = {w_nohint}")
    print(f"wall_s hint    = {w_hint}")
    if w_nohint and w_hint:
        print(f"ratio hint/no-hint = {w_hint / w_nohint:.4f}")


if __name__ == "__main__":
    main()
