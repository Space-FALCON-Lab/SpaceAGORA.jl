#!/usr/bin/env python3
"""Bar charts of the adaptive routing profiles (R4, R5, R6) relative to the best
static route, one figure per paper-benchmark phase.

    python3 scripts/plot_adaptive_vs_best_static.py                       # every outdir
    python3 scripts/plot_adaptive_vs_best_static.py output/performance/paper_benchmarks/20260903_2148*
    python3 scripts/plot_adaptive_vs_best_static.py --group-by=route --phases=B7,B13

Reads the paper harness's raw CSV from each outdir given (default: every
directory under output/performance/paper_benchmarks) AND the per-point worker
row files under <outdir>/<phase>/**/worker_rows/perf_*.csv, so a phase that is
still running, or a run that died before the raw file was written, still
plots whatever points have completed. Where the same point appears in more
than one outdir the newest outdir wins, so a phase re-run after a fix
supersedes the earlier attempt without the earlier outdir having to be moved.

Vertical axis: speedup of the adaptive profile relative to the best static
route at the SAME launch point, i.e. best_static_median / adaptive_median.
1.0 is parity, above 1.0 the router beat every pinned route, below 1.0 it
paid regret. A launch point is (case, thread_count, process_workers,
mc_samples): every mode in a phase is launched on that same grid, and a
pinned process route present at a point means the same process pool was
available to the adaptive route there, so this is the "given what it had,
did the router choose well" comparison. (reporting.jl's regret table instead
keys on the cores a route actually spent; the two agree wherever the
thread-backed and process-backed routes are launched at the same size and
differ on the worker-ladder phases, where the launch-point comparison is
the stricter one for the router.)

Horizontal axis: one group of bars per value of the phase's primary
independent variable (thread count; the workers x threads split on the
fixed-budget phases; process workers on the worker-ladder phases), with one
bar per adaptive profile inside each group. --group-by=route flips that: one
group per profile, one bar per variable value. A phase with a second
independent variable (sample count on the Monte Carlo phases) gets one panel
per (case, value of that variable). The route that was best static at each
point is written under its group, and points whose serial baseline is under
the noise floor are hatched.

Writes <plot-dir>/<phase>_relative_speedup.png and a CSV of the plotted
values. Nothing here needs the run to be complete.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_ROOT = os.path.join(REPO, "output", "performance", "paper_benchmarks")

STATIC_MODES = ("serial", "outer_threads", "outer_process", "inner_only", "outer_inner_static")
ADAPTIVE = [  # mode -> (label, colour)
    ("outer_inner_adaptive", "R4", "#6b7fa3"),
    ("full_smart", "R5", "#d98c3f"),
    ("policy_v2", "R6", "#3f9a6e"),
    ("policy_v2_learn", "R6L", "#8ccf9f"),
]
STATIC_SHORT = {
    "serial": "serial", "outer_threads": "threads", "outer_process": "process",
    "inner_only": "inner", "outer_inner_static": "static",
}
NOISE_FLOOR_SERIAL_S = 3.0  # reporting.jl's PPB_NOISE_FLOOR_SERIAL_S

PHASE_LABELS = {
    "B1": "Constellation scaling, gravity-only baseline",
    "B2": "Constellation scaling, L20 harmonics",
    "B3": "Constellation scaling, atmosphere and GRAM surrogate",
    "B4": "Monte Carlo throughput scaling",
    "B5": "Routing profile comparison at max threads",
    "B6": "Small-workload control (below the floor)",
    "B7": "Heavy constellation thread ladder",
    "B8": "Heavy Monte Carlo process throughput",
    "B9": "Router evaluation: spacecraft count",
    "B10": "Router evaluation: atmosphere and GRAM usage",
    "B11": "Router evaluation: force and actuator model count",
    "B12": "Router evaluation: interacting vs. independent",
    "B13": "Router evaluation: thread vs. process budget split",
    "B14": "Router evaluation: mission duration and output cadence",
    "B15": "Router evaluation: nested campaign aspect ratio",
}


def phase_sort_key(p: str):
    return (0, int(p[1:])) if p[1:].isdigit() else (1, p)


# ── Loading ──────────────────────────────────────────────────────────────────

def _to_int(v, default=1) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _load_csv_rows(path: str):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _phase_from_worker_path(outdir: str, path: str):
    rel = os.path.relpath(path, outdir).split(os.sep)
    return rel[0] if rel and rel[0].startswith("B") else None


def load_outdirs(outdirs: list[str]):
    """Returns {point_key: {"outdir", "phase", "times": [...], "mode", ...}}.

    point_key = (phase, case, mode, thread_count, process_workers, mc_samples).
    Rows from the raw CSV and from worker_rows are merged per outdir (worker
    rows only add repeats the raw file does not carry yet); across outdirs
    the newest one holding a point wins outright.
    """
    points = {}
    for outdir in sorted(outdirs, key=os.path.basename):
        per_point = defaultdict(lambda: {"raw": [], "worker": []})
        for raw in glob.glob(os.path.join(outdir, "paper_benchmarks_raw_*.csv")):
            for r in _load_csv_rows(raw):
                phase = r.get("phase_id", "")
                if not phase:
                    continue
                per_point[_key(phase, r)]["raw"].append(r)
        for wf in glob.glob(os.path.join(outdir, "B*", "**", "worker_rows", "perf_*.csv"), recursive=True):
            phase = _phase_from_worker_path(outdir, wf)
            if phase is None:
                continue
            for r in _load_csv_rows(wf):
                per_point[_key(phase, r)]["worker"].append(r)
        for key, src in per_point.items():
            rows = src["raw"] if len(src["raw"]) >= len(src["worker"]) else src["worker"]
            times = [float(r["wall_time_s"]) for r in rows
                     if r.get("success", "true").lower() in ("true", "1", "") and r.get("wall_time_s")]
            if not times:
                continue
            points[key] = {
                "outdir": outdir, "phase": key[0], "case": key[1], "mode": key[2],
                "threads": key[3], "workers": key[4], "samples": key[5],
                "times": times, "median": statistics.median(times), "n": len(times),
            }
    return points


def _key(phase: str, r: dict):
    return (
        phase, r["case"], r["mode"],
        _to_int(r.get("thread_count") or r.get("julia_threads"), 1),
        _to_int(r.get("process_workers"), 1),
        _to_int(r.get("mc_samples"), 1),
    )


# ── Relative speedups ────────────────────────────────────────────────────────

def compute_relative(points: dict):
    """Per (phase, case, threads, workers, samples): best static median, its
    mode, the serial baseline, and each adaptive mode's speedup against it."""
    by_launch = defaultdict(dict)
    serial_at = {}
    for p in points.values():
        lk = (p["phase"], p["case"], p["threads"], p["workers"], p["samples"])
        by_launch[lk][p["mode"]] = p
        if p["mode"] == "serial":
            serial_at[(p["phase"], p["case"], p["samples"])] = p["median"]
    out = []
    for lk, modes in by_launch.items():
        statics = {m: p for m, p in modes.items() if m in STATIC_MODES}
        if not statics:
            continue
        best_mode = min(statics, key=lambda m: statics[m]["median"])
        best = statics[best_mode]["median"]
        serial = serial_at.get((lk[0], lk[1], lk[4]))
        for mode, label, _ in ADAPTIVE:
            if mode not in modes:
                continue
            out.append({
                "phase": lk[0], "case": lk[1], "threads": lk[2], "workers": lk[3], "samples": lk[4],
                "mode": mode, "profile": label, "median_s": modes[mode]["median"],
                "n": modes[mode]["n"], "best_static_mode": best_mode, "best_static_s": best,
                "speedup": best / modes[mode]["median"],
                "serial_s": serial if serial is not None else float("nan"),
                "below_floor": (serial is not None and serial < NOISE_FLOOR_SERIAL_S),
            })
    return out


# ── Layout ───────────────────────────────────────────────────────────────────

def _independent_variables(rows: list[dict]):
    """Which of threads / workers / samples vary within this phase, and how to
    split them into the primary (x groups) and secondary (panels) variable."""
    varies = {v: len({r[v] for r in rows}) > 1 for v in ("threads", "workers", "samples")}
    if varies["threads"] and varies["workers"]:
        # Fixed-budget split (B13, B15): one label carries both.
        primary = ("split", lambda r: f"{r['workers']}w×{r['threads']}t",
                   lambda r: (r["workers"], r["threads"]))
        secondary = ["samples"] if varies["samples"] else []
    elif varies["threads"]:
        primary = ("threads", lambda r: f"{r['threads']} thr", lambda r: r["threads"])
        secondary = [v for v in ("workers", "samples") if varies[v]]
    elif varies["workers"]:
        primary = ("workers", lambda r: f"{r['workers']} wkr", lambda r: r["workers"])
        secondary = ["samples"] if varies["samples"] else []
    elif varies["samples"]:
        primary = ("samples", lambda r: f"{r['samples']} smp", lambda r: r["samples"])
        secondary = []
    else:
        primary = ("threads", lambda r: f"{r['threads']} thr", lambda r: r["threads"])
        secondary = []
    return primary, secondary


def _secondary_label(r: dict, secondary: list[str]) -> str:
    names = {"workers": "workers", "samples": "samples", "threads": "threads"}
    return ", ".join(f"{r[v]} {names[v]}" for v in secondary)


def plot_phase(phase: str, rows: list[dict], plot_dir: str, group_by: str) -> str:
    primary, secondary = _independent_variables(rows)
    _, plabel, psort = primary
    panels = defaultdict(list)
    for r in rows:
        panels[(r["case"], tuple(r[v] for v in secondary))].append(r)
    panel_keys = sorted(panels, key=lambda k: (k[0], k[1]))
    n = len(panel_keys)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 4.2 * nrows), squeeze=False)
    present = [(m, l, c) for m, l, c in ADAPTIVE if any(r["mode"] == m for r in rows)]

    for ax, pk in zip(axes.flat, panel_keys):
        prow = panels[pk]
        xvals = sorted({psort(r): plabel(r) for r in prow}.items())
        xlabels = [lbl for _, lbl in xvals]
        if group_by == "route":
            groups = [lbl for _, lbl, _ in present]
            series = [(lbl, None) for lbl in xlabels]
        else:
            groups = xlabels
            series = [(lbl, col) for _, lbl, col in present]
        nser = max(1, len(series))
        width = 0.8 / nser
        best_under = defaultdict(list)
        cmap = plt.get_cmap("Blues")
        for si, (slabel, colour) in enumerate(series):
            xs, ys, hatches = [], [], []
            for gi, glabel in enumerate(groups):
                if group_by == "route":
                    prof, xl = glabel, slabel
                else:
                    prof, xl = slabel, glabel
                match = [r for r in prow if r["profile"] == prof and plabel(r) == xl]
                if not match:
                    continue
                r = match[0]
                xs.append(gi + (si - (nser - 1) / 2) * width)
                ys.append(r["speedup"])
                hatches.append("//" if r["below_floor"] else "")
                best_under[gi].append(STATIC_SHORT.get(r["best_static_mode"], r["best_static_mode"]))
            if group_by == "route":
                colour = cmap(0.35 + 0.6 * si / max(1, nser - 1))
            bars = ax.bar(xs, ys, width=width * 0.95, color=colour, edgecolor="#333", linewidth=0.5, label=slabel)
            for b, h in zip(bars, hatches):
                if h:
                    b.set_hatch(h)
        ax.axhline(1.0, color="#333", linewidth=0.9, linestyle="--")
        ticks = []
        for gi, glabel in enumerate(groups):
            best = sorted(set(best_under.get(gi, [])))
            ticks.append(glabel + ("\nbest: " + "/".join(best) if best else ""))
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(ticks, fontsize=8)
        ax.set_ylabel("speedup vs. best static route")
        title = pk[0] + (f"  ({_secondary_label(prow[0], secondary)})" if secondary else "")
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8, title="profile" if group_by != "route" else primary[0], loc="best")
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)

    note = "hatched: serial baseline under the %.0f s noise floor" % NOISE_FLOOR_SERIAL_S
    fig.suptitle(f"{phase}: {PHASE_LABELS.get(phase, '')}\nadaptive profiles relative to the best static route at the same launch point",
                 fontsize=12)
    fig.text(0.99, 0.01, note, ha="right", va="bottom", fontsize=8, color="#555")
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    os.makedirs(plot_dir, exist_ok=True)
    path = os.path.join(plot_dir, f"{phase}_relative_speedup.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


# ── Entry point ──────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdirs", nargs="*", help="paper benchmark outdirs (default: all under output/performance/paper_benchmarks)")
    ap.add_argument("--plot-dir", default=None, help="where to write PNGs (default: <root>/plots_relative)")
    ap.add_argument("--phases", default="", help="comma-separated phase ids to plot (default: every phase with data)")
    ap.add_argument("--group-by", choices=("budget", "route"), default="budget",
                    help="x groups are budget values with one bar per profile (default), or profiles with one bar per budget value")
    args = ap.parse_args(argv)

    outdirs = args.outdirs or sorted(d for d in glob.glob(os.path.join(DEFAULT_ROOT, "*")) if os.path.isdir(d))
    outdirs = [d for d in outdirs if os.path.isdir(d)]
    if not outdirs:
        print("no outdirs found", file=sys.stderr)
        return 1
    plot_dir = args.plot_dir or os.path.join(DEFAULT_ROOT, "plots_relative")

    points = load_outdirs(outdirs)
    rel = compute_relative(points)
    if not rel:
        print("no adaptive-vs-static points found in", ", ".join(outdirs), file=sys.stderr)
        return 1
    wanted = {p.strip() for p in args.phases.split(",") if p.strip()}
    phases = sorted({r["phase"] for r in rel if not wanted or r["phase"] in wanted}, key=phase_sort_key)

    os.makedirs(plot_dir, exist_ok=True)
    csv_path = os.path.join(plot_dir, "relative_speedups.csv")
    fields = ["phase", "case", "threads", "workers", "samples", "mode", "profile", "median_s", "n",
              "best_static_mode", "best_static_s", "speedup", "serial_s", "below_floor"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in sorted(rel, key=lambda r: (phase_sort_key(r["phase"]), r["case"], r["samples"], r["workers"], r["threads"], r["profile"])):
            w.writerow({k: r[k] for k in fields})

    for phase in phases:
        rows = [r for r in rel if r["phase"] == phase]
        path = plot_phase(phase, rows, plot_dir, args.group_by)
        worst = min(rows, key=lambda r: r["speedup"])
        print(f"{phase:4s} {len(rows):3d} points -> {path}   worst: {worst['profile']} {worst['case']} "
              f"t{worst['threads']} w{worst['workers']} s{worst['samples']} = {worst['speedup']:.2f}x")
    print(f"values: {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
