#!/usr/bin/env python3
"""Plot the paper's routing comparisons: serial vs best static route vs R6.

Reads the same aggregated CSVs as make_paper_routing_tables.py and reuses its
row builder, so a figure and the table beside it cannot disagree about which
route won a point or what the serial baseline was.

Each phase gets one figure: rows are machines, columns are the two things the
comparison is made of -- raw median wall time on the left, speedup against that
point's serial baseline on the right. Every pinned static route is drawn faintly
behind the best-static line, because at several points their near-coincidence is
itself the result (see finding 2).

Usage:
    python3 scripts/make_paper_routing_plots.py RUN_DIR_OR_CSV [...] \
        [--out DIR] [--label NAME] [--format png,pdf]
"""
from __future__ import annotations

import argparse
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_paper_routing_tables import (  # noqa: E402
    ADAPTIVE_LABEL, NOISE_FLOOR_S, PHASE_AXIS, PHASE_TITLE, STATIC_PARALLEL,
    load, phase_rows,
)

C_SERIAL, C_STATIC, C_R6, C_FAINT = "#6b6b6b", "#1f6fb4", "#c23b22", "#b8c6d4"
MARK_STATIC, MARK_R6 = "s", "o"

ROUTE_LABEL = {
    "outer_threads": "outer threads", "outer_process": "outer process",
    "inner_only": "inner only", "outer_inner_static": "outer+inner static",
}

# Phases whose axis is a resource count, where perfect scaling is a meaningful
# reference. P1's axis is problem size, where it is not.
IDEAL_REF = {"P2", "P3", "P4"}


def machine_of(df) -> str:
    if "machine" in df.columns:
        vals = df.machine.dropna().unique()
        if len(vals):
            return str(vals[0])
    return "unknown machine"


def _xmap(rows, phase):
    """Order value -> x coordinate. Splits are a categorical sequence, not a
    numeric one: 1x32/2x16/4x8 are equally spaced choices of how to divide one
    budget, so placing them at their worker counts both crowds the left end and
    implies a magnitude relation the axis does not have."""
    if phase == "P5":
        return {r["order"]: i for i, r in enumerate(rows)}
    return {r["order"]: r["order"] for r in rows}


def _panel(ax, rows, *, value_fn, ylabel, logy, phase, ideal=False):
    """Draw one panel. value_fn maps a seconds value to what is plotted."""
    xm = _xmap(rows, phase)
    xs = [xm[r["order"]] for r in rows]
    labels = [str(r["axis"]) for r in rows]

    for mode in STATIC_PARALLEL:
        pts = [(xm[r["order"]], value_fn(r["all_static"][mode], r))
               for r in rows if r["all_static"].get(mode)]
        pts = [(x, y) for x, y in pts if y is not None]
        if len(pts) > 1:
            ax.plot([x for x, _ in pts], [y for _, y in pts],
                    color=C_FAINT, lw=0.9, zorder=1)

    def series(key):
        out = [(xm[r["order"]], r[key], r) for r in rows if r.get(key)]
        return ([x for x, _, _ in out], [value_fn(v, r) for _, v, r in out])

    if value_fn is not _speedup:  # serial is the flat reference; as a ratio it is 1.0
        sx, sy = series("serial_s")
        ax.plot(sx, sy, color=C_SERIAL, lw=1.6, ls="--", marker="^", ms=4,
                label="serial (R0)", zorder=3)
    else:
        ax.axhline(1.0, color=C_SERIAL, lw=1.4, ls="--", zorder=3)

    bx, by = series("best_static_s")
    ax.plot(bx, by, color=C_STATIC, lw=1.9, marker=MARK_STATIC, ms=5,
            label="best static route", zorder=4)
    ax_, ay = series("adaptive_s")
    ax.plot(ax_, ay, color=C_R6, lw=1.9, marker=MARK_R6, ms=5,
            label=f"{ADAPTIVE_LABEL} (adaptive)", zorder=5)

    if ideal and value_fn is _speedup and xs:
        counts = [r["order"] for r in rows]
        lo = min(counts)
        ax.plot(xs, [c / lo for c in counts], color="#999999", lw=1.0, ls=":",
                label="perfect scaling", zorder=2)

    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
        # A decade-only log axis labels one tick on these ranges; label the
        # 1-2-3-5 subdivisions in plain seconds instead of 10^n.
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 3.0, 5.0),
                                              numticks=14))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", labelsize=8)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo * 0.82, hi * 1.18)
    if phase == "P5":
        ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
        ax.set_xlim(min(xs) - 0.35, max(xs) + 0.35)
    else:
        ax.set_xscale("log", base=2)
        ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
        ax.minorticks_off()
    ax.grid(alpha=0.25, lw=0.6)

    # Mark points the harness cannot resolve rather than dropping them.
    for r in rows:
        if r["serial_s"] and r["serial_s"] < NOISE_FLOOR_S:
            x = xm[r["order"]]
            half = 0.3 if phase == "P5" else x * 0.08
            ax.axvspan(x - half, x + half, color="#f0d0d0", alpha=0.35, zorder=0)


def _raw(v, _rec):
    return v


def _speedup(v, rec):
    s = rec.get("serial_s")
    return s / v if s and v else None


def _verify_panels(axes, per_machine, phase):
    """Read each drawn line back off the axes and check it against the records
    it was built from. A figure that disagrees with its own source data is worse
    than no figure, so this runs on every generation rather than on request."""
    bad = []
    for i, (machine, rows) in enumerate(per_machine):
        xm = _xmap(rows, phase)
        byx = {xm[r["order"]]: r for r in rows}
        for col, kind in ((0, "raw"), (1, "speedup")):
            named = {l.get_label(): l for l in axes[i][col].get_lines()
                     if not l.get_label().startswith("_")}
            for lbl, key in (("best static route", "best_static_s"),
                             (f"{ADAPTIVE_LABEL} (adaptive)", "adaptive_s"),
                             ("serial (R0)", "serial_s")):
                if lbl not in named:
                    continue
                for x, y in zip(*named[lbl].get_data()):
                    rec = byx[x]
                    want = rec[key]
                    if kind == "speedup":
                        want = rec["serial_s"] / rec[key]
                    if abs(y - want) > 1e-9 * max(1.0, abs(want)):
                        bad.append(f"{phase}/{machine}/{kind}/{lbl} x={x}: "
                                   f"plotted {y!r}, source {want!r}")
    if bad:
        raise SystemExit("plot does not match its source data:\n  " + "\n  ".join(bad))


def figure(phase, per_machine, out_dir, formats, case_label=""):
    n = len(per_machine)
    fig, axes = plt.subplots(n, 2, figsize=(10.2, 3.5 * n), squeeze=False)
    ideal = phase in IDEAL_REF
    _, axis_name = PHASE_AXIS[phase]

    for i, (machine, rows) in enumerate(per_machine):
        _panel(axes[i][0], rows, value_fn=_raw, ylabel="median wall time (s)",
               logy=True, phase=phase)
        _panel(axes[i][1], rows, value_fn=_speedup, ylabel="speedup vs serial",
               logy=False, phase=phase, ideal=ideal)
        for ax in axes[i]:
            ax.set_xlabel(axis_name)
        axes[i][0].set_title(machine, loc="left", fontsize=10, fontweight="bold")
        axes[i][1].set_title("ratio to that point's serial baseline", loc="left",
                             fontsize=9, color="#444444")

    _verify_panels(axes, per_machine, phase)

    handles, labels = axes[0][0].get_legend_handles_labels()
    handles.append(Line2D([], [], color=C_FAINT, lw=0.9))
    labels.append("other pinned static routes")
    if ideal:
        handles.append(Line2D([], [], color="#999999", lw=1.0, ls=":"))
        labels.append("perfect scaling")
    fig.legend(handles, labels, loc="lower center", ncol=min(5, len(labels)),
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.005))

    title = f"{phase} — {PHASE_TITLE[phase]}"
    if case_label:
        title += f"\n{case_label}"
    fig.suptitle(title, fontsize=11.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0.045, 1, 0.97))

    stem = f"fig_{phase.lower()}" + (f"_{case_label}" if case_label else "")
    stem = stem.replace(" ", "_")
    written = []
    for fmt in formats:
        p = os.path.join(out_dir, f"{stem}.{fmt}")
        fig.savefig(p, dpi=200, bbox_inches="tight")
        written.append(p)
    plt.close(fig)
    return written


def summary_figure(frames, out_dir, formats):
    """One panel: R6 divided by the best pinned static route at every measured
    point. 1.0 is parity with an oracle that already knew the right route; below
    it R6 is faster than any pinned route measured there."""
    phases = ["P1", "P2", "P3", "P4", "P5"]
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    machine_style = [(C_R6, "o"), (C_STATIC, "^")]

    for pi, phase in enumerate(phases):
        for mi, (name, df) in enumerate(frames):
            color, marker = machine_style[mi % len(machine_style)]
            ys = []
            for r in phase_rows(df, phase):
                if r.get("adaptive_s") and r.get("best_static_s"):
                    ys.append(r["adaptive_s"] / r["best_static_s"])
            if not ys:
                continue
            off = (mi - (len(frames) - 1) / 2) * 0.20
            xs = [pi + off + (k - len(ys) / 2) * 0.012 for k in range(len(ys))]
            ax.scatter(xs, ys, s=26, color=color, marker=marker, alpha=0.8,
                       edgecolors="white", linewidths=0.5, zorder=3,
                       label=name if pi == 0 else None)
            med = statistics.median(ys)
            ax.plot([pi + off - 0.085, pi + off + 0.085], [med, med],
                    color=color, lw=2.4, zorder=4)

    ax.axhline(1.0, color="#444444", lw=1.3, ls="--", zorder=2)
    ax.text(0.995, 1.03, "parity with the best pinned route", transform=
            ax.get_yaxis_transform(), fontsize=8, color="#444444", va="bottom",
            ha="right")
    # The two deepest P1 points are the batched-RHS heuristic failing on every
    # pinned route (finding 2), not the router out-routing them. Say so on the
    # figure rather than letting a 10x bar read as a routing win.
    ax.annotate("N=64 and N=256: every pinned route hits the\n"
                "batched-RHS heuristic (finding 2), not a routing loss",
                xy=(0.12, 0.10), xytext=(0.52, 0.16), fontsize=8, color="#555555",
                arrowprops=dict(arrowstyle="->", color="#999999", lw=0.9))
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(0.5, 0.7, 1.0, 1.5, 2.0, 3.0),
                                          numticks=14))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(range(len(phases)))
    short = {"P1": "constellation size", "P2": "thread budget",
             "P3": "MC ladder (cheap)", "P4": "MC ladder (heavy)",
             "P5": "worker/thread split"}
    ax.set_xticklabels([f"{p}\n{short[p]}" for p in phases], fontsize=8.5)
    ax.set_ylabel(f"{ADAPTIVE_LABEL} time / best static route time")
    ax.grid(alpha=0.25, lw=0.6, axis="y")
    ax.set_title("How close the adaptive router gets to the best pinned route\n"
                 "every measured point; thick bar is the per-phase median",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=9, loc="lower right", bbox_to_anchor=(1.0, 0.02))
    fig.tight_layout()

    written = []
    for fmt in formats:
        q = os.path.join(out_dir, f"fig_summary_regret.{fmt}")
        fig.savefig(q, dpi=200, bbox_inches="tight")
        written.append(q)
    plt.close(fig)
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sources", nargs="+")
    ap.add_argument("--out", default="output/paper_routing_plots")
    ap.add_argument("--label", action="append", default=[])
    ap.add_argument("--format", default="png,pdf")
    args = ap.parse_args()

    formats = [f.strip() for f in args.format.split(",") if f.strip()]
    os.makedirs(args.out, exist_ok=True)

    frames = []
    for i, src in enumerate(args.sources):
        df = load(src)
        name = args.label[i] if i < len(args.label) else machine_of(df)
        frames.append((name, df))

    written = []
    for phase in ["P1", "P2", "P3", "P4", "P5"]:
        cases = set()
        for _, df in frames:
            cases |= {r["case"] for r in phase_rows(df, phase)}
        # P1 puts every case on one axis; the others are one figure per workload.
        groups = [None] if phase == "P1" else sorted(cases)
        for case in groups:
            per_machine = []
            for name, df in frames:
                rows = phase_rows(df, phase)
                if case is not None:
                    rows = [r for r in rows if r["case"] == case]
                rows = sorted(rows, key=lambda r: r["order"])
                if rows:
                    per_machine.append((name, rows))
            if per_machine:
                written += figure(phase, per_machine, args.out, formats,
                                  case_label=case or "")

    written += summary_figure(frames, args.out, formats)

    for p in written:
        print("wrote", p)


if __name__ == "__main__":
    main()
