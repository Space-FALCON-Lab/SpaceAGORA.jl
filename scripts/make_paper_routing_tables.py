#!/usr/bin/env python3
"""Build the paper's routing comparison tables from a P-series benchmark run.

Reads one or more `paper_benchmarks_aggregated_*.csv` files (the P1-P5 phases in
benchmarks/studies/paper_parallelization_benchmarks) and writes, per machine and
per phase, a table of

    serial | best static route | R6 (policy_v2)

with the raw median wall times and the ratio of each against that point's serial
baseline. The best static route is the fastest *parallel* static route measured
at the same point -- serial is reported in its own column, so including it in
the minimum would hide the comparison the table exists to make.

Usage:
    python3 scripts/make_paper_routing_tables.py RUN_DIR_OR_CSV [...] \
        [--out DIR] [--label NAME]

Each positional argument is either an aggregated CSV or a run directory holding
one. --label renames the machine for that source (repeatable, positional order).
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
ADAPTIVE_LABEL = "R6"

PHASE_TITLE = {
    "P1": "Constellation size scaling (fixed thread budget)",
    "P2": "Thread scaling at a fixed spacecraft count",
    "P3": "Monte Carlo resource ladder, one spacecraft per sample",
    "P4": "Monte Carlo resource ladder, compute-bound samples",
    "P5": "Monte Carlo over constellations, worker/thread split at a fixed budget",
}

# Phases whose axis is derived from the case name itself (P1's spacecraft count)
# put every case in ONE table, one row per axis value; phases that measure
# several distinct workloads over the same axis get one table each.
GROUP_BY_CASE = {"P1": False, "P2": True, "P3": True, "P4": True, "P5": True}

# Below this serial baseline the harness treats a point as unreportable router
# performance -- dispatch overhead and machine noise are the same size as the
# difference being measured. Marked, not dropped: the small rungs are still the
# evidence that the router does not *lose* where there is nothing to win.
NOISE_FLOOR_S = 3.0

# Which column is the x-axis of each phase's table, and how to label it.
PHASE_AXIS = {
    "P1": ("n_sat", "spacecraft"),
    "P2": ("thread_count", "threads"),
    "P3": ("process_workers", "budget"),
    "P4": ("process_workers", "budget"),
    "P5": ("split", "workers x threads"),
}


def _mission_s_from_case(case: str):
    """Simulated mission length when the case name carries one (the iso-work ladder)."""
    m = re.search(r"_(\d+)s$", case)
    return int(m.group(1)) if m else None


def _n_sat_from_case(case: str) -> int:
    m = re.search(r"(\d+)sat", case)
    return int(m.group(1)) if m else 1


def load(path: str) -> pd.DataFrame:
    if os.path.isdir(path):
        hits = sorted(glob.glob(os.path.join(path, "paper_benchmarks_aggregated_*.csv")))
        if not hits:
            raise SystemExit(f"no aggregated CSV under {path}")
        path = hits[-1]
    df = pd.read_csv(path)
    df["_source"] = path
    if "machine" not in df.columns:
        # The aggregated CSV groups the machine column away; the raw CSV beside
        # it keeps it, and it is what names the table.
        raw = path.replace("_aggregated_", "_raw_")
        if os.path.exists(raw):
            try:
                machines = pd.read_csv(raw, usecols=["machine"]).machine.dropna().unique()
                if len(machines):
                    df["machine"] = machines[0]
            except (ValueError, KeyError):
                pass
    return df


def axis_values(df: pd.DataFrame, phase: str) -> pd.DataFrame:
    col, _ = PHASE_AXIS[phase]
    if col == "split":
        df = df.copy()
        df["split"] = [
            f"{int(w)}x{int(t)}" for w, t in zip(df.process_workers, df.thread_count)
        ]
        df["_order"] = df.process_workers.astype(int)
    elif col == "n_sat":
        df = df.copy()
        df["n_sat"] = [_n_sat_from_case(c) for c in df.case]
        df["_order"] = df.n_sat
    else:
        df = df.copy()
        df["_order"] = df[col].astype(int)
    return df


def phase_rows(df: pd.DataFrame, phase: str) -> list[dict]:
    """One record per (case, axis point): serial, best static, adaptive."""
    sub = df[df.phase_id == phase]
    if sub.empty:
        return []
    sub = axis_values(sub, phase)
    axis_col, _ = PHASE_AXIS[phase]
    out = []
    for (case, axis, order), grp in sub.groupby(["case", axis_col, "_order"]):
        times = {}
        for mode, mgrp in grp.groupby("mode"):
            t = mgrp.wall_time_median_s.dropna()
            if len(t):
                times[mode] = float(t.min())
        # The harness runs serial ONCE per phase where it is thread-independent
        # (a single-simulation phase's thread ladder), and propagates that median
        # into every row's serial_median_s. Reading the column rather than
        # looking for a serial row at this axis point is what gives a thread
        # ladder a ratio at every rung instead of only at the rung serial ran on.
        serial = times.get("serial")
        if serial is None and "serial_median_s" in grp.columns:
            baseline = grp.serial_median_s.dropna()
            baseline = baseline[baseline > 0.0]
            if len(baseline):
                serial = float(baseline.median())
        statics = {m: t for m, t in times.items() if m in STATIC_PARALLEL}
        best_mode = min(statics, key=statics.get) if statics else None
        rec = {
            "case": case,
            "mission_s": _mission_s_from_case(case),
            "axis": axis,
            "order": order,
            "serial_s": serial,
            "best_static_mode": best_mode,
            "best_static_s": statics.get(best_mode) if best_mode else None,
            "adaptive_s": times.get(ADAPTIVE),
            "all_static": statics,
        }
        out.append(rec)
    out.sort(key=lambda r: (r["case"], r["order"]))
    return out


def _case_groups(phase: str, rows: list[dict]) -> list[tuple[str, list[dict]]]:
    if GROUP_BY_CASE.get(phase, True):
        return [(case, [r for r in rows if r["case"] == case])
                for case in sorted({r["case"] for r in rows})]
    return [("", sorted(rows, key=lambda r: r["order"]))]


def _floor_mark(rec: dict) -> str:
    s = rec.get("serial_s")
    return "*" if s is not None and s < NOISE_FLOOR_S else ""


def _ratio(num, den):
    if num is None or den is None or den <= 0:
        return None
    return num / den


def _fmt(x, digits=3):
    return "—" if x is None else f"{x:.{digits}f}"


def markdown_table(phase: str, rows: list[dict]) -> str:
    _, axis_label = PHASE_AXIS[phase]
    lines = []
    marked = False
    for case, crows in _case_groups(phase, rows):
        if case:
            lines.append(f"\n**{case}**\n")
        show_mission = any(r.get("mission_s") for r in crows)
        mission_head = " mission (h) |" if show_mission else ""
        mission_rule = "---:|" if show_mission else ""
        lines.append(
            f"| {axis_label} |{mission_head} serial (s) | best static (s) | route | {ADAPTIVE_LABEL} (s) "
            f"| serial/best static | serial/{ADAPTIVE_LABEL} | {ADAPTIVE_LABEL}/best static |"
        )
        lines.append(f"|---|{mission_rule}---:|---:|---|---:|---:|---:|---:|")
        for r in crows:
            marked = marked or bool(_floor_mark(r))
            sp_static = _ratio(r["serial_s"], r["best_static_s"])
            sp_adapt = _ratio(r["serial_s"], r["adaptive_s"])
            vs_static = _ratio(r["adaptive_s"], r["best_static_s"])
            mission_cell = (
                f" {r['mission_s'] / 3600:.2f} |" if show_mission and r.get("mission_s") else
                (" — |" if show_mission else "")
            )
            lines.append(
                f"| {r['axis']}{_floor_mark(r)} |{mission_cell} {_fmt(r['serial_s'])} | {_fmt(r['best_static_s'])} "
                f"| {r['best_static_mode'] or '—'} | {_fmt(r['adaptive_s'])} "
                f"| {_fmt(sp_static, 2)}x | {_fmt(sp_adapt, 2)}x | {_fmt(vs_static, 2)} |"
            )
    if marked:
        lines.append(
            f"\n`*` serial baseline under {NOISE_FLOOR_S:.0f} s — at this size the "
            "point measures dispatch overhead and machine noise, not routing."
        )
    return "\n".join(lines)


def full_markdown_table(phase: str, rows: list[dict]) -> str:
    """Every measured mode, raw medians and speedup against serial."""
    _, axis_label = PHASE_AXIS[phase]
    modes = ["serial"] + STATIC_PARALLEL + [ADAPTIVE]
    lines = []
    for case, crows in _case_groups(phase, rows):
        present = [
            m
            for m in modes
            if any(
                (m == "serial" and r["serial_s"] is not None)
                or (m == ADAPTIVE and r["adaptive_s"] is not None)
                or (m in r["all_static"])
                for r in crows
            )
        ]
        heading = f"**{case}** — " if case else ""
        lines.append(f"\n{heading}median wall time in seconds (speedup vs serial)\n")
        lines.append("| " + axis_label + " | " + " | ".join(present) + " |")
        lines.append("|---" * (len(present) + 1) + "|")
        for r in crows:
            cells = []
            for m in present:
                if m == "serial":
                    t = r["serial_s"]
                elif m == ADAPTIVE:
                    t = r["adaptive_s"]
                else:
                    t = r["all_static"].get(m)
                if t is None:
                    cells.append("—")
                else:
                    sp = _ratio(r["serial_s"], t)
                    cells.append(f"{t:.3f}" + (f" ({sp:.2f}x)" if sp else ""))
            lines.append(f"| {r['axis']} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def latex_table(phase: str, rows: list[dict], machine: str) -> str:
    _, axis_label = PHASE_AXIS[phase]
    show_mission = any(r.get("mission_s") for r in rows)
    cols = "lrrrlrrrr" if show_mission else "lrrlrrrr"
    ncols = 9 if show_mission else 8
    mission_head = "mission [h] & " if show_mission else ""
    caption = (
        f"{PHASE_TITLE.get(phase, phase)} on {machine}. Median wall time over the "
        "timed repeats; ratios are against the serial baseline at the same point. "
        "The best static route is the fastest pinned parallel route measured at "
        "that point, named in the route column."
    )
    if show_mission:
        caption += (
            " The ladder holds the total work fixed rather than the mission "
            "length, so that every rung's serial baseline is measurable; the "
            "simulated mission length of each rung is given alongside it."
        )
    head = (
        "\\begin{table}[htbp]\n\\centering\n"
        f"\\caption{{{caption}}}\n"
        f"\\label{{tab:routing_{phase.lower()}_{re.sub(r'[^a-zA-Z0-9]', '', machine)}}}\n"
        f"\\begin{{tabular}}{{{cols}}}\n\\toprule\n"
        f"{axis_label} & {mission_head}serial [s] & best static [s] & route & {ADAPTIVE_LABEL} [s] & "
        f"$T_1/T_\\text{{static}}$ & $T_1/T_\\text{{{ADAPTIVE_LABEL}}}$ & "
        f"$T_\\text{{{ADAPTIVE_LABEL}}}/T_\\text{{static}}$ \\\\\n\\midrule\n"
    )
    body = []
    for case, crows in _case_groups(phase, rows):
        if case and len({r["case"] for r in rows}) > 1:
            label = case.replace("_", "\\_")
            body.append(f"\\multicolumn{{{ncols}}}{{l}}{{\\textit{{{label}}}}} \\\\\n")
        for r in crows:
            sp_static = _ratio(r["serial_s"], r["best_static_s"])
            sp_adapt = _ratio(r["serial_s"], r["adaptive_s"])
            vs_static = _ratio(r["adaptive_s"], r["best_static_s"])
            route = (r["best_static_mode"] or "--").replace("_", "\\_")
            mission_cell = ""
            if show_mission:
                mission_cell = (
                    f"{r['mission_s'] / 3600:.2f} & " if r.get("mission_s") else "-- & "
                )
            body.append(
                f"{r['axis']} & {mission_cell}{_fmt(r['serial_s'])} & {_fmt(r['best_static_s'])} & {route} & "
                f"{_fmt(r['adaptive_s'])} & {_fmt(sp_static, 2)} & {_fmt(sp_adapt, 2)} & "
                f"{_fmt(vs_static, 2)} \\\\\n"
            )
    tail = "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    return head + "".join(body) + tail


def warmth_table(cold: pd.DataFrame, warm: pd.DataFrame, phase: str) -> str:
    """Cold-store against warm-store for one phase, mode by mode.

    R6 carries state across campaigns -- the persisted RHS-calibration verdicts
    and inner-policy hints -- so the same point measured on a machine that has
    never run the workload and on one that has is not the same measurement. The
    static routes form no such state and are the control: if they move between
    the two columns, the difference is machine noise rather than warmth.
    """
    _, axis_label = PHASE_AXIS[phase]
    crows = {(r["case"], r["axis"]): r for r in phase_rows(cold, phase)}
    wrows = {(r["case"], r["axis"]): r for r in phase_rows(warm, phase)}
    keys = [k for k in wrows if k in crows]
    if not keys:
        return ""
    # Order by the axis value, not by the case name: a phase whose axis lives in
    # the case name (P1's spacecraft count) would otherwise sort 1, 1024, 16, 256.
    keys.sort(key=lambda k: (wrows[k]["order"], k[0]))
    # A phase present in both runs but measured only once (P1/P2/P5 are not
    # re-run for the warm pass) would print a column of exact 1.00s, which reads
    # as a result rather than as the same numbers twice.
    if all(
        crows[k]["adaptive_s"] == wrows[k]["adaptive_s"]
        and crows[k]["serial_s"] == wrows[k]["serial_s"]
        for k in keys
    ):
        return ""
    lines = [
        f"\n**{phase}** — {PHASE_TITLE[phase]}\n",
        f"| {axis_label} | serial cold | serial warm | best static cold | best static warm "
        f"| {ADAPTIVE_LABEL} cold | {ADAPTIVE_LABEL} warm | {ADAPTIVE_LABEL} warm/cold |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for k in keys:
        c, w = crows[k], wrows[k]
        gain = _ratio(w["adaptive_s"], c["adaptive_s"])
        lines.append(
            f"| {w['axis']} | {_fmt(c['serial_s'])} | {_fmt(w['serial_s'])} "
            f"| {_fmt(c['best_static_s'])} | {_fmt(w['best_static_s'])} "
            f"| {_fmt(c['adaptive_s'])} | {_fmt(w['adaptive_s'])} | {_fmt(gain, 2)} |"
        )
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sources", nargs="+")
    ap.add_argument("--out", default=None, help="directory for the generated files")
    ap.add_argument("--label", action="append", default=[], help="machine label per source")
    ap.add_argument(
        "--cold",
        default=None,
        help="a run of the same phases taken from an empty calibration store; adds a "
             "cold-vs-warm section for every phase both runs contain",
    )
    args = ap.parse_args()

    md_parts = ["# Routing comparison tables (R6 vs serial vs best static route)\n"]
    tex_parts = []
    for i, src in enumerate(args.sources):
        df = load(src)
        machine = (
            args.label[i]
            if i < len(args.label)
            else (str(df.machine.iloc[0]) if "machine" in df.columns and len(df) else f"source{i}")
        )
        md_parts.append(f"\n## {machine}\n")
        md_parts.append(f"\nSource: `{df._source.iloc[0]}`\n")
        for phase in ["P1", "P2", "P3", "P4", "P5"]:
            rows = phase_rows(df, phase)
            if not rows:
                continue
            md_parts.append(f"\n### {phase} — {PHASE_TITLE[phase]}\n")
            md_parts.append(markdown_table(phase, rows))
            md_parts.append("\n\nAll measured routes:\n")
            md_parts.append(full_markdown_table(phase, rows))
            md_parts.append("\n")
            tex_parts.append(latex_table(phase, rows, machine))

    if args.cold:
        cold_df = load(args.cold)
        warm_df = load(args.sources[0])
        md_parts.append("\n## Calibration-store warmth\n")
        md_parts.append(
            "\nR6 persists what it learns: RHS-calibration verdicts and inner-policy "
            "hints outlive the campaign that produced them. A point measured on a "
            "machine that has never run the workload is therefore not the same "
            "measurement as one taken after the store has converged. The static "
            "routes hold no such state and act as the control.\n"
        )
        for phase in ["P1", "P2", "P3", "P4", "P5"]:
            section = warmth_table(cold_df, warm_df, phase)
            section and md_parts.append(section + "\n")

    md = "\n".join(md_parts)
    tex = "\n".join(tex_parts)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "routing_tables.md"), "w") as fh:
            fh.write(md)
        with open(os.path.join(args.out, "routing_tables.tex"), "w") as fh:
            fh.write(tex)
        print(f"wrote {args.out}/routing_tables.md and routing_tables.tex")
    else:
        print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
