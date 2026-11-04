#!/usr/bin/env python3
"""Verify every table in the paper-routing findings doc against its source CSVs.

Parses the markdown tables out of docs/architecture/paper_routing_figures_20260915.md
and recomputes each cell from the aggregated benchmark CSVs, so the doc cannot
drift from the measurements it cites. Exits non-zero on any mismatch.

Usage: python3 scripts/check_paper_routing_doc.py [DOC]
"""
from __future__ import annotations
import csv, re, sys

DOC = sys.argv[1] if len(sys.argv) > 1 else "docs/architecture/paper_routing_figures_20260915.md"
STATIC = ["outer_threads", "outer_process", "inner_only", "outer_inner_static"]

SOURCES = {
    "warm": "output/performance/paper_benchmarks/20260915_181642/paper_benchmarks_aggregated_20260915_181642.csv",
    "cold": "output/performance/paper_benchmarks/20260915_181642_cold_store/cold_store_aggregated.csv",
    "trx":  "output/performance/paper_benchmarks_trx50/20260916_143516/paper_benchmarks_aggregated_20260916_143516.csv",
}

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None

data = {k: list(csv.DictReader(open(v))) for k, v in SOURCES.items()}
fails: list[str] = []

def median(src, phase, axis_col, axis_val, mode, case=None):
    for r in data[src]:
        if (r["phase_id"] == phase and r[axis_col] == str(axis_val) and r["mode"] == mode
                and (case is None or r["case"] == case)):
            return num(r["wall_time_median_s"])
    return None

def serial_of(src, phase, axis_col, axis_val, case=None):
    """Serial at this point: its own row if measured there, else the serial_median_s column."""
    direct = median(src, phase, axis_col, axis_val, "serial", case)
    if direct is not None:
        return direct
    for r in data[src]:
        if (r["phase_id"] == phase and r[axis_col] == str(axis_val)
                and (case is None or r["case"] == case)):
            return num(r["serial_median_s"])
    return None

def best_static(src, phase, axis_col, axis_val, case=None):
    c = [(num(r["wall_time_median_s"]), r["mode"]) for r in data[src]
         if r["phase_id"] == phase and r[axis_col] == str(axis_val)
         and r["mode"] in STATIC and (case is None or r["case"] == case)]
    return min(c) if c else (None, None)

def chk(label, claimed, actual, tol=6e-4):
    if actual is None:
        fails.append(f"{label}: doc says {claimed}, no data"); return
    if abs(claimed - actual) > tol * max(1.0, abs(actual)):
        fails.append(f"{label}: doc says {claimed}, data says {actual:.4f}")

def rows_of(table_text):
    out = []
    for line in table_text.strip().splitlines():
        if not line.startswith("|"): continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if all(set(c) <= set("-: ") for c in cells): continue
        out.append(cells)
    return out[1:]  # drop header

doc = open(DOC).read()

def table_after(heading, index=0):
    i = doc.index(heading)
    chunk = doc[i:]
    tables = re.findall(r"(?:^\|.*\n)+", chunk, re.M)
    return rows_of(tables[index])

# --- P1 / P2 / TRX50 P1: serial | best static | route | R6 | serial/static | serial/R6 ---
SPECS = [
    ("### P1 — constellation size, 12 threads", 0, "warm", "P1", "n_sat",        None, 1),
    ("### P2 — thread budget at 4096 spacecraft", 0, "warm", "P2", "thread_count", None, 0),
    ("P1, constellation size at 32 threads:",   0, "trx",  "P1", "n_sat",        None, 0),
]
for heading, ti, src, phase, axis, case, mission_col in SPECS:
    for cells in table_after(heading, ti):
        axis_val = int(cells[0])
        off = 1 + mission_col          # skip the mission column where present
        ser, bs, route, r6 = float(cells[off]), float(cells[off+1]), cells[off+2], float(cells[off+3])
        rs, rr = float(cells[off+4]), float(cells[off+5])
        a_ser = serial_of(src, phase, axis, axis_val, case)
        a_bs, a_rt = best_static(src, phase, axis, axis_val, case)
        a_r6 = median(src, phase, axis, axis_val, "policy_v2", case)
        tag = f"{src}/{phase} {axis}={axis_val}"
        chk(f"{tag} serial", ser, a_ser); chk(f"{tag} best static", bs, a_bs); chk(f"{tag} R6", r6, a_r6)
        if route != a_rt: fails.append(f"{tag} route: doc {route}, data {a_rt}")
        if a_ser and a_bs: chk(f"{tag} serial/static", rs, round(a_ser/a_bs, 2), 6e-3)
        if a_ser and a_r6: chk(f"{tag} serial/R6", rr, round(a_ser/a_r6, 2), 6e-3)

# --- P3/P4 warm: budget | serial | best static | route | R6 | serial/R6 | R6/best static ---
for heading, phase, case in [
    ("Warm, `independent_1sat_1hr`, 256 samples:", "P3", "independent_1sat_1hr"),
    ("Warm, `montecarlo_heavy_aerobraking`, 32 samples:", "P4", "montecarlo_heavy_aerobraking"),
]:
    for cells in table_after(heading):
        b = int(cells[0]); ser, bs, route, r6 = float(cells[1]), float(cells[2]), cells[3], float(cells[4])
        sr6, r6bs = float(cells[5]), float(cells[6])
        a_ser = serial_of("warm", phase, "process_workers", b, case)
        a_bs, a_rt = best_static("warm", phase, "process_workers", b, case)
        a_r6 = median("warm", phase, "process_workers", b, "policy_v2", case)
        tag = f"warm/{phase} budget={b}"
        chk(f"{tag} serial", ser, a_ser); chk(f"{tag} best static", bs, a_bs); chk(f"{tag} R6", r6, a_r6)
        if route != a_rt: fails.append(f"{tag} route: doc {route}, data {a_rt}")
        if a_ser and a_r6: chk(f"{tag} serial/R6", sr6, round(a_ser/a_r6, 2), 6e-3)
        if a_bs and a_r6: chk(f"{tag} R6/best static", r6bs, round(a_r6/a_bs, 2), 6e-3)

# --- P3/P4 cold: budget | serial | outer_threads | outer_process | R6 | R6/best static ---
for heading, phase, case in [
    ("`independent_1sat_1hr`, 256 samples:", "P3", "independent_1sat_1hr"),
    ("`montecarlo_heavy_aerobraking`, 32 samples:", "P4", "montecarlo_heavy_aerobraking"),
]:
    i = doc.index("Cold-store versions of the same two tables (superseded):")
    j = doc.index(heading, i)
    tables = re.findall(r"(?:^\|.*\n)+", doc[j:], re.M)
    for cells in rows_of(tables[0]):
        b = int(cells[0]); ser, ot, op, r6, r6bs = (float(cells[1]), float(cells[2]),
                                                    float(cells[3]), float(cells[4]), float(cells[5]))
        tag = f"cold/{phase} budget={b}"
        chk(f"{tag} serial", ser, serial_of("cold", phase, "process_workers", b, case))
        chk(f"{tag} outer_threads", ot, median("cold", phase, "process_workers", b, "outer_threads", case))
        chk(f"{tag} outer_process", op, median("cold", phase, "process_workers", b, "outer_process", case))
        a_r6 = median("cold", phase, "process_workers", b, "policy_v2", case)
        chk(f"{tag} R6", r6, a_r6)
        a_bs, _ = best_static("cold", phase, "process_workers", b, case)
        if a_bs and a_r6: chk(f"{tag} R6/best static", r6bs, round(a_r6/a_bs, 2), 6e-3)

# --- P5: split | 16x8 best static | 16x8 R6 | ratio | 8x16 best static | 8x16 R6 | ratio ---
for cells in table_after("### P5 — worker/thread split of a fixed budget of 12"):
    w = int(cells[0].split("x")[0])
    for case, o, lbl in [("mcgrid_16sat_8mc", 1, "16x8"), ("mcgrid_8sat_16mc", 4, "8x16")]:
        bs_d, r6_d, rt_d = float(cells[o]), float(cells[o+1]), float(cells[o+2])
        a_bs, _ = best_static("warm", "P5", "process_workers", w, case)
        a_r6 = median("warm", "P5", "process_workers", w, "policy_v2", case)
        tag = f"warm/P5 {lbl} {w}x{12//w}"
        chk(f"{tag} best static", bs_d, a_bs); chk(f"{tag} R6", r6_d, a_r6)
        if a_bs and a_r6: chk(f"{tag} ratio", rt_d, round(a_r6/a_bs, 2), 6e-3)

print(f"checked {len(SOURCES)} sources against {DOC}")
if fails:
    print(f"\n!! {len(fails)} MISMATCH(ES)")
    for x in fails: print("  -", x)
    sys.exit(1)
print("every table cell in the doc matches its source CSV")
