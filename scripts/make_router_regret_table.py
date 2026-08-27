#!/usr/bin/env python3
r"""
Emit the consolidated router-regret table from paper-benchmark raw CSVs.

One row per (workload, thread count): the best static profile, its wall time,
and how far R4 and R5 land from it. Regret is measured against the best static
route actually observed rather than against a nominated baseline, because the
claim under test is that the adaptive profiles match or beat whichever fixed
route would have been the right choice -- which is only meaningful if the
comparison is to the best one, not to a convenient one.

Usage:
  python3 scripts/make_router_regret_table.py [--label=NAME] OUT.tex RAW1.csv [RAW2.csv ...]

--label overrides the LaTeX label (default tab:router_regret_consolidated). It
exists because this script now emits more than one table into the same document
-- a consolidated one and a Monte Carlo one -- and two label commands with the
same name make every \cref to either of them resolve to whichever came last,
silently. Pass a distinct label per table.
"""
import csv, sys, collections, statistics, os

STATIC = ["serial", "outer_threads", "outer_process", "inner_only", "outer_inner_static"]
PROFILE = {"serial": r"\texttt{R0}", "outer_threads": r"\texttt{R1\_a}",
           "outer_process": r"\texttt{R1\_b}", "inner_only": r"\texttt{R2}",
           "outer_inner_static": r"\texttt{R3}"}
# Readable workload names; anything unmapped falls back to the raw case name.
LABEL = {
    "gravity_4096sat_l50_vacuum_1hr": "Gravity only, 4096 spacecraft",
    "interact_64sat_1hr": "Interacting, 64 spacecraft",
    "interact_256sat_1hr": "Interacting, 256 spacecraft",
    "independent_1sat_1hr": "Independent, 256 samples",
    "stack256_e1_harm": "256 spacecraft, harmonics only",
    "stack256_e3_aero": "256 spacecraft, plus aerodynamics",
    "stack256_e4_nbody": "256 spacecraft, plus third body",
    "stack256_e5_6dof": "256 spacecraft, coupled 6-DOF",
    "stack32_e6_actuated": "32 spacecraft, actuated",
    "gravity_1024sat_l50_vacuum_15min": "1024 spacecraft, 15-minute arc",
    "gravity_1024sat_l50_vacuum_1hr": "1024 spacecraft, one-hour arc",
    "cadence_1024sat_none": "1024 spacecraft, no output",
    "cadence_1024sat_60s": "1024 spacecraft, 60\\,s output",
    "heavy_1024sat_l50_6hr": "1024 spacecraft, six-hour arc",
    "cadence_1024sat_10s": "1024 spacecraft, 10\\,s output",
    "cadence_1024sat_1s": "1024 spacecraft, 1\\,s output",
    "montecarlo_high_accuracy": "Monte Carlo, high-accuracy gravity",
    "montecarlo_multi_sat": "Monte Carlo, four spacecraft",
    "montecarlo_mars_aerobraking": "Monte Carlo, Mars aerobraking",
    "montecarlo_heavy_aerobraking": "Monte Carlo, long-arc aerobraking",
    "independent_1sat_1hr": "Independent samples, one spacecraft",
    "multi_128_high_fidelity": "128 spacecraft, full force stack",
}

def load(paths):
    runs = collections.defaultdict(list)
    for path in paths:
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if r.get("success", "true").lower() != "true":
                    continue
                try:
                    w = float(r["wall_time_s"])
                except (KeyError, ValueError):
                    continue
                # mc_samples is part of the key, not a field to average over.
                # A Monte Carlo case run at 16 and 64 samples is two different
                # workloads; keying only on (case, threads, mode) silently took
                # the median ACROSS the sample ladder and reported one number
                # for both. Non-Monte-Carlo rows all carry mc_samples=1, so this
                # is a no-op for them.
                try:
                    mc = int(r.get("mc_samples", 1) or 1)
                except ValueError:
                    mc = 1
                runs[(r["case"], int(r["thread_count"]), mc, r["mode"])].append(w)
    # Median across repeats: the harness already discards warm-up solves, and the
    # median is insensitive to a single contended repeat in a way the mean is not.
    return {k: statistics.median(v) for k, v in runs.items()}

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    args = sys.argv[1:]
    label = "tab:router_regret_consolidated"
    rest = []
    for a in args:
        if a.startswith("--label="):
            label = a.split("=", 1)[1]
        else:
            rest.append(a)
    if len(rest) < 2:
        sys.exit(__doc__)
    out_path, raw_paths = rest[0], rest[1:]
    agg = load(raw_paths)
    cases = sorted({k[0] for k in agg}, key=lambda c: list(LABEL).index(c) if c in LABEL else 99)

    rows, wins, ties, losses = [], 0, 0, 0
    best_counts = collections.Counter()
    fixed_penalty = collections.defaultdict(list)
    for case in cases:
        for t in sorted({k[1] for k in agg if k[0] == case}):
          for mc in sorted({k[2] for k in agg if k[0] == case and k[1] == t}):
            modes = {k[3]: v for k, v in agg.items()
                     if k[0] == case and k[1] == t and k[2] == mc}
            static = {m: v for m, v in modes.items() if m in STATIC}
            r4, r5 = modes.get("outer_inner_adaptive"), modes.get("full_smart")
            if not static or r5 is None:
                continue
            bm = min(static, key=static.get)
            bv = static[bm]
            best_counts[PROFILE[bm]] += 1
            for m, v in static.items():
                fixed_penalty[PROFILE[m]].append(100.0 * (v / bv - 1.0))
            g4 = 100.0 * (r4 / bv - 1.0) if r4 else None
            g5 = 100.0 * (r5 / bv - 1.0)
            if g5 < -2.0:   wins += 1
            elif g5 <= 2.0: ties += 1
            else:           losses += 1
            best5 = r"\textbf{%.3f}" % r5 if r5 <= bv else "%.3f" % r5
            rows.append((LABEL.get(case, case.replace("_", r"\_")), t, mc, PROFILE[bm], bv,
                         ("%.3f" % r4) if r4 else "--", best5,
                         ("$%+.0f\\%%$" % g4) if g4 is not None else "--",
                         "$%+.0f\\%%$" % g5))

    with open(out_path, "w") as fh:
        fh.write("\\begin{table}[htbp]\n\\centering\n")
        fh.write("\\caption{Adaptive routing against an ORACLE fixed route, across workloads. "
                 "Regret is relative to the fastest static profile measured for that workload and "
                 "thread count, selected after the fact; negative values indicate the adaptive "
                 "profile was faster than any fixed route. This baseline is deliberately "
                 "unattainable: the fastest fixed route changes with the workload and the thread "
                 "budget, so no strategy a practitioner could commit to in advance achieves the "
                 "Time column. Read this table for how close adaptive routing comes to perfect "
                 "hindsight, and the accompanying aggregate table for what it is worth against a "
                 "choice that could actually be made.}\n")
        fh.write("\\label{%s}\n\\small\n" % label)
        fh.write("\\setlength{\\tabcolsep}{4pt}\n\\renewcommand{\\arraystretch}{0.92}\n")
        fh.write("\\begin{tabular}{lrrlrrrrr}\n\\toprule\n")
        fh.write("& & & \\multicolumn{2}{c}{Best static route} & "
                 "\\multicolumn{2}{c}{Adaptive (s)} & \\multicolumn{2}{c}{Regret} \\\\\n")
        fh.write("\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n")
        fh.write("Workload & Threads & Samples & Profile & Time (s) & \\texttt{R4} & \\texttt{R5} & "
                 "\\texttt{R4} & \\texttt{R5} \\\\\n\\midrule\n")
        prev = None
        for lbl, t, mc, prof, bv, r4, r5, g4, g5 in rows:
            shown = lbl if lbl != prev else ""
            prev = lbl
            fh.write("%s & %d & %d & %s & %.3f & %s & %s & %s & %s \\\\\n"
                     % (shown, t, mc, prof, bv, r4, r5, g4, g5))
        fh.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    n = wins + ties + losses
    print("wrote %s: %d rows" % (out_path, len(rows)))
    print("R5 vs best static -- faster: %d, within 2%%: %d, slower: %d (of %d)" % (wins, ties, losses, n))
    # The cost of committing to one fixed route. This is the claim the table
    # actually supports: not that the adaptive router wins everywhere, but that
    # no fixed route does, and that guessing wrong is expensive.
    print("best-static profile by workload: %s" % dict(best_counts))
    for prof, pen in sorted(fixed_penalty.items(), key=lambda kv: -max(kv[1])):
        print("  committing to %-16s worst-case regret %+.0f%% over %d workloads"
              % (prof, max(pen), len(pen)))


if __name__ == "__main__":
    main()
