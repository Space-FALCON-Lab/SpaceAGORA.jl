# Paper routing figures — R6 against serial and the best static route

Record of the P1–P5 benchmark phases built for the parallelization paper's four
routing comparisons, the results measured on `space-falcon-1`, and the findings
that need stating alongside them.

- **Code:** `origin/main` 4b42b544, run from a dedicated worktree so the numbers
  cite one hash; harness changes on `bench/paper-routing-figures`.
- **Machine:** `space-falcon-1`, Ryzen 9 9900X, 12 physical cores / 24 threads,
  60 GB. Thread ladder `1,2,4,8,12`, process-worker cap 12.
- **Run:** `output/performance/paper_benchmarks/20260915_181642` (gitignored).
  Cold-store P3/P4 preserved in the `_cold_store` sibling directory.
- **Tables:** regenerate with
  `python3 scripts/make_paper_routing_tables.py <run dir> --out <dir>`;
  markdown for reading, booktabs LaTeX for the manuscript.

## What the four comparisons are

| Phase | Question | Axis |
|---|---|---|
| P1 | Does the router track the best route as a constellation grows? | 1 → 4096 spacecraft at a fixed budget |
| P2 | Does it track it as the budget grows? | 1 → 12 threads at 4096 spacecraft |
| P3 | Does it scale a Monte Carlo campaign of single spacecraft? | budget 1 → 12, 256 samples |
| P4 | The same, on compute-bound samples | budget 1 → 12, 32 aerobraking samples |
| P5 | Does it split one budget between samples and satellites? | six worker×thread splits of 12 |

Each phase measures `serial`, every applicable pinned static route, and
`policy_v2` (R6) at the same point. "Best static" in every table below is the
fastest *parallel* pinned route at that point; serial has its own column.

## Results

### P1 — constellation size, 12 threads

| N | mission | serial | best static | route | R6 | serial/static | serial/R6 |
|---:|---:|---:|---:|---|---:|---:|---:|
| 1 | 1153 h | 8.309 | 8.429 | inner_only | 8.630 | 0.99 | 0.96 |
| 16 | 143 h | 10.943 | 2.561 | outer_inner_static | 2.672 | 4.27 | 4.10 |
| 64 | 115 h | 11.630 | 14.081 | outer_inner_static | 5.725 | 0.83 | 2.03 |
| 256 | 34 h | 11.670 | 7.146 | outer_threads | 7.379 | 1.63 | 1.58 |
| 1024 | 6.8 h | 11.890 | 3.219 | outer_threads | 3.187 | 3.69 | 3.73 |
| 4096 | 1.6 h | 11.997 | 2.626 | outer_inner_static | 2.649 | 4.57 | 4.53 |

### P2 — thread budget at 4096 spacecraft

| threads | serial | best static | route | R6 | serial/static | serial/R6 |
|---:|---:|---:|---|---:|---:|---:|
| 1 | 12.117 | 11.778 | outer_threads | 12.052 | 1.03 | 1.01 |
| 2 | 12.117 | 6.649 | outer_threads | 6.624 | 1.82 | 1.83 |
| 4 | 12.117 | 3.836 | inner_only | 3.922 | 3.16 | 3.09 |
| 8 | 12.117 | 2.721 | inner_only | 2.878 | 4.45 | 4.21 |
| 12 | 12.117 | 2.663 | outer_inner_static | 2.652 | 4.55 | 4.57 |

### P3 / P4 — Monte Carlo resource ladders

Both ladders were measured twice: once from the fresh worktree's empty
calibration store, and once after it had converged. The warm columns are the
ones to quote; the cold ones are in the `_cold_store` directory and in finding 6.

Warm, `independent_1sat_1hr`, 256 samples:

| budget | serial | best static | route | R6 | serial/R6 | R6/best static |
|---:|---:|---:|---|---:|---:|---:|
| 1 | 8.770 | 8.964 | outer_threads | 9.046 | 0.97 | 1.01 |
| 2 | 9.188 | 4.645 | outer_threads | 3.460 | 2.66 | 0.75 |
| 4 | 9.006 | 2.589 | outer_threads | 1.760 | 5.12 | 0.68 |
| 8 | 9.274 | 1.555 | outer_threads | 1.470 | 6.31 | 0.95 |
| 12 | 9.368 | 1.414 | outer_process | 1.531 | 6.12 | 1.08 |

Warm, `montecarlo_heavy_aerobraking`, 32 samples:

| budget | serial | best static | route | R6 | serial/R6 | R6/best static |
|---:|---:|---:|---|---:|---:|---:|
| 1 | 23.581 | 24.255 | outer_process | 22.712 | 1.04 | 0.94 |
| 2 | 22.206 | 12.085 | outer_threads | 8.211 | 2.70 | 0.68 |
| 4 | 22.123 | 6.358 | outer_process | 4.138 | 5.35 | 0.65 |
| 8 | 22.259 | 3.592 | outer_process | 2.959 | 7.52 | 0.82 |
| 12 | 22.602 | 2.882 | outer_process | 3.486 | 6.48 | 1.21 |

Cold-store versions of the same two tables (superseded):

`independent_1sat_1hr`, 256 samples:

| budget | serial | outer_threads | outer_process | R6 | R6/best static |
|---:|---:|---:|---:|---:|---:|
| 1 | 9.081 | 9.130 | 9.208 | 9.051 | 0.99 |
| 2 | 9.261 | 4.787 | 4.945 | 3.446 | 0.72 |
| 4 | 9.306 | 2.639 | 2.786 | 1.801 | 0.68 |
| 8 | 9.185 | 1.628 | 1.754 | 2.453 | 1.51 |
| 12 | 9.600 | 1.524 | 1.491 | 1.567 | 1.05 |

`montecarlo_heavy_aerobraking`, 32 samples:

| budget | serial | outer_threads | outer_process | R6 | R6/best static |
|---:|---:|---:|---:|---:|---:|
| 1 | 23.283 | 24.502 | 23.635 | 23.484 | 0.99 |
| 2 | 22.430 | 12.135 | 12.878 | 8.312 | 0.68 |
| 4 | 22.858 | 6.818 | 6.288 | 4.215 | 0.67 |
| 8 | 22.776 | 4.284 | 3.624 | 3.805 | 1.05 |
| 12 | 23.982 | 3.378 | 2.984 | 3.538 | 1.19 |

### P5 — worker/thread split of a fixed budget of 12

`mcgrid_16sat_8mc` (serial ≈ 10.6 s) and `mcgrid_8sat_16mc` (serial ≈ 9.2 s):

| split | 16×8 best static | 16×8 R6 | ratio | 8×16 best static | 8×16 R6 | ratio |
|---|---:|---:|---:|---:|---:|---:|
| 1×12 | 1.594 | 1.546 | 0.97 | 1.491 | 1.360 | 0.91 |
| 2×6 | 2.758 | 2.730 | 0.99 | 1.906 | 1.920 | 1.01 |
| 3×4 | 2.795 | 2.760 | 0.99 | 2.464 | 2.012 | 0.82 |
| 4×3 | 3.047 | 2.792 | 0.92 | 2.650 | 2.010 | 0.76 |
| 6×2 | 2.760 | 1.835 | 0.66 | 1.924 | 1.988 | 1.03 |
| 12×1 | 1.492 | 1.395 | 0.94 | 1.370 | 1.289 | 0.94 |

## Findings

**1. The best static route changes identity along every axis.** It moves three
times down P1's column, twice down P2's, and three times across each of P5's.
No single pinned route is the right answer at every size, budget or split, which
is the case for routing at all — and it is measured here rather than asserted.

**2. At 64 spacecraft every static route is slower than serial.** 14.08–14.25 s
against serial's 11.63 s, while R6 reaches 5.73 s. The raw telemetry attributes
it: the static routes run with `rhs_plan_source=none` and let the auto heuristic
choose the batched-RHS configuration, which is a net loss at that width, while
R6 runs `plan=cache/satellite_batch, allotment=1, static`. R6's win here is
calibration, not route selection, and the paper should say so.

**3. R6 beats every pure route at small Monte Carlo budgets by being impure.**
The dispatch trace at budget 8 reads `workers=8 local_slots=4 pool_size=8`: it
dispatches to eight worker processes *and* runs four samples on the coordinator.
Neither `outer_threads` nor `outer_process` can express that, which is why R6 is
~1.5x faster than both at budgets 2 and 4 on both Monte Carlo workloads.

**4a. The cold-store result at budget 8 was an artifact, and is withdrawn.**
Measured warm, R6 runs budget 8 in 1.470 s against the best static route's
1.555 s on P3 (cold: 2.453 against 1.628) and 2.959 s against 3.592 s on P4
(cold: 3.805 against 3.624). The controls move 1-5% between the two runs, which
is machine noise; R6 moves 40% and 22% at that one rung and is unchanged
elsewhere. The earlier reading -- that the bandit had discarded a measured-better
arm -- was a converging calibration store and nothing more.

**4b. At the full budget R6 pays for exploration, and a five-repeat median
charges it.** The per-repeat routes at budget 12 are identical on both
workloads: process, process, *threads*, *threads*, process. Its exploiting
campaigns are the fastest measurements at that point on either workload --
1.057 s against the best static route's 1.414 s on P3, 2.524 s against 2.832 s
on P4 -- but repeats 3 and 4 cost 2.5 s and 4.1 s, and the median lands on one
of them. At budget 8 there is no flapping and R6 wins outright. Quote both: the
median over a short campaign sequence is what a user with five campaigns to run
experiences, and the steady state is what the routing itself achieves.

**5. R6's variance is far higher than the static routes'.** At P4's budget 12 its
repeats span 2.53–3.97 s while `outer_process` holds 2.96–3.07 s. Its median
over five repeats is correspondingly unstable, and the apparent regressions at
budgets 8 and 12 did not reproduce when those points were measured in isolation
(R6 1.48 s vs threads 1.68 s at P3 budget 8; 2.96 s vs process 3.05 s at P4
budget 12). Any single-median comparison of an adaptive route against a pinned
one should be read with that spread in view.

**6. Cold and warm calibration stores are different measurements.** This run
started from a fresh worktree, so the RHS-calibration and inner-policy stores
began empty and were still converging while P3 and P4 ran; every prior B- and
L-series number used the repository's accumulated store. A warm re-run of P3/P4
is recorded in the "Calibration-store warmth" section of the generated tables.

## Methodology notes

- **The measurability floor drove the workload design.** At one fixed mission
  length five of P1's six rungs sit under the harness's 3 s floor, so the size
  ladder holds the *work* fixed instead and moves the mission length per rung
  (`PPC_L50_ISO_MISSION_S`). P3 clears the floor with 256 samples rather than 64.
- **The quiet guard is load-bearing.** Three of P4's five rungs were rejected
  outright on the first attempt — "Machine has 57.0% of its CPU already in use"
  — because a previous phase's process pool was still winding down and a global
  OOM had just kicked the desktop. Those points were re-run on an idle machine
  rather than kept.
- **Twelve process workers can exhaust this box.** A 12-worker aerobraking rung
  plus a browser and an editor drove a global OOM that killed a desktop process.
  The benchmark survived it; the desktop did not.

## The S1 speedups rest on a different serial baseline

`paper_scenarios`' S1 reports up to 21.6x at L50/4096 on twelve threads where P2
measures 4.57x. Both are arithmetically correct; they divide by different
denominators, and the difference is 4.5x of baseline, not of parallel
performance.

Two differences were measured, and only one of them matters:

- **Integrator step ceiling, which cancels.** S1 builds its workload with
  `dt_max_orbit=2.0` against the ppb catalog's `20.0` at identical tolerances,
  so S1 does ten times the RHS evaluations per simulated second. That scales
  both sides of a ratio equally and drops out of it.
- **The RHS execution mode of the baseline, which does not.** S1 pins
  `rhs_mode=serial` for its serial rows and `flat` for its parallel rows; ppb's
  `serial` mode turns batched RHS off under profile R0 but leaves the execution
  plan on `auto`, which still takes the batched kernel on one thread. Measured
  on `gravity_4096sat_l50_vacuum_5800s`, one thread, three repeats: `auto`
  11.84-12.13 s against `serial` 53.90-53.93 s -- **4.5x**.

The two reconcile completely. S1's 160.77 s serial baseline, divided by ten for
the step ceiling and scaled to P2's mission, predicts 51.8 s; ppb measures
53.9 s with the serial RHS pinned. Re-basing S1's 4096 point on the `auto`
baseline gives 35.7 / 7.45 = 4.8x, against P2's 4.57x.

**Consequence for the paper.** Every S1 speedup -- the 21.6x, and the
10.23x/11.80x/10.17x of the CPU-only scaling work -- divides by a serial run
4.5x slower than the fastest single-threaded code the library runs on its own
defaults. Both definitions are defensible ("all parallelism off" against "best
serial implementation"), but a headline speedup invites the question of whether
the baseline was optimised, and only one of the two answers it well. Pick one,
state it, and do not print numbers from both conventions in the same table.

## Open

- **TRX50.** The same five phases at budget 32 are still outstanding; the box has
  been occupied by another user's job throughout.
