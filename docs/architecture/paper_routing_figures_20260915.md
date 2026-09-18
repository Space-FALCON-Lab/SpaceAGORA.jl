# Paper routing figures — R6 against serial and the best static route

Record of the P1–P5 benchmark phases built for the parallelization paper's four
routing comparisons, the results measured on `space-falcon-1` (12 cores) and
`space-falcon-lab-TRX50-AERO-D` (64 cores), and the findings that need stating
alongside them. Both machines ran the same code at the same commit.

- **Code the numbers were measured against:** `origin/main` 4b42b544, run from a
  dedicated worktree so every figure in this document cites one hash. Harness
  changes on `parallelization-paper-results`.
- **Base the branch now sits on:** `origin/main` 31e04bc8, which moved 21 `src/`
  files after the numbers were taken. P1 and P5 were re-measured against it,
  with the calibration store controlled, and the difference is smaller than
  run-to-run noise -- see "Re-measured against the branch's own base" below.
  4b42b544 remains the hash to quote, since it is what P2/P3/P4 were measured
  on and what the tables report.
- **Machines:** `space-falcon-1`, Ryzen 9 9900X, 12 physical cores / 24 threads,
  60 GB, thread ladder `1,2,4,8,12`, process-worker cap 12; and
  `space-falcon-lab-TRX50-AERO-D`, 64 physical cores, thread ladder
  `1,2,4,8,16,32`, process-worker cap 32.
- **Run:** `output/performance/paper_benchmarks/20260915_181642` (gitignored).
  Cold-store P3/P4 preserved in the `_cold_store` sibling directory.
- **Tables:** regenerate with
  `python3 scripts/make_paper_routing_tables.py <run dir> --out <dir>`;
  markdown for reading, booktabs LaTeX for the manuscript.
- **Figures:** regenerate with
  `python3 scripts/make_paper_routing_plots.py <run dir> [<run dir> ...] --out <dir>`;
  PNG for reading, PDF for the manuscript. See "Figures" below.

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
| 1 | 1153 h | 8.463 | 8.207 | outer_inner_static | 8.294 | 1.03 | 1.02 |
| 16 | 143 h | 10.985 | 2.500 | outer_threads | 2.514 | 4.39 | 4.37 |
| 64 | 115 h | 11.287 | 13.943 | outer_threads | 5.710 | 0.81 | 1.98 |
| 256 | 34 h | 11.812 | 7.117 | outer_inner_static | 7.082 | 1.66 | 1.67 |
| 1024 | 6.8 h | 11.773 | 3.149 | inner_only | 3.207 | 3.74 | 3.67 |
| 4096 | 1.6 h | 11.651 | 2.630 | inner_only | 2.692 | 4.43 | 4.33 |

### P2 — thread budget at 4096 spacecraft

| threads | serial | best static | route | R6 | serial/static | serial/R6 |
|---:|---:|---:|---|---:|---:|---:|
| 1 | 11.613 | 11.670 | inner_only | 11.766 | 1.00 | 0.99 |
| 2 | 11.613 | 6.616 | inner_only | 6.683 | 1.76 | 1.74 |
| 4 | 11.613 | 3.780 | outer_inner_static | 4.151 | 3.07 | 2.80 |
| 8 | 11.613 | 2.730 | inner_only | 2.616 | 4.25 | 4.44 |
| 12 | 11.613 | 2.606 | inner_only | 2.679 | 4.46 | 4.34 |

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
| 1 | 23.581 | 24.255 | outer_threads | 22.712 | 1.04 | 0.94 |
| 2 | 22.206 | 12.085 | outer_process | 8.211 | 2.70 | 0.68 |
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

`mcgrid_16sat_8mc` (serial 10.12–10.88 s across the six splits) and
`mcgrid_8sat_16mc` (serial 9.06–9.59 s):

| split | 16x8 best static | 16x8 R6 | ratio | 8x16 best static | 8x16 R6 | ratio |
|---|---:|---:|---:|---:|---:|---:|
| 1x12 | 1.608 | 1.563 | 0.97 | 1.472 | 1.372 | 0.93 |
| 2x6 | 2.734 | 2.697 | 0.99 | 1.928 | 1.929 | 1.00 |
| 3x4 | 2.781 | 2.756 | 0.99 | 2.480 | 2.065 | 0.83 |
| 4x3 | 2.981 | 2.787 | 0.93 | 2.788 | 2.003 | 0.72 |
| 6x2 | 2.831 | 1.740 | 0.61 | 1.976 | 1.979 | 1.00 |
| 12x1 | 1.465 | 1.411 | 0.96 | 1.357 | 1.254 | 0.92 |

## TRX50 — the same five phases at budget 32

- **Machine:** `space-falcon-lab-TRX50-AERO-D`, 64 physical cores, thread ladder
  `1,2,4,8,16,32`, process-worker cap 32. Idle and quiet-gated throughout.
- **Run:** `20260916_143516`, 163 aggregated rows, none below the 3 s floor.
  Pulled to `output/performance/paper_benchmarks_trx50/` (gitignored).
- **Same code, same cases, same iso-work mission lengths.** The constellation
  phases' serial baselines land at 8.6–12.2 s, matching `space-falcon-1`'s
  8.5–11.8 s, so the two machines' columns are directly comparable. (The Monte
  Carlo phases carry their own baselines: 9.2–10.4 s on P3, 23.6–26.6 s on P4.)

P1, constellation size at 32 threads:

| N | serial | best static | route | R6 | serial/static | serial/R6 |
|---:|---:|---:|---|---:|---:|---:|
| 1 | 8.607 | 8.567 | inner_only | 8.701 | 1.00 | 0.99 |
| 16 | 10.635 | 2.199 | outer_threads | 2.221 | 4.84 | 4.79 |
| 64 | 11.389 | 34.457 | outer_inner_static | 3.740 | 0.33 | 3.05 |
| 256 | 11.554 | 20.907 | inner_only | 2.833 | 0.55 | 4.08 |
| 1024 | 12.195 | 4.963 | outer_inner_static | 5.053 | 2.46 | 2.41 |
| 4096 | 11.532 | 2.213 | outer_threads | 2.457 | 5.21 | 4.69 |

P2, thread budget at 4096 spacecraft, reaches 5.83x at 16 threads and falls
back to 5.44x at 32 — the ceiling is between 8 and 16 threads on this box, and
R6 tracks the best static route to within 4% at every rung.

P3/P4 are the clearest scaling result in the set: on the cheap ladder R6 runs
9.70x at budget 32 against serial and matches the best static route exactly;
on the heavy ladder the best static route reaches 16.51x while R6 reaches
12.71x. P5's two workloads peak at 7.48x and 12.58x, both at the 32x1 split,
with R6 the fastest measured route at that point on both.

Full tables for both machines: `output/paper_routing_tables/routing_tables.md`
and `.tex`, regenerated with

```
python3 scripts/make_paper_routing_tables.py \
    output/performance/paper_benchmarks/20260915_181642 \
    output/performance/paper_benchmarks_trx50/20260916_143516 \
    --cold output/performance/paper_benchmarks/20260915_181642_cold_store/cold_store_aggregated.csv \
    --out output/paper_routing_tables
```

## Figures

`output/paper_routing_plots/` (gitignored), one PDF and one PNG each. Every
phase figure has machines as rows and, as columns, the two halves of the
comparison: raw median wall time on the left, ratio to that point's serial
baseline on the right. Each pinned static route is drawn faintly behind the
best-static line, because at several points their near-coincidence is itself
the result.

| File | Shows |
|---|---|
| `fig_summary_regret` | R6 against the best pinned route at all 69 measured points, by phase |
| `fig_regret_distribution` | the same regret as a distribution, drawn against the measured noise floor |
| `fig_repeat_variance` | how reproducible one measured point is, R6 against pinned routes and serial |
| `fig_p1` | Constellation size, 1 -> 4096 spacecraft, both machines |
| `fig_p2_gravity_4096sat_l50_vacuum_5800s` | Thread budget at 4096 spacecraft |
| `fig_p3_independent_1sat_1hr` | Monte Carlo resource ladder, cheap samples |
| `fig_p4_montecarlo_heavy_aerobraking` | Monte Carlo resource ladder, compute-bound samples |
| `fig_p5_mcgrid_16sat_8mc`, `fig_p5_mcgrid_8sat_16mc` | Worker/thread split of one budget |

Three things the figures carry that the tables do not. P1's raw panel shows the
static routes rising *above* the serial line at N=64 and N=256 while R6 stays
below it -- finding 2 as a picture. P5's faint lines cross: `outer_threads`
falls and `outer_process` rises across the splits, and the best-static and R6
lines ride the upper envelope of the two, which is the whole argument for
routing in one panel. And P2/P3/P4 carry a perfect-scaling reference, so the
distance between measured and ideal is legible rather than inferred.

The summary figure is the one to lead with, but it needs its caveat: R6's two
deepest P1 points are the batched-RHS heuristic failing on every pinned route,
not the router beating them, so the figure annotates them rather than banking
a 10x win. Read with those two excluded, R6's per-phase median ratio runs
0.82-1.02 over the ten machine-phase combinations: at or just above parity on
the constellation and thread axes (1.006-1.017), and below it on all four
Monte Carlo ladders (0.824-1.001), where the mixed dispatch of finding 3 is
something no pinned route can express.

`fig_regret_distribution` is the one to lead an adaptive-policy section with.
A single median per point cannot carry that claim -- `policy_v2` is the least
reproducible mode in the set -- so it plots the whole distribution of
R6 / best-static with the identical-code noise floor behind it, and states how
many points clear that floor in each direction: **20 of 66 faster than
measurement noise explains, 8 slower**. The three finding-2 calibration
artifacts are excluded and counted in the annotation rather than banked, since
they are the batched-RHS heuristic failing on every pinned route and not the
router out-routing them. Its right panel is a per-phase ECDF rather than more
histograms, because 9-24 points per phase is too few to bin honestly; it is
what shows that the wins are concentrated in the Monte Carlo phases while P1
and P2 sit on parity.

`fig_repeat_variance` is finding 5 as a distribution rather than an anecdote,
and is the honest counterweight to put beside any R6 win. Pooled over every run
on disk, with each point's warm-up repeat dropped, a single point's spread
separates the three route classes by roughly a decade, visible as three
separated ECDF curves. Dropping the warm-up repeat is not cosmetic: repeat 1
runs 2.0x the steady-state time for `policy_v2` and 1.4x for `outer_process`,
so including it would measure process-pool startup rather than routing.

The spread statistic is **IQR/median, not (max-min)/median**, and that matters
once runs with different repeat counts are pooled. The range of *k* draws grows
with *k*, so a max-min spread compares the repeat count as much as the
variance: the same P1 points read 3.3% for `policy_v2` at 3 repeats and 9.2% at
11, which is an artifact and not a regression. On IQR the same comparison is
3.3% against 3.2%, and on CV 1.6% against 2.7% -- i.e. unchanged, which is the
right answer for identical code on one machine. More repeats do not reduce a
point's spread; they make its *median* precise from the same spread, which is
what the repeat-count bootstrap in Open measures.

R6 does not converge over a five-campaign sequence -- repeats 2 through 5 sit
flat at 1.005-1.010 of each other -- so this spread is the steady-state cost of
adaptivity, not a transient that more campaigns would settle.

Regenerate the regret distribution by naming two runs of identical code with
`--noise-pair`, which is what supplies the floor:

```
python3 scripts/make_paper_routing_plots.py <run> <run> \
    --noise-pair output/performance/paper_benchmarks_trx50/20260916_143516 \
    --noise-pair output/performance/paper_benchmarks_trx50_control/20260917_141035 \
    --noise-pair output/performance/paper_benchmarks_trx50_p234_control/20260917_193014 \
    --out output/paper_routing_plots
```

Both scripts share one loader, so a figure and the table beside it cannot
disagree about which route won a point. `make_paper_routing_plots.py` also
reads each drawn line back off the axes and checks it against the records it
came from, failing the run rather than writing a figure that disagrees with
its own source data; `scripts/check_paper_routing_doc.py` does the same for
every table cell in this document.

## Re-measured against the branch's own base

The branch was rebased onto `origin/main` 31e04bc8, which had moved 21 `src/`
files since the numbers were taken -- `simulation/campaigns/constellation_ensemble.jl`
and `simulation/engine/setup.jl` among them, both on the path the constellation
phases exercise. All five phases were re-run on TRX50 to find out whether that
mattered.

A straight re-run cannot answer it. The original run's calibration store was
still converging while P1 ran, a fresh run starts warm, and finding 4a already
put store warmth at 40% on one point. So each phase group was measured twice,
back to back on an idle machine, from one snapshot of the store restored
byte-for-byte between them (md5 checked each time), with byte-identical harness
files and `src/` differing by exactly the 21 files. That gives two comparisons:

| comparison | `src/` | store + day | n | median | median abs dev | p90 | outside +/-8% |
|---|---|---|---:|---:|---:|---:|---:|
| **A: new base vs old base** | **differs** | same | 163 | 1.001 | 1.5% | 4.5% | **6** |
| **B: old base vs original run** | same | **differs** | 163 | 1.002 | 1.8% | 6.4% | **11** |

**The `src/` change is performance-neutral.** It is smaller than run-to-run
variance on every statistic: fewer points outside the band (6 against 11), a
smaller median deviation, and a smaller p90. Per phase, A puts 1 of 30 outside
on P1, 1 of 25 on P2, **0 of 24 on P3**, 2 of 24 on P4 and 2 of 60 on P5.

Half of A's six movers are `policy_v2`, all three at the top budget, and two
more are barely over the line (`outer_process` 0.909 and 0.912). The sixth is
`serial` at P5's 8x4 split (1.085) -- the unrouted baseline, where a move is
machine variance by construction. P2's `policy_v2` at 32 threads looks like the
largest effect in the table at 0.405, and is the clearest argument for running
the control: comparison B moves that same point by 2.489 with identical code.
The three measurements are 2.176 s, 5.418 s and 2.194 s; the control simply
caught R6 in a bad campaign sequence, and a two-run comparison would have
banked a 2.5x code win that does not exist.

Two things fell out of this that matter more than the rebase question.

**The N=4096 rung is where P1's static routes are least reproducible.** The two
static points that looked like a code regression on the naive comparison --
`outer_inner_static` +41% and `inner_only` -17% -- reproduce with *identical*
code, at 1.465 and 0.816 in comparison B.

**P1's and P2's best-static route labels are noise**, which is finding 1's
correction below. Three runs of P1, two of them the same binary, name a
different winning route at most rungs, because the routes sit within 0.1-1.9%
of each other. The dip is untouched by this: N=64's best static is 34.46 s,
34.74 s and 34.46 s across the three runs, and N=256's is 20.91 s, 21.25 s and
20.72 s. Finding 2 stands; finding 1 did not.

### What run-to-run noise actually looks like

Comparison B is 163 points of identical code re-measured on a quiet machine,
which is the first proper noise characterization this harness has. The 8%
band used by `compare_paper_routing_runs.py` sits near the p90; the tail is
much longer than that, and it is not evenly distributed.

| phase | n | median abs dev | p90 | max |
|---|---:|---:|---:|---:|
| P1 | 30 | 1.8% | 14.8% | 46.5% |
| P2 | 25 | 2.5% | 7.6% | 148.9% |
| P3 | 24 | 1.1% | 3.3% | 11.2% |
| P4 | 24 | 2.3% | 5.2% | 7.6% |
| P5 | 60 | 1.7% | 5.2% | 11.5% |

| mode | n | median abs dev | p90 | max |
|---|---:|---:|---:|---:|
| `policy_v2` | 36 | 2.2% | 10.2% | 148.9% |
| `outer_process` | 24 | 2.1% | 6.4% | 11.2% |
| `inner_only` | 12 | 2.0% | 3.8% | 18.4% |
| `serial` | 31 | 1.8% | 3.8% | 6.0% |
| `outer_threads` | 36 | 1.7% | 5.2% | 10.9% |
| `outer_inner_static` | 24 | 1.5% | 3.6% | 46.5% |

A typical point is reproducible to about 2%, which is why the sub-1% gaps
between P1's static routes cannot support an argmin. The tails belong to
`policy_v2` -- consistent with finding 5, and now quantified against a control
rather than inferred from one run's repeat spread.

Runs, all gitignored: `paper_benchmarks_trx50_rebased/20260917_112754` and
`paper_benchmarks_trx50_p234_new/20260917_171125` (new base);
`paper_benchmarks_trx50_control/20260917_141035` and
`paper_benchmarks_trx50_p234_control/20260917_193014` (old base). Compare any
two with `scripts/compare_paper_routing_runs.py`.

## Findings

**1. The best static route changes identity with the shape of the budget, and
only there. The earlier version of this finding, which also claimed it changes
down P1's and P2's columns, was reading noise.** What decides whether an argmin
is a result is the gap to the second-best route, and that gap differs by orders
of magnitude between the phases:

| phase | gap, best to second-best static | argmin reproducible across runs? |
|---|---|---|
| P1 | 0.1-1.9% | no |
| P2 | 0.4-9.1% (median ~2%) | no |
| P3 | 0.9-5.9% | marginal |
| P4 | 3.6-19.3%, growing with budget | yes above budget 4 |
| P5 | 0.1-13.7% on the thread splits; **47-1188%** on the process splits | yes on the process splits |

P5 is where the claim holds outright. `outer_process` wins the 8x4, 16x2 and
32x1 splits by 47% to 1188%, and three separate TRX50 runs -- two of them
running identical code -- name it every time. Threads win the narrow splits and
processes win the wide ones, by margins nothing in this harness could
manufacture. P4 supports a weaker version: `outer_threads` at budget 1,
`outer_process` from budget 2 on, with the margin growing to 19.3% at 12.

P1 and P2 do not support it. There the three in-process routes sit within a few
percent of each other at every point, which is below run-to-run variance, so
which one is "best" is decided by noise. Three TRX50 runs of P1 name a different
winner at most rungs, and the two 12-core runs disagree at all six. The route
column in P1's and P2's tables should be read as "any of the three", not as a
result; the *times* in those tables are reproducible, and the route labels are
not. Quoting them as a changing argmin would be claiming a signal the data does
not contain.

**2. At mid constellation sizes every static route is slower than serial, on
both machines, and the cause is the RHS plan rather than the route.** On
`space-falcon-1` at 64 spacecraft the three static routes take 13.94–13.99 s
against serial's 11.29 s while R6 reaches 5.71 s. On TRX50 the same dip is
deeper and spans two rungs: 34.46–34.89 s at N=64 and 20.91–21.31 s at N=256,
against serial's 11.39 s and 11.55 s, while R6 reaches 3.74 s and 2.83 s. (All
figures here are the warm run tabulated above; the first run's N=64 numbers are
in finding 7.) It is not noise — repeats of a given static route agree to within
2.5%, and the three routes agree with *each other* to within 1.3% at N=64 and
1.9% at N=256, which is the first clue: the choice of outer route is irrelevant
to the cost.

The raw telemetry attributes it. The static routes run
`rhs_plan_source=none, rhs_plan_mode=none, rhs_batch_parallel=auto` and let the
auto heuristic pick the batched-RHS configuration; R6's fast repeats run
`rhs_plan_source=sweep, rhs_plan_mode=satellite_batch, allotment=1,
scheduler=static`. TRX50 supplies the control the local run could not: R6's
second repeat at N=64 fell back to `rhs_plan_source=cache,
rhs_plan_mode=heuristic, allotment=0` and took 34.46 s — the static routes'
time to within 0.1% — then returned to 3.33 s on the third repeat under the
swept plan. The same binary, the same route, the same width; only the RHS plan
differs, and it is worth 10x. Serial escapes because it runs
`rhs_batch_parallel=off` entirely.

R6's win at these sizes is calibration, not route selection, and the paper
should say so. It also means the static-route columns at N=64 and N=256 are
measuring the auto heuristic's failure, not a routing result — quote them as
such or the comparison flatters R6 for the wrong reason.

**3. R6 beats every pure route at small Monte Carlo budgets by being impure.**
The dispatch trace at budget 8 reads `workers=8 local_slots=4 pool_size=8`: it
dispatches to eight worker processes *and* runs four samples on the coordinator.
Neither `outer_threads` nor `outer_process` can express that, which is why R6 is
1.3-1.7x faster than both at budgets 2 and 4 on both Monte Carlo workloads.

**4a. The cold-store result at budget 8 was an artifact, and is withdrawn.**
Measured warm, R6 runs budget 8 in 1.470 s against the best static route's
1.555 s on P3 (cold: 2.453 against 1.628) and 2.959 s against 3.592 s on P4
(cold: 3.805 against 3.624). The controls move 1-6% between the two runs, which
is machine noise; R6 moves 40% and 22% at that one rung and is unchanged
elsewhere. The earlier reading -- that the bandit had discarded a measured-better
arm -- was a converging calibration store and nothing more.

**4b. At the full budget R6 pays for exploration, and a five-repeat median
charges it.** The per-repeat routes at budget 12 are identical on both
workloads: process, process, *threads*, *threads*, process. Its exploiting
campaigns are the fastest measurements at that point on either workload --
1.057 s against `outer_process`'s fastest repeat of 1.335 s on P3, and 2.524 s
against 2.832 s on P4 -- but the two exploring repeats cost 2.52 s and 1.53 s
on P3 and 4.14 s and 3.49 s on P4, and in both cases the five-repeat median
lands on the second of them (1.531 s and 3.486 s). At budget 8 there is no
flapping and R6 wins outright. Quote both: the median over a short campaign
sequence is what a user with five campaigns to run experiences, and the steady
state is what the routing itself achieves.

**5. R6's variance is far higher than the static routes'.** Every mode pays a
warm-up cost on its first repeat; what separates them is what happens after it.
At P4's budget 12, repeats 2-5 of `outer_process` span 2.832-2.970 s -- a 5%
spread -- while R6's span 2.524-4.136 s, a 64% one. P3's budget 12 is the same
shape: `outer_process` 1.335-1.492 s against R6's 1.057-2.521 s. The pinned
routes converge once warm and the adaptive one keeps exploring, so R6's
five-repeat median is the least stable number in the set. Any single-median
comparison of an adaptive route against a pinned one should be read with that
spread in view, and finding 4b's two readings quoted together.

**6. Cold and warm calibration stores are different measurements.** This run
started from a fresh worktree, so the RHS-calibration and inner-policy stores
began empty and were still converging while P3 and P4 ran; every prior B- and
L-series number used the repository's accumulated store. A warm re-run of P3/P4
is recorded in the "Calibration-store warmth" section of the generated tables.

**8. R6 stops exploring once the calibration store converges, and can be locked
onto a plan 2.4x slower than one it had already found.** This is the most
consequential thing the 11-repeat work turned up, and it was invisible at 3 and
5 repeats.

At P1's N=1024 rung, two independent runs from a partly-converged store explore
and improve monotonically:

```
ctrl: r1 cache/heuristic@0 = 5.243 | r2 sweep/satellite_batch@1 = 3.236 | r3 = 2.075
new : r1 cache/heuristic@0 = 5.213 | r2 sweep/satellite_batch@1 = 3.325 | r3 = 2.065
```

The two runs agree on the final value to 0.5%, and both are *still improving*
when the repeats run out. For scale, every pinned static route at that rung
takes 4.86-5.32 s and serial takes 12.3 s, so `sweep/satellite_batch@1` at
~2.07 s is **2.4x faster than the best pinned route and 6x faster than
serial** -- by a wide margin the largest routing win anywhere in the P-series.

A third run, started from a store those runs had left converged, never found it:
all 11 repeats ran `cache/heuristic@0` at 4.85-5.23 s, which is parity with the
pinned routes. The exploration rate across the whole P1/P2 set tells the same
story -- share of `policy_v2` campaigns that ran a `sweep`: 33% and 22% from the
partly-converged store, **5%** from the converged one, with 4 of 12 points never
exploring at all.

**It is localized, not pervasive, and that matters for how it is quoted.** A
cleared-store 11-repeat run of P1 settles the frequency. It finds
`sweep/satellite_batch@1` on repeat 1 at N=1024 and converges to a median of
1.987 s against the converged store's 5.01 s -- 2.5x, slightly better than the
2.07 s the 3-repeat runs had reached, as expected from runs that were still
improving when they stopped. But at the other five rungs the two stores agree
within 3%:

| N | serial | best static | R6 converged | R6 cold | cold/converged |
|---:|---:|---:|---:|---:|---:|
| 1 | 8.65 | 8.54 | 8.93 | 8.64 | 0.97 |
| 16 | 10.70 | 2.11 | 2.11 | 2.16 | 1.02 |
| 64 | 11.23 | 34.67 | 3.25 | 3.29 | 1.01 |
| 256 | 11.67 | 21.02 | 2.80 | 2.85 | 1.02 |
| **1024** | 11.84 | 4.96 | **5.01** | **1.99** | **0.40** |
| 4096 | 11.38 | 2.17 | 2.23 | 2.16 | 0.97 |

So the lock-in is severe where it lands and absent elsewhere: one of six P1
rungs, with exploration over the whole phase running 26% of campaigns from cold
against 6% from converged. The claim to make is "a converged store can
foreclose a large win at some points", not "R6 is 2.4x slower with a warm
store". The N=64 and N=256 columns are finding 2's calibration artifacts and
both stores handle them identically.

**A limit of this experiment, which bounds what the later phases can show.**
The store is shared across a run and fills as the phases progress, so a
"cleared-store run" is only genuinely cold for whichever phase meets a case
first. P2's case *is* P1's top rung (`gravity_4096sat_l50_vacuum_5800s`), so by
the time P2 runs, P1 has calibrated it in both runs and the two are comparing
identical store states. P2 duly shows nothing -- every rung within 3%,
exploration 4 sweeps against 3 of 66 campaigns -- and that is a null by
construction rather than evidence of no effect. It should not pad the
denominator. P3, P4 and P5 carry cases no earlier phase touches, so they remain
real tests.

Three things follow, and they matter more than the routing numbers themselves.

*The cache suppresses the search that would fix it.* R6 exploits a cached
verdict without re-testing, so a store that converged early on a mediocre plan
keeps R6 on it indefinitely. This is finding 4a's mirror image: there a cold
store made R6 look worse than it is, here a converged store does the same by a
different mechanism, and the second is worse because it is stable and therefore
looks like a real measurement.

*"Warm store" is not one operating point.* The methodology used throughout this
document, and the warm/cold contrast in finding 6, treats store warmth as a
transient to be got past. It is not: the converged state is a persistent input
that selects which plan R6 runs, and two warm stores can differ by 2.4x.
Anything quoting a single R6 number has to say which store produced it.

*The headline R6 result is a floor, not a ceiling.* The P-series medians were
all measured from stores in some partly-converged state. At least at this rung
R6's reachable performance is far better than its reported median, and the
limiter is the policy's exploration schedule rather than its routing.

Worth an issue against the inner-policy cache independently of the paper: a
cached verdict should carry an expiry or a re-test probability, so a converged
store cannot permanently foreclose an arm the policy has never compared against
on this machine.

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

**7. Every phase was measured twice, and they agree.** P1, P2 and P5 were
re-run warm after P3/P4 showed that store warmth matters, so the whole set now
shares one convention -- and the two runs, a day apart, are a reproducibility
check the single-run B- and L-series never had. R6 moves 0.91-1.06 across all of
P1 and P2 and 0.95-1.03 across all twelve P5 points; the static routes move by
a median of 1.6% over the 69 points measured in both runs, worst case 7.8%. The N=64 result reproduces at 13.94 s against 5.71 s (first run: 14.08 s
against 5.73 s), and P5's 6x2 win at 1.740 s against 2.831 s (first run: 1.835 s
against 2.760 s). Only P3 and P4 moved materially between cold and warm, and
only at one rung each -- which is the finding, not noise.

## The S1 speedups rest on a different serial baseline

`paper_scenarios`' S1 reports up to 21.6x at L50/4096 on twelve threads where P2
measures 4.46x. Both are arithmetically correct; they divide by different
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
baseline gives 35.7 / 7.45 = 4.8x, against P2's 4.46x.

**Consequence for the paper.** Every S1 speedup -- the 21.6x, and the
10.23x/11.80x/10.17x of the CPU-only scaling work -- divides by a serial run
4.5x slower than the fastest single-threaded code the library runs on its own
defaults. Both definitions are defensible ("all parallelism off" against "best
serial implementation"), but a headline speedup invites the question of whether
the baseline was optimised, and only one of the two answers it well. Pick one,
state it, and do not print numbers from both conventions in the same table.

## Open

- **Which serial baseline the paper adopts.** S1 divides by an
  `rhs_mode=serial` baseline (53.90–53.93 s) and the P-series divides by
  `rhs_mode=auto` (11.84–12.13 s) on the same case. Every S1 speedup is
  therefore ~4.5x larger than the corresponding P-series figure for arithmetic
  reasons alone; re-based, S1 gives 4.8x against P2's 4.46x. The manuscript has
  to pick one convention and restate the other set. This is the one decision
  still needed from the author.

- **The auto heuristic's mid-size failure is a live defect, not just a
  measurement.** Finding 2 shows the batched-RHS heuristic costing 3x at N=64
  and 2x at N=256 on a 64-core box, on every static route. R6 routes around it;
  nothing else does. Worth an issue against the RHS calibration path
  independently of the paper.

- **Warm P3/P4 on TRX50.** The local methodology re-ran the two Monte Carlo
  ladders against a converged calibration store; TRX50's ran once from the
  release's store. The TRX50 P3/P4 columns are therefore cold-ish and may
  understate R6 at the low budgets, the same way finding 4a's withdrawn
  regression did locally. Optional — the qualitative result is unchanged.

- **Every P-series R6 number predates finding 8 and should be relabelled by
  store state.** The tables above quote R6 medians measured from stores in
  various partly-converged states, which finding 8 shows is a selector for which
  plan R6 runs rather than a transient. Two 11-repeat runs are being measured to
  bracket this properly: `paper_benchmarks_trx50_r11/20260918_140559` is the
  converged-store end (P1 and P2 complete, P3-P5 abandoned once the lock-in was
  found), and a cleared-store pass of all five phases is in flight. Until both
  land, quote R6 with the store state named.

- **The default 8% band in `compare_paper_routing_runs.py` is now known to be a
  p90, not a ceiling.** The measured noise section above gives per-phase and
  per-mode figures; a `policy_v2` point can move 149% between identical-code
  runs and a P1 static point 46%. The script still applies one flat band, which
  is the right default for a first look but will keep flagging `policy_v2` as a
  change when it is variance. Teaching it the per-mode bands measured here
  would remove most of that.
