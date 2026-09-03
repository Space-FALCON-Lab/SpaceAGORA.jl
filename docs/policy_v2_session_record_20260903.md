# Policy V2 (profile R6): Session Record, 2026-09-02/03

Branch `policy-v2` (worktree `.claude/worktrees/policy-v2`, cut from
`3be5b0ec` on `router-eval-expanded-b6`). Reference machine: AMD Ryzen 9
9900X, 12 physical cores / 24 SMT, 60 GB, Julia 1.12.1. Nine commits on the
branch (9a9022ac .. the commit carrying this section); not pushed.

Companion: `docs/adaptive_policy_session_record_20260902.md` (the shipped
policy's state and methodology), and the review findings in the memory note
`project_adaptive_policy_review_20260902`.

---

## 1. Outcome in one paragraph

One switch, `SPACEAGORA_PARALLEL_POLICY_V2`, bundled as profile **R6**
(= R5 + the switch; harness and probe mode `policy_v2`). Off is byte-identical
to the shipped algorithm. On, four changes driven by the 2026-09-02 paper run:
Monte Carlo routes to whichever outer axis has more cores; the route selector
measures both parallel arms before exploiting any; a cached "heuristic won"
calibration verdict is honoured on long solves; callback and effector widths
come from the static rule rather than the per-call hint store. Paired probes
show R6 at parity with the best static route on every constellation and Monte
Carlo point measured, with no significant loss against R5 anywhere and wins of
21% (constellation, t8) to 71% (Monte Carlo, narrow threads / wide pool).

---

## 2. What the paper run showed (commit 3be5b0ec, outdir 20260902_153738)

**Monte Carlo.** The shipped default (`_priority_outer_route_montecarlo`)
returns threads whenever `nthreads > 1`. Pinned routes on the B12/B13 grids:

| case, (threads, workers) | pinned threads | pinned process |
|---|---|---|
| heavy_aerobraking (12, 1) | 10.2 s | 59.2 s |
| heavy_aerobraking ( 6, 2) | 13.9 s | 31.5 s |
| heavy_aerobraking ( 4, 3) | 19.3 s | 20.8 s |
| heavy_aerobraking ( 3, 4) | 24.1 s | 16.2 s |
| heavy_aerobraking ( 2, 6) | 34.7 s | 11.9 s |
| heavy_aerobraking ( 1,12) | 57.5 s |  7.4 s |
| independent_1sat (12,12)  |  1.65 s | 1.01 s |
| independent_1sat ( 8,12)  |  1.88 s | 1.03 s |

"Process when workers >= threads" reproduces every row. The adaptive profiles
chose threads at every point with threads available. Exploration could not
rescue it: `_route_is_proven` declares the default proven once it beats serial,
before process is tried, and the harness routes each point cold.

**Constellations.** R5 regret against best static was within +-3% or a large
win everywhere except the eight-thread rung (gravity_4096 +28%, surrogate +25%,
heavy_1024_6hr +8%, cadence_1024 +8%). Every losing t8 point carried
`rhs_plan_source=sweep` with verdict heuristic on all repeats -- the
always-resweep rule for solves over 1 s -- while the same cases at t12 hit a
cached heuristic verdict. The large wins (stack256_e3..e5, -26% to -46%) are
credited by telemetry to the sweep-pinned satellite_batch plan, not to the
per-call width layer; on the callback-dominated surrogate case the static
inner_only route was best at both thread counts.

---

## 3. What R6 changes

| Mechanism | File | Under V2 |
|---|---|---|
| Monte Carlo default | `routing/outer_route_selection.jl` `_priority_outer_route_montecarlo` | process when `process_max_workers >= outer_thread_budget` and the campaign clears the size threshold; else threads |
| Process worker cap | `routing/outer_route_state.jl` `_outer_process_worker_cap` | `usable_core_budget()` capped by `SPACEAGORA_PERF_PROCS` when set |
| Exploration guard | `_any_candidate_proven` + `must_measure` | forced exploration stops when ANY arm is proven, but threads and process must both have been measured; serial never required |
| Calibration cache | `engine/rhs_calibration.jl` `_rhs_calib_cached_verdict` | a cached heuristic verdict is honoured regardless of solve length; a cached PLAN still re-sweeps on long solves |
| Per-call width | `policy/adaptive_decision.jl` | `adaptive_active = false` -> static `min(items, budget)`; hint store and AIMD untouched but not taken |
| Hint gate | `policy/policy_telemetry.jl` `_hint_layer_pays` | fails closed on an unmeasured source |
| Hint confidence | `policy/persistent_hints.jl` `_hint_mean_and_width` | width scaled by the arm's std (was raw ns; chooser was greedy) |
| Observation context | `policy/context.jl` `policy_context_hint`, `SharedBuffers.policy_context` | worker-task observations reach the scoped context (TLS is not inherited by `@spawn`/Polyester) |

Harness/probe fixes that were needed to measure this at all:

- `benchmarks/studies/parallelization_performance/execution.jl`
  `ppc_resolve_outer_backend` now resolves the route UNDER the mode's env.
  Before, adaptive profiles were routed with shipped defaults whatever they
  declared, so a paired R6-vs-R5 campaign probe compared two identical routes.
- `scripts/paired_campaign_probe.jl`: per-arm route history (the two arms
  shared one in-process bandit state), per-pair timing output, and
  `--profile=` (default `smoke`, a 600 s heavy arc; the paper phases use
  `full`, 21600 s -- the probe's historic numbers were smoke-profile).
- `scripts/paired_profile_probe.jl`: `policy_v2 -> R6` in its mode table.

---

## 4. Measurements (paired, order-alternating, two-sided sign test)

Constellations, 15 pairs, `--isolate-calib`:

| case | threads | R6 vs R5 | R6 vs inner_only (best static at t8) |
|---|---|---|---|
| gravity_4096_l50 | 12 | -4.4% (10-5, n.s.) | +3.6% (5-10, n.s.) |
| heavy_1024_l50_6hr | 12 | -0.2% (8-7, n.s.) | -1.0% (9-6, n.s.) |
| interact_256_aero | 12 | +6.2% (5-10, n.s.) | +2.1% (5-10, n.s.) |
| gravity_4096_l50 | 8 | **-20.7% (15-0, p=0.0001)** | -5.2% (10-5, n.s.) |
| heavy_1024_l50_6hr | 8 | **-3.8% (12-3, p=0.035)** | +1.8% (4-11, n.s.) |
| interact_256_aero | 8 | +4.0% (6-9, n.s.) | **-17.3% (15-0, p=0.0001)** |

Monte Carlo, full profile, 64 samples, 9 pairs, per-arm route history:

| (threads, workers) | R6 vs R5 heavy | R6 vs R5 indep. | R6 vs pinned process heavy | indep. |
|---|---|---|---|---|
| (12, 12) | **-27.7% (9-0)** | -37% median, 6-3 (see note) | -1.1% (5-4, n.s.) | -1.0% (5-4, n.s.) |
| ( 6,  2) | **-1.6% (8-1, p=0.039)** | **-1.7% (8-1, p=0.039)** | -- | -- |
| ( 2,  6) | **-71.2% (9-0)** | **-67.0% (9-0)** | +0.1% (4-5, n.s.) | +1.1% (4-5, n.s.) |

Note: at (12,12) R6's first three timed campaigns on `independent_1sat_1hr`
ran 0.55-0.75 s before settling at 0.26-0.29 s against R5's 0.40-0.47 s; the
probe warms each arm with one campaign and that is not enough for the process
pool on a freshly built case. Steady state is -40%, matching B12. The probe
should warm process arms twice.

The smoke-profile pass (600 s arc) gave the same signs at every point.

---

## 5. Tests

`test/unit/parallel/policy_v2_tests.jl`: 90 assertions, all pass. PR contract
gates all pass (the mode-profile parity gate now checks `policy_v2` vs R6).

Two PRE-EXISTING failures on the shipped path, unchanged by this branch:

- `test/unit/parallel/outer_route_persistence_tests.jl:95` -- the beaten-default
  case; the shipped selector picks serial. Passes under R6.
- `test/suites/02_callbacks_parallel_and_smoke_tests.jl:1043,1071` -- expects an
  `efficient_satisfied` AIMD window. Since `433b9e39` seeds desire at the cap
  (> budget), every window is `efficient_deprived`; that commit did not touch
  the test. Fixable by starving the controller once before asserting.

---

## 6. Not done (as of the first pass; see section 7 for what the second pass added)

- Deleting the AIMD path and inner-hint cold eviction.
- Knob pruning; user docs for R6.
- Thermal/control callback width calibration (state-advancing / integrator-bound).
- The two pre-existing failing tests on the shipped path (section 5).

---

## 7. Second pass (2026-09-03, commits fc44c0cd .. f53b821a)

Added under the same switch: a pre-solve density-callback width sweep with the
no-regret floor; an in-campaign split race (three samples per worker per
width, else the widest split; widths are never explored one per campaign);
no forced route trials (the core-budget default stands unless held history
says otherwise); a full collection on the coordinator before dispatch. The
campaign probe gained `--src-runner` (dispatch through the shipped
`run_monte_carlo(threads=:auto)`, pinned arms through the shipped pool),
`--warm-campaigns`, `--trace-route`.

Constellations, 15 pairs, with the density sweep (interact_256_aero threads
its density callback):

| case | threads | R6 vs R5 | R6 vs inner_only |
|---|---|---|---|
| interact_256_aero | 12 | +0.4% (7-8, n.s.) | +2.0% (6-9, n.s.) |
| gravity_4096_l50 | 12 | +4.1% (5-10, n.s.) | -3.3% (9-6, n.s.) |
| heavy_1024_l50_6hr | 12 | -0.1% (8-7, n.s.) | +0.9% (7-8, n.s.) |
| interact_256_aero | 8 | -5.4% (8-7, n.s.) | -- |
| gravity_4096_l50 | 8 | **-22.1% (13-2, p=0.007)** | -- |
| heavy_1024_l50_6hr | 8 | **-4.4% (15-0, p=0.0001)** | -- |

Monte Carlo through the PRODUCTION runner (`--src-runner --profile=full`,
64 samples, 9 pairs, 2 warm campaigns), final code:

| (threads, workers) | R6 vs R5 heavy | R6 vs R5 indep. | R6 vs pinned process (same pool) |
|---|---|---|---|
| (12, 12) | **-29.1% (9-0)** | -38.0% median (7-2; first 4 pairs are pool warm-up, steady 0.27 vs 0.44 s) | heavy -2.2%, indep. +3.4%, both n.s. |
| ( 2,  6) | **-71.0% (9-0)** | **-67.0% (9-0)** | -- |
| ( 6,  2) | **-2.7% (9-0)** | +0.7% (4-5, n.s.) | -- |

R6 chose process on every campaign at the two wide-pool points and threads at
the narrow-pool point; no route trial was spent.

What the production-runner probe found on the way, each fixed in f53b821a:

- With one serial warm-up sample, a race on a pre-existing pool timed each
  worker's JIT into its batches (0.64 s steady against 0.27 s). Process-route
  warm-up is now one sample per worker of the widest split.
- One sample per worker per width is not a measurement (two runs of the same
  code ranked [4, 8, 12] differently). Three rounds per worker, else no race.
- Without a race the shipped split selector explored w4 and w8 for two
  campaigns each: eight consecutive campaigns at 0.65 s. Under V2 the split
  selector never force-explores.
- A process campaign following a threaded campaign in the same coordinator
  ran 2.6x slow (19.3 s vs 7.5 s): coordinator full collections stalling the
  async dispatch loop. Collect before dispatch; the harness always did.
- The must-measure threads trial cost a 26-45 s campaign on the heavy case
  and, on the cheap case, landed while the pool was still warming; the
  bandit then exploited threads for the rest of the run. Route trials are
  off under V2.
- R5's shipped guard spends a campaign on SERIAL for Monte Carlo (55 s
  against 10 s threaded on the heavy case): its exploration order tries
  `:none` before the default is proven. Previously misread as a JIT cost.
- The process pool warms over several campaigns on a fresh case, not one:
  the pinned-process arm read 0.66, 0.60, 0.28 s over its first three
  campaigns on independent_1sat_1hr.
