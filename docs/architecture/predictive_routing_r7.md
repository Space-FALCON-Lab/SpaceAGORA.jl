# R7: the predictive campaign planner

The original investigation, including its measurements and limitations, is
preserved at this same relative path in the private archive described by the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).
This page keeps current methods and configuration; historical results do not
qualify performance or scientific equivalence at a later revision.

R7 computes a campaign plan from its shape, available workers and calibrated
machine constants before dispatch. It defaults to the best predicted
static-equivalent plan and requires a margin before selecting a deviation. A
guard can narrow a mixed plan while it runs. These are safeguards around an
approximate cost model, not a bound on actual runtime regret.

The planner selects consumers and inner thread budgets. For a plan with one
consumer and at most one effective inner thread, dispatch also selects the
static serial sample environment: adaptive policy and the pre-solve RHS
calibration sweep are disabled, and density, control and thermal callbacks are
not threaded. The environment is restored after dispatch. Larger effective
budgets retain the existing inner policy. An explicit inner budget remains a
ceiling; suppression does not make an undeclared budget explicit.

These scheduling settings do not change configured physical models or solver
settings. Scheduling independence and trajectory reproducibility still require
validation for the workload and native-model history involved.

## The switch

`SPACEAGORA_CAMPAIGN_PLANNER` selects `bandit` (the default R6 path) or
`predictive`. `campaign_planner_mode()` reads it on each call, so a `withenv`
around one campaign can select the planner. Unrecognized values throw an
`ArgumentError`. Profile `R7` selects predictive planning with the R6 inner
policy settings.

## The plan space

For `n` samples, `T` coordinator threads and `W_p` affordable pool workers:

| Plan | Consumers | Static-equivalent form |
|---|---|---|
| `:none` | One, with the available inner budget | Serial campaign |
| `:threads` | Up to `min(n, T)` tasks | Full outer width, inner budget 1 |
| `:process` | `W_p` workers plus `L` local slots | `L = 0`, inner budget 1 |

Local slots are bounded by the configured cap, remaining samples and
`ParallelProfiles.mixed_local_slots`. Route feasibility still comes from
`outer_route_candidates` and `effective_process_workers`, including memory and
native-GRAM constraints. The planner replaces the cost decision after those
checks, not the checks themselves. Without a measured inner-speedup curve,
narrower budget-1 widths offer no modeled benefit, and `:none` is limited to
the single-sample/single-thread or no-feasible-parallel-route cases.

### Inner thread budgets

`SimulationEngine.rhs_inner_speedup_curve(stem)` reads candidate timings from
the RHS calibration store and forms
`speedup[b] = t(1) / min(t(w) for w <= b)`. No one-thread timing means no curve.
The configuration's RHS signature stem is associated with the campaign workload
signature by `campaign_route_features`; features without that association
cannot supply a curve.

With a curve, candidate plans can give a sample `b > 1` threads: narrower
threads-route splits, serial samples at budget `T`, or local slots in a mixed
plan. Process workers remain single-threaded. Every split respects its thread
share; mixed local slots leave a coordinator thread for dispatch. A flat curve
keeps the static budget-1 plan through tie-breaking.

Applying an RHS-only speedup to a complete sample is an optimistic assumption:
solver, callback and other work may not scale with it. A sweep also need not
represent contention from concurrent samples. The margin does not prove either
assumption valid. `SPACEAGORA_PREDICTIVE_INNER_CURVE=0` disables this input.

## The cost model

The planner uses a deterministic greedy list schedule: the next free consumer
takes the next queued sample, with ties resolved by index. Costs are in units
of one uncontended sample time. Samples are modeled as identical; systematic
sample-to-sample workload variation is outside that assumption.

A pool worker costs `1 + remote_overhead` per sample. Heap-sharing consumers
can carry a slowdown `s_heap(k)`, applied according to
`SPACEAGORA_PREDICTIVE_HEAP_MODEL`:

- `locals` (default): charge coordinator local slots in mixed plans; leave the
  threads-route static plan uncharged.
- `none`: do not charge heap contention.
- `usl`: apply the calibrated USL term to both local slots and threads-route
  tasks.

An uncalibrated machine has no measured contention term. Machine constants are
loaded and traced even when a selected model does not use them. The historical
study explains the choice of `locals` and the limitations of transferring an
allocation-kernel fit to whole simulation samples; it does not establish that
this mapping is accurate for another workload or host.

### Fixed cost of a short process campaign

The common sample-time scale cancels out of the core schedule comparison.
Seconds-valued costs need a positive sample-time estimate for conversion; they
are omitted when that estimate is unavailable.

`PREDICTIVE_POOL_FIXED_S` adds a 25 ms process-dispatch cost in sample units,
only when the sample count fits in one pool round and the cost is at least
`PREDICTIVE_POOL_FIXED_MIN_SHARE` (5%) of the predicted makespan before this
addition. Multi-round campaigns and longer single-round campaigns are unchanged
by this term. It is distinct from cold-pool `pool_startup_s` and can apply to an
already warm pool.

Conversion uses correction `sample_time_s` when available, otherwise the last
valid pool-worker sample mean for that workload signature in this process. The
fallback is locked, remains active with `SPACEAGORA_CAMPAIGN_CORRECTIONS=off`,
and neither modifies nor persists corrections. Forced process measurements also
populate it. Without a prior observation the fixed term is zero, so an untimed
warm-up can affect later timed route choices.

The 25 ms value and 5% cutoff are heuristics informed by historical measurements.
Unit tests establish selection rules, including protected long single rounds
and multi-round pool choices; they do not establish a wall-time gain. Further
TRX50 measurements are deferred for software delivery and remain unverified.

### The local-slot heap term

With `local_heap_slope = c` in the constants file's `[campaign]` table, `locals`
uses `s(k) = 1 + c k (k - 1)`. The coefficient is fitted from sample-level work
ratios, not chosen to make the planner rank particular plans correctly. Without
that field, local slots fall back to the calibrated USL term; without machine
constants they have no contention term.

The heap correction scales only the excess:
`s = 1 + heap_scale * (s_model - 1)`. Its observation is
`(s_observed - 1) / (s_model - 1)`. A scalar correction changes a curve's level,
not its shape. One machine-wide slope does not establish equivalent contention
across workloads. The archived coordinator-load experiments did not qualify a
replacement; no coordinator-load term ships.

## The final round and the cold pool

The optional `[campaign]` constants are shared model inputs:

| Field | Meaning | When absent |
|---|---|---|
| `round_tail` | Add `round_tail * (H_k - 1)` to each consumer class's final-round finish, where `H_k` is the harmonic number of that class's final-round width. | No tail term |
| `pool_startup_s` | Delay cold pool workers' first samples; local slots start immediately. Convert seconds to sample units only when the workload's sample time is known. | No startup term |
| `local_heap_slope` | Pairwise local-slot contention coefficient described above. | Calibrated USL fallback |

The tail form models exponential-tailed sample durations; a final round of one
has no excess. Applying it to heap-side consumers as well as workers is an
assumption. A cold-pool charge applies only to workers not yet used by predictive
process dispatch in this process. It is omitted when sample time is unknown.

`scripts/calibration/extract_campaign_cost_terms.py` derives these terms from archived runs
and `--write-constants` writes the fingerprinted constants file's `[campaign]`
table, with source provenance. `save_machine_constants` retains that table
through recalibration. Historical fitted values are not portable defaults;
measure and retain provenance for the machine being used.

## The decision rule

A non-static plan must beat the best static-equivalent plan by the configured
relative margin. Ties prefer static-equivalent plans, then `:process`, then
fewer local slots. The process preference is only a tie-break: a faster
predicted threads plan wins on makespan. All comparisons remain predictions,
subject to the model's assumptions and measurement error.

## The guard

A mixed plan (`L > 0`) is observed inside a single dispatch. There is no barrier
or second dispatch. `DispatchAdmission` flags can close consumers before their
next job, never interrupt in-flight work or drop a sample already taken. Flags
cannot reopen. Trimming local slots closes the highest-numbered slots, making
the narrowing deterministic.

The guard decides once every consumer has reported at least one sample and at
least a full consumer count has completed. It compares consumers, not merely
class totals, because local completions can arrive before a worker's first
result. A short campaign may drain before that condition; its verdict is then
traced as `never_decided`, with nothing left to change.

It distinguishes work (`elapsed_s`, timed inside a sample) from occupancy
(the interval between that consumer taking a job and receiving its result).
Occupancy includes communication and scheduling that sample work cannot see.
Steady worker occupancy excludes a worker's first sample to avoid treating
one-time dispatch costs as persistent starvation.

| Observation | Action |
|---|---|
| A sample failed | Close local slots, retaining the process route. |
| Worker occupancy exceeds the guard factor relative to locals, and the reachable local-only plan is predicted to finish sooner | Close pool workers only if `SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH` is enabled; otherwise report `route_switch_disabled`. |
| That occupancy condition holds but the alternative is no faster | Keep the plan (`threads_no_better`). |
| Steady worker occupancy indicates coordinator starvation, or local slowdown exceeds the thrashing threshold | Trim local slots toward the contention curve's throughput peak. |
| No qualifying observation | Keep the plan. |

A merely slower local slot may still add throughput. Trimming is based on harm
to workers or local thrashing. The USL peak is `sqrt((1-alpha)/beta)`; an absent
curve, or one whose peak contradicts the observed need to narrow, falls back to
halving. These fallback decisions are assumptions. If no steady worker
observation exists, the worker-starvation condition cannot fire.

### Route-switch limits

The route switch defaults to `false`; its verdict is still visible in tracing.
Closing pool workers leaves the existing `L` local slots consuming the queue.
The reachable threads width is therefore `L`, not a newly widened
`min(remaining, T)`. Its estimate uses the observed local occupancy and any
applicable contention scaling. Missing contention information is optimistic,
not evidence that widening or switching is free.

The guard cannot select a route excluded by feasibility checks. When a route
change and a local trim both qualify, the route change takes precedence.
`predictive_replan` accepts narrowing on the process route or a move to threads
with no local slots; it rejects unsupported replans. Static-equivalent plans
have no local slots to close and dispatch without this guard.

`MonteCarloResult` reports the consumers admitted at completion. The dispatch
trace retains the initial choice, observation, verdict and final plan. These
must be distinguished when interpreting a benchmark.

## Online corrections

The predictive path neither loads nor writes the bandit's route-state file.
It persists bounded corrections to shared parameters, not per-plan rewards:

| Parameter | Prior | Observation |
|---|---|---|
| `heap_scale` | `1` | Local-over-worker work slowdown from a guarded mixed plan that charged contention. |
| `round_tail` | The constants file's value | A successful pure-process campaign on a warm pool; no correction when the prior is absent/zero. |
| `sample_time_s` | Unknown, per workload | Workers' mean sample work, used to convert pool startup and the short-campaign fixed cost to sample units. |

For a parameter with a prior, evidence moves toward an observation by at most
`step * abs(prior)` per campaign, remaining nonnegative. The applied value
retains `prior_share` of the prior and blends the remaining share with the
evidence. With defaults `0.05` and `0.25`, one update moves the applied value by
at most `0.0375` of the prior. Sample time has no prior: its first observation
sets the estimate, and later changes are bounded by `step` of its current value.
Evidence expires after more than 20 campaigns without an observation; sample
time then becomes unknown. These values are assumed controls, not measured
universal convergence rates. A machine/code/schema mismatch makes the persisted
state inapplicable.

### The leash

A non-static choice moves at most one local slot from the last plan for the
same shape, in a direction the corrected model prices no worse. A static plan
can always be selected. The first campaign is governed by the margin alone.
Inner-budget plans step along the shape's available budget ladder; changing a
mixed plan's budget starts from one local slot at that budget. The recorded plan
is the one remaining after any guard trim. This bounds the size of a plan
change, not measured regret on an unseen workload.

### The file, and how to reset it

Corrections live in
`output/parallel_policy_state/campaign_corrections_<fingerprint>.toml`;
`SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH` overrides the path.
`SPACEAGORA_CAMPAIGN_CORRECTIONS` accepts `on` (default), `read` or `off`.
`read` uses existing evidence without writing; `off` neither reads nor writes.
Writes are atomic and occur for shapes with a pool. Delete the corrections file
to reset that evidence.

For a cold-store benchmark, exclude the file from the starting snapshot and set
corrections to `off`; otherwise the campaign can rebuild correction state while
the benchmark is meant to remain cold. This does not disable the transient
sample-time fallback for the short-campaign fixed cost. Use a fresh process to
start without that history, and record warm-ups and prior pool dispatches with
the results.

## The calibration store: per-candidate timings and a code token

The RHS store is
`output/parallel_policy_state/rhs_calibration_<machine>.toml`, schema 2.
Besides the verdict, its candidate table records mode, allotment, scheduler,
width, score, repetitions and round. The final paired round records the winner
and heuristic default. Schema-1 rows remain readable but carry no candidate
curve. A verdict stored without a new sweep retains previous candidate timings.

Signatures include `_RHS_CALIB_CODE_TOKEN` (currently `2026-09-23`). Change this
version label when RHS execution changes invalidate calibration: old rows then
remain in the file but do not match the new signatures. The campaign corrections
use the same code token. Store reuse must be qualified against both machine and
code, rather than assumed valid because the workload name is unchanged.

## Constants

These defaults come from `PredictivePlannerConfig` and `CampaignCorrectionRules`.
They preserve the current implementation; none is a new performance claim.

| Setting | Default | Role and standing |
|---|---|---|
| `SPACEAGORA_PREDICTIVE_MARGIN` | `0.15` | Assumed minimum predicted gain for a deviation; not a measured error bound. |
| `SPACEAGORA_PREDICTIVE_GUARD_FACTOR` | `1.5` | Assumed guard threshold. |
| `SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX` | `max(0, T - 1)` | Cap leaving a coordinator thread for pool dispatch, further bounded by feasibility. |
| `SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD` | `0.0` | Assumed relative remote cost; short samples can violate it. |
| `SPACEAGORA_PREDICTIVE_LOCAL_THRASH` | `3.0` | Assumed local-over-worker slowdown threshold. |
| `SPACEAGORA_PREDICTIVE_HEAP_MODEL` | `locals` | Scope of the contention model; measured fits need workload qualification. |
| `SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH` | `false` | Observe the switch verdict without acting on it by default. |
| `SPACEAGORA_PREDICTIVE_INNER_CURVE` | `1` | Use a suitable measured RHS curve when available. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STEP` | `0.05` | Maximum evidence update relative to its prior. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_PRIOR_SHARE` | `0.25` | Prior's retained share. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STALE_CAMPAIGNS` | `20` | Expiry interval without an observation. |
| `SPACEAGORA_CAMPAIGN_CORRECTIONS` | `on` | Persist bounded shared-parameter corrections. |

The USL `alpha` and `beta` come from `MachineConstants`: `alpha` is fitted to the
arithmetic calibration ladder and `beta` to the allocation ladder's excess over
it. They are calibration-kernel parameters, not proven whole-sample slowdown
coefficients.
The optional `[campaign]` parameters are measured inputs, absent by default as
described above. A single scale cannot correct a wrong contention-curve shape.

## Measurement hooks and validation

`SPACEAGORA_PREDICTIVE_FORCE_PLAN` is unset by default. It accepts a plan key
such as `process@w32+l4` or `threads@w8+l0+b4` for controlled measurement. The
route, worker count, local slots and budget must be feasible; invalid keys throw.
A mixed forced plan still uses the guard, but forced campaigns do not fold or
save corrections or update the shape's last plan. Unset leaves normal selection
unchanged.

`SPACEAGORA_POOL_DISPATCH_PROBE=1` attributes dispatch costs before the timed
region. `SPACEAGORA_POOL_WORKER_GC` selects `incremental` (default), `full` or
`off` collection after dispatch. Incremental collection does not guarantee that
a later full collection cannot interrupt a campaign. `SPACEAGORA_POOL_DISPATCH_CACHE`
is off by default; enabling it requires that captured campaign-function state
remain unchanged between reused dispatches. The harness
`SPACEAGORA_PPC_SAMPLE_FN_REUSE` defaults to `1`; `0` rebuilds its sample closure
per batch. Measure these mechanisms separately from per-sample work.

`benchmarks/studies/predictive_cost_terms/plan_parity.jl` writes each sample's
saved times and full state history in sample order through the planner's actual
dispatch. Use it to compare feasible plans, inner budgets and revisions with the
same inputs and seeds. The archived successful cases do not qualify arbitrary
workloads or history-dependent native models. Unit tests in
`test/unit/parallel/predictive_planner_tests.jl` and
`predictive_corrections_tests.jl` cover model decisions and correction rules;
their historical fixtures are not fresh timing results.

The paper harness's optional X1/X3/X4/X5a/X5b explorations require their
explicit phase selection and compatible budgets. Their historical plan keys
require `SPACEAGORA_PPB_PAPER_BUDGET=32`, no preview, and sufficient thread or
worker budgets. Defaults omit these experiments. P6ps uses each matching
thread/worker sub-run's serial baseline; X5a/X5b use their shared phase baseline.
Ambiguous baselines remain missing rather than being filled with unrelated rows.

## Tracing

`SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1` reports candidates, constants and cost
terms, the chosen plan and reason, any leash adjustment, observed work and
occupancy, the guard verdict and the final admitted consumers. It also reports
whether corrections were read or updated. Distinguish `static_equivalent_best`,
`predicted_gain` and `margin_not_met`; a selected plan is not proof of a speedup.

Use repeated comparisons against the same static baselines and report failures,
noise and unsupported cases. Remaining model limitations include sample-to-sample
cost spread, concurrent-sample effects on inner speedups, workload-dependent
contention, and cold remote costs when no sample-time evidence exists. A new
fit or enabled route-switch policy requires its own validation.
