# The predictive campaign planner

The predictive planner (`src/orchestration/simulation/campaigns/predictive_planner.jl`,
with its dispatch in `adaptive_routing.jl`) computes a Monte Carlo campaign's
plan before dispatch from the campaign's shape, the available workers and the
machine's calibrated constants. It defaults to the best predicted
static-equivalent plan and leaves it only when another plan clears a margin. A
guard can narrow a mixed plan while it runs. These are safeguards around an
approximate cost model; they do not bound how much slower than the best plan a
campaign can actually run.

This page describes the current model and its controls. The measurements behind
the defaults are kept in the lab's research records (see
[Record a result](../src/user/studies_benchmarks.md#record-a-result)); a
historical result does not establish performance or scientific equivalence at a
later revision.

The planner chooses consumers and inner thread budgets. When a plan has one
consumer and at most one effective inner thread, dispatch also runs the sample
in the static serial environment (`_SERIAL_SAMPLE_ENV` in
`adaptive_routing.jl`): the adaptive policy and the pre-solve RHS calibration
sweep are off, and the density, control and thermal callbacks are not threaded.
The environment is restored after dispatch. Larger effective budgets keep the
existing inner policy. An explicit inner budget remains a ceiling; switching the
inner machinery off does not make an undeclared budget explicit.

None of these scheduling settings changes the configured physical models or
solver settings. Whether scheduling leaves trajectories unchanged still has to
be checked for the workload and the native-model history involved.

## The switch

`SPACEAGORA_CAMPAIGN_PLANNER` selects `bandit` (the default, the learning route
planner) or `predictive`. `campaign_planner_mode()` reads it on every call, so a
`withenv` block around one campaign selects the planner for that campaign. An
unrecognized value throws an `ArgumentError`. The parallel profile `R7` is
profile `R6` with `SPACEAGORA_CAMPAIGN_PLANNER=predictive` and every other
setting unchanged.

## The plan space

For `n` samples, `T` coordinator threads and `W_p` affordable pool workers:

| Plan | Consumers | Static-equivalent form |
|---|---|---|
| `:none` | One, with the available inner budget | Serial campaign |
| `:threads` | Up to `min(n, T)` tasks | Full outer width, inner budget 1 |
| `:process` | `W_p` workers plus `L` local slots | `L = 0`, inner budget 1 |

Local slots are bounded by the configured cap, the remaining samples and
`ParallelProfiles.mixed_local_slots`. Which routes are feasible still comes from
`outer_route_candidates` and `effective_process_workers`, including the memory
and native-GRAM checks; the planner replaces the cost decision after those
checks, not the checks themselves. Without a measured inner-speedup curve,
narrower budget-1 widths offer no modeled benefit, and `:none` is chosen only
for a single sample or a single thread, or when no parallel route is feasible.

### Inner thread budgets

`SimulationEngine.rhs_inner_speedup_curve(stem)` reads candidate timings from
the RHS calibration store and forms
`speedup[b] = t(1) / min(t(w) for w <= b)`. Without a one-thread timing there is
no curve. `campaign_route_features` associates the configuration's RHS
signature stem with the campaign's workload signature; features without that
association cannot supply a curve.

With a curve, candidate plans can give a sample `b > 1` threads: narrower
threads-route splits, serial samples at budget `T`, or local slots in a mixed
plan. Process workers stay single-threaded. Every split respects its thread
share, and mixed local slots leave a coordinator thread for dispatch. A flat
curve keeps the static budget-1 plan through tie-breaking.

Applying an RHS-only speedup to a whole sample is optimistic: solver, callback
and other work may not scale with it, and a calibration sweep need not reflect
contention from concurrent samples. The margin does not make either assumption
valid. `SPACEAGORA_PREDICTIVE_INNER_CURVE=0` turns this input off.

## The cost model

The planner uses a deterministic greedy list schedule: the next free consumer
takes the next queued sample, with ties broken by index. Costs are in units of
one uncontended sample time, and samples are modeled as identical, so a
systematic difference in work between samples is outside the model.

A pool worker costs `1 + remote_overhead` per sample. Consumers that share the
coordinator's heap can carry a slowdown `s_heap(k)`, applied according to
`SPACEAGORA_PREDICTIVE_HEAP_MODEL`:

- `locals` (default): charge coordinator local slots in mixed plans; leave the
  threads-route static plan uncharged.
- `none`: charge no heap contention.
- `usl`: apply the calibrated USL term to both local slots and threads-route
  tasks.

The calibrated term is fitted to allocation kernels, not to whole samples.
`locals` is the default because, on the workloads and hosts measured, that
term described local slots well and the threads-route plan badly; this does not
show the mapping holds for another workload or host. An uncalibrated machine
has no contention term. Machine constants are loaded and traced even when the
selected model does not use them.

### Fixed cost of a short process campaign

The common sample-time scale cancels out of the core schedule comparison.
Costs given in seconds need a positive sample-time estimate to convert to
sample units, and are left out when none is available.

`PREDICTIVE_POOL_FIXED_S` (0.025 s) is the fixed cost of dispatching to a warm
pool. It is added to a process plan, in sample units, only when the sample
count fits in one pool round and the cost is at least
`PREDICTIVE_POOL_FIXED_MIN_SHARE` (5%) of that plan's predicted makespan before
the addition. Multi-round campaigns and longer single-round campaigns are
unaffected. It is separate from the cold-pool `pool_startup_s` and can apply to
a pool that is already warm.

The conversion uses the correction `sample_time_s` when available, otherwise
the last valid pool-worker sample mean for that workload signature in this
process. That fallback is locked, stays active with
`SPACEAGORA_CAMPAIGN_CORRECTIONS=off`, and neither changes nor persists
corrections; forced process measurements also populate it. Without a prior
observation the fixed term is zero, so an untimed warm-up can change later
timed route choices.

Both constants are heuristics informed by measurements of warm-pool dispatch.
The unit tests check the selection rules (long single rounds and multi-round
campaigns keep the pool); they do not show a wall-time gain.

### The local-slot heap term

With `local_heap_slope = c` in the constants file's `[campaign]` table, the
`locals` model uses `s(k) = 1 + c k (k - 1)` (`PairwiseHeapTerm`). The
coefficient is fitted from sample-level work ratios, not chosen to make the
planner rank particular plans correctly. Without that field, local slots fall
back to the calibrated USL term; without machine constants they have no
contention term.

The heap correction scales only the excess,
`s = 1 + heap_scale * (s_model - 1)`, and its observation is
`(s_observed - 1) / (s_model - 1)`. A scalar correction changes a curve's level,
not its shape, and one machine-wide slope does not mean contention is the same
across workloads. No coordinator-load term ships: none has been qualified.

## The final round and the cold pool

The optional `[campaign]` constants are model inputs shared by every plan:

| Field | Meaning | When absent |
|---|---|---|
| `round_tail` | Add `round_tail * (H_k - 1)` to each consumer class's final-round finish, where `H_k` is the harmonic number of that class's final-round width. | No tail term |
| `pool_startup_s` | Delay cold pool workers' first samples; local slots start at once. Seconds convert to sample units only when the workload's sample time is known. | No startup term |
| `local_heap_slope` | The pairwise local-slot coefficient above. | Calibrated USL fallback |

The tail form models sample durations with an exponential tail; a final round of
one sample has no excess. Applying it to heap-side consumers as well as to pool
workers is an assumption. The cold-pool charge applies only to workers that
predictive process dispatch has not yet used in this process, and is left out
when the sample time is unknown.

`scripts/calibration/extract_campaign_cost_terms.py` derives these terms from
archived traced runs, and `--write-constants` writes them into the `[campaign]`
table of the fingerprinted machine constants file, with their source.
`save_machine_constants` carries that table across a recalibration. Fitted
values do not transfer between machines: measure them on the machine in use and
keep their provenance.

## The decision rule

A non-static plan must beat the best static-equivalent plan by the configured
relative margin. Ties go to static-equivalent plans, then to `:process`, then to
fewer local slots. The process preference is only a tie-break: a threads plan
with a shorter predicted makespan wins. Every comparison is a prediction,
subject to the model's assumptions and measurement error.

## The guard

A mixed plan (`L > 0`) is watched inside a single dispatch; there is no barrier
or second dispatch. `DispatchAdmission` flags can close consumers before their
next job, but never interrupt work in flight or drop a sample already taken, and
a closed consumer cannot reopen. Trimming local slots closes the
highest-numbered slots first, so the narrowing is deterministic.

The guard decides once every consumer has reported at least one sample and at
least as many samples as there are consumers have completed. It compares
consumers rather than class totals, because local completions can arrive before
a worker's first result. A short campaign can drain before that; its verdict is
then traced as `never_decided`, with nothing left to change.

It distinguishes work (`elapsed_s`, timed inside a sample) from occupancy (the
interval between a consumer taking a job and receiving its result). Occupancy
includes communication and scheduling that sample work cannot see. Steady worker
occupancy leaves out a worker's first sample so that one-time dispatch costs are
not read as persistent starvation.

| Observation | Action |
|---|---|
| A sample failed | Close the local slots and keep the process route. |
| Worker occupancy exceeds the guard factor relative to local slots, and the reachable local-only plan is predicted to finish sooner | Close the pool workers only if `SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH` is on; otherwise report `route_switch_disabled`. |
| That occupancy condition holds but the alternative is no faster | Keep the plan (`threads_no_better`). |
| Steady worker occupancy shows the coordinator starving the workers, or local slowdown exceeds the thrashing threshold | Trim local slots toward the contention curve's throughput peak. |
| No qualifying observation | Keep the plan. |

A local slot that is merely slower than a worker can still add throughput, so
trimming responds only to harm to the workers or to local thrashing. The USL
peak is `sqrt((1-alpha)/beta)` (for the pairwise term, `1/sqrt(c)`). With no
curve, or one whose peak contradicts the observed need to narrow, the guard
halves the local slots instead; both fallbacks are assumptions. Without a
steady worker observation, the worker-starvation condition cannot fire.

### Route-switch limits

The route switch defaults to `false`, and its verdict is still traced. It was
left off because in testing a switch left the pool idle, and the next
campaign's first round then read the cost of waking the pool as evidence to
switch again. Closing pool workers leaves the existing `L` local slots consuming
the queue, so the reachable threads width is `L`, not a newly widened
`min(remaining, T)`. Its estimate uses the observed local occupancy and any
applicable contention scaling. Missing contention information makes it
optimistic; it is not evidence that widening or switching is free.

The guard cannot select a route that the feasibility checks excluded. When a
route change and a local trim both qualify, the route change wins.
`predictive_replan` states the two legal moves (narrowing on the process route,
or moving to threads with no local slots) and refuses anything else; the guard
itself closes consumers inside the running dispatch. Static-equivalent plans
have no local slots to close and dispatch without the guard.

`MonteCarloResult` reports the consumers admitted at completion. The dispatch
trace keeps the initial choice, the observation, the verdict and the final
plan; keep them apart when interpreting a benchmark.

## Online corrections

The predictive path neither reads nor writes the bandit's route-state file. It
persists bounded corrections to shared model parameters, not rewards for
particular plans:

| Parameter | Prior | Observation |
|---|---|---|
| `heap_scale` | `1` | Local-over-worker work slowdown from a guarded mixed plan that charged contention. |
| `round_tail` | The constants file's value | A successful pure-process campaign on a warm pool; no correction when the prior is absent or zero. |
| `sample_time_s` | Unknown, per workload | The workers' mean sample work, used to convert pool startup and the short-campaign fixed cost to sample units. |

For a parameter with a prior, evidence moves toward an observation by at most
`step * abs(prior)` per campaign and never goes negative. The value the planner
uses keeps `prior_share` of the prior and takes the rest from the evidence.
With the defaults `0.05` and `0.25`, one update moves the applied value by at
most `0.0375` of the prior. Sample time has no prior: its first observation sets
the estimate, and later changes are bounded by `step` of its current value. A
parameter with no observation for more than `stale_campaigns` (default 20)
campaigns reverts to its prior, and sample time becomes unknown. These values
are chosen controls, not measured convergence rates. A file whose schema
version, machine fingerprint or code token differs is ignored as if absent.

Corrections are kept only for shapes with a pool of at least two workers; other
shapes are planned from the constants file alone.

### The leash

A non-static choice moves at most one local slot from the last plan for the same
shape, in a direction the corrected model prices no worse. A static plan can
always be selected. The first campaign is governed by the margin alone.
Inner-budget plans step along the shape's available budget ladder; changing a
mixed plan's budget starts from one local slot at that budget. The recorded plan
is the one left after any guard trim. This bounds how far a plan can change, not
how much slower it can be on an unseen workload.

### The file, and how to reset it

Corrections live in `campaign_corrections_<fingerprint>.toml` in the directory
of the machine constants file, by default `output/parallel_policy_state/` under
the working directory; `SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH` overrides the
path. `SPACEAGORA_CAMPAIGN_CORRECTIONS` accepts `on` (default), `read` or
`off`: `read` uses existing evidence without writing, and `off` neither reads
nor writes. Writes are atomic. Delete the file to reset the evidence.

For a cold-store benchmark, leave the file out of the starting snapshot and set
corrections to `off`; otherwise the campaign can rebuild correction state while
the benchmark is meant to stay cold. This does not turn off the in-process
sample-time fallback for the short-campaign fixed cost, so start a fresh process
to begin without that history, and record warm-ups and earlier pool dispatches
with the results.

## The calibration store: per-candidate timings and a code token

The RHS store is `output/parallel_policy_state/rhs_calibration_<machine>.toml`
(`SPACEAGORA_RHS_CALIBRATION_PATH` overrides it), schema 2. Besides the verdict,
each row's `candidates` table records every candidate's mode, allotment,
scheduler, thread width, time (`ns`), repetition count and sweep round (0 for
the final paired round), and marks the plan the heuristic would have run
(`default`). Schema-1 rows still load but carry no candidate curve. A verdict
stored without a new sweep keeps the previous candidate timings.

Signatures end with `_RHS_CALIB_CODE_TOKEN` (currently `2026-09-23`, in
`src/orchestration/simulation/engine/rhs_calibration.jl`). Change this label
when a change to RHS execution invalidates calibration: old rows then stay in
the file but no longer match. The campaign corrections use the same token. A
store is valid only for its machine and code; an unchanged workload name does
not make it reusable.

## Constants

These defaults come from `PredictivePlannerConfig` and `CampaignCorrectionRules`.

| Setting | Default | Role |
|---|---|---|
| `SPACEAGORA_PREDICTIVE_MARGIN` | `0.15` | Minimum predicted relative gain for a deviation; an assumption, not a measured error bound. |
| `SPACEAGORA_PREDICTIVE_GUARD_FACTOR` | `1.5` | Guard threshold (assumed). |
| `SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX` | `max(0, T - 1)` | Cap that leaves a coordinator thread for pool dispatch; feasibility can lower it further. |
| `SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD` | `0.0` | Relative remote cost per sample; assumed zero, which short samples can violate. |
| `SPACEAGORA_PREDICTIVE_LOCAL_THRASH` | `3.0` | Local-over-worker slowdown past which the guard trims (assumed). |
| `SPACEAGORA_PREDICTIVE_HEAP_MODEL` | `locals` | Scope of the contention model; a fit needs qualifying for each workload. |
| `SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH` | `false` | Trace the switch verdict without acting on it. |
| `SPACEAGORA_PREDICTIVE_INNER_CURVE` | `1` | Use a matching measured RHS curve when one exists. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STEP` | `0.05` | Largest evidence move per campaign, relative to the prior. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_PRIOR_SHARE` | `0.25` | Share of the prior kept in the applied value. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STALE_CAMPAIGNS` | `20` | Campaigns without an observation before a correction expires. |
| `SPACEAGORA_CAMPAIGN_CORRECTIONS` | `on` | Persist bounded shared-parameter corrections. |

The USL `alpha` and `beta` come from `MachineConstants`: `alpha` is fitted to the
arithmetic calibration ladder and `beta` to the allocation ladder's excess over
it. They describe the calibration kernels, not proven whole-sample slowdowns.
The optional `[campaign]` parameters are measured inputs, absent by default. A
single scale cannot correct a contention curve of the wrong shape.

## Measurement hooks and validation

`SPACEAGORA_PREDICTIVE_FORCE_PLAN` is unset by default. For controlled
measurement it accepts a plan key such as `process@w32+l4` or
`threads@w8+l0+b4`; the route, worker count, local slots and budget must be
feasible, and an invalid key throws. A forced mixed plan still runs under the
guard, but forced campaigns neither fold nor save corrections nor update the
shape's last plan.

`SPACEAGORA_POOL_DISPATCH_PROBE=1` attributes dispatch costs before the timed
region. `SPACEAGORA_POOL_WORKER_GC` selects the collection workers run after a
dispatch: `incremental` (default), `full` or `off`. An incremental collection
does not prevent a later full collection from interrupting a campaign.
`SPACEAGORA_POOL_DISPATCH_CACHE` is off by default; turning it on requires that
the campaign function's captured state stays unchanged between reused
dispatches. In the parallelization benchmark harness,
`SPACEAGORA_PPC_SAMPLE_FN_REUSE` defaults to `1`, and `0` rebuilds the sample
closure for every batch. Measure these mechanisms separately from per-sample
work.

`benchmarks/studies/predictive_cost_terms/plan_parity.jl` writes each sample's
saved times and full state history, in sample order, through the planner's
actual dispatch. Use it to compare feasible plans, inner budgets and revisions
on the same inputs and seeds; earlier successful comparisons do not qualify
arbitrary workloads or history-dependent native models. The unit tests in
`test/unit/parallel/predictive_planner_tests.jl` and
`test/unit/parallel/predictive_corrections_tests.jl` cover the model's decisions
and the correction rules; their fixtures are not timing results.

## Tracing

`SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1` reports the candidates, constants and
cost terms, the chosen plan and the reason for it, any leash adjustment, the
observed work and occupancy, the guard verdict and the final admitted
consumers, and whether corrections were read or updated. The selection reasons
are `static_equivalent_best`, `predicted_gain` and `margin_not_met`; a selected
plan is not proof of a speedup.

Compare repeatedly against the same static baselines, and report failures,
noise and unsupported cases. Known limits of the model include sample-to-sample
cost spread, the effect of concurrent samples on inner speedups,
workload-dependent contention, and cold remote costs when there is no
sample-time evidence. A new fit, or turning on the route switch, needs its own
validation.
