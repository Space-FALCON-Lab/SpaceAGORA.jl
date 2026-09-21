# R7: the predictive campaign planner

R7 replaces the R6 campaign-level learner with a planner. The campaign plan is
computed once, before any sample runs, from the campaign shape and the
machine's calibrated constants; it defaults to the plan a pinned static route
would run; it deviates only when the model predicts a gain larger than a
margin; and a guard after the first round of samples re-plans the remainder if
the observation contradicts the prediction.

R7 changes nothing inside a sample. The engine-level inner policy (the RHS
calibration sweep, the callback width rules) is exactly what R6 has. R7 decides
the campaign-level plan in `_run_campaign_adaptive` and nothing else.

## Why

R6 learns the plan: a bandit over `:none | :threads | :process`, a split race
over widths inside the first campaign of an unseen shape, and persisted state
so the next session exploits what this one paid for. Three things follow from
that, and the third is why R7 exists:

1. Exploration costs whole campaigns. An arm's first campaign is its cold one
   (pool spin-up, worker JIT, the collection that follows), so a rounds-tie
   trial is charged twice before it is judged once.
2. What it learns is keyed by a workload signature. A campaign shape seen once
   is never exploited.
3. The learner's own mistakes are unbounded. On the TRX50's P1-P5 grid, three
   points ended more than 10% behind the best pinned static route: P4 at budget
   32 (ratio 1.146), P5 `mcgrid_16sat_8mc` at 1x32 (1.109), and P5
   `mcgrid_8sat_16mc` at 1x32, where 10 of 11 campaigns threw.

A planner has the opposite failure mode. It can be wrong about the cost of a
plan, but the margin rule bounds what being wrong costs: inside the margin it
runs the static plan, which is the baseline it is being compared against.

## The switch

`SPACEAGORA_CAMPAIGN_PLANNER`:

- `bandit` (the default, and the value of an unset variable) -- R6 unchanged,
  bit for bit. Nothing on that path calls anything in this document.
- `predictive` -- the planner described here.

Read on every call through `campaign_planner_mode()`, so a `withenv` around one
campaign selects the planner for that campaign. An unrecognized value throws an
`ArgumentError` naming the two valid values rather than falling back, because a
typo in a profile's environment must not silently run the other planner.

Profile `R7` emits this variable together with everything `R6` emits.

## The plan space (v1)

With `n` samples, `T` coordinator threads and `W_p` affordable pool workers:

| Plan | Consumers | Static-equivalent |
|---|---|---|
| `:none` | 1, with the whole thread budget | yes |
| `:threads` at `W = min(n, T)` | `W` tasks, inner budget 1 | yes (`outer_threads`) |
| `:process` at `W_p` workers plus `L` local slots | `W_p + L`, inner budget 1 | when `L = 0` (`outer_process`) |

`L` runs from `0` to `min(SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX, n - W_p,
ParallelProfiles.mixed_local_slots(features, tuning, W_p))`. `:none` is
enumerated only for a single sample, or a single-threaded coordinator with no
pool, or as the fallback when no parallel route is feasible at all.

Route *candidacy* is still the shipped rule (`outer_route_candidates`,
`effective_process_workers`, `mixed_local_slots`). Whether a pool may be used
at all on a given shape and machine is a feasibility question -- memory, native
GRAM's thread-safety, the machine class -- and R7 replaces only the cost
decision that follows it.

Narrower widths are not enumerated, in either route. For the process route the
workers are one thread each, so a narrower pool only idles cores. For the
threads route, v1 hands every concurrent sample an inner budget of 1 regardless
of width, so a narrower split buys nothing to trade against the concurrency it
gives up.

### Limitation: every concurrent sample gets one thread

R7 v1 never assigns an inner thread budget above 1 under an outer split. The
budget is declared as an environment ceiling around the whole dispatch, which
`capped_inner_thread_budget` then leaves alone.

This is the rule "unknown means no gain" applied honestly: how much faster one
sample runs on `k` threads is a property of the workload that no a-priori
model here measures, and assuming a speedup the planner cannot predict is
exactly the kind of invention this planner exists to avoid. The visible
consequence is at `n < T` on the threads route, where a pinned `outer_threads`
run gives each of `n` concurrent samples `fld(T, n)` threads and R7 gives it
one. Making the inner budget a real decision requires an inner-speedup
measurement and is the first item of a v2.

## The cost model

Plans are ranked by predicted makespan in units of one uncontended sample time
`t_1`. `t_1` is unknown a priori and cancels out of the ranking, which is why
the planner never needs to know how long a sample takes -- only how the
campaign's time scales with how the samples are spread.

The schedule is a deterministic greedy list schedule: every consumer takes the
next queued sample the moment it is free, ties to the lower index. That is
exactly what the one-queue dispatchers in `monte_carlo.jl` do, so the model and
the dispatcher agree by construction. Samples are modeled as identical, because
a Monte Carlo campaign's seeds are draws from one distribution and the planner
has no per-sample information that would justify a spread.

Per-class sample times:

- **Pool worker**: `s_w = 1 + o_r`, with `o_r` the `remotecall_fetch` round
  trip relative to `t_1`.
- **Coordinator local slot, or threads-route task**: `s_heap(k) = k /
  usl_speedup(alpha, beta, k)`, the inverse parallel efficiency of
  allocation-heavy work on one heap, at `k = L` local slots or `k = W`
  threads-route tasks. Both use the same term because they are the same thing:
  samples sharing one process's heap and allocator.

With no machine constants, `s_heap(k) = 1` at every width. That is not a claim
that contention is absent; it is the planner declining to model what it has not
measured, and the margin rule is what carries the safety in that state. An
uncalibrated machine therefore ranks plans purely by round count.

## The decision rule

The chosen plan is the best static-equivalent plan, unless some other plan is
predicted to beat it by at least `margin` relative. Among plans that tie on
predicted makespan the ranking prefers, in order: static-equivalent plans, the
`:process` route, then fewer local slots. So the planner ends as close to the
static plan as the predicted gain allows -- on 64 samples over 8 workers it
takes the five local slots that buy the fifth round, not the seven it is
allowed.

The preference for `:process` over `:threads` among static-equivalents is a
tie-break, not an override. When the threads route is predicted faster it wins
on the makespan and the preference never fires.

## The guard

A plan that deviates (`L > 0`) and has more samples than consumers is
dispatched in two rounds. The first round is one sample per consumer; the mixed
dispatcher tags each sample with the class that ran it (`class_sink`, added to
`_run_monte_carlo_mixed` for this purpose -- nothing in a `MonteCarloResult`
otherwise records it), and the means of the two classes are compared against
the ratio the plan predicted, `s_heap / s_w`. The pool workers are the
uncontended reference, so the observed ratio is directly comparable.

| Observation | Verdict |
|---|---|
| ratio <= `guard_factor` | keep the plan |
| ratio > `guard_factor` | halve the local slots |
| ratio > 2 x `guard_factor` | drop the local slots to zero |
| any sample failed | drop the local slots to zero |

The guard only ever moves toward a static-equivalent plan: fewer local slots,
never more, never a different route. A guard that could widen would be an
optimizer running inside the campaign, which is what R7 exists not to be.

A static-equivalent plan is dispatched in ONE round, not two. The only
reduction the guard can make is to the local slots, which are already zero, so
a second round would buy no information and cost a synchronization barrier --
every consumer waiting on the round's straggler. This also keeps the static
plans dispatch-for-dispatch comparable with the pinned static routes they are
measured against.

`MonteCarloResult` reports the plan the remainder ran under: `route`,
`threads = W_p + L` (or `W`, or 1), and `local_slots = L`. The dispatch trace
carries the first round's plan and the change.

## No learning state

The predictive path records no bandit feedback and neither loads nor saves
`campaign_route_state_path()`. Two campaigns of the same shape in the same
process get the same plan for the same reason, and a campaign shape seen once
is planned as well as one seen a hundred times.

## Constants

Every number the planner uses, and what kind of number it is.

| Name | Default | Class | Source |
|---|---|---|---|
| `SPACEAGORA_PREDICTIVE_MARGIN` | `0.15` | ASSUMED | No measurement of this model's error exists, so no margin can be derived from one. 0.15 is a round number chosen to sit above the ~6% p90 identical-code noise floor of the reduced-scale harness and below the 25-40% mixed-dispatch wins R6 measured at the P3/P4 mid budgets, i.e. large enough to refuse noise and small enough to keep the wins that motivated mixed dispatch. Tune it, do not trust it. |
| `SPACEAGORA_PREDICTIVE_GUARD_FACTOR` | `1.5` | ASSUMED | Same standing. A local slot running 1.5x slower than predicted relative to a pool worker is outside anything the round-count model explains; the second threshold at 2x that is the "stop entirely" case. |
| `SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX` | `Threads.nthreads() - 1` | DERIVED | `ParallelProfiles.mixed_local_slots` keeps thread 1 free for the `@async` feeders that keep the pool supplied, so `T - 1` is the most slots the coordinator can offer under R6's own measured practice. |
| `PREDICTIVE_REMOTE_OVERHEAD` (`o_r`) | `0.0` | ASSUMED | A `remotecall_fetch` round trip is milliseconds; the campaign samples in the P3-P5 measurements are 0.03-1 s. The term is below the noise floor of anything the planner could measure at that ratio. Named rather than inlined so a workload with sub-millisecond samples has one place to change. |
| `alpha` in `s_heap` | `MachineConstants.usl_alpha_base` | SOURCED (fit), ASSUMED (mapping) | Fitted per machine by `scripts/calibrate_machine.jl` on the allocation kernel. Applying it to campaign samples claims only that a sample's contention with its neighbors on one heap has the same SHAPE as that kernel's, not the same magnitude. |
| `beta` in `s_heap` | `MachineConstants.usl_beta_alloc` | SOURCED (fit), ASSUMED (mapping) | As above. |
| `s_heap` with no constants | `1.0` | DERIVED | "Unknown means no gain" applied to a cost: an unmeasured contention term is not modeled, and the margin rule carries the safety. |
| Static tie-break prefers `:process` | -- | SOURCED | TRX50 cold-11 (`paper_benchmarks_trx50_cold11`): the one-heap threads route never beat the pool by more than 11% and lost to it by 5-60% at W >= 16 on P3/P4. |

## The measurements the unit tests encode

`test/unit/parallel/predictive_planner_tests.jl` checks the planner's decision
at seven points from the TRX50 cold-11 run. They are used as the expected
DIRECTION of each decision -- which plan must win, which must not be chosen --
and never as a target the model is fitted to. Medians in seconds.

| Point | Shape | Measured | Required decision |
|---|---|---|---|
| P3 `independent_1sat_1hr` | n=256, W_p=2, T=2 | R6 mixed 3.576, threads 4.735, process 5.072 | mixed, `L=1` |
| P3 `independent_1sat_1hr` | n=256, (32, 32) | R6 mixed 1.069, process 1.124, threads 2.743 | `:process`; never `:threads` |
| P4 `montecarlo_heavy_aerobraking` | n=32, (32, 32) | process 1.678, threads 1.856, R6 1.924 | `:process` at `L=0` |
| P4 `montecarlo_heavy_aerobraking` | n=32, (8, 8) | R6 mixed 2.757, process 3.781 | mixed, `L>0` |
| P5 `mcgrid_8sat_16mc` | n=16, (4, 8) | threads 1.384, R6 mixed 4+7 1.500 | the static-equivalent plan |
| P5 `mcgrid_16sat_8mc` | n=8, (1, 32) | threads 1.896, R6 2.104 | `:threads` at `W=8`, budget 1 |
| P5 `mcgrid_8sat_16mc` | n=16, (1, 32) | threads 1.016; R6 threw in 10 of 11 | `:threads` at `W=16`, budget 1 |
| P5 `mcgrid_8sat_16mc` | n=16, (32, 1) | process 0.728, threads 1.467 | `:process` at `L=0` |

## Tracing

`SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1` prints, per campaign:

```
[predictive] shape n=... threads=... pool=... local_cap=... candidates=[...] constants=loaded|absent margin=... guard_factor=... local_slots_max=...
[predictive]   candidate <route>@w<W>+l<L> makespan=... consumers=... s_worker=... s_heap=... [static]
[predictive] chosen <plan> reason=... gain=...
[predictive] guard round1=...s worker_mean=...ms/<n> local_mean=...ms/<n> predicted_ratio=... observed/predicted=... verdict=... replan=... local_slots=A->B
```

`reason` is one of `static_equivalent_best` (the static plan was best
outright), `predicted_gain` (a deviation cleared the margin),
`margin_not_met` (a deviation was predicted but too small to take), or
`no_static_candidate`.

## What a v2 would add

- An inner-speedup measurement, so the inner thread budget becomes a decision
  and the narrower widths become worth enumerating.
- A remote round-trip term with a measured value, for workloads whose samples
  are short enough for it to matter.
- A per-sample cost spread, for campaigns whose samples are known to differ
  (an aerobraking grid whose corners run far longer than its center).
