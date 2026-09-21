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

A plan that deviates (`L > 0`) is watched from inside its own dispatch. There
is no second round and no barrier: the guard rides along on the dispatcher's
completion hook, decides once, and acts by closing consumers that then drop out
at their next take.

That is a change from the first version, which dispatched one round, stopped,
decided, and dispatched the remainder. The barrier was not the sample time --
it was a second process dispatch, and a dispatch's fixed cost (pool wake-up,
the `CachingPool` closure re-serialization that `Distributed.clear!` forces at
the end of every dispatch, the workers' collection) is paid per dispatch, not
per sample. Measured on this repo's workstation: a guard round of eleven
samples cost 0.24-0.36 s against 0.24-0.26 s for the remaining fifty-three.

**Admission.** `DispatchAdmission` carries one flag per pool worker and one per
local slot, which that consumer checks BEFORE taking its next job -- before,
because a consumer that took a job and then withdrew would drop that sample on
the floor. Closing is one-way and never interrupts work in flight: the consumer
finishes the sample it holds, declines the next, and leaves. Closing every
local slot is "reduce L to zero"; reducing L to `k > 0` closes the `L - k`
highest-numbered slots, so which consumers stop is deterministic; closing the
pool class leaves the local slots to finish the campaign. Flags can only be
closed, never reopened, which is what keeps the guard from becoming an
optimizer that widens a running campaign.

**When it decides.** On the first completion at which EVERY consumer has
reported at least one sample, and at least `consumers` samples are done. Per
consumer, not per class and not a count of completions: the local slots run in
this process and return far sooner than a pool worker, and on the workstation's
P3 shape twenty-four local samples landed before the first worker sample did.
A rule counting completions decided there, comparing a twenty-four-sample mean
against a one-sample mean and calling it a round. The price of the stricter
rule is that the decision lands later, and on a short campaign may not land at
all; that is the honest answer, since a guard with no fair observation of a
class has nothing to say about it. If the campaign drains before that, the verdict
is computed anyway and traced as `never_decided`; there is nothing left to act
on. The observation runs under a lock on the completing consumer's task, and
after the decision every later completion early-outs on one atomic load.

Two quantities come out of the observation, and they answer two different
questions.

**Work (`mean_s`)** is the class's mean `elapsed_s`: what happened INSIDE the
sample, timed where the sample ran. Its ratio between the classes tests the
plan's own `s_heap / s_w` -- the contention the model predicted. The pool
workers are the uncontended reference, so observed over predicted is directly
comparable and `1.0` is the model being right.

**Occupancy (`occupancy_s`)** is the class's mean wall from that consumer
TAKING the job (`take_sink`) to that sample's result being in hand
(`finished_ns`). That brackets one consumer's cost for one
sample: the work plus that sample's own share of the `remotecall_fetch` round
trip, the closure's serialization, and the scheduling around it. Measured from
a dispatch-wide start instead it carries the pool's shared setup and the
queueing of whoever was served earlier, and reads 965 ms for a 38 ms sample. This is what the `remote_overhead`
constant stands in for a priori, and it is the measurement that says which
SIDE of the machine the campaign belongs on.

The second quantity exists because the first cannot see the thing that
actually decided the P3 point. Measured on this repo's workstation,
`independent_1sat_1hr` at 64 samples of ~38 ms: the two classes' `elapsed_s`
differ by 21%, well inside any sane guard factor, while the campaign on the
pool takes 1.7x the campaign on threads. Nothing inside a sample's own timing
can reveal a cost paid around it.

So the guard has two directions.

| Observation | Verdict |
|---|---|
| any sample failed | drop the local slots to zero (`:process@L=0`) |
| worker occupancy > `guard_factor` x local occupancy, **and** the pinned threads plan predicted (from the observed occupancies) to finish the remainder sooner than continuing | move the remainder to `:threads` at `min(remaining, T)`, budget 1, pool idle -- **only when `route_switch` is on, which it is not by default**; otherwise `route_switch_disabled` and the plan stands |
| worker occupancy past the factor but threads no faster | keep the plan (`threads_no_better`) |
| work ratio > `guard_factor` | halve the local slots |
| work ratio > 2 x `guard_factor` | drop the local slots to zero |
| otherwise | keep the plan |

### The route switch is built, measured, and off

`SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH` defaults to `false`. The verdict is
still computed and traced with it off, so the evidence is readable without
being acted on; only the move is suppressed.

It is off because it measured worse, and the reason is instructive. On this
repo's workstation, `independent_1sat_1hr` at 64 samples over 8 workers plus 3
local slots, five repeats with the switch ON: the campaigns that switched ran
0.853 and 0.983 s, and the campaigns that did not ran 0.398 and 0.514 s. The
verdict alternated between repeats, and the first round's process dispatch
alternated with it -- 0.736, 0.052, 0.499, 0.052, 0.376 s for the same eleven
samples.

The cheap first rounds are exactly the campaigns whose remainder stayed on the
pool. Switching to the threads route leaves the pool idle for the rest of the
campaign, so the NEXT campaign's first round pays to wake it again -- the
`CachingPool` closure re-serialization that `Distributed.clear!` forces at the
end of every dispatch, and the workers' post-campaign collection landing at the
start of the next one -- and the guard reads that cost as evidence that it
should switch again. The mechanism manufactures its own evidence.

The underlying premise did not survive the measurement either. What is slow at
this point is `:process@L=0` (0.722 s against the threads route's 0.427 s); the
MIXED plan the planner actually chose runs 0.398-0.514 s when left alone, which
is competitive with the pinned threads route. The campaign was never on the
wrong side of the machine; only the pure-pool static plan was.

What a working version needs is a measurement that is not taken in the round
that pays the dispatch's one-time costs. That is the same lesson
`steady_per_sample_s` already encodes for the bandit, and it is the v2 item
below.

Both directions end on a static-equivalent plan -- `:process@L=0` on one side,
the pinned threads plan on the other -- and `predictive_replan` enforces that:
it accepts a reduced local-slot count on the same route, or a move to
`:threads` with no local slots, and throws on anything else, so a future caller
cannot quietly widen through it. The guard never widens, never invents a plan
the planner would not have enumerated, and never moves to a route
`outer_route_candidates` did not offer (native GRAM point density withdraws the
threads route; a single-threaded coordinator has none to move to).

When both directions qualify, the route change wins: being on the wrong side of
the machine costs more than holding too many local slots on the right one.

The threads plan is priced from the observed LOCAL class, because a local slot
and a threads-route task are the same thing -- a sample on this process's heap.
Widening from `L` slots to `W` tasks is charged the USL ratio
`s_heap(W) / s_heap(L)` when the machine has constants, and nothing when it does
not; that flat case is ASSUMED and optimistic for the threads plan, and the
ratio is clamped at 1 so widening is never predicted to make a sample cheaper.

This is the direction the machine-specific part of the decision lives in. The
worker-to-local occupancy ratio is about 1.6 on this workstation at 8 workers,
1.07-1.11 on the TRX50 at W <= 8, and inverts there at W >= 16 -- one number,
measured in the campaign's own first round, in place of a constant that would
have to be right on every machine.

A static-equivalent plan is dispatched without a guard at all: the only thing
it could close is a local slot, and there are none. It is a plain dispatch,
which is what keeps it comparable with the pinned static route it is equivalent
to.

`MonteCarloResult` reports what the campaign ended up running -- the consumers
still admitted when the queue drained: `route`, `threads` and `local_slots`.
The trace carries the chosen plan, the observation, the verdict, how many
consumers of each class were closed, and the plan it ended on.

Because there is only ever one dispatch, its end IS the end of the campaign,
so the process dispatch's end-of-campaign hooks (the fire-and-forget `GC.gc` on
each pool worker and the coordinator GC debt the local slots leave) fire where
they always did. The `mid_campaign` suppression the two-round version needed is
gone with the round it existed for.

**The width the switch can reach.** Closing the pool class leaves the `L` local
slots consuming the queue, so the threads plan reachable without a barrier is
at width `L` -- not `min(remaining, T)`, which would mean starting consumers
mid-dispatch, i.e. the widening the guard is not allowed to do. The comparison
is priced at `L` for that reason. It makes the switch much harder to trigger,
and correctly so: three local slots do not beat eight pool workers plus those
same three slots, and the earlier version fired only because it compared
against a plan it could not run.

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
| `SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD` (`o_r`) | `0.0` | ASSUMED, and measurably wrong on short samples | The assumption was that a `remotecall_fetch` round trip is milliseconds against samples of 0.03-1 s, so the term sits below any noise floor the planner could measure. The reduced-scale run below refutes that at the bottom of the stated domain: on `independent_1sat_1hr`, 64 samples of ~38 ms, 8 pool workers against 8 coordinator threads, the pool campaign ran 0.722 s against the threads route's 0.427 s for the same per-sample work -- a worker class roughly 1.6x a coordinator thread, not 1.0x. The default stays zero because one point on one machine is not a value to hard-code, and it now has a better answer than a constant: the guard measures the same quantity in the campaign's own first round (see The guard) and moves the remainder to the threads route when the pool is not paying for itself. The field remains for a machine that has measured its own and wants the FIRST round planned correctly too; declaring 0.6 on that point moves the planner to the threads route from the start and to parity with it (0.409 s). |
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

## Reduced-scale measurements on this repo's workstation

12 physical / 24 logical cores, Julia at `--threads=8`, the
`parallelization_performance` single-case worker, 3 repeats, `--warmup=1`,
uncalibrated (no machine constants present). Medians in seconds; repeat 1 is
the cold one in every column. Identical-code noise on this harness is about 2%
on the median and 6% at p90, so differences below those are not differences.

| Case | Mode | Repeats | Median | Plan chosen |
|---|---|---|---|---|
| `independent_1sat_1hr`, n=64, 8 workers | `predictive` | 3.454, 0.628, 0.767 | 0.767 | `process@w8+l3` |
| `independent_1sat_1hr`, n=64, 8 workers | `outer_process` | 1.876, 0.626, 0.722 | 0.722 | pinned pool |
| `independent_1sat_1hr`, n=64, 8 workers | `outer_threads` | 0.559, 0.427, 0.407 | 0.427 | pinned threads |
| `independent_1sat_1hr`, n=64, 8 workers | `policy_v2` (R6 bandit) | 3.206, 0.395, 0.647 | 0.647 | `process` + 4 local slots, one dispatch |
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, `local_slots_max=0` | 3.715, 0.518, 0.717 | 0.717 | `process@w8+l0`, one dispatch |
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, `o_r=0.6` | 2.864, 0.409, 0.404 | 0.409 | `threads@w8` |
| `mcgrid_8sat_16mc`, n=16, 1 worker | `predictive` | 3.734, 1.323, 1.341 | 1.341 | `threads@w8` |
| `mcgrid_8sat_16mc`, n=16, 1 worker | `outer_threads` | 1.590, 1.427, 1.427 | 1.427 | pinned threads |

With the guard's second direction added (5 repeats, `independent_1sat_1hr`;
3 repeats, `mcgrid_8sat_16mc`):

| Case | Mode | Repeats | Median |
|---|---|---|---|
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, route switch OFF (shipped) | 3.513, 0.631, 0.741, 0.703, 0.715 | 0.715 |
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, route switch ON | 3.674, 0.514, 0.983, 0.398, 0.853 | 0.853 |
| `mcgrid_8sat_16mc`, n=16, 1 worker | `predictive` | 3.822, 1.381, 1.395 | 1.395 |

The shipped default is stable across repeats (0.631-0.741 after the cold one)
and lands on the pinned pool's 0.722. P5 is unchanged, as it must be: no pool
is affordable there, so the plan is static and the guard never arms.

Both tables above were measured with the two-round guard, and are superseded by
the single-dispatch numbers below; they are kept because they are what the
barrier cost.

With the barrier-free guard (5 repeats on `independent_1sat_1hr`, 3 on
`mcgrid_8sat_16mc`; repeat 1 is cold in every row):

| Case | Mode | Repeats | Median |
|---|---|---|---|
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, defaults | 3.293, 0.589, 0.642, 0.663, 0.644 | **0.644** |
| `independent_1sat_1hr`, n=64, 8 workers | `predictive`, route switch ON | 3.629, 0.538, 0.649, 0.675, 0.890 | 0.675 |
| `mcgrid_8sat_16mc`, n=16, 1 worker | `predictive` | 3.870, 1.407, 1.380 | 1.407 |

Removing the barrier took the P3 point from 0.767 (two rounds, WS7a-1) through
0.715 (two rounds without the mid-campaign GC hooks) to 0.644, which is where
R6's bandit sits (0.647) and inside the spread of the pinned pool (0.722). The
whole campaign is now one dispatch of 0.46-0.55 s after the cold repeat.

The route switch, on, fired once in five repeats, and that repeat was the
slowest of them (0.890 against 0.538-0.675). It remains off.

**A caveat that matters more than the ordering.** Repeat-to-repeat spread on
this case is 0.398-0.756 for identical plans and identical code -- far wider
than the ~2% median / ~6% p90 identical-code noise quoted for this harness
elsewhere. Nothing here separates 0.644 from 0.647 or from 0.675, and no
single-repeat comparison on this case means anything.

Three things to read off it.

The static-equivalent plan really is the pinned static route: forced to it
(`local_slots_max=0`) the planner runs 0.717 against the pinned pool's 0.722,
which is parity to within the median noise. That is the property the margin
rule exists to protect, and it holds. On the P5 shape, where no pool is affordable, the
planner returns the static plan and runs at parity with it (1.341 against
1.427, inside the noise floor). On the P3 shape it matches the pinned pool
(0.767 against 0.722, at the edge of the noise floor) but both are behind the
pinned threads route, because with `o_r = 0` the model cannot see that a pool
worker costs more than a thread on a 38 ms sample. That is a failure of ONE
ASSUMED CONSTANT, not of the model: declared at the measured value the same
planner picks the threads route and reaches 0.409.

Third, the guard round is not free. R6's bandit picks the same route with four
local slots and dispatches it once (0.647); the planner picks three local slots
and dispatches twice (0.767). The gap is partly the slot count and partly the
barrier -- every consumer waits on the round's straggler before the remainder
goes out -- and on a campaign whose whole wall is under a second the barrier is
a visible fraction of it. The clean fix is not to split at all: the class means
the guard needs are available from the `finished_ns` stamps and the class tags
of the first completions of a SINGLE dispatch, and a re-plan could then change
what the still-idle consumers do rather than requiring a synchronization point.
That is a v2 item, not a v1 tuning knob.

## The pool's per-dispatch fixed cost

A pool worker's first sample of a campaign comes home long after a local slot
has finished the same work, and on a campaign whose whole wall is half a second
that is most of it. The P3 trace above shows it plainly: on
`independent_1sat_1hr` at 64 samples, 8 workers and 3 local slots, each pool
worker's FIRST sample occupied 313-370 ms while reporting 39-47 ms of work,
and the median occupancy of every later sample on the same workers was 43.8 ms.
A local slot's first sample occupied what it worked, 43-47 ms. The cost is paid
once per worker per dispatch, which is why it is invisible on a long campaign
and why the same campaign on the TRX50, where it runs 1-9 s, shows the pool
winning.

`SPACEAGORA_POOL_DISPATCH_PROBE=1` attributes it. It runs once per process,
before a dispatch's timed region, on the campaign's own closure and the pool's
own workers (see `_probe_pool_dispatch_cost`). On this repo's workstation
(12 physical / 24 logical cores, `--threads=8`, 8 workers, the P3 point):

| Part | Measured |
|---|---|
| The campaign closure, serialized | 2357 bytes, 0.1 ms |
| Round trip of a named function, idle warm worker | 0.02 ms, worst of five 0.11-0.12 ms |
| Round trip of a fresh anonymous closure | 0.02-0.04 ms |
| First call of the campaign's closure on a worker | 145-291 ms above the sample's own work |
| Second call through the same `CachingPool` | 33-60 ms |
| Call after `clear!` -- a cache miss with the code already compiled | 0.47-0.64 ms |
| Call through a fresh `CachingPool` and a fresh closure | 0.33-0.56 ms |
| Round trip issued right after `remote_do(GC.gc, w)` on all eight | 587-714 ms |
| The same after `GC.gc(false)` | 7.8-12.7 ms |
| Round trip while three CPU-bound tasks run on the coordinator | 0.02-0.06 ms, against 0.02-0.09 ms idle |

The first three rows of the middle block are one-time compilation per worker
per process, not per dispatch: the miss and hit paths through
`Distributed.exec_from_cache` are different methods, so each is JIT-compiled at
its own first use. What a dispatch actually pays to build a fresh
`CachingPool` and ship the closure again is the two rows at 0.3-0.6 ms -- about
4 ms across eight workers, against a 500 ms campaign.

So the per-dispatch pool churn is not what the first sample waits for. The
previous campaign's collection is. `_run_campaign_with_route_env` fires
`remote_do(GC.gc, w)` at every worker when a dispatch ends, to keep a worker's
collection out of the middle of the next round; a full collection on these
workers takes long enough that a round trip issued while it runs waits 587-714
ms, and the next campaign starts well inside that window. The hook moved the
stall from the middle of a round to the front of the next one, where it is
still in the way. A young-generation collection does the same job in 8-13 ms.

`SPACEAGORA_POOL_WORKER_GC` selects which one: `incremental` (the default,
`GC.gc(false)`), `full` (`GC.gc()`, what the hook did before), or `off`. Five
repeats each, `predictive`, the P3 point, repeat 1 cold in every row:

| Collection | Repeats | Median |
|---|---|---|
| `full` | 4.086, 0.487, 0.678, 0.656, 0.661 | 0.661 |
| `incremental` | 4.173, 0.838, 0.395, 0.401, 0.387 | **0.401** |
| `off` | 4.506, 0.806, 0.692, 0.398, 0.420 | 0.692 |

Under `incremental` the per-consumer trace shows no fixed cost left at all: on
the warm repeats every pool worker's first sample occupies what it works
(40-67 ms against a 43-45 ms median for the rest), where under `full` it
occupied 299-407 ms for 37-54 ms of work.

`off` is not the answer, and for exactly the reason the hook was written: two
of its five repeats ran at 0.398 and 0.420 s and two at 0.692 and 0.806 s, and
the trace caught the cause in one of them -- a single sample reporting 509 ms
of work on a worker that had stopped to collect in the middle of the round.
The collection has to happen; it has to be one that ends before the next
campaign starts.

Two limits on that default. It is measured over five repeats of one 64-sample
campaign, not over a long session: a worker that only ever collects its young
generation between campaigns will eventually need a full collection, and when
Julia's own heuristic fires it, it fires mid-campaign -- which is the stall the
`off` row prices. And it is measured on this shape, where each worker allocates
on the order of 80 MB per campaign.

### What the P3 point costs now

Five repeats each, the same session, the same machine, repeat 1 cold in every
row. "Before" is this branch with `SPACEAGORA_POOL_WORKER_GC=full`,
`SPACEAGORA_POOL_DISPATCH_CACHE=0` and `SPACEAGORA_PPC_SAMPLE_FN_REUSE=0`,
which is what the code did before this change; "after" is the shipped
defaults.

| Mode | Before | After |
|---|---|---|
| `outer_process` (pinned pool) | 1.801, 0.649, 0.715, 0.739, 0.728 -> 0.728 | 1.591, 0.358, 0.322, 0.323, 0.327 -> **0.327** |
| `policy_v2` (R6 bandit) | 4.078, 0.415, 0.638, 0.635, 0.637 -> 0.637 | 4.112, 0.396, 0.451, 0.386, 0.389 -> **0.396** |
| `predictive` (R7) | 4.086, 0.487, 0.678, 0.656, 0.661 -> 0.661 | 4.188, 0.459, 0.423, 0.398, 0.401 -> **0.423** |
| `outer_threads` (pinned threads) | not affected: no pool, no collection hook | 0.591, 0.457, 0.450, 0.430, 0.422 -> 0.450 |

Those before medians reproduce the ones measured on this point before the
change (0.722, 0.647, 0.644) to within 3%, and the `outer_threads` column,
which none of the three switches can touch, is the control.

The caveat stated elsewhere in this document -- that repeat-to-repeat spread on
this case ran 0.40-0.76 s for identical plans and identical code, so no
difference inside it means anything -- was itself a symptom. That spread is
what a campaign costs depending on how much of the previous campaign's full
collection it ran inside. With the young-generation collection the warm
repeats of a row span 0.322-0.327 (`outer_process`) and 0.386-0.451
(`policy_v2`), and every before/after difference above is far outside both the
old spread and the new one.

What it changes for the router: on this shape the pool now wins outright. The
pinned pool is 0.327 against the pinned threads route's 0.450, where before it
was 0.728 against 0.427, and both adaptive profiles land between them rather
than behind both.

`SPACEAGORA_POOL_DISPATCH_CACHE` is the other half. `CachingPool` keys its
worker-side cache on the identity of the function it is handed, so a dispatch
that built a fresh pool and a fresh wrapper closure missed on every worker; the
pool now holds both between campaigns, keyed on the campaign function and the
worker set, and drops them when either changes or when the pool is shut down.
On this shape it is worth the 0.3-0.6 ms per worker above and nothing more, and
it is kept because what it removes scales with the size of the closure a
campaign captures, which 2357 bytes is the small end of. It costs one
assumption: that a campaign function's captured state is not mutated between
campaigns that dispatch it, where before the assumption only had to hold within
one campaign. `SPACEAGORA_POOL_DISPATCH_CACHE=0` restores the per-dispatch
pool.

## Tracing

`SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1` prints, per campaign:

```
[predictive] shape n=... threads=... pool=... local_cap=... candidates=[...] constants=loaded|absent margin=... guard_factor=... local_slots_max=...
[predictive]   candidate <route>@w<W>+l<L> makespan=... consumers=... s_worker=... s_heap=... [static]
[predictive] chosen <plan> reason=... gain=...
[predictive] dispatch=...s n=... failures=... decided_after=<k>/<n> samples seen=<a>/<W>w,<b>/<L>l
[predictive] guard worker_mean=...ms/<n> local_mean=...ms/<n> predicted_ratio=... observed/predicted=...
[predictive] guard occupancy worker=...ms local=...ms worker/local=... rest continue=...s threads=...s
[predictive] guard verdict=... closed=<w>w/<l>l plan=<route>@w<W>+l<L>-><route>@w<W>+l<L>
```

An unguarded (static-equivalent) plan prints one `dispatch=...s ... unguarded`
line instead.

The guard's `verdict` is one of `observation_matches`, `local_slots_slower`,
`local_slots_far_slower`, `workers_occupying_more_than_threads`,
`route_switch_disabled`, `threads_no_better`, `sample_failure`,
`not_observed`, `nothing_to_reduce`, or `never_decided` when the campaign
drained before both classes had been seen.

`reason` is one of `static_equivalent_best` (the static plan was best
outright), `predicted_gain` (a deviation cleared the margin),
`margin_not_met` (a deviation was predicted but too small to take), or
`no_static_candidate`.

## What a v2 would add

- An inner-speedup measurement, so the inner thread budget becomes a decision
  and the narrower widths become worth enumerating.
- A remote round-trip term measured BEFORE the dispatch rather than during it,
  so the FIRST plan is right and not only the plan the guard corrects to. A
  cheap pool probe at plan time would supply it.
- Re-measuring the route switch now that its evidence comes from inside the
  dispatch and its width is the one it can reach. Both defects the first
  version had are gone; whether anything is left to gain is unmeasured, which
  is why it is still off by default.
- A per-sample cost spread, for campaigns whose samples are known to differ
  (an aerobraking grid whose corners run far longer than its center).
