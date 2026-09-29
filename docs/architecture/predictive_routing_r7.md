# R7: the predictive campaign planner

R7 replaces the R6 campaign-level learner with a planner. The campaign plan is
computed once, before any sample runs, from the campaign shape and the
machine's calibrated constants; it defaults to the plan a pinned static route
would run; it deviates only when the model predicts a gain larger than a
margin; and a guard after the first round of samples re-plans the remainder if
the observation contradicts the prediction.

R7 changes nothing inside a sample. The engine-level inner policy (the RHS
calibration sweep, the callback width rules) is exactly what R6 has. R7 decides
the campaign-level plan in `_run_campaign_adaptive` -- including, when the
calibration store has measured it, how many threads each sample gets -- and
nothing else.

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

Narrower widths at budget 1 are not enumerated, in either route. For the
process route the workers are one thread each, so a narrower pool only idles
cores. For the threads route a narrower split at budget 1 buys nothing to trade
against the concurrency it gives up; narrower threads-route widths appear only
with a larger inner budget, and only when a measured curve says what that
budget buys (next section).

### Inner thread budgets, when a curve is measured

v1 never assigned an inner thread budget above 1 under an outer split: how much
faster one sample runs on `k` threads is a property of the workload that no
a-priori model here measures, and "unknown means no gain". The visible cost was
at `n < T`, where four 1024-spacecraft samples on 32 threads used four.

The RHS calibration sweep does measure it, for every shape it sweeps: it times
satellite_batch at the full budget and at a half and a quarter of it, and the
flat route at a geometric ladder of widths including 1. The store used to keep
only the winner. It now keeps every candidate's timing (see The calibration
store below), and `SimulationEngine.rhs_inner_speedup_curve(stem)` turns them
into `speedup[b] = t(1) / min(t(w) for w <= b)` for `b` up to the widest width
measured, over every row of the workload's signature stem whatever budget or
outer-split flag it was measured under. No one-thread timing, no curve.

When the workload has a curve, the planner adds plans whose samples run on
`b > 1` threads, every one of them a deviation that has to clear the margin:

| Plan | Consumers | Inner budget |
|---|---|---|
| `:threads` at `W = fld(T, b)`, `b = 2, 4, 8, ...`, and at `W = min(n, T)` | `W` tasks | `fld(T, W)` |
| `:none` | one sample at a time | `T` |
| `:process` at `W_p` workers plus `L` local slots, `L * b <= T - 1` | `W_p + L` | `b` for the local slots, 1 for the workers |

A sample at budget `b` is priced at `1 / speedup(b)` of a one-thread sample,
times whatever heap term the plan carries. The pool's workers are one-thread
processes and are never given more. Without a curve the plan space is exactly
v1's, and the static plan stays at budget 1 whether or not a curve exists, so a
flat curve changes nothing: every `b > 1` plan ties with its budget-1 twin and
the tie goes to the static plan. The budget is declared as the environment
ceiling around the dispatch, as budget 1 always was, and never exceeds the
split's own share (`W b <= T`), which is the cap `capped_inner_thread_budget`
applies.

How the planner finds the curve: it sees an opaque sample function and the
campaign's features, never a sample's configuration. `campaign_route_features(args;
samples)` is where a configuration and the features meet, so it records the
configuration's RHS signature stem under the features' workload signature, and
the planner reads it back. A campaign whose features were not built from a
configuration has no curve.

Two things are ASSUMED and named. The curve is the RHS evaluation's speedup
applied to the whole sample; a sample also spends time outside the RHS (the
implicit solver's linear algebra, callbacks), so the curve is an upper bound.
And it is measured one sample at a time (or under whatever split its row was
swept under), not with `W` samples beside it. The margin carries both.

`SPACEAGORA_PREDICTIVE_INNER_CURVE=0` plans without a curve.

The leash covers these plans too: a plan without local slots steps one rung at
a time along the shape's own ladder of budgets (the distinct budgets its
threads-route and serial plans offer, 1 included); a mixed plan with local
slots at budget `b` steps one slot at a time within that budget, starting from
one slot when the previous plan was at another budget.

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
- **Coordinator local slot, or threads-route task**: `s_heap(k)`, the inverse
  parallel efficiency of work on one heap, at `k = L` local slots or `k = W`
  threads-route tasks. Both use the same term because they are the same thing:
  samples sharing one process's heap and allocator.

`s_heap` has three models, and the question they answer is WHERE the calibrated
term applies, not whether the machine has one.

- `:locals` (default): the term applies to the coordinator's local slots in a
  mixed plan, `s_heap(L)`, and nowhere else. The threads-route static plan is
  priced at `s_heap = 1`.
- `:none`: nowhere. Plans are ranked purely by round count.
- `:usl`: everywhere -- local slots and threads-route tasks alike.

`SPACEAGORA_PREDICTIVE_HEAP_MODEL` selects between them. On an uncalibrated
machine all three coincide. Machine constants are loaded and traced under every
model; one that declines to use them does not pretend they are absent.

The default is `:locals` because that is the split two TRX50 runs measured, and
the next two sections are those measurements: the term is right about
coordinator local slots and wrong about the pinned-threads plan.

### Why the term applies to local slots and not to the threads plan

Two cold 11-repeat TRX50 runs of the same P3-P5 grid, one with the term applied
everywhere and one with it applied nowhere, disagree in opposite directions at
opposite ends of the grid. The rule that survives both is the `:locals` split.

**It is right about local slots.** P3 `independent_1sat_1hr` at 32 workers x 32
threads (job `20260922-063133-1886099`, 11 repeats, 3 warm-ups):

| | Median |
|---|---|
| pinned process | 0.430 s |
| R6, mixed `w32+l31` | 0.552 s |
| R7 with `heap_model=none`, chose `w32+l31`, guard closed every local slot | 0.599 s |
| R7 with the term applied, chose pure process a priori | 0.433 s |

Thirty-one local slots beside thirty-two `@async` feeders thrash the
coordinator, and the guard cannot recover what round one already spent: it
closed every slot and still finished at 1.39x the pinned route, a failure.
Charged in the PLAN, the same shape never opens them. The term has to be in the
plan, not in the recovery.

**It is wrong about the threads plan.** All twelve P5 splits pass with
`heap_model=none`, including the two the everywhere-model lost outright
(2x16 at 0.92 and 4x8 at 0.96 of the pinned route). The reason is in the next
section.

**What the split costs, stated plainly.** Charging the local slots makes the
planner more cautious than R6 at the P4 mid budgets, and two points are worse
for it while still passing:

| Point | pinned process | R6 | R7 |
|---|---|---|---|
| P4 at 8 | 3.800 | 2.563 (`w8+l7`) | 3.071 (`w8` guard-trimmed to `l3`); 3.144 a priori `w8+l4` |
| P4 at 16 | 2.274 | 1.831 (`w16+l15`) | 2.264 (pure process) |
| P4 at 32 | 1.008 | 1.896 | 1.098 (pure process) |
| P3 at 16 | 0.762 | 0.793 (`w16+l15`) | 0.779 (`w16` trimmed to `l6`) |

At P4 at 16 the term forecloses R6's mixed win -- 1.83 against 2.26 -- and at
P4 at 8 it gives four local slots where seven was faster. Both still pass the
criterion against the pinned routes, and P4 at 32 is where the caution pays
(1.098 against R6's 1.896). The trade is deliberate: the two-model disagreement
above is a failure on one side and a 20% shortfall on the other, and a failure
is the thing the margin rule exists to prevent.

### The USL heap model applied to the threads plan, measured and refuted

TRX50 cold 11-repeat run, job `20260921-202331-989854`, tree `197c53c39`,
constants loaded (`usl_alpha_base = 0.156`, `usl_beta_alloc = 0.00605`,
fingerprint `f48f5e44b83aeac4`). R7 passed 34 of 36 points and failed two, both
P5 `mcgrid_8sat_16mc`:

| Point | R7 chose | R7 | Pinned threads |
|---|---|---|---|
| 2 workers x 16 threads | mixed, ended `w2+l10` (guard-trimmed) | 1.524 s | 0.957 s |
| 4 workers x 8 threads | mixed, ended `w4+l2` | 1.846 s | 1.382 s |

At those constants `s_heap(16) = 4.8` and `s_heap(8) = 3.0`. That prices the
pinned-threads static plan at several times one round, so any plan carrying
pool workers wins the ranking however small the pool -- a two-worker pool beat
a sixteen-wide threads plan -- and the local-slot count is suppressed along
with it. The measured threads route at those points is not several times a
round; it is the fastest thing there.

The same run shows the milder form of it everywhere else, as excess caution
rather than a failure:

| Point | R7 | R6 |
|---|---|---|
| P4 at 8 | mixed `w8+l4`, 3.14 s | mixed `w8+l7`, 2.58 s |
| P4 at 16 | pure process, 2.29 s | mixed `w16+l15`, 1.87 s |
| P5 16sat at 4x8 | mixed `w4+l4`, 1.55 s | threads, 1.47 s |

What is refuted is the MAPPING, not the fit. `usl_alpha_base` and
`usl_beta_alloc` are measured honestly by `scripts/calibrate_machine.jl` on the
allocation kernel. Applying them to whole campaign samples claimed that a
sample's contention with its neighbors has the same shape AND a usable
magnitude -- the contract flagged the shape claim as ASSUMED, and the magnitude
is what these numbers refute. With `:none` the same shapes rank as
`threads@16` at 2x16, `threads@8` at 4x8, and mixed with full local slots at
the P3/P4 mid budgets, which is what was measured to be fastest at each.

Turning `:usl` back on wants a **sample-level contention measurement, not a
kernel fit**: how much slower one campaign sample runs with `k` of its own kind
beside it on this heap. That is measurable -- a short probe of the same sample
at two widths would do it -- and it is a different number from the one
`calibrate_machine.jl` produces. Until it exists, the planner ranks by round
count and the margin rule carries the safety.

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
| workers' STEADY occupancy inflated past `guard_factor` against the work they report (the pool is being starved through the coordinator) | trim the local slots to the contention curve's peak |
| local slowdown against a worker past `local_thrash` | trim the local slots to the contention curve's peak |
| otherwise | keep the plan |

Being slower than a pool worker is not by itself a reason to close a local
slot. A slot at 1.5x a worker still adds two thirds of a worker's throughput,
and the rule this replaces -- trim whenever the observed/predicted ratio passed
`guard_factor` -- closed slots that were paying their way: on the TRX50's P4 at
8 threads it cut seven to three, where R6's seven ran 2.563 s against the
pinned pool's 3.800 s. Two things make a slot worth closing, and neither is
"slower": the pool being harmed through the coordinator, and the slots
thrashing among themselves.

The workers' STEADY occupancy is their cost per sample after their first,
which is free of the one-time per-worker dispatch cost that would otherwise
make every dispatch look like a starved one. On a short campaign no steady
worker observation exists when the guard decides -- the decision fires once
every consumer has reported ONCE, and a worker's one sample is its first -- so
the starvation rule simply cannot fire there, which is the conservative way
round.

A trim is **sized from the curve, not by halving**. Throughput with `W` workers
and `L` local slots is `W/c_w + L/c_l(L)`; modeling `c_l(L)` as the observation
scaled along the calibrated shape makes the second term proportional to `L /
s_heap(L)`, which is `usl_speedup` itself, so the maximizing width is
`usl_peak_workers = sqrt((1-alpha)/beta)`. The observation sets the level, the
curve sets the shape, and the argmax depends only on the shape. Two fallbacks,
both ASSUMED and both narrower than the running plan: a curve whose peak is at
or beyond the running width contradicts the observation that fired the guard,
so the width is halved; a machine with no curve has nothing better than halving
either.

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

## No learning state, and what is persisted instead

The predictive path records no bandit feedback and neither loads nor saves
`campaign_route_state_path()`. Nothing anywhere records how well a plan or a
route did. What it does persist is a small file of bounded corrections to the
cost model's shared parameters (see Online corrections below), and it can be
switched off with `SPACEAGORA_CAMPAIGN_CORRECTIONS=off`, which gives back the
original property exactly: two campaigns of the same shape get the same plan
for the same reason, and a campaign shape seen once is planned as well as one
seen a hundred times.

## The final round and the cold pool

### What the archived rows show

The TRX50 gate run `trx50_targeted_cold_20260923_140152` (P3 at 32, P4 at 8
and 16; `outer_process`, R6 and R7; 11 timed repeats after 3 warm-ups, cold
store) has 55 timed campaigns that ran the pure process route: both modes at P3
at 32 and P4 at 16 (R7 chose `w+l0` in every campaign), and `outer_process` at
P4 at 8. For each, the campaign wall in sample units minus the round count the
list schedule predicts is the cost the model does not see. Divided by the
round count it is not a constant:

| Rows | Rounds | Wall excess per round (sample units, median) | Final-round tail (median) |
|---|---|---|---|
| P4 at 16, `outer_process` | 2 | 0.262 | 0.2204 |
| P4 at 16, R7 (`w16+l0`) | 2 | 0.268 | 0.2250 |
| P4 at 8, `outer_process` | 4 | 0.097 | 0.2261 |
| P3 at 32, `outer_process` | 8 | 0.085 | 0.2217 |
| P3 at 32, R7 (`w32+l0`) | 8 | 0.087 | 0.2278 |

The per-round figure falls roughly as one over the round count, which says the
cost is paid once per campaign, not once per round. It is also not a dispatch
or collection cost: on a warm pool a worker's first sample occupies its
consumer for its own work plus a median of 1.7 ms (64 timed campaigns across
the six traced logs, p90 3.6 ms). What it is, is the last round's straggler.
In earlier rounds a consumer that drew a slow sample simply takes fewer of the
rest; in the last round everyone waits for the slowest of the `k` samples in
flight. For sample times with an exponential tail the expected excess of the
slowest of `k` over the mean is `sigma (H_k - 1)`, `H_k` the `k`-th harmonic
number, and dividing each row's excess by `H_k - 1` gives the last column:
0.220 to 0.228 across two workloads whose samples differ seventeen-fold in
length (P3 about 50 ms, P4 about 0.85 s) and across final rounds of 8, 16 and
32. The pooled median over the 55 rows is **0.2261**, and that is the
`round_tail` constant.

The model therefore charges it once, per consumer class (the pool's workers;
the consumers on this process's heap, whether local slots or threads-route
tasks), at `round_tail * (H_k - 1)` on top of the class's list-schedule finish,
`k` the number of its consumers that ran its final round. A final round of one
sample costs nothing extra. It is charged to both classes because it is the
samples' spread, not the pool's; charging only the pool would price every
threads plan a quarter of a sample ahead of the equivalent pool plan, which the
cold-11 run refutes (the pool matched or beat the threads route at W >= 16).
That transfer to the heap side is ASSUMED: the archived rows measure it on the
pool only.

The first campaign of a process meets a cold pool. In each of the six logs'
first campaign, the median over workers of first take plus first occupancy
minus first work was 712.7, 786.5, 793.2, 835.9, 948.7 and 955.6 ms; the
median, **0.8146 s**, is the `pool_startup_s` constant. The planner charges it
as a delay before each pool worker's first sample (the local slots start at
once, so it shifts the first samples onto them), only for workers that have
not run a predictive process dispatch in this process, and only when the
signature's sample time is known (from the corrections file) so that seconds
can be turned into sample units. Otherwise it is not charged.

`scripts/extract_campaign_cost_terms.py` computes both from an archived run
and, with `--write-constants`, writes them into the machine's fingerprinted
constants file as a `[campaign]` table (`round_tail`, `pool_startup_s`,
`source`). `save_machine_constants` carries that table across a
re-calibration. A machine without the table is priced exactly as before: no
tail, no start-up.

### The ranking on the archived rows

Every plan the run measured, priced by the planner's model (the Python mirror
in the script, with the TRX50's `usl_alpha_base = 0.156`, `usl_beta_alloc =
0.00605`) in sample units. Rank 1 is fastest. "Heap observed" replaces the heap
term by the local slots' slowdown the run's own traces measured for that plan
(mean local first work over mean worker first work, median over campaigns).

| Point | Plan | Measured median (s) | Rank | Model v1 | Rank | + tail | Rank | + tail, heap observed | Rank |
|---|---|---|---|---|---|---|---|---|---|
| P3 at 32 | `w32+l0` | 0.438 | 1 | 8.000 | 1 | 8.692 | 1 | 8.692 | 1 |
| P3 at 32 | `w32+l31` | 0.552 | 2 | 11.306 | 2 | 11.991 | 2 | 13.543 (x6.77) | 2 |
| P4 at 16 | `w16+l0` | 2.292 | 2 | 2.000 | 1 | 2.538 | 1 | 2.538 | 2 |
| P4 at 16 | `w16+l15` | 1.843 | 1 | 4.455 | 2 | 4.979 | 2 | 2.057 (x1.53) | 1 |
| P4 at 8 | `w8+l0` | 3.742 | 3 | 4.000 | 2 | 4.388 | 3 | 4.388 | 3 |
| P4 at 8 | `w8+l7` | 2.620 | 1 | 4.380 | 3 | 4.380 | 2 | 3.113 (x1.28) | 1 |
| P4 at 8 | `w8+l4` | 3.002 | 2 | 3.081 | 1 | 3.388 | 1 | 3.388 (x1.29) | 2 |

What this says, plainly:

- The tail term fixes what it is responsible for. The pure pool moves to last
  at P4 at 8, as measured, and at P4 at 16 the tie among `w16+l0` through
  `w16+l6` (all exactly two rounds in v1, which is why R7 took the static
  plan) is broken in favor of the mixed plans: `w16+l6` ranks 4% ahead.
- It does not by itself reproduce the measured order at P4 at 16 or between
  `w8+l4` and `w8+l7`. What stands in the way is the heap term's SHAPE: the
  calibrated curve charges 4.45x at fifteen local slots where the traces
  measured 1.53x, and 2.19x at seven where they measured 1.28x.
- With the heap term at its measured value, the model with the tail ranks all
  three points exactly as measured, and its predicted times land on the
  measurements: `w16+l15` at 2.057 x 0.896 s = 1.84 s against 1.843 s, `w8+l7`
  at 3.113 x 0.85 s = 2.65 s against 2.620 s.
- At P4 at 16 the decision is still the static plan (the best mixed plan's 4%
  is inside the 15% margin). A heap slowdown measured per width, rather than
  one scale on the calibrated curve, is what would move it; see Online
  corrections for what the scale can and cannot do.

## The local-slot heap term

Under `:locals` the local slots were priced with the USL fit that
`calibrate_machine.jl` measures on the allocation kernel. That reuse was
ASSUMED, and the archived traces refute it: at fifteen local slots it charges
4.45x where the traces measured 1.52x, and no correction under the rules below
can close that gap (the prior's share alone keeps it near 2.2x).

The term the planner now charges local slots is fitted from sample-level
measurements: `s(k) = 1 + c k (k - 1)`, a cost per PAIR of samples sharing the
heap, with `c` the `local_heap_slope` field of the constants file's
`[campaign]` table. A linear `1 + c k` was the first candidate and fails: fitted
to the P4 points it prices 31 slots at 2.1x where P3 measured 6.8x, and the
planner would put P3 at 32 on 26 local slots, the failure the P3 row of the
cold P3-P5 run records (1.39x the pinned pool).

The points, all from `trx50_targeted_cold_20260923_140152`, each a first-round
ratio of work measured inside the samples (local-slot mean over pool-worker
mean), median over the 11 timed campaigns:

| k | Point | Measured `s(k)` | How | Fitted |
|---|---|---|---|---|
| 4 | P4 at 8, R7 `w8+l4` | 1.289 | the guard's own `local_slowdown` | 1.072 |
| 7 | P4 at 8, R6 `w8+l7` | 1.247 | R6 dispatch trace, mean local `first_work` over mean worker `first_work` | 1.252 |
| 15 | P4 at 16, R6 `w16+l15` | 1.523 | as above | 2.262 |
| 31 | P3 at 32, R6 `w32+l31` | 6.754 | as above | 6.590 |

`c` = **0.006011**, least squares of `s - 1` on `k (k - 1)` through the origin
over the four points (`scripts/extract_campaign_cost_terms.py`). The script
also inverts each row's blended mean sample time against the pure pool's
(1.179, 1.384, 1.654 and 9.391 at the four points) as a cross-check; the
inversion is ill-posed where few samples reach the local slots (P3 at 31 gives
5.2 or 9.4 depending on the search range) and is not fitted. The two archived
P3-P5 runs named for this (`trx50_ppb_cold_20260922_103231`,
`trx50_ppb_converged_20260922_170434`) carry no dispatch traces, only CSV rows,
and are not used for the fit. Their P4 inversions (1.00 at k = 1 and 3, 1.36 at
7, 1.72 at 15) agree with the traced points; their P3 inversions (1.18 at 3, 2.0
at 7, 7.2 at 15) are several times steeper, which is the limit below.

Unmeasured: every `k` other than 4, 7, 15 and 31 on this machine, and every `k`
on any other machine. A machine whose constants file has no `local_heap_slope`
falls back to the calibrated USL term (`usl_alpha_base`, `usl_beta_alloc`),
and one with no constants at all to no term: both are what the planner did
before, the first of them the conservative side to be wrong on.

Ranking under the refit prior with no corrections (round tail included):

| Point | Plan | Measured (s) | Model | Rank measured / model |
|---|---|---|---|---|
| P3 at 32 | `w32+l0` / `w32+l31` | 0.438 / 0.552 | 8.692 / 13.180 | 1, 2 / 1, 2 |
| P4 at 8 | `w8+l7` / `w8+l4` / `w8+l0` | 2.620 / 3.002 / 3.742 | 3.113 / 3.388 / 4.388 | 1, 2, 3 / 1, 2, 3 |
| P4 at 16 | `w16+l15` / `w16+l0` | 1.843 / 2.292 | 2.786 / 2.538 | 1, 2 / **2, 1** |

R7 now chooses `w8+l7` at P4 at 8 (the measured best; predicted gain 29%),
`w32+l0` at P3 at 32 (the best mixed plan, `l15`, is 12.8% ahead, inside the
margin) and still `w16+l0` at P4 at 16. That last point is what one slope
cannot do: the two workloads' local slots contend differently at the same
width, and the least-squares slope sits between them. A slope fitted to the P4
points alone (0.002688) ranks both P4 points as measured and sends P3 at 32 to
21 local slots; every slope that ranks all three (0.00323 to 0.00483) is
chosen by the ranking, not fitted to a measurement, and is not used. A
workload-dependent term -- the obvious candidate is the coordinator's dispatch
load, which is seventeen times higher for P3's 50 ms samples than for P4's
0.85 s ones -- is the next measurement, not an assumption to make here.

### A coordinator-load term, fitted and not shipped

The coordinator thread that runs the local slots also dispatches and collects
every sample, and that load per unit time is the warm-pool per-dispatch cost
(1.7 ms median, see The final round) over the sample time: `r = 0.0017 /
t_sample`, about 0.002 for P4 and 0.034 for P3. Four two-parameter forms
carrying `r` were fitted on the four traced points above, by least squares on
`s - 1`, with `t_sample` the pure pool's mean sample time at the point's width:

| Form | a | b | Fit at k = 4, 7, 15, 31 (measured 1.289, 1.247, 1.523, 6.754) |
|---|---|---|---|
| `1 + a k(k-1) + b k r` | 0.002133 | 3.544 | 1.054, 1.139, 1.547, 6.757 |
| `1 + a k(k-1) + b k (W+k) r` | 0.002468 | 0.05157 | 1.035, 1.115, 1.563, 6.754 |
| `1 + a k(k-1) + b k W r` | 0.002461 | 0.1017 | 1.036, 1.115, 1.562, 6.754 |
| `1 + a k(k-1) (1 + b r)` | 0.002487 | 0.1078 | 1.032, 1.114, 1.564, 6.754 |

Validated without refitting against the inversions from
`trx50_ppb_cold_20260922_103231` (blended means against the pure pool at the
same width; ill-posed where few samples reach the local slots, which applies
least at P3 k = 15, where 30 of 256 samples did): P4 k = 7 measured 1.36,
predicted 1.11-1.14; P4 k = 15 measured 1.72, predicted 1.55-1.57; P3 k = 3
measured 1.18, predicted 1.04-1.49; P3 k = 7 measured 2.0, predicted
1.30-2.17; P3 k = 15 measured 7.2, predicted 2.37-3.43. Every form
under-predicts P3 at mid widths by two to three times.

With no corrections, every form ranks all three archived points as measured
and chooses `w8+l7` at P4 at 8 and `w16+l15` at P4 at 16 (predicted gain
17.7-18.4%). None keeps P3 at 32 static outside the margin: each prices a
mid-width mixed plan (`l12` to `l17`) 11.5-12.8% ahead of `w32+l0`, so the
static plan survives only because the margin is 15%, and the cold P3-P5 run
measured the mixed plans R7 took at that point at up to twice the pure pool's
time. That fails the acceptance test set for the term, so no coordinator-load
form ships and the one-slope fit above stands.

`heap_scale` now scales the term's EXCESS, `s = 1 + scale (s_model - 1)`, so
for the pairwise term it is a scale on the fitted slope; the guard's
observation maps to it as `(s_observed - 1) / (s_model - 1)`. The P4-at-8
guard observations sit on the fitted curve at seven slots (1.247 against
1.252), so they move nothing, and P4 at 16 never runs a mixed plan to observe:
the corrections cannot repair that point either, by design.

## Online corrections

Every campaign the guard already measures what the model predicted. The
corrections file keeps a bounded running correction to the three parameters
those measurements bear on, and the next campaign is priced with them:

| Parameter | Prior | Observed from |
|---|---|---|
| `heap_scale` (machine-wide; multiplies the heap term wherever the heap model charges one) | 1 | a guarded mixed plan whose heap model charged a term: the scale that would have predicted the guard's observed local-over-worker slowdown |
| `round_tail` (machine-wide) | `[campaign] round_tail` in the constants file; not corrected, and zero, when the file has none | a pure process campaign on a warm pool with no failures: `(wall / sample - R) / (H_k - 1)` |
| `sample_time_s` (per workload signature) | none | the pool workers' mean sample work; used only to turn `pool_startup_s` into sample units |

Rules, all enforced in `predictive_planner.jl` and all ASSUMED values (named
fields, environment-configurable, see the constants table):

- **Bounded step.** A parameter's evidence starts at the prior and moves toward
  each observation by at most `step x |prior|` per campaign (default 5%).
- **The prior keeps a share.** The value the planner uses is
  `prior_share x prior + (1 - prior_share) x evidence` (default share 0.25), so
  however many campaigns agree, the parameter moves at most three quarters of
  the way from its prior. One campaign moves it by at most 3.75% of the prior.
- **Staleness.** A parameter with no observation in `stale_campaigns` campaigns
  (default 20) reverts to its prior and is dropped from the file. The sample
  time has no prior and simply becomes unknown again.
- **Parameters, never arms.** Every correction is to a parameter all candidates
  share. No per-plan or per-route reward is stored anywhere, and no campaign is
  ever run for the purpose of learning.
- **The decision rule is unchanged.** The margin and the static-equivalent
  default apply to the corrected model exactly as to the uncorrected one.
- **Keyed by machine and code.** The file carries the machine fingerprint and a
  code token; a file written for another machine, another code token or
  another schema is treated as empty. The code token is the RHS calibration
  store's (`_RHS_CALIB_CODE_TOKEN`) when the engine defines one, so a change to
  the RHS execution that invalidates the store invalidates the corrections too.

What one machine-wide heap scale can do is bounded by the measurement above.
Fed the P4-at-8 guard observations (observed over predicted 0.71 to 0.85 on all
eleven repeats), it settles at `0.25 + 0.75 x 0.8 = 0.85` of the calibrated
curve, which is right at four local slots and still far too high at fifteen
(the curve's shape, not its level, is what is wrong there). It does not undo
the P3-at-32 outcome: with the corrections P4 would induce, P3 at 32 still
resolves to `w32+l0` (the best mixed plan is 5.5% ahead, inside the margin),
which is what `test/unit/parallel/predictive_corrections_tests.jl` checks.

### The leash

The only exploration. A NON-static plan may differ from the plan the same
campaign shape (workload signature, sample count, thread count, pool size)
ran last time by at most one local slot, and only in the direction the
corrected model favors: the step is taken only if the model prices it no worse
than the previous plan, otherwise the previous plan is kept. A previous plan on
another route counts as zero local slots, so every walk away from the static
plan starts at one slot. A static-equivalent plan always stands; retreating to
it is never exploration. The first campaign of a shape (no previous plan) is
bounded by the margin rule alone, as before.

The regret this admits is one local slot's share of the campaign per campaign:
the plan that runs is at most one slot away from a plan that has already run
on this shape, and it moves only toward the model's choice. The plan recorded
for a shape is the one the campaign ended on, after any guard trim.

### The file, and how to reset it

`output/parallel_policy_state/campaign_corrections_<fingerprint>.toml`, beside
the constants file (`SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH` overrides). Written
atomically after each campaign whose shape has a pool; shapes without a pool
never touch it. `SPACEAGORA_CAMPAIGN_CORRECTIONS` is `on` (the default), `read`
(use the file, never write it) or `off` (neither). To reset: delete the file.
To run without it: `SPACEAGORA_CAMPAIGN_CORRECTIONS=off`.

Cold-store benchmark runs must not write it, and the snapshot a cold run
starts from must not contain it. Emptying the calibration store (the harness's
`cold_store`) moves the file away with the rest of
`output/parallel_policy_state`; a cold arm then also needs
`SPACEAGORA_CAMPAIGN_CORRECTIONS=off` in its environment, so that its own
campaigns do not write one as they go. A run that wants the corrections sets
`on` explicitly.

### The dynamics are untouched

The corrections, the tail and the leash change only which plan runs. A plan
decides which consumer runs which sample and nothing inside the sample, so it
cannot move a bit of any trajectory. Measured rather than asserted:
`benchmarks/studies/predictive_cost_terms/plan_parity.jl` runs one
`mcgrid_8sat_16mc` campaign (16 samples of 8 spacecraft, 1 h, the harness's
predictive-mode environment) through the planner's own dispatch under several
plans and writes every sample's full state history, in sample order, raw
Float64. This repo's workstation, `--threads=8`, 2 pool workers at one thread:

| Tree | Plan | Saved steps | Bytes | `cmp` against tip `threads@w8` |
|---|---|---|---|---|
| tip | `threads@w8` | 7936 | 4126720 | -- |
| tip | `process@w2+l3` (two worker processes, three local slots) | 7936 | 4126720 | identical |
| tip | `none` (one consumer) | 7936 | 4126720 | identical |
| base `adb343566` | `threads@w8` | 7936 | 4126720 | identical |
| base `adb343566` | `none` | 7936 | 4126720 | identical |

## The calibration store: per-candidate timings and a code token

The RHS calibration store (`output/parallel_policy_state/rhs_calibration_<machine>.toml`)
is at schema 2. Each `[[calibrations]]` row is what it was -- `mode`,
`allotment`, `scheduler`, `elapsed_mean_ns`, `solve_ns`, `sweep_ns`,
`honoured_ns` and the vote counts -- plus a `candidates` table: one entry per
candidate the sweep measured, with `mode`, `allotment`, `scheduler`, `width`
(the threads the plan used at the budget it was measured under), `ns` (its
score from the last round it took part in, its most precise reading), `reps`
and `round`, and the two entries of the final paired round (`round = 0`): the
swept winner and the plan the runtime heuristic would have run, marked
`default = true`, so a reader can see the margin the verdict was formed
against. Schema-1 rows load as before and simply carry no timings; readers of
schema 1 ignore the new key. A verdict stored without a sweep keeps the row's
previous timings, since they describe the shape rather than the verdict.

Every signature now ends in `code=<_RHS_CALIB_CODE_TOKEN>`, a constant in
`rhs_calibration.jl` bumped by hand whenever the RHS execution changes. A store
written by older code no longer matches anything and every shape is treated
as cold: swept again rather than trusted. The converged P1/P5 pass on
2026-09-22 is why: it replayed a 4096-spacecraft verdict from a 2026-09-18
store, formed before the harmonics kernel changes, and ran 18-51% behind. The
token is set to `2026-09-23`, which no store written before it carries --
including the TRX50 snapshots of 2026-09-18 and 2026-09-22. Old rows stay in
the file, unmatched. The campaign corrections file is keyed by the same token.

### Trajectories do not depend on the inner budget

The same `plan_parity.jl` dump, `mcgrid_8sat_16mc`, 16 samples, `--threads=8`,
the predictive-mode environment (RHS calibration on, so samples at `b > 1`
sweep and run the plans their sweeps choose):

| Plan | Inner budget | Bytes | `cmp` against `threads@w8`, budget 1 |
|---|---|---|---|
| `threads@w8` | 1 | 4126720 | -- |
| `threads@w4` | 2 | 4126720 | identical |
| `threads@w2` | 4 | 4126720 | identical |
| `none` (`threads@w1`) | 8 | 4126720 | identical |
| `process@w2+l3`, local slots at 2 | 2 | 4126720 | identical |

The store that run left behind shows what the sweeps chose between (the
satellite_batch route at widths 2-8 and the flat route at widths 1-2) and gives
this shape a flat curve (`speedup = 1` at every width to 8), so the planner
keeps budget 1 here; the dumps at the wider budgets were run by forcing the
plan.

## Constants

Every number the planner uses, and what kind of number it is.

| Name | Default | Class | Source |
|---|---|---|---|
| `SPACEAGORA_PREDICTIVE_MARGIN` | `0.15` | ASSUMED | No measurement of this model's error exists, so no margin can be derived from one. 0.15 is a round number chosen to sit above the ~6% p90 identical-code noise floor of the reduced-scale harness and below the 25-40% mixed-dispatch wins R6 measured at the P3/P4 mid budgets, i.e. large enough to refuse noise and small enough to keep the wins that motivated mixed dispatch. Tune it, do not trust it. |
| `SPACEAGORA_PREDICTIVE_GUARD_FACTOR` | `1.5` | ASSUMED | Same standing. A local slot running 1.5x slower than predicted relative to a pool worker is outside anything the round-count model explains; the second threshold at 2x that is the "stop entirely" case. |
| `SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX` | `Threads.nthreads() - 1` | DERIVED | `ParallelProfiles.mixed_local_slots` keeps thread 1 free for the `@async` feeders that keep the pool supplied, so `T - 1` is the most slots the coordinator can offer under R6's own measured practice. |
| `SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD` (`o_r`) | `0.0` | ASSUMED, and measurably wrong on short samples | The assumption was that a `remotecall_fetch` round trip is milliseconds against samples of 0.03-1 s, so the term sits below any noise floor the planner could measure. The reduced-scale run below refutes that at the bottom of the stated domain: on `independent_1sat_1hr`, 64 samples of ~38 ms, 8 pool workers against 8 coordinator threads, the pool campaign ran 0.722 s against the threads route's 0.427 s for the same per-sample work -- a worker class roughly 1.6x a coordinator thread, not 1.0x. The default stays zero because one point on one machine is not a value to hard-code, and it now has a better answer than a constant: the guard measures the same quantity in the campaign's own first round (see The guard) and moves the remainder to the threads route when the pool is not paying for itself. The field remains for a machine that has measured its own and wants the FIRST round planned correctly too; declaring 0.6 on that point moves the planner to the threads route from the start and to parity with it (0.409 s). |
| `alpha` in `s_heap` | `MachineConstants.usl_alpha_base` | SOURCED (fit), **REFUTED (mapping)** | Fitted per machine by `scripts/calibrate_machine.jl` on the allocation kernel. Applying it to campaign samples claimed the same shape and a usable magnitude; the magnitude is refuted by two TRX50 runs. Consulted under `heap_model = :usl`, and under `:locals` only as the fallback when the constants file has no `local_heap_slope`. |
| `beta` in `s_heap` | `MachineConstants.usl_beta_alloc` | SOURCED (fit), **REFUTED (mapping)** | As above. |
| `s_heap` with no constants | `1.0` | DERIVED | "Unknown means no gain" applied to a cost: an unmeasured contention term is not modeled, and the margin rule carries the safety. |
| `SPACEAGORA_PREDICTIVE_LOCAL_THRASH` | `3.0` | ASSUMED | The observed local-slot slowdown, relative to a pool worker, past which the guard trims even though the pool is unharmed. Well above `guard_factor` because a slot at 1.5x a worker still adds two thirds of a worker's throughput; the old rule trimmed there and cost the P4-at-8 point. Nothing measures where the real threshold is. |
| `SPACEAGORA_PREDICTIVE_HEAP_MODEL` | `locals` | SOURCED (the split), MEASURED-OFF (`:usl`) | Where the calibrated term applies. `:locals` is the split the two TRX50 runs measured -- right about local slots at P3 at 32, wrong about the threads plan at every P5 split. `:usl` | `:usl` charges the alloc-kernel USL fit as whole-sample contention. Measured on the TRX50 with its own constants it cost two of thirty-six points outright and made the planner systematically over-cautious elsewhere; see The USL heap model, measured and refuted. `:none` reproduces the measured winners at those points. |
| Static tie-break prefers `:process` | -- | SOURCED | TRX50 cold-11 (`paper_benchmarks_trx50_cold11`): the one-heap threads route never beat the pool by more than 11% and lost to it by 5-60% at W >= 16 on P3/P4. |
| `round_tail` (`[campaign]` table of the constants file) | absent: no term | SOURCED | `trx50_targeted_cold_20260923_140152`: median over the 55 timed pure-process rows (`p3_32_outer_process`, `p3_32_predictive`, `p4_16_outer_process`, `p4_16_predictive`, `p4_8_outer_process`; 11 each) of `(wall / mean_sample - R) / (H_k - 1)` = 0.2261 sample units; per-file medians 0.2204-0.2278. Extracted by `scripts/extract_campaign_cost_terms.py`. |
| Tail shape `H_k - 1`, final round only | -- | DERIVED | Expected excess of the slowest of `k` exponential-tailed samples over their mean. Supported, not assumed blindly: it collapses the per-file excess to 0.220-0.228 across final rounds of 8, 16 and 32, where a per-round charge spreads 0.085-0.268. |
| Tail charged to heap-side classes too | -- | ASSUMED | The rows measure it on the pool only. Charging the pool alone would price every threads plan ahead of the equivalent pool plan, which the cold-11 run refutes at W >= 16. |
| `pool_startup_s` (`[campaign]` table) | absent: no term | SOURCED | Same run: for the first campaign of each of the six traced logs, the median over pool workers of first take + first occupancy - first work (712.7-955.6 ms); median of the six = 0.8146 s. Warm-pool campaigns: 1.7 ms median, so it is charged to cold pools only. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STEP` | `0.05` | ASSUMED | Largest per-campaign move of a correction's evidence, as a fraction of the prior. Nothing measures what it should be; small enough that one campaign moves a parameter by at most 3.75% of its prior after the prior's share. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_PRIOR_SHARE` | `0.25` | ASSUMED | The prior's minimum weight in the value used. At 0.25 a parameter can travel at most three quarters of the way to what campaigns observe. |
| `SPACEAGORA_CAMPAIGN_CORRECTION_STALE_CAMPAIGNS` | `20` | ASSUMED | Campaigns without an observation after which a correction reverts to its prior. Longer than a harness point's 14 campaigns, so one benchmark point does not expire what the previous one learned. |
| Leash of one local slot | -- | ASSUMED | The exploration bound: the regret it admits is one slot's share of a campaign per campaign. |
| `SPACEAGORA_CAMPAIGN_CORRECTIONS` | `on` | -- | `on`, `read` or `off`; see Online corrections. |
| `_RHS_CALIB_CODE_TOKEN` | `2026-09-23` | -- | Not a measurement: a version label. Bumped by hand on any change to how an RHS evaluation executes; its only requirement is that no earlier store carries it. |
| Inner-speedup curve applied to the whole sample | -- | ASSUMED | The sweep times the RHS alone; the solver's own work does not speed up with the inner budget. An upper bound, carried by the margin. |
| Curve measured without `W` samples beside it | -- | ASSUMED | A row's timings come from one solve's sweep; `W` concurrent samples at `b` threads share memory bandwidth the sweep did not. |
| `SPACEAGORA_PREDICTIVE_INNER_CURVE` | `1` | -- | `0` plans without an inner-speedup curve. |
| `local_heap_slope` (`[campaign]` table) | absent: the USL term | SOURCED (fit), form DERIVED | 0.006011, least squares of `s - 1` on `k (k - 1)` over the four first-round per-class work ratios of `trx50_targeted_cold_20260923_140152` (k = 4, 7, 15, 31; see The local-slot heap term). The pairwise form is chosen because the linear one fails P3 at 32. Ranks P3 at 32 and P4 at 8 as measured, not P4 at 16. |

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
On this shape it is worth the 0.3-0.6 ms per worker above and nothing more.
What it removes scales with the size of the closure a campaign captures, which
2357 bytes is the small end of, but it costs one assumption: that a campaign
function's captured state is not mutated between campaigns that dispatch it,
where before the assumption only had to hold within one campaign. A caller
cannot be expected to know that, so the cache ships off:
`SPACEAGORA_POOL_DISPATCH_CACHE=1` enables it, and the default keeps the
per-dispatch pool.

## Tracing

`SPACEAGORA_CAMPAIGN_DISPATCH_TRACE=1` prints, per campaign:

```
[predictive] shape n=... threads=... pool=... local_cap=... candidates=[...] constants=loaded|absent heap_model=none|usl margin=... guard_factor=... local_slots_max=...
[predictive] terms heap_scale=... round_tail=... startup=... pool_cold=... sample_time_s=...|unknown corrections=on|read|none campaigns=... inner_curve=none|<speedup at 1,2,...>
[predictive]   candidate <route>@w<W>+l<L> [budget=<b>] makespan=... consumers=... s_worker=... s_heap=... [static]
[predictive] chosen <plan> reason=... gain=...
[predictive] leash leashed|leash_held -> <plan>          (only when the leash changed the plan)
[predictive] dispatch=...s n=... failures=... decided_after=<k>/<n> samples seen=<a>/<W>w,<b>/<L>l
[predictive] guard worker_mean=...ms/<n> local_mean=...ms/<n> predicted_ratio=... observed/predicted=...
[predictive] guard occupancy worker=...ms steady=...ms/<n> local=...ms local_slowdown=... worker_degradation=... worker/local=... rest continue=...s threads=...s
[predictive] guard verdict=... closed=<w>w/<l>l plan=<route>@w<W>+l<L>-><route>@w<W>+l<L>
```

An unguarded (static-equivalent) plan prints one `dispatch=...s ... unguarded`
line instead.

A campaign whose shape has a pool, with corrections on or read, ends with
`[predictive] corrections campaigns=... heap_scale=<evidence>|prior
round_tail=<evidence>|prior mode=...`.

The guard's `verdict` is one of `observation_matches`, `local_slots_slower`,
`local_slots_thrashing`, `workers_degraded`,
`workers_occupying_more_than_threads`,
`route_switch_disabled`, `threads_no_better`, `sample_failure`,
`not_observed`, `nothing_to_reduce`, or `never_decided` when the campaign
drained before both classes had been seen.

`reason` is one of `static_equivalent_best` (the static plan was best
outright), `predicted_gain` (a deviation cleared the margin),
`margin_not_met` (a deviation was predicted but too small to take), or
`no_static_candidate`.

## What a v2 would add

- An inner-speedup curve measured with the samples that will run beside it,
  rather than one sweep's view of one sample (see Inner thread budgets).
- A remote round-trip term measured BEFORE the dispatch rather than during it,
  so the FIRST plan is right and not only the plan the guard corrects to. A
  cheap pool probe at plan time would supply it.
- A sample-level contention measurement, which is what would justify turning
  the `:usl` heap model back on.
- Re-measuring the route switch now that its evidence comes from inside the
  dispatch and its width is the one it can reach. Both defects the first
  version had are gone; whether anything is left to gain is unmeasured, which
  is why it is still off by default.
- A per-sample cost spread, for campaigns whose samples are known to differ
  (an aerobraking grid whose corners run far longer than its center). The
  final-round tail covers the random spread of identically distributed
  samples; a systematic spread is a different thing.
- A workload-dependent local-slot term: one pairwise slope ranks P3 at 32 and
  P4 at 8 as measured but not P4 at 16 (see The local-slot heap term).
