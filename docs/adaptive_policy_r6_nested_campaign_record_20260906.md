# R6 on nested Monte Carlo constellations: record of the 2026-09-06 work

Branch `router-eval-expanded-b6`, continuing from
`docs/policy_v2_session_record_20260903.md` (profile R6 / policy V2) and
`docs/adaptive_policy_session_record_20260902.md`. Reference machine: AMD
Ryzen 9 9900X, 12 physical cores / 24 SMT, 60 GB, Julia 1.12.1. Pool workers
are `--threads=1`.

Goal set for this work: R6 at least as fast as the best static route at every
launch point of B15 (within the harness's ~8% noise band, and strictly ahead
where only R6 can use both the process pool and the coordinator's threads),
with no loss on the other router phases (B8–B14). Regression set: B8–B15 under
`--lean-modes`.

> Status: IN PROGRESS. Sections marked TODO are filled as measurements land.

---

## 1. The defect that started this

B15 (`mcgrid_16sat_8mc`, 12 threads, 1 worker, 8 samples): R6 ran **10.33 s
per sample** against **2.70–2.96 s** for the best static route and **4.56 s**
for R5. Same case, same launch point, same outer route (`threads`).

What it was not:

- Not the process route. `outer_backend_actual = threads` for R6 and for the
  winning static route alike.
- Not the cached calibration verdict. The dead run of 2026-09-06 04:5x shows
  R4 and R5 honouring the identical `cache/heuristic` verdict at the identical
  point at 4.6–4.9 s. Re-running R6 with the verdict retained and only the
  budget corrected gave 3.04 s.

What it was — a lie about the inner thread budget:

1. The harness dispatches the 8 samples with `Threads.@threads` under an env
   that sets `SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` and never sets
   `SPACEAGORA_INNER_THREAD_BUDGET` (`execution.jl`, `ppc_run_sample_batch`).
2. Unset, `effective_inner_thread_budget()` returns the whole pool — 12
   (`env_config.jl:153`). The calibration store corroborated it: every
   16-sat `outer=1` signature sat at budget = the process's full thread count
   at each split (12, 6, 4, 3, 2) while 8 samples ran concurrently.
3. R6's `allow_inner_with_outer=true` lifts the outer clamp for the callback
   layer, and V2's width is the static `min(items, budget)` = 12
   (`adaptive_decision.jl:82`), held for the whole solve — the AIMD that would
   have shrunk it was retired in `67c5d7da`.
4. Eight samples × 12-wide callbacks on 12 cores. Superlinear in width.

Measured ladder, advertised budget → per-sample wall (all else identical):

| advertised budget | per-sample | callback threads enabled | RHS verdict |
|---|---|---|---|
| 12 (unset → pool) | 10.33 s | 550 | cache/heuristic |
| 4 | 3.77 s | 550 | cache/satellite_batch (pinned) |
| 2 | 3.04 s | 0 | cache/heuristic |
| 1 | 3.07 s | 0 | none (calibration skipped) |
| R5 at 12 | 4.56 s | 608 | cache/heuristic |

Budget 2 with the same verdict ran 3.04 s: the verdict is irrelevant. Budget 4
with a pinned plan still threaded 550 times: the RHS plan and the callback
width are independent controls, which is why the verdict only ever correlated
by accident (both key on the budget).

Two dispatchers were lying. The production `threads=:auto` path already set
the share (`adaptive_routing.jl`, "instead of oversubscribing"); the
integer-threads `run_monte_carlo(f, seeds; threads=N)` did not, and the paper
harness bypasses `run_monte_carlo` entirely.

## 2. What was attempted, in order

| # | change | result | kept |
|---|---|---|---|
| 1 | `outer_split_env_pairs(worker_count)` around the threaded dispatch in `run_monte_carlo(f, spec)` | integer-threads API now advertises `fld(nthreads, workers)`; user budget wins | yes (`f1b94ba0`) |
| 2 | harness `ppc_mode_env_pairs` emits the same share for the threads backend | harness worker, unmodified env: 3.06 s/sample, 0 callback threads, wall 3.25 s vs static 3.19 s | yes (`f1b94ba0`) |
| 3 | row env recorded from the RESOLVED mode | audit string stops saying `backend=auto` and hiding backend-conditional pairs | yes (`f1b94ba0`) |
| 4 | purge the 15 `outer=1` calibration entries formed under the overstatement | backup `rhs_calibration_*.toml.bak_20260906_outer1_purge` | yes |
| 5 | mixed dispatch: process route fills the coordinator's spare threads | end-to-end at (2 workers, 6 threads), 16sat_8mc: 3.48 s wall vs static 4.05 / 5.85 (R6 before: 10.36) | yes |
| 6 | harness measures adaptive profiles through `run_monte_carlo(threads=:auto)` | required to see 5 at all: the harness's own `pmap` over W one-thread workers can only show W-way concurrency | yes |
| 7 | B15 gains the 1024 spacecraft-hour grid and an `inner_only` control | pending the B15 re-run | yes |
| 8 | `satellite_batch` RHS: Polyester width follows the inner budget (fix A) | see §2.1 | yes |
| 9 | V2 records one policy observation per RHS call, not one per satellite (fix B) | see §2.1 | yes |

### 2.1 The residual after the budget fix, and how it was attributed

With the budget advertised, R6 at (1 worker, 12 threads) was still +17% at 16
spacecraft and +48% at 32 (`mcgrid_32sat_4mc`: 5.5 s wall vs 3.7 s), with
callbacks at width 1 on both sides. Every hypothesis below was measured on that
point through the unmodified harness worker (4 samples, per-sample budget 3):

| ladder | variant | wall | per-sample | verdict |
|---|---|---|---|---|
| A | R6 | 5.25 | 4.86 | baseline |
| A | R6, `RHS_CALIBRATE=off` | 5.34 | 4.93 | calibration is not it |
| A | R6, all five inner layers off | 5.44 | 5.13 | inner layers are not it |
| A | static `outer_threads` | 3.25 | 2.81 | — |
| B | R6 / inner-off / static: GC and allocation | 2731 / 2738 / 2598 MB per sample, GC 5–10% | allocation is not it |
| C | R6, ONE sample, budget 12 | 0.63 | 0.63 | 6.8× faster than static alone (4.30): satellite_batch 12-wide vs a serial flat/1 RHS |
| D | RHS route forced to `satellite_batch`: R6 vs static | 5.06 vs 3.20 | | the ~1.9 s R6-only cost is route-independent |

Two CPU profiles of the 4-sample batch (12 threads, all threads sampled)
separated the two effects:

- Both modes run `_spacecraft_dynamics_dispatch!`'s **Polyester `@batch`** for
  the `satellite_batch` RHS, and that loop sized `minbatch` from
  `Polyester.num_cores()` -- never from the plan's allotment (always 1 for
  that route) or the advertised inner budget. Four concurrent samples each ran
  a 12-wide batch on 12 cores: 4.9 s per sample against 0.63 s for one sample
  alone. Any threaded outer campaign with >= 16 spacecraft per sample was
  oversubscribing by construction, under every profile; this is also the
  ~1.4x "adaptive tax" R4/R5 paid at 16 spacecraft with a pinned
  `satellite_batch/1`. **Fix A**: `_rhs_batch_workers(p) = min(num_cores,
  inner_thread_budget)`, `minbatch = cld(n, workers)` at all seven `@batch`
  sites, and the batch gate refuses a 1-wide batch. With no budget advertised
  the bound is the pool, so a single simulation is unchanged.
- Only R6 showed `record_policy_observation!` (~550 samples) and `locks-mt.jl
  trylock` with **1181 self samples** -- threads spinning on a lock. The
  effector accumulation records a policy observation **per satellite per RHS
  call from inside the Polyester body**, on every worker thread; under V2 the
  solve's context is captured at setup and shared by all of them. The EMA that
  observation maintains is read only by the shipped hint layer, which V2's
  static width rule never consults. **Fix B**: under V2 observe once per RHS
  call (satellite 1, which already feeds the effector cost model). R4/R5 keep
  the per-satellite observation.

Also learned on the way, and worth its own line: the static `outer_threads`
route's single-sample RHS at 32 spacecraft is a serial `flat/1` at 4.30 s
where `satellite_batch` runs 0.63 s -- the heuristic's route choice for the
pinned profiles is leaving a 7x on the table in single-simulation constellation
runs of this size. Not R6's problem; recorded for whoever owns the heuristic.

Instrumentation added to the harness row for this: `sample_gc_time_sum_s`,
`sample_alloc_mb_sum`, and `execution_scope` now names the mixed split
(`adaptive_mixed_w2_l5_sample_batch`).

Considered and not done:

- Defaulting `effective_inner_thread_budget()` to 1 under an active outer
  split with no advertised budget. Safe by construction, but it changes the
  shipped R4/R5 behaviour in any such context and would need its own paired
  measurement; with 1 and 2 in place the only remaining path to it is a
  caller who sets `OUTER_PARALLEL_ACTIVE=1` by hand.
- Threaded sub-campaigns inside pool workers. Workers are `--threads=1` by
  design (the native GRAM lock); the spare capacity at a (W, T) split is the
  coordinator's T threads, not the workers'. Mixed dispatch uses that instead.

## 3. Mixed dispatch — design

TODO: final form after measurement. Design as implemented:

Pool workers are `--threads=1` and run one sample each, so at a split of W
workers on a T-thread coordinator the process route alone uses W cores and
the threads route alone uses T. On this 12-core box the (2, 6) split leaves
either 10 or 6 cores idle, and B15 measured the two static routes at 5.85 s
and 4.05 s there for eight samples that fit in one round at 8 slots.

The process route now runs samples on the coordinator's spare threads from
the same job queue the pool workers consume:

    L = min(T − 1, usable_cores − W)

Thread 1 is kept for the `@async` feeders (sticky to it; a compute-bound
sample there starves them). W + L never exceeds the usable core budget, so
(12, 1) stays pure process and (1, 12) is already threads (one worker
withdraws the process route). Native GRAM point density: at most one local
slot (one lock per process). Under memory-aware routing the local working
sets are charged against the budget (`memory_local_slot_cap`).

The V2 Monte Carlo rule compares `mixed_capacity` = W + L against the thread
budget instead of W alone, so a pool of any size wins the comparison whenever
the route is affordable — the "outer processes for Monte Carlo" prior,
extended to fill the rest of the machine.

Expected slots per B15 split on this machine: (1,12)→threads 8; (2,6)→2+5;
(3,4)→3+3; (4,3)→4+2; (6,2)→6+1; (12,1)→12+0.

### 3.1 Fix C: the residual overheads, and the machine-independent form of the rule

Three residuals were left after §4.2/§4.3 (all measured through the harness):

1. **(12,1) and B12-independent, B8**: R6 ran the same per-sample times as the
   static process route but lost 12–36 % of wall to the race × pre-dispatch
   GC: five dispatches per campaign, each preceded by an unconditional
   `GC.gc()`, on the one route whose per-sample cost does not depend on the
   split width (pool workers are one thread each), so the race could not find
   anything. Fix: the collection is owed only after a dispatch that ran
   samples on the coordinator's heap (`_GC_DEBT`, set by a threaded or mixed
   dispatch, cleared by the collection); the race is gated to the threads
   route; the race result reports `W + L` in-flight so the harness label reads
   `w3_l3` rather than `w0_l3`.
2. **Small-sample ties on the (2,6)/(3,4) splits**: 1.05–1.11 at 4–8 samples,
   every one a launch point where mixed dispatch and the threads route need
   the same number of rounds, `cld(n, W+L) == cld(n, T)`.
3. **The static routes' `satellite_batch` width** (§2.1): the batch kernel ran
   at the whole inner budget and nothing else.

The first version of the tie rule was a threshold fitted to the 36 B15 points
of this box — "tie → threads at 1–2 rounds, → pool at 3" — and it agreed with
the faster measured route at every tie point here. It was dropped before it
ran, on the question of whether the policy would hold on another machine. It
would not: the rule assumes a thread-slot sample and a pool-slot sample cost
the same, which is true at ≤ 12 threads after the budget fix and false on the
64-core TRX50, whose own thread ladder (paper_scenarios S1, 2048 spacecraft)
degrades from 2.69 s at 8–16 threads to 3.30 s at 32 and 6.54 s at 64. On that
box a 2-round tie at (32, 32) would hand 64 samples to 32 shared-heap threads
instead of 63 mostly-isolated slots. The same assumption sat under two more
places: local slots filled every spare thread (`T − 1`, validated only to 5),
and the batch kernel had no width search.

So fix C makes each of those decisions a measurement on the host rather than
a constant from this one:

- **Rounds rule** (`mc_route_rounds`, `default_outer_route`): fewer rounds
  wins outright. At equal rounds — equal capacities included, since W
  isolated workers against T threads on one heap are different per-sample
  costs and not the same route twice — the pool is the *cold* answer, chosen
  because its downside is bounded (a pool sample's cost is fixed by its
  isolation; a thread sample's grows with the coordinator's thread count),
  and the tie is then
  **measured**: `select_outer_route!` spends one campaign on the other
  parallel arm (`OuterRouteTuning.explore_route_ties`,
  `tie_explore_min_campaigns = 1`, reason `explore_tie`) and the existing
  route bandit exploits the faster one from then on. This is the only place
  V2 forces a trial, because a rounds tie bounds what the trial can cost;
  `explore_routes` stays off everywhere else for the reason recorded on it.
- **Local slots** (`_campaign_route_plan`): a slot is a threads-route split
  of the coordinator (its inner budget is its share, exactly as under the
  threads route at that width), so once the threads route has split history
  on the shape — which the tie's threads campaign supplies through the split
  race — the slots are capped at the width the split selector measured best.
  Cold, the cap is absent and the slots fill the spare threads as before.
- **`satellite_batch` width** (`_rhs_plan_candidates`, `_rhs_batch_workers`):
  the calibration trial now carries batch rungs at budget/2 and budget/4
  beside the full-budget arm; a pinned batch plan's `allotment > 1` bounds the
  Polyester width (1 keeps the legacy meaning, the whole budget, so every
  cached verdict written before this stays valid). Two rungs, not the flat
  ladder's seven, for the arm-count reason recorded on the trial.

What this gives up: on the first campaign of a tie shape the policy is not
measured, it is defaulted (to the pool), and on this box that default costs
5–11 % at 1–2-round ties once — the measured B15 points at (2,6)/(3,4) with
4–8 samples. From the second campaign the choice is the measured one. The
harness now reflects that use: one `OuterRouteState` per worker process,
shared by a point's repeats (`_ppc_adaptive_route_state`), so repeat 1 is the
cold answer, repeat 2 the exploration, repeat 3 the exploitation, and the
row's median is the majority. Points remain cold with respect to each other,
and the state path is pinned to a file no run writes.

Unverifiable here: whether the defaults are right on a box that is not this
one. The light set (§4.4) and the full B8–B15 set are what to run on the
TRX50 when it is reachable.

### 3.2 The policy at this point in time, top down (fix C `883ede03`, fixes D and E in the commit that follows)

**0. What R6 is.** `SPACEAGORA_PARALLEL_PROFILE=R6` is R5's adaptive
machinery plus `SPACEAGORA_PARALLEL_POLICY_V2=1`, which flips a set of
`OuterRouteTuning` fields (`mc_route_by_core_budget`, `split_race`,
`gc_before_dispatch`, `memory_aware`, `mixed_dispatch`, `explore_route_ties`)
and the inner policy's V2 behaviours (static callback widths, cached
calibration verdicts honoured, class-gate bypass). Everything below is what
those switches turn on.

**1. The outer route of a Monte Carlo campaign** (`_campaign_route_plan`,
once per campaign). Input: `OuterRouteFeatures` — a probe simulation's shape
(satellite count, effectors, density family, mission time, …) plus
`montecarlo_samples = n`.

- *Capacities.* `budget = T`, the coordinator's threads. Pool workers
  `W = effective_process_workers` = min(usable physical cores, memory cap for
  this workload, `--process-workers`). Local slots
  `L = min(T − 1, usable_cores − W)`, capped by memory (and to 1 under native
  GRAM). `cap = W + L`.
- *Default rule* (`default_outer_route`). If the pool is affordable (≥ 2
  workers, enough samples or mission time), count rounds:
  `rp = cld(n, cap)`, `rt = cld(n, T)`. `rp < rt` → `:process` (mixed
  dispatch); `rp > rt` → `:threads`; a tie — equal capacities included —
  → `:process` as the cold answer, because a pool sample's cost is fixed by
  its isolation while a thread sample's grows with T.
- *Bandit on top* (`select_outer_route!`). Per shape signature, per-route
  statistics (per-sample cost, success rate) persisted in `OuterRouteState`;
  with history the route is picked by lower-confidence-bound UCB. V2 never
  forces exploration, with one exception: on a rounds tie it spends one
  campaign on the other parallel arm (reason `explore_tie`), then exploits
  whichever measured cheaper. Fix E is what makes that comparison fair: the
  arm is credited with its steady per-sample cost, not its cold first
  campaign.

**2. Split width and local slots.**

- *Threads route.* Candidate widths are a geometric ladder from `T/4` to
  `T`. Unseen shape → the in-campaign split race: one warm-up sample, a batch
  at each width, the rest at the fastest per-sample width; every width's
  result is recorded so the next campaign exploits without racing. Seen shape
  → the split bandit picks.
- *Process route.* No race (a pool worker's per-sample cost does not depend
  on W), widest split. Local slots `L` as above, capped at the width the
  threads route measured best on this shape once such history exists — a
  local slot is a threads split of the coordinator.

**3. Dispatch.**

- `:threads` — `run_monte_carlo` with W concurrent samples; each advertises
  `SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` and an inner budget of `fld(T, W)` so
  its inner layer cannot oversubscribe.
- `:process` (mixed) — one job channel; W `@async` feeders each blocking on
  `remotecall_fetch` to a `--threads=1` worker, plus L `Threads.@spawn` local
  consumers under an inner budget of `fld(T, L)`. Pool workers are adopted
  from the harness or spawned and warmed through the sample closure.
- GC: a full collection before a process dispatch only when a threaded or
  mixed dispatch has left heap behind since the last one (`_GC_DEBT`).

**4. Inside a sample: the inner (RHS) layer.** Budget = the advertised share,
or the whole pool for a lone simulation.

- *RHS plan*, i.e. how a step's per-satellite work is parallelised:
  `satellite_batch` (Polyester over satellites, width ≤ budget and ≤ a pinned
  plan's allotment), a `flat` queue at some allotment with a static or
  dynamic scheduler, or serial. Chosen, in order, by: a cached *verdict* for
  the shape's signature
  (`v6|machine|budget|sats-bucket|effectors|harmonics|density|outer`); else
  the pre-solve *calibration sweep* (every candidate timed for 5 warm-up + 10
  calls; ladder `[1, 2, 4, …, budget] × {static, dynamic}` plus batch at
  `budget, budget/2, budget/4`); else the in-run *width trial* when the solve
  is long enough to amortise it; else the per-call *heuristic*. Verdict rules
  as shipped: solve < 1 s → the cached verdict is honoured; solve ≥ 1 s → a
  pinned plan is re-swept on every solve and a heuristic verdict honoured only
  after three consecutive votes. Fix D replaces "every solve" with "once the
  solves run on the verdict have cost 20× the sweep".
- *Callbacks* (density and the rest) run at the static width
  `min(items, budget)` under V2 — no per-call AIMD — with the V2 pre-solve
  density-callback width calibration on GRAM-surrogate shapes.
- *Observation* for the policy telemetry: one per RHS call (satellite 1),
  not one per satellite.

**5. What persists, per machine.** `output/parallel_policy_state/`: the
outer-route bandit state (routes and split widths, keyed by shape signature,
profile, machine and T), the RHS calibration store (verdicts, solve length,
votes; fix D adds the sweep's cost and the solve time run on the verdict),
the cost-model constants, and the inner policy's hints. A new machine starts
cold and pays: one split race per threads shape, one tie exploration per tie
shape, one sweep per RHS signature bucket (three for a heuristic verdict).
Every later campaign replays.

**6. Measured versus constant.** Measured on the host: core and memory
capacities, rounds, the route at a tie, split widths, the local-slot cap, the
RHS plan and width, the sweep's cost and the solve length (D), the steady
per-sample route cost (E). The remaining constants are dimensionless or
structural: the 1 s long-solve threshold, three heuristic votes, the 5 %
re-verification share (D), the `T/4` ladder floor, "second half" for the
steady estimate (E), and the direction of the tie default (the pool). None of
them encodes a core count.

**7. Where the light set puts it** (L8–L14 on both machines, §4.4): parity
or better wherever the store was warm; ahead by 1.4–10× where the static
heuristics mis-thread (the 256-satellite stacks, 64-body interaction, the
exponential atmosphere); 25 % ahead on the mixed split. The losses are the
two cold-store mechanisms D and E address: the per-solve sweep on long
solves (L9/L10 on the TRX50, 1.35–1.61×) and the pool judged on its cold
campaign at a tie (L12 on the TRX50, 1.61×).

## 4. Measurements

### 4.1 B15, 128 spacecraft-hour grid, fixes 1–4 only (run `20260906_205919`, split 1 of 6 only)

| case (1 worker, 12 threads) | best static wall | R6 wall | R6 / best |
|---|---|---|---|
| mcgrid_16sat_8mc | 2.87 (`outer_threads`) | 3.37 | 1.17 |
| mcgrid_32sat_4mc | 3.73 (`outer_threads`) | 5.52 | 1.48 |
| mcgrid_8sat_16mc | 1.55 (`outer_threads`) | 1.59 | 1.03 |

The run died at split 2 (see §6). The two lagging points are what §2.1 traced.

### 4.1a Point validation of fixes A+B (harness worker, mcgrid_32sat_4mc, 1 worker × 12 threads, 3 repeats)

| | before A+B | after A+B |
|---|---|---|
| R6, 4 concurrent samples, wall | 5.25 s | **2.09 s** |
| R6, per-sample | 4.86 s | 1.97 s |
| static `outer_threads`, 4 samples, wall | 3.25 s | 2.24 s (fix A alone) |
| R6, one sample | 0.63 s | 0.56 s |

R6 goes from +62% to −7% against the static route at the point that was
worst; the static route itself gains from fix A; a single simulation is
unchanged.

### 4.2 B15 with mixed dispatch and fixes A+B, 128 + 1024 grids (run `20260906_220830`, 6h56m, `78bf070f`)

Lean ladder: `outer_process`, `outer_threads`, `inner_only`, `policy_v2`. Cell = R6 wall / best static wall at that launch point (< 1 is R6 ahead).

| case | (1,12) | (2,6) | (3,4) | (4,3) | (6,2) | (12,1) |
|---|---|---|---|---|---|---|
| mcgrid_32sat_4mc  | 1.03 | 1.11 | 1.06 | 0.99 | 1.06 | 0.97 |
| mcgrid_16sat_8mc  | 1.01 | 1.05 | 1.09 | 0.89 | 0.68 | 1.20 |
| mcgrid_8sat_16mc  | 0.99 | 0.91 | 0.87 | 0.82 | 1.06 | 1.25 |
| mcgrid_64sat_16mc | 0.92 | 0.93 | 0.70 | 0.72 | 0.67 | 1.00 |
| mcgrid_32sat_32mc | 0.96 | 0.80 | 0.84 | 0.92 | 0.69 | 1.22 |
| mcgrid_16sat_64mc | 0.95 | 0.84 | 0.78 | 0.79 | 0.81 | 1.12 |

36 launch points: min 0.67, median 0.94, max 1.25; 23 strictly ahead, 30 within
the 8% band. Best static is `outer_threads` at (1,12)–(3,4) and `outer_process`
from (4,3) on. `inner_only` is never within 3× of the outer routes, at 64
spacecraft per sample included: inner scaling on this workload is ~1.5× at 12
threads, so outer-first holds across the whole grid on this machine.

Misses: two small-sample ties at (2,6)/(3,4) where mixed dispatch needs the
same number of rounds as the threads route (§6, rounds rule), and four at
(12,1) -- the pure-process split, where R6 has no extra slots and should equal
`outer_process` -- at 1.12–1.25. See §4.2a.

### 4.2a The (12,1) misses are wall, not work

At the pure-process split R6 has no local slot and dispatches exactly as
`outer_process` does, yet reads 1.12–1.25× on four cases. Its **per-sample**
times are lower at every case (0.61 vs 0.71, 1.36 vs 1.49, 2.70 vs 2.84,
1.53 vs 1.63, 5.10 vs 5.41 s); the **wall** carries a fixed ~0.35 s on the
8–16-sample campaigns and intermittent 1–2 s stalls on the 32/64-sample ones
(32sat×32mc repeats: R6 11.04 / 8.72 / 10.56 vs process 9.88 / 8.67 / 8.56 --
repeat 2 matched to 0.05 s). Two things sit inside the harness's adaptive
timing and nowhere else:

- `ppc_run_adaptive_batch` builds the probe config and route features after
  `batch_started`; the pinned path resolves its route before the clock.
- V2's `gc_before_dispatch` runs a full `GC.gc()` before a process dispatch,
  inside the campaign. Its rationale is a threaded campaign leaving a
  multi-GiB heap that the process route's async feeders then stall on
  (measured 19.3 s vs 7.5 s in production). Every harness point is a fresh
  process whose heap the harness already collected, so here the collection is
  pure overhead; with a single coordinator thread it is also serial.

Planned: move the probe outside the clock; make the pre-dispatch collection
conditional on a "GC debt" flag that a threaded dispatch sets and the
collection clears, so the production case keeps its measured benefit and a
clean heap pays nothing. Both after the regression run.

### 4.3 Regression: B8–B14 (+B12) lean, run `20260907_050609`, `78bf070f`

Cell = R6 wall / best static wall (< 1 is R6 ahead); "before" = run
`20260902_153738` (pre-fix R4/R5/R6 ladder) for the same launch point.

**B13** `montecarlo_heavy_aerobraking`, 64 samples, fixed 12-core budget split six ways (1h05m):

| split | (1,12) | (2,6) | (3,4) | (4,3) | (6,2) | (12,1) |
|---|---|---|---|---|---|---|
| R6 / best static | 0.99 | **0.80** | **0.82** | **0.87** | **0.90** | 1.11 |
| R6 dispatch | threads | mixed 2+5 | mixed 3+3 | mixed 4+2 | mixed 6+1 | process |
| before (R6) | 1.02 | 1.05 | 0.98 | 1.52 | 2.94 | 1.01 |

The two splits where the shipped profiles lost by 1.5× and 2.9× are now
ahead of both static routes. The one miss is the (12,1) wall overhead of
§4.2a (7.81 vs 7.05 s), not a routing error.

**B12** interacting vs. independent (1h16m):

| case | t | R6 / best | before | note |
|---|---|---|---|---|
| interact_64sat_1hr | 8 / 12 | 0.24 / 0.18 | 0.24 / 0.18 | unchanged (§ B12 record, calibration effect) |
| interact_256sat_1hr | 8 / 12 | 0.89 / 0.73 | 0.89 / 0.70 | unchanged |
| independent_1sat_1hr | 1 / 8 / 12 | **3.40 / 1.72 / 1.65** | 0.95 / 1.05 / 0.94 | REGRESSION, route correct (`process`) |

The independent campaign is 256 one-satellite samples of ~37 ms each. New
for it in this run is the runner dispatch: at 256 samples the split race runs
(warm 12, then 36 at each of widths 4/8/12, then the rest), and every raced
batch goes through `_run_campaign_with_route_env`, which performs V2's full
`GC.gc()` before each dispatch -- five collections inside one timed campaign,
single-threaded at t=1. The previous run dispatched with the harness's own
`pmap`: no race, no collection. Same mechanism, one collection, is the (12,1)
overhead of §4.2a.

Two changes follow, both `adaptive_routing.jl`: collect before a dispatch
only when a threaded dispatch has run since the last collection (a "GC debt"
flag -- the production case the collection was measured for keeps it; a
clean heap pays nothing, and a race pays it at most once); and race split
widths only for the `:threads` route, where a narrower width buys each sample
inner budget -- for the process route (and mixed, whose local slots do not
depend on W) a narrower width only idles workers, so the widest split is the
answer and the race is pure cost.

**B8** `montecarlo_heavy_aerobraking`, worker ladder [1,2,4,8] × {16, 64} samples (1h14m): R6 / process route

| workers | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| 16 samples | 1.03 (route `none` = serial) | 1.04 | 1.08 | 1.18 |
| 64 samples | 1.01 | 1.01 | **1.36** | 1.09 |

Every point is the runner's pre-dispatch collection (+0.4–0.8 s per
campaign, once per dispatch, including the `none` route); the 1.36 is the one
point where the race is possible (candidates [2,4]: 12 samples with half the
pool idle, four dispatches, four collections). Same mechanism as B12; fix C.

**B9** `gravity_4096sat_l50_vacuum_1hr`, single simulation, thread ladder (16m): R6 / best static 1.00 / 1.09 / 1.04 at t = 1 / 8 / 12 (before: 1.09 / 1.11 / 1.06). No outer split, nothing here changed by this work; the t8 rung is the pre-existing open item of the 2026-09-02 record.

**B14** 1024-satellite duration/cadence cases, single simulation, thread ladder (1h25m): R6 / best static

| case | t=1 | t=8 | t=12 |
|---|---|---|---|
| `cadence_1024sat_10s` | 1.02 | 0.99 | 1.01 |
| `cadence_1024sat_1s` | 1.02 | 1.01 | 1.00 |
| `heavy_1024sat_l50_6hr` | 1.03 | 1.01 | 1.07 |

(before: 0.98–1.07). The best static route is `inner_only` at t ≤ 8 and
`outer_threads` / `outer_inner_static` at t = 12; R6 picks the same layer at
every rung. The 1.07 at heavy t = 12 is the fix-C class-gate/race overhead on
a 2.3 s run, inside the 8 % band.

**B10** 256-satellite atmosphere / GRAM-usage cases, single simulation, thread ladder (2h22m): R6 / best static

| case | t=1 | t=8 | t=12 |
|---|---|---|---|
| `atmo256_exponential_10min` | 0.96 | 0.86 | 0.75 |
| `atmo256_gram_live_10min` | 1.00 | 1.00 | 1.01 |
| `atmo256_gram_live_nbody_10min` | 1.00 | 1.00 | 1.01 |
| `atmo256_gram_surrogate_10min` | 0.99 | 1.03 | 1.04 |

(before: 0.68–1.10). The GRAM-live cases are pinned to the GRAM lock, so every
route is the same 19–41 s and R6 matches. The exponential case is where R6 is
clearly ahead of every static route: the static heuristics choose a narrower
batch width for a 256-satellite 10-minute run than the calibrated
`satellite_batch` width R6 keeps from its cached verdict.

TODO: B11 as it lands.

## 5. Changes to SpaceAGORA itself (`src/`)

| file | change |
|---|---|
| `simulation/campaigns/monte_carlo.jl` | `outer_split_env_pairs`; threaded dispatch advertises the share; `MonteCarloResult` gains `route`, `local_slots` (3-arg constructor unchanged); `_run_monte_carlo_mixed`; `_run_monte_carlo_process` is now the no-local-slots case |
| `simulation/campaigns/adaptive_routing.jl` | plan carries `local_slots`/`local_slots_at`; process branch dispatches mixed under `outer_split_env_pairs(local_slots)`; race result carries the route; trace prints local slots |
| `simulation/campaigns/adaptive_routing.jl` (fix C) | `_GC_DEBT` (collect before a process dispatch only after a threaded/mixed one); split race gated to `:threads`; race result reports `W + L`; local slots capped by the threads route's measured best width |
| `simulation/engine/setup.jl` | `_rhs_batch_workers`, `_rhs_batch_minbatch`; `_rhs_batch_parallel_enabled(p, n)` refuses a 1-wide batch |
| `simulation/engine/setup.jl` (fix C) | `_rhs_batch_workers` honours a pinned `satellite_batch` plan's `allotment > 1` |
| `simulation/engine/rhs_calibration.jl` (fix C) | `_make_calib_satellite_batch_plan(allotment)`; batch rungs at budget/2, budget/4 in `_rhs_plan_candidates`; `_rhs_plan_width` and the store round-trip carry the batch width |
| `simulation/engine/dynamics_rhs.jl` | the seven `satellite_batch` `@batch` loops sized by the inner budget; V2 records one effector observation per RHS call |
| `SpaceAGORA.jl` | imports and exports `adopt_process_workers!` |
| `parallel/routing/outer_route_selection.jl` | `mixed_local_slots`, `mixed_capacity`; V2 Monte Carlo rule compares mixed capacity |
| `parallel/routing/outer_route_selection.jl` (fix C) | `mc_route_rounds`, `mc_route_tie`; the V2 Monte Carlo rule counts rounds, ties default to the pool; `select_outer_route!` explores a tie's other arm once (`explore_tie`) |
| `parallel/routing/outer_route_state.jl` | `outer_route_mixed_dispatch()` (`SPACEAGORA_PARALLEL_MIXED_DISPATCH`, V2 only); `OuterRouteTuning.mixed_dispatch` |
| `parallel/routing/outer_route_state.jl` (fix C) | `OuterRouteTuning.explore_route_ties`, `tie_explore_min_campaigns` |
| `parallel/routing/machine_topology.jl` | `memory_local_slot_cap` |
| `parallel/process/worker_pool.jl` | `adopt_process_workers!` (exported) |

Harness (`benchmarks/`): `ppc_mode_env_pairs` share for the threads backend;
`ppc_run_adaptive_batch` (adaptive modes through the runner,
`SPACEAGORA_PPC_ADAPTIVE_VIA_RUNNER=0` restores pinned dispatch); row env from
the resolved mode plus the runner's additions; B15 cases and ladder;
`--lean-modes` (`81075084`). Fix C: the probe/features are built before the clock; one `OuterRouteState` per worker process across a point's repeats; `--light` = phases L8–L15 (`PPB_LIGHT_PHASES`), the B8–B15 axes at ~100 points, in `PPB_ROUTER_PHASES` for the regret summary.

Tests: `test/unit/parallel/outer_split_budget_tests.jl`,
`test/unit/parallel/mixed_dispatch_tests.jl`,
`test/unit/parallel/rhs_batch_budget_tests.jl`,
`test/unit/parallel/mc_route_tie_tests.jl` (fix C: rounds rule, tie exploration and exploitation both ways, local-slot cap, race gating, GC debt, batch rungs).

## 6. Open items

- Lesson, recorded so it is not repeated: the paper harness spawns a fresh
  Julia process per point that recompiles `SpaceAGORA` from the working tree.
  Editing `src/` under a live run fed a half-edited tree to a worker and killed
  the fixes-1–4 baseline run at split 2 (`UndefVarError: adopt_process_workers!`
  — an export appended before the matching `using .ParallelProcess:` import).
  The first split's rows stand; nothing else from that run does.
- The static routes' heuristic picks a serial `flat/1` RHS for a 32-spacecraft
  single simulation where `satellite_batch` is 6.8× faster (§2.1). Not touched.
- `_multibody_thread_decision` reads `ENV` once per satellite per RHS call, on
  every route. Not touched.
- The tie default, the local-slot cap and the batch rungs are measured on
  this box only (§3.1). The TRX50 run is the check that matters and needs
  the machine reachable; the `--light` set is sized for it.
- B11 of run `20260907_050609` (pre-fix-C code) was stopped at 33 of 60 rows
  to start fix C; its axis is covered by L11 on the final code.
