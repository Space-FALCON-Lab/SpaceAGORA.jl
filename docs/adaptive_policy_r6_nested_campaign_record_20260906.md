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

### 4.3 Regression: B8–B14 lean

TODO

## 5. Changes to SpaceAGORA itself (`src/`)

| file | change |
|---|---|
| `simulation/campaigns/monte_carlo.jl` | `outer_split_env_pairs`; threaded dispatch advertises the share; `MonteCarloResult` gains `route`, `local_slots` (3-arg constructor unchanged); `_run_monte_carlo_mixed`; `_run_monte_carlo_process` is now the no-local-slots case |
| `simulation/campaigns/adaptive_routing.jl` | plan carries `local_slots`/`local_slots_at`; process branch dispatches mixed under `outer_split_env_pairs(local_slots)`; race result carries the route; trace prints local slots |
| `simulation/engine/setup.jl` | `_rhs_batch_workers`, `_rhs_batch_minbatch`; `_rhs_batch_parallel_enabled(p, n)` refuses a 1-wide batch |
| `simulation/engine/dynamics_rhs.jl` | the seven `satellite_batch` `@batch` loops sized by the inner budget; V2 records one effector observation per RHS call |
| `SpaceAGORA.jl` | imports and exports `adopt_process_workers!` |
| `parallel/routing/outer_route_selection.jl` | `mixed_local_slots`, `mixed_capacity`; V2 Monte Carlo rule compares mixed capacity |
| `parallel/routing/outer_route_state.jl` | `outer_route_mixed_dispatch()` (`SPACEAGORA_PARALLEL_MIXED_DISPATCH`, V2 only); `OuterRouteTuning.mixed_dispatch` |
| `parallel/routing/machine_topology.jl` | `memory_local_slot_cap` |
| `parallel/process/worker_pool.jl` | `adopt_process_workers!` (exported) |

Harness (`benchmarks/`): `ppc_mode_env_pairs` share for the threads backend;
`ppc_run_adaptive_batch` (adaptive modes through the runner,
`SPACEAGORA_PPC_ADAPTIVE_VIA_RUNNER=0` restores pinned dispatch); row env from
the resolved mode plus the runner's additions; B15 cases and ladder;
`--lean-modes` (`81075084`).

Tests: `test/unit/parallel/outer_split_budget_tests.jl`,
`test/unit/parallel/mixed_dispatch_tests.jl`,
`test/unit/parallel/rhs_batch_budget_tests.jl`.

## 6. Open items

- Lesson, recorded so it is not repeated: the paper harness spawns a fresh
  Julia process per point that recompiles `SpaceAGORA` from the working tree.
  Editing `src/` under a live run fed a half-edited tree to a worker and killed
  the fixes-1–4 baseline run at split 2 (`UndefVarError: adopt_process_workers!`
  — an export appended before the matching `using .ParallelProcess:` import).
  The first split's rows stand; nothing else from that run does.
- The static routes' heuristic picks a serial `flat/1` RHS for a 32-spacecraft
  single simulation where `satellite_batch` is 6.8× faster (§2.1). Not touched.
- Found in the B15 re-run, to fix after it: (a) a race result reports
  `threads = best_w` (workers only), so `execution_scope` reads
  `adaptive_mixed_w0_l3` where it should read `w3_l3` -- label only, one line
  in `_run_campaign_split_race`; (b) every small-sample loss so far (1.05–1.11
  at 4–8 samples on the (2,6)/(3,4) splits) is a launch point where mixed
  dispatch and the threads route need the same number of rounds,
  `cld(n, W+L) == cld(n, T)`, so the pool's dispatch round-trip buys nothing.
  Candidate rule for the V2 Monte Carlo default: take the process (mixed)
  route only when it saves a round; fits 12 of the first 13 mixed points,
  costs a 9% win on 16sat×64mc at (2,6). The regression set is insensitive
  to it (many-sample or single-sample campaigns), so it can land after that
  run without invalidating it.
- The 1024 spacecraft-hour B15 rungs are running (`20260906_220830`);
  `inner_only` loses by 3–8× at every point so far, including 64
  spacecraft per sample -- outer-first holds across this grid on this box.
