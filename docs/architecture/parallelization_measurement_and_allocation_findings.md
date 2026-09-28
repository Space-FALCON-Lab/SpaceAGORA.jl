# Parallelization: Measurement Validity, Scaling Limits, and Allocation Findings

Investigation record, 2026-08-21, branch `router-eval-expanded-b6`.
Reference machine unless stated otherwise: AMD Ryzen 9 9900X, **12 physical
cores / 24 SMT threads**, 60 GB RAM, Julia 1.12.1. Final benchmark target is the
TRX50 box (**64 physical cores**), where several conclusions below will differ —
see [Transferability](#9-transferability-and-what-must-be-re-measured).

---

## 1. The headline: every published parallelization number measured JIT compilation

`ppc_worker_cmd` (`benchmarks/studies/parallelization_performance/execution.jl`)
built the worker argv without `--warmup`. The worker is a fresh subprocess that
re-parses its own CLI, where `warmup` defaulted to `0` from
`SPACEAGORA_PPC_WARMUP`, and the profile-default fallback was guarded by
`warmup < 0` — which `0` never satisfies. So the `full` profile's `warmup = 1`
was dead code and **every worker ran with zero warm-up**, making the timed run
the process's first run.

| case | first solve (what the harness timed) | steady-state solve |
|---|---|---|
| `gravity_16sat_l20_vacuum` | 21.4 s | **0.017 s** |
| `multi_64_high_fidelity` | 38.2 s | **0.18 s** |
| `multi_256_high_fidelity` | 47.5 s | 2.87 s |

The B6 CSV reports 43.6–47.9 s for `gravity_16sat_l20_vacuum` across all five
modes and both thread counts — i.e. ~99.9 % compilation. This single defect
explains why every mode landed within noise of every other, why router regret
was ±5 % noise, and Falcone review point 7's "R5 is ~25 % slower than the best
route and slower than serial".

**Nothing in B1–B6 as previously run measured parallelism.** All prior
parallelization results should be regarded as void and re-run.

---

## 2. Harness corrections

| fix | file | effect |
|---|---|---|
| Forward `--warmup` to the worker; `-1` unset sentinel so profile defaults apply | `parallelization_performance/{cli,execution}.jl` | wall times became real solve time |
| Provision + warm the `Distributed` pool before the clock starts | `parallelization_performance/execution.jl` | `addprocs` + worker `include`/`using SpaceAGORA` + first-solve JIT were inside the timed region, which made *more workers cost more* |
| `outer_tasks` → honest `SPACEAGORA_OUTER_PARALLEL_ACTIVE` | `parallelization_performance/{modes,execution}.jl` | see §3 |
| Thread ladder in **physical** cores (`SPACEAGORA_PPC_PHYSICAL_CORES` to override) | `parallelization_performance/cli.jl` | detects 12 here → `[1,2,4,6,9,12]`; TRX50's 64 → `[1,2,4,8,16,32,64]` |
| Calibration gated on `mode.policy_adaptive` (`auto` for R4/R5, `off` for R0-R3) | `parallelization_performance/modes.jl` | **superseded an earlier blanket `off`, which disabled the adaptive profiles' own mechanism and understated them by 5.2x — see §8b** |

`repeats` do not amortize compilation: each `(case, mode, threads, repeat)` is a
separate subprocess.

---

## 3. `OUTER_PARALLEL_ACTIVE` was a label, not a fact

The RHS router reads this flag as the factual claim *"an enclosing outer split
already owns the thread pool, do not start a nested one"*, and responds by
serialising the RHS. The harness asserted it for **every** threads/process-backed
mode — including single-simulation constellation cases, which dispatch exactly
one outer task and have no outer split at all. That pinned the RHS to one worker
for most of B1/B2/B3/B5/B6.

It now tracks the actual concurrent outer-unit count.

---

## 4. Source fixes

### 4.1 Router: `outer_active` must clamp width, not switch route

`setup.jl` had three branches that responded to the nested-split hazard by
routing to `:satellite_batch, allotment=1`. But the constraint is "do not start a
nested *thread* split", which clamping the allotment to 1 already satisfies —
both the harmonics batch kernel and the flat queue's `threaded_foreach`
degenerate to inline serial calls at allotment 1.

Routing away instead gave up the batched kernel's load-coefficients-once sweep:

| 1024 sat / L50 / 1 thread | µs per RHS call |
|---|---|
| `:flat_constellation_effector_queue` | **1216** |
| `:satellite_batch` | **4918** (≈4×) |

Same ≈4× ratio at 64 and 256 satellites. Fixed in all three branches.

### 4.2 R0 was not actually serial

`profile_forces_serial` reached only `_rhs_batch_parallel_enabled`, so an R0
"true serial baseline" launched with >1 thread still threaded the flat queue:
**1272 µs/RHS-call at 1 thread vs 387 at 12**, same profile. That silently
deflates every speedup measured against it. The budget is now clamped to 1.

### 4.3 O(N²) Jacobian prototype → O(N)

`_build_block_diagonal_jac_prototype` allocated `zero(u0)` plus a full-length
`Vector{Float64}` copy plus a full-length `findall` **per satellite**. Replaced
with a single stamped pass. Structurally identical output on every layout
`ComponentArrays` can construct (verified against the original).

| satellites | before | after |
|---|---|---|
| 1024 | 8.7 ms | 0.3 ms |
| 4096 | 62.5 ms | 2.2 ms |

### 4.4 Vacuum heat-rate short-circuit

`_compute_stage_heat_rates!` ran a full `sample_planet_frame` (`rtolatlong`:
atan/asin/sqrt) plus a density sample per satellite per RHS stage even under
`NoAtmosphereModel`, only to fall out of the `rho <= 0` guard with a zeroed
buffer. ~41 % of the solve in a 1024-satellite L20 vacuum profile. Now a static
type check. Full 6 h solve: **10.22 s → 7.88 s**.

---

## 5. Real scaling, measured post-warm-up

### 5.1 Vacuum (harmonics only) — scales, saturates by workload size

`heavy_1024sat_l50_6hr`: 7.75 s serial → **2.00 s @ 4t (3.87×)** → 2.28 s @ 12t
(*worse*).
`heavy_4096sat_l50_1hr`: 9.68 s serial → 3.62 @ 4t → 2.59 @ 8t → **1.88 @ 12t
(5.14×)**, still climbing.

Isolated RHS, best speedup and where it occurs:

| N | best | at threads |
|---|---|---|
| 256 | 2.3× | 8 |
| 1024 | 4.2× | 4 |
| 4096 | 4.9× | 16 |

**Interpretation.** What moves with N is *where the curve flattens*, not how high
it gets — a ~4–5× asymptote reached at a wider thread count, not exceeded.
Caveat: per-satellite serial cost rises with N (0.98 → 1.23 → 1.59 µs) while the
parallel per-satellite cost is flat (~0.3 µs), so part of the apparent gain is
the *baseline degrading*, not the parallel path improving. Do not present these
three rows as a table without re-running them under one code vintage — they mix
pre- and post-fix measurements.

### 5.2 CORRECTED: adaptive routing is what makes constellations scale

> **The tables previously in §5.2/§5.3 were invalid.** They were taken with
> `SPACEAGORA_RHS_CALIBRATE=off`, a setting introduced during this investigation
> (§2) and applied to *every* benchmark mode. That disabled the RHS plan
> calibration which makes R4/R5 adaptive — i.e. it switched off the mechanism the
> adaptive profiles exist to evaluate. Every "multi-effector constellations do
> not scale / invert below serial" conclusion came from that artifact. See §8b.

Corrected profile ladder, medians of 2 reps, `serial` at 1 thread as baseline:

| case | serial | R1a static @4 | R4 adaptive @4 | R1a static @12 | **R4 adaptive @12** |
|---|---|---|---|---|---|
| `multi_256_high_fidelity` | 2.674 s | 1.585 (1.69x) | 0.887 (3.02x) | 4.331 (**0.62x**) | **0.681 (3.93x)** |
| `heavy_1024sat_fullstack_1hr` | 23.985 s | 14.353 (1.67x) | 7.227 (3.32x) | 30.504 (**0.79x**) | **5.603 (4.28x)** |
| `heavy_256sat_coupled6dof_2hr` | 79.625 s | 48.177 (1.65x) | 26.118 (3.05x) | 88.995 (**0.89x**) | **19.276 (4.13x)** |
| `heavy_1024sat_l50_6hr` (vacuum) | 7.779 s | **1.976 (3.94x)** | 2.103 (3.70x) | 2.173 (3.58x) | 2.439 (3.19x) |

R5 (`full_smart`) tracks R4 within noise throughout (0.689 / 5.811 / 19.822 /
2.336 at 12 threads).

Three results, all consistent across a 4x range of satellite count and across
3-DOF and 6-DOF:

1. **Static outer-threading inverts at 12 threads on every multi-effector case**
   (0.62x, 0.79x, 0.89x) — worse than running serial.
2. **Adaptive routing reaches 3.9-4.3x on the same cases at the same thread
   count**, i.e. **4.6x to 6.4x faster than static outer-threading**
   (30.504/5.603 = 5.4x; 4.331/0.681 = 6.4x; 88.995/19.276 = 4.6x).
3. **On the simple vacuum case adaptive is marginally *worse* than static**
   (3.19x vs 3.58x at 12 threads): calibration costs a few percent where there is
   nothing to adapt. Worth reporting honestly — it is the price of the mechanism.

The distinction that matters is therefore not "does SpaceAGORA scale" but "does
the execution strategy match the workload". A single heavy batchable effector
needs no routing help; a heterogeneous effector set with atmosphere and callbacks
cannot be threaded naively at all.

---

## 6. Thread vs process: the strongest result in this investigation

`montecarlo_heavy_aerobraking`, 24 samples, matched core budgets:

| cores | outer_threads | speedup | outer_process | speedup | process ÷ thread |
|---|---|---|---|---|---|
| 1 | 26.34 s | 1.00× | — | — | — |
| 2 | 15.46 s | 1.70× | 14.17 s | 1.86× | 1.09× |
| 4 | 9.82 s | 2.68× | 7.81 s | 3.37× | 1.26× |
| 8 | 12.63 s | 2.09× ↓ | 4.46 s | 5.91× | 2.83× |
| 12 | 13.38 s | 1.97× ↓ | **3.09 s** | **8.53×** | **4.33×** |

Process efficiency: 93 / 84 / 74 / 71 %. Threads peak at 4 cores and invert.

**This case uses `ExponentialAtmosphereModel` — no native library, no lock.** So
process isolation wins for reasons unrelated to GRAM: separate heaps, separate
GC. That generalizes further than the native-library-contention framing in the
current manuscript.

**Structural constraint for any routing policy:** a single interacting
constellation is one ODE over a shared state vector, so the process route has
nothing to split. Thread-vs-process is only a live question where the workload
has independent outer units.

---

## 7. Allocation: large, reducible — but NOT the thread ceiling

> **Scope caveat.** Everything in §7 was measured under `outer_threads` (R1a), a
> *static* profile, which is the configuration that inverts at 12 threads
> (§5.2). The allocation reductions are real and profile-independent, but the
> "thread ceiling" this section chases is a property of static outer-threading,
> **not of the shipped adaptive configuration**, which reaches 3.9-4.3x on the
> same cases. Read §7.3/§7.4 as an analysis of why the static route degrades.

`heavy_1024sat_fullstack_1hr`, identical physics, only thread count varies:

| threads | wall | GC | GC % | **allocated** | collections |
|---|---|---|---|---|---|
| 1 | 23.73 s | 2.69 s | 11.3 % | **29,483 MiB** | 474 |
| 4 | 16.54 s | 2.71 s | 16.4 % | **41,300 MiB** | 846 |
| 12 | 31.65 s | 5.27 s | 16.7 % | **61,252 MiB** | 1402 |

`heavy_1024sat_l50_6hr` (vacuum) for contrast: 1,638 / 1,856 / 2,259 MiB — ~110×
less per unit of simulated time, over a mission 6× longer.

**+31.8 GB purely from adding threads**, on unchanged physics. GC alone does not
explain the regression (16.7 % at 12 threads); subtracting it, non-GC time still
rose 21.0 → 26.4 s, implicating allocator contention and page faults. All of
these are process-wide, which is exactly why separate processes escape the
ceiling and threads do not.

### 7.1 Root cause: heterogeneous effector tuple

Measured directly (not inferred from the profiler):

| pattern | allocation |
|---|---|
| `for x in` heterogeneous 3-tuple | **240 bytes** |
| `for x in` homogeneous 1-tuple | 0 bytes |
| runtime index `e[i]` into het. 3-tuple | **144 bytes** |
| `any(generator)` over het. 3-tuple | 0 bytes |

`eltype(dynamic_effectors)` is the abstract `AbstractForceTorqueModel`, and
`e[i]` infers to a 3-way `Union`. So iterating or runtime-indexing the effector
tuple boxes, scaling as satellites × effectors × RHS calls. The vacuum cases have
a **1-tuple, which infers concretely and allocates 0** — precisely why they scale
and the atmosphere cases do not.

> **`Profile.Allocs` mis-attributes under inlining.** It blamed 47 % of serial
> allocation on `_uses_j2_gravity_effector`, a trivial `isa` loop that cannot
> allocate (it is inlined into the per-satellite density sampling in
> `effector_sampling.jl`), and 3.27 GiB on an `any(generator)` that measures
> 0 bytes. Treat its line numbers as leads and confirm with `@allocated`.

### 7.2 Allocation fixes applied

All four are validated (full suite `EXIT=0`, serial-vs-threaded parity
bit-identical at fixed thread count). Their measured effect on
`heavy_1024sat_fullstack_1hr` is at the end of this section.

1. **Env-snapshot hoist.** `_rhs_effector_static_cost_ns` called
   `_effector_cost_ns_per_item_default()` -> `_parse_positive_float_env`, which
   does `string(default)` (unconditional String allocation), then `strip`, then a
   parse — **432-448 bytes per call, from inside the RHS**. The value is already
   snapshotted per run in `RhsPlanEnvConfig.effector_cost_ns_per_item_default`; it
   now reads that field. **448 -> 0 bytes/call.** Same hoist commit `d4f5cce4`
   performed for other knobs; these cost-model readers were missed. The other two
   432-byte readers are called only from the snapshot builder (once per run) and
   needed no change.
2. **Serial tuple peel.** `_accumulate_effector_chain!` replaces
   `for effector in dynamic_effectors` with a recursive `Base.tail` peel, giving
   each level a concretely typed `first(effs)`.
3. **Flat-queue selection mask.** `_prepare_rhs_flat_work_items!` indexed the
   effector tuple `num_sats x n_effectors` times per RHS call for a decision that
   depends only on `eff_idx`. Hoisted to a homogeneous `NTuple{N,Bool}` built once
   per call (homogeneous tuples index allocation-free).
4. **Flat-queue item peel — implemented, measured, then REVERTED.** Both queue
   variants were routed through a shared `_flat_item_accumulate!` peel. It was
   correct (parity bit-identical, suite green) and cut allocation ~3 GB at 12
   threads, but produced **no wall-time change** on any multi-effector case, and
   consolidating the two bodies removed the literal
   `... && continue` form of the double-counting skip guard that
   `ci_rhs_parallel_route_gate` asserts textually in the packet worker. With no
   performance justification, relaxing an architecture gate was not warranted, so
   the change was reverted. Anyone retrying this must keep both the
   `&& return nothing` and `&& continue` guard spellings intact, or move the gate
   deliberately.

**Measured effect** (`heavy_1024sat_fullstack_1hr`, allocation in MiB):

| threads | original | shipped (fixes 1-3) | with reverted fix 4 |
|---|---|---|---|
| 1 | 29,483 | **20,199** | 20,199 |
| 4 | 41,300 | **25,928** | 22,722 |
| 12 | 61,252 | **45,880** | 42,674 |

The third column is recorded only to document what fix 4 would have bought; it is
not the shipped state.

Serial/low-thread wall time improved ~10-20 % (fullstack best case 16.50 s ->
12.94 s at 4 threads, peak speedup 1.56x -> 1.90x).



### 7.3 The thread ceiling is NOT (mainly) allocation

The four fixes above cut allocation by 25-45 % and left the scaling curve's
**shape completely unchanged** — see §5.3. GC is only 13-14 % of wall at 12
threads, and wall grows 12.9 -> 28.2 s from 4 to 12 threads, so roughly 86 % of
the growth is non-GC, non-allocation time.

**Allocation was a symptom that was over-weighted, not the cause** of the static
route's ceiling. Note also that the ceiling itself only binds the static profile:
the adaptive profile reaches 3.9-4.3x on these same cases (§5.2), so this is a
characterisation of one route's degradation rather than a limit on the simulator.

> **`Profile.Allocs` was wrong three times in this investigation**, each time by
> over-attributing through inlining: 47 % of serial allocation blamed on
> `_uses_j2_gravity_effector` (a trivial `isa` loop that provably cannot
> allocate); 3.27 GiB on an `any(generator)` measuring 0 bytes; and ~25 GiB on six
> `partials` accumulate lines whose removal freed 3 GB and no time. **Confirm
> every attribution with `@allocated` before acting on it.**

One structural candidate identified but **not yet measured**: the flat queue's
reduction is `for worker_id in 1:workers, for sat_idx in 1:num_sats` over six
components — an O(workers x satellites) *serial* pass every RHS call, plus a
matching O(workers x satellites) zeroing of `partials`. Both grow linearly in
worker count, which is the right shape. A back-of-envelope puts it at ~0.75 s
over a full solve, far short of the ~15 s to be explained, so it is at best
partial.

### 7.4 Hunting the thread ceiling: what has been ruled out

The aero RHS **in isolation** (N=1024, harmonics + SRP + aero + density, µs per
RHS call) inverts, so this is not a solver-spine or callback effect:

| threads | 1 | 2 | 4 | 8 | 12 |
|---|---|---|---|---|---|
| aero RHS | 4957 | 2678 | **1637** | 3431 | 4755 |
| vacuum L50 RHS | 969 | 396 | **229** | 249 | 247 |

The vacuum path plateaus; the aero path *inverts back to serial speed*. Ruled out
by direct measurement:

| hypothesis | test | result |
|---|---|---|
| GC thread count | `--gcthreads` 1/2/4/6 vs default 12 | 24.5-31.4 s, all within noise |
| False sharing / scheduler layout | static stride-W vs dynamic contiguous chunk 32/128 | 1911 -> 1599 µs @4t, 4580 -> 4083 @12t; marginal |
| Dispatch mechanism | persistent channel pool vs plain `@spawn` | no consistent win (spawn better @8, worse @12) |
| Allocation | four allocation fixes, -25..45 % | curve shape unchanged (§7.3) |
| Effector boxing | flat-queue item peel | no wall-time change (§7.2) |
| Route choice | forced `flat` vs `satellite_batch` | matters (1.7x) but does not explain inversion |

Aggregate CPU from the phase profile (`heavy_1024sat_fullstack_1hr`, all threads
summed, core-seconds): **43.3 @1t -> 194.6 @4t -> 973.5 @12t** against wall times
of 22.3 / 13.6 / 28.3 s. On a 12-core box that is ~34 threads' worth of CPU being
consumed — Julia runs 12 GC threads alongside 12 compute threads, and Polyester
maintains a further pool. GC is 24 % of aggregate CPU at 12 threads and idle
`poptask` another 24 %. A 12-thread profile of the aero RHS is dominated by
channel-wait frames (`wait`, `condition.jl:wait`, `channels.jl:take!`,
`_persistent_foreach_worker_loop`).

**Root cause not established.** The phenomenon is well characterised and the
above are excluded; the remaining suspects are oversubscription between Julia's
threads, its GC threads and Polyester's pool, and memory-bandwidth saturation.
Note run-to-run variance on these micro-measurements is ~16 %, so future work
needs more repeats than used here.

---

## 8. GRAM

### 8.1 Two documented blockers no longer reproduce

- `multi_4_gram_live`'s comment documents 14,540 accepted steps for a **10 s**
  window and no completion inside a 10-minute timeout. Measured: **1621 steps for
  the full 1200 s mission, solving in 1.23 s.**
- `multi_16_gram_live` is excluded from B6 for an unbounded native leak climbing
  to 30+ GB. Measured: **`Success` in 5.36 s, RSS flat at 3.8–4.2 GB.**

This branch's ephemeris-sharing work (`92142ade`, `7ae8c891`) appears to have
fixed both. **B6's exclusion comment is now stale**, and the interacting-GRAM
axis review point 8 asks for may be recoverable. One clean run is
non-reproduction, not proof — confirm at full duration under a memory watchdog
before re-enabling.

### 8.2 The isolated pool is unreachable for constellation workloads

`SPACEAGORA_GRAM_ISOLATED_POOL=on` produced **`pool_models len = 0` at N=4, N=16
and N=64** — it never engaged. Any pool timing comparison taken this way is a
null *measurement*, not a null result. Two independent causes:

1. Benchmark mode `outer_threads` (R1a) sets
   `SPACEAGORA_DENSITY_CALLBACK_PARALLEL => "off"`, which `withenv` applies over
   any value passed alongside it.
2. **Structural:** `gram_runtime_stats` reports `density_calls = 0` even under
   `full_smart`/`inner_only`. The constellation aero path samples density via
   `effector_sampling.jl` → `_density_model_for_sat(p, sat_idx)` →
   `_density_state_from_kinematics!`, which is **per-satellite**.
   `_gram_isolated_pool_batch_eval!` lives in the density *callback's batched*
   branch, which this workload never enters. A secondary gate also applies:
   density-callback auto-threading requires a ≥16 thread budget.

The crash fix on this branch (`_warm_gram_pool_model!`) made the pool *safe*, not
*reachable*. It cannot be offered as a routing strategy on current code.

### 8.3 Per-satellite instances: now safe, but slower

`_initialize_density_model_instances!` deep-copied GRAM per satellite with **no
warm-up**, carrying the same CSPICE `trcpkg` hazard the pool fix addressed (a
fresh clone reaches CSPICE on its first `update()`; CSPICE shares a global
call-trace stack). The same serial warm-up is now applied, guarded to the native
GRAM type.

`heavy_16sat_gram_nbody_l50`, RSS sampled *during* the run:

| config | peak RSS | solve |
|---|---|---|
| per_sat=0, 1t | 4.4 GB | 2.853 s |
| per_sat=0, 4t | 4.3 GB | 2.574 s |
| per_sat=1, 1t | 5.9 GB | 3.551 s |
| per_sat=1, 4t | 5.8 GB | 3.233 s |

No hang, no crash — where this path was previously documented as hanging. But
isolation is **24–26 % slower** and costs +1.5 GB (~95 MB per clone), and thread
scaling is flat (~1.10×) either way. **The GRAM lock is not the bottleneck at
N=16.** Whether it becomes one at higher N is untested.

### 8.4 New case family

`heavy_<N>sat_gram_nbody_l50` (N ∈ 16/32/64/128/256): live GRAM Earth atmosphere
+ L50 harmonics + Sun/Moon third-body gravity (SPICE calls under `SPICE_LOCK`) +
aero. Unlike `multi_4`/`multi_16_gram_live`, it produces real thread scaling —
at N=64, **11.26 s → 5.75 s at 4 threads (1.96×)**, turning over to 6.24 s at 8 —
so it is a valid test bed where the smaller GRAM cases were flat.

---

## 8b. Adaptive routing: it works, and the benchmark was hiding it

Investigated 2026-08-22 on `heavy_1024sat_fullstack_1hr` and others.

### 8b.1 The defect was in this investigation's own harness change

§2 records adding `SPACEAGORA_RHS_CALIBRATE => "off"` to
`ppc_mode_env_pairs`, on the reasoning that pre-solve calibration installs an
`rhs_plan_override` which short-circuits the routing policy and so contaminates a
profile-ladder comparison.

That reasoning is sound for the **static** profiles R0-R3, which exist to measure
one fixed route. Applying it to **R4/R5 as well disabled the mechanism those
profiles are defined by** — calibration is part of "the full adaptive routing
policy SpaceAGORA ships with" (`modes.jl`). Measured cost of the mistake, at 12
threads:

| `heavy_1024sat_fullstack_1hr` | wall | vs serial |
|---|---|---|
| calibration off (as benchmarked) | 29.90 s | 0.80x |
| calibration on (as shipped) | **5.71 s** | **4.20x** |

**Fixed**: `"SPACEAGORA_RHS_CALIBRATE" => (mode.policy_adaptive ? "auto" : "off")`.
Static profiles keep a fixed route; adaptive profiles get the adaptive mechanism.

### 8b.2 Two attempted router improvements, both negative

Recorded so they are not retried blindly.

**Widening the calibration ladder.** `_rhs_plan_candidates` built its allotment
ladder from `viable_workers = active_sats / min_sats_per_worker` (256 at 1024
satellites) and then clamped with `min(a, budget)`, so at a 12-thread budget
`[2, 128, 256]` collapsed to `[2, 12]` — intermediate widths were unreachable.
Rebuilt geometrically from the budget (`[1, 2, 4, 8, 12]`; `[1..64]` on the
TRX50). Correct in principle, but **no measurable effect**: 5.63 s (original
ladder) vs 5.96 s (new) at 12 threads, within noise, because the winning
candidate is `:satellite_batch` whose selection does not depend on the flat
allotment ladder. Retained for the wider search space, not claimed as a gain.

**Calibrating the global inner width.** `SPACEAGORA_INNER_THREAD_BUDGET=4` gives
16.35 s where the uncapped static run gives 29.31 s, so a global width cap looked
promising as a new calibration dimension. Implemented, measured **worse** (18.46
s vs 5.71 s) and **removed**. Two reasons: applying a pinned width requires
suppressing the route sweep, because `_rhs_execution_plan` returns
`rhs_plan_override` before consulting the budget; and the route the sweep prefers
(`:satellite_batch`) takes its width from Polyester's own pool and honours
neither `allotment` nor `INNER_THREAD_BUDGET`. The whole avenue is dominated by
simply letting calibration pick the route.

### 8b.3 Probe fidelity: two hypotheses tested and rejected

Before the harness defect was found, the calibration probe was suspected of
mis-ranking routes. Both plausible explanations were tested and **rejected**:

- *Unrepresentative probe state* — sweeping at mission fraction 0.0 / 0.5 / 0.9
  gives the same ranking.
- *Tight probe loop keeping worker threads hot* — inserting 0 / 200 / 2000 us
  gaps between probe calls does not change the ranking.

The probe ranks `:satellite_batch` fastest and the full-solve evidence agrees.
**The calibrator's judgement was never the problem.**

---

## 9. Transferability and what must be re-measured

- **SMT hurts.** Every workload regressed 16 → 24 threads on this 12-physical-core
  box. A ladder topped out at `Sys.CPU_THREADS` puts its worst point last.
- **Nothing tuned here transfers to 64 cores.** Fork/join overhead, the
  density-callback ≥16-thread gate (which *passes* on the TRX50 and fails here),
  and every saturation point are machine-dependent. Policy evaluation must run on
  the target machine. `persistent_hints.jl` / `output/parallel_policy_state` must
  be keyed by machine.
- **Known traps.** `SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=1`, documented as
  enabling scaling to 32–128 threads, measured **2× slower** (2.67 → 5.31 s @12t)
  — spinning workers starve the other parallel regions.
  `SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER` 4→256: no effect. Persistent
  pool vs plain `@spawn`: 56 vs 25 µs/dispatch in isolation, a wash end-to-end.

---

## 10. Other defects found, not yet fixed

- **`num_steps_to_save` is inert.** Only
  `analysis/verification/telemetry_verification/` reads it; the solver uses
  `save_everystep` with no `saveat`. Measured `nsave = 2173` for
  `num_steps_to_save` ∈ {10, 300, 3000}. **B6's output-cadence axis
  (`gravity_16sat_l20_vacuum_sparse_output`) therefore measures nothing** — it is
  byte-identical work to its baseline.
- **`reset_fsal!` costs an extra full RHS per accepted step** (~15 % of thread-1
  samples): an always-true callback marks the state modified. If the
  thermal/event callbacks do not actually mutate `u`, `u_modified!(integrator,
  false)` removes ~1/7 of all RHS evaluations. Flagged from profiling; untested,
  and it needs event-handling validation.

---

## 11. New benchmark phases

- **B7 — Heavy Constellation Scaling.** `heavy_1024sat_l50_6hr`,
  `heavy_4096sat_l50_1hr`, `heavy_1024sat_fullstack_1hr`,
  `heavy_256sat_coupled6dof_2hr`. Serial baselines 7.9 / 8.0 / 28.5 / 88.7 s.
  The 1024-vs-4096 pair is the load-bearing contrast: it localises *where*
  per-RHS fork/join overhead stops being amortised.
- **B8 — Heavy Monte Carlo Process Throughput.** `montecarlo_heavy_aerobraking`
  (~1.05 s per sample), worker ladder 1–64.

Both are excluded from `--preview` (its 64-satellite cap would drop every case
and then run the two largest anyway).

### B9–B14 — Expanded router evaluation (review point 8)

Added 2026-08-22. One phase per workload axis the review names, replacing B6 —
whose numbers are not reportable, since every case it draws on is 0.017–2.9 s
serial and it measured 0.0–3.3% regret across every route, inside the ~16%
variance §7.4 records. B6 is retained as a below-the-floor control. Three
differences from B6 beyond sizing: the full static ladder R0–R3 runs at every
point (B6 omitted R2/R3, which biases `best_static` *downward* and flatters the
router); the thread budget has four points rather than two; and the router's
actual selection is now recorded rather than inferred.

- **B9** spacecraft count, `gravity_{256,1024,4096}sat_l50_vacuum_1hr`.
- **B10** atmosphere fidelity at N=256/600 s: vacuum → exponential → GRAM
  surrogate → live GRAM → live GRAM + SPICE third-body. The load-bearing phase
  for review point 7's "R5 is ~25% slower than the best static route on GRAM".
- **B11** force/actuator model count, one model added per rung.
- **B12** interacting vs. independent propagation at matched satellite-hours.
- **B13** one fixed 64-core budget split seven ways between processes and
  threads — the hybrid nesting `OUTER_PARALLEL_ACTIVE`/`INNER_THREAD_BUDGET`
  arbitrate, which no previous phase exercised.
- **B14** mission duration and output cadence.

Three things the sizing calibration established, all of which changed the design:

1. **The actuator path is O(N^2.2).** `stack{32,64,256}_e6_actuated` measures
   0.583 / 2.728 / 33.4 s per 10 s of simulated mission, against a near-linear
   `e5`. Cause: `get_control_callbacks` builds one `PeriodicCallback` per control
   effector and each loops over *all* satellites, so N per-satellite models cost
   N² `calcControlEffect!` calls per tick — N intended, N(N−1) returning
   immediately on a `sat_idx` the model does not own. Invisible at the 8
   satellites `multi_8sat_magnetorquer_attitude` uses, dominant at 256. The
   shared-model-with-per-spacecraft-vectors shape `BaseThrusterModel` enforces
   avoids it. **Not fixed** — B11 runs two sub-ladders instead.
2. **Output cadence is flat; output on/off is not.** At N=1024 over a 10 s
   mission: 0.051 s with output off, 3.13 / 3.09 / 3.17 s at `saveat` = 60 / 10 /
   1 s. Pinning `SPACEAGORA_SOLVER_SAVE_EVERYSTEP=0` *and*
   `SPACEAGORA_SOLVER_SAVE_ON=0` does not move it (3.25 / 3.25 / 3.24 s), so the
   62× is neither per-step solution storage nor the solver's save path but a
   fixed cost of the saving callback being in the `CallbackSet` at all. §10's
   note that `num_steps_to_save` is inert stands, and is now moot for the
   harness: the live knob is `simulation_settings.results`, which this harness
   had hardcoded off, so no benchmark case before B14 produced any output.
3. **Live GRAM is no longer the blocker §8.1 hoped.** `atmo256_gram_live_10min`
   at N=256 runs 0.623 s per 10 s of mission with no leak and no step-size
   collapse — 2.6× the precomputed surrogate and comfortably measurable. The
   interacting-GRAM axis point 8 asked for is recoverable, and B10 takes it.

---

## 12. Paper framing

Rewritten 2026-08-22 after §8b corrected the calibration-gating artifact. The
earlier version of this section advised narrowing the contribution and dropping
the adaptive-routing claim; **that advice was based on invalid measurements and
is withdrawn.**

### 12.1 The result

At a matched 12-core budget, on identical physics:

| workload | serial | static outer-threads | **adaptive (R4/R5)** |
|---|---|---|---|
| `multi_256_high_fidelity` | 1.00x | 0.62x | **3.93x** |
| `heavy_1024sat_fullstack_1hr` | 1.00x | 0.79x | **4.28x** |
| `heavy_256sat_coupled6dof_2hr` | 1.00x | 0.89x | **4.13x** |
| `heavy_1024sat_l50_6hr` (vacuum) | 1.00x | **3.58x** | 3.19x |
| MC campaign, process-isolated (§6) | 1.00x | 1.97x (threads) | **8.53x** (process) |

Adaptive routing is worth **4.6x-6.4x over static outer-threading** on
multi-effector constellations, and static outer-threading is *slower than serial*
on all three. On the single-effector vacuum case the ordering reverses and
adaptive costs ~10 %.

### 12.2 Why this is the strongest available framing

**Lead with the contrast, not a single speedup.** The claim is that execution
strategy must match workload coupling, and the evidence is that the *same*
strategy which delivers 3.58x on a vacuum constellation delivers **0.79x** on a
multi-effector one, while routing recovers 4.28x. A naive threaded simulator is
not merely leaving performance on the table — it is slower than not threading at
all. That makes contribution 3 demonstrated rather than asserted.

**The negative result strengthens it.** Static outer-threading inverting below
serial, consistently across 256-1024 satellites and 3-DOF/6-DOF, is the evidence
that routing is load-bearing. Report it.

**Report the cost of the mechanism.** Adaptive is ~10 % *worse* on the workload
that needs no adaptation (3.19x vs 3.58x). Reporting that is what makes the rest
credible.

**Report capability alongside ratios.** "1024 satellites with harmonics + SRP +
aero + density, one hour of mission time in 5.6 s" is ~640x real time and is
directly comparable against Basilisk/Tudat. Ratios belong in the scaling
analysis; throughput belongs in the capability claim.

**Add a methods-validity subsection.** Two defects in this investigation each
changed conclusions by more than any physics or algorithm change: the warm-up
omission (§1, ~99.9 % of every prior number was JIT) and the calibration gating
(§8b, 5.2x). Both generalise to any JIT-compiled simulator benchmark with a
configurable execution policy. Including them demonstrates rigour and explains
why the numbers moved.

### 12.3 What to claim

- **Workload-aware routing works and is necessary**: 4.6-6.4x over static
  threading on coupled multi-effector constellations; static threading is
  net-negative there.
- **Process isolation for independent campaigns**: 8.53x vs 1.97x at 12 cores on
  identical physics, mechanism = heap/GC isolation, so it generalises beyond
  native-library contention (§6).
- **Honest limits**: adaptive costs ~10 % where unnecessary; the static route's
  thread ceiling is characterised but its root cause is unresolved (§7.4).

### 12.4 Benchmark hygiene this implies

Any phase intended to measure the *shipped* configuration must let R4/R5 keep
calibration. Only the static profiles should pin a route. B1-B3 and B7 use
`serial`/`outer_threads`/`full_smart`, so their static columns are a lower bound
on shipped performance, not a measurement of it — state this where those tables
appear.

## 13. Reproduction

```bash
# post-warm-up steady-state timing for one case/mode/thread count
ST_CASE=heavy_1024sat_l50_6hr ST_MODE=outer_threads ST_REPS=2 \
  julia --threads=4 --project=. <scratch>/steady.jl

# allocation + GC accounting
AP_CASE=heavy_1024sat_fullstack_1hr julia --threads=12 --project=. <scratch>/allocprobe.jl

# full harness, warm-up now forwarded
SPACEAGORA_PPC_PROFILE=full SPACEAGORA_PPC_CASES=<case> \
SPACEAGORA_PPC_MODES=serial,outer_threads SPACEAGORA_PPC_THREADS=1,4,12 \
  julia --project=. benchmarks/studies/parallelization_performance.jl full

# phase dry-run (verifies case resolution + thread ladder without executing)
SPACEAGORA_PPB_PHASES=B7,B8 SPACEAGORA_PPB_DRY_RUN=1 \
  julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl
```

## 14. Validation status

Every change below is in the working tree of `router-eval-expanded-b6` and has
passed: full `test/runtests.jl` (`EXIT=0`, 128 test summaries), the PR-tier
contract gates (including `ci_rhs_parallel_route_gate` and
`ci_rhs_calibration_gate`), and serial-vs-threaded trajectory parity
(`pos_rel_max = 0.0`) at a fixed thread count.

| change | file | performance outcome |
|---|---|---|
| `--warmup` forwarded to worker; `-1` sentinel | `parallelization_performance/{cli,execution}.jl` | measurements became valid at all |
| Distributed pool provisioned + warmed pre-clock | `parallelization_performance/execution.jl` | removed startup bias from MC throughput |
| `outer_tasks` → honest `OUTER_PARALLEL_ACTIVE` | `parallelization_performance/{modes,execution}.jl` | unpinned the RHS for single-simulation cases |
| Physical-core thread ladder | `parallelization_performance/cli.jl` | TRX50 ladder `[1,2,4,8,16,32,64]` |
| Calibration gated on `mode.policy_adaptive` | `parallelization_performance/modes.jl` | static profiles pin a route; adaptive profiles keep calibration (worth 4.6-6.4x, §8b) |
| Calibration allotment ladder geometric in the budget | `engine/rhs_calibration.jl` | wider search space; **no measurable gain** (§8b.2) |
| Router: `outer_active` clamps width, not route | `engine/setup.jl` (3 branches) | 4918 → 1216 µs/RHS-call |
| R0 budget clamped to 1 | `engine/setup.jl` | R0 is a real serial baseline (see caveat below) |
| O(N) Jacobian prototype | `engine/execution.jl` | 62.5 → 2.2 ms at 4096 sat |
| Vacuum heat-rate short-circuit | `callbacks/thermal_callbacks.jl` | 6 h solve 10.22 → 7.88 s |
| Env-snapshot hoist | `engine/dynamics_rhs.jl`, `engine/setup.jl` | 448 → 0 bytes/call |
| Serial tuple peel | `engine/dynamics_rhs.jl` | part of the 25-45 % allocation cut |
| Flat-queue selection mask | `engine/dynamics_rhs.jl` | part of the same cut |
| ~~Flat-queue item peel~~ | — | implemented, measured, **reverted**: no wall-time effect and it erased a guard spelling `ci_rhs_parallel_route_gate` asserts |
| GRAM per-satellite warm-up | `engine/setup.jl` | path no longer hangs; 24-26 % slower, not recommended on |
| B7/B8 phases, heavy case family, 2 plots | `paper_parallelization_benchmarks/*`, `parallelization_performance/cases.jl` | workloads large enough to measure |

### Parity caveat: the gate was previously vacuous

Serial-vs-threaded parity is bit-identical **at a fixed thread count**. It is
*not* bit-identical when the harness runs its serial reference at
`maximum(threads)`, which is what `ppc_run_worker_parity` does.

That is a consequence of fixing R0, not a regression. Before the fix, R0 at
8 threads still threaded the flat queue, so the "serial" reference and the
threaded candidate were executing the same code and matched trivially. Now the
reference is genuinely serial (effectors summed in tuple order per satellite) and
the candidate is genuinely threaded (per-worker partials, then a worker-major
reduction). Both are deterministic; they group the sums differently, so they
round differently — measured `pos_rel_max = 1.02e-2` on
`multi_64_high_fidelity` over a 1200 s drag propagation, the amplification of
~1e-16 rounding.

Bisect confirming this (all at 8 threads, `multi_64_high_fidelity`): baseline
0.0; all changes 1.0198549452768975e-2; mask reverted *identical*; peel reverted
*identical*; all changes at 1 thread 0.0.

**Action required:** the parity gate needs a physics-based tolerance, or it
should test determinism within a mode and use a tolerance across modes. This is
directly review point 7's complaint that the manuscript "refers to tolerance
gates without defining their values or scientific basis" — there is now a
concrete case that forces the question.

### Not resolved

- **The thread ceiling.** Peak 1.6-1.8x at 4 threads on multi-effector
  workloads, inverting by 12. Allocation is not the cause (§7.3). Next step is
  wall-clock instrumentation of RHS phases.
- **Flat-queue reduction/zeroing** are O(workers x satellites) serial passes per
  RHS call — right shape, unmeasured magnitude.
- **`num_steps_to_save` is inert** (§10), so B6's output-cadence axis measures
  nothing.
- **`reset_fsal!`** costs an extra full RHS per accepted step (§10), untested.
- **GRAM isolated pool** is unreachable from the constellation RHS (§8.2).
- **B6's stale exclusions**: `multi_4_gram_live` and `multi_16_gram_live` no
  longer reproduce their documented failures (§8.1); confirm at full duration
  before re-enabling.
