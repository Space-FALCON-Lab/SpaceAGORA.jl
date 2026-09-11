# LEO GRAM Constellation Benchmark — Handoff

## Status

Built a serial-vs-threaded benchmark for a 2048-satellite LEO constellation (150 km
altitude, harmonics gravity + real native GRAM atmosphere + drag, 1-hour mission).
Along the way, found and fixed three real bugs in this codebase (not the benchmark
script), root-caused (but did not fix) one bug in the vendored GRAMSuite native
library, and fixed a routing inefficiency plus two smaller hot-path costs. Measured
result after those fixes: serial 434.6s, 8-thread parallel 395.2s (1.10x).

A follow-up re-profile of that 1.10x result (see "Bottleneck Profile" below) found the
real explanation wasn't residual serial overhead — it was a **second, independent
policy gate**: `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto` hardcodes a 16-thread
minimum before GRAM atmosphere sampling will ever thread, so at `--threads=8` the most
expensive part of the whole RHS ran on one thread the entire time. Setting it to `on`
instead validated at **361.7s** (8-thread parallel, down from 395.2s) — a further
~8.4% improvement — before hitting a *third*, structural ceiling: native GRAM's own
required global lock, which caps how much of that work can ever run concurrently
regardless of Julia-level thread budget.

## Entry Points

- Serial baseline: `julia --project=. --threads=1 benchmarks/studies/leo_2048_constellation_gram_scaling.jl serial`
- Parallel: `julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_scaling.jl parallel`
- Thread-count sweep (1/2/4/8/16/32/64), **not yet run** — several hours:
  `julia --project=. benchmarks/studies/leo_2048_constellation_gram_thread_scaling.jl`
- Re-profile with per-thread utilization breakdown (used for the "Bottleneck Profile"
  section below; also validates the density-callback fix via its `<auto|on>` arg):
  `julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_profile.jl <auto|on>`

## Source Fixes Made (all validated against the full test suite, 4039/4039 clean each time)

1. **`in_atmosphere` never initialized from actual state** — `src/simulation/engine/setup.jl`
   (`_initialize_in_atmosphere_flags!`) / `src/simulation/engine/execution.jl`. The flag
   defaulted to `false` for every satellite and only flipped on an altitude-crossing
   event; an orbit that starts (and stays) below `EI` never crosses, so it silently
   never got the vacuum-predicted GRAM cache or the finer `dt_max_atmosphere` step size.
   Now set from each satellite's starting altitude at solve setup.

2. **Real data race in aero force caching, unrelated to GRAM** —
   `src/simulation/engine/setup.jl` (`_initialize_save_cache_buffers!`) /
   `src/simulation/engine/execution.jl`. `SaveCache.drag_cache`/`lift_cache`/
   `cross_cache` (`src/core/types/runtime_types.jl`) started as empty `Vector`s and
   grew lazily via `resize!` in `_store_vector_cache!`
   (`src/dynamics/coupled/aerodynamic_wrench_models.jl:139`). `AerodynamicCoefficientfM`'s
   wrench is dispatched per-satellite across persistent worker threads whenever gravity
   is threaded across 2+ satellites (**any** density model, not just GRAM), so two
   satellites on different workers could both see the cache as too short and both
   `resize!` the same shared `Vector` concurrently. Manifested as a catchable
   `ConcurrencyViolationError` at small N and a GC segfault (real memory corruption) at
   N=2048. **This likely affects the existing paper benchmark cases too**
   (`multi_16_aero_surrogate_cached`, `multi_64_high_fidelity`, which use
   `ExponentialAtmosphereModel` and would have hit the same race) — worth an audit/rerun.
   Fixed by pre-sizing all three caches to `num_sats` before any threaded dispatch begins.

3. **Routing bug: 2-effector constellations never reach the flat queue** —
   `src/simulation/engine/setup.jl` (`_rhs_flat_has_batch_privileged_effector`,
   wired into `_rhs_execution_plan`'s flat-queue eligibility check). The
   `n_effectors >= _rhs_flat_min_effectors()` (default 3) gate silently excluded any
   2-effector constellation (e.g. harmonics + drag — a very common shape) from ever
   reaching `flat_constellation_effector_queue`, falling through to
   `_satellite_batch_saturates_pool` and forcing `:satellite_batch` once
   `active_sats >= budget`. That heuristic is blind to the fact that harmonics (and the
   inverse-square/J2 batch kernels added earlier this branch) get an algorithmic win
   under the flat route that `:satellite_batch` structurally cannot replicate at any
   thread count. Directly measured: ~1.7x faster per RHS call once routed correctly
   (26.75ms serial-mode vs 15.57ms flat-mode, same scenario). Fix bypasses the
   effector-count gate specifically when a harmonics or non-gravity-gradient
   inverse-square/J2 effector is present.

4. **Redundant per-call sort in the flat queue** —
   `src/simulation/engine/dynamics_rhs.jl` (`_prepare_rhs_flat_work_items!`). Sorted
   every RHS call even when only one effector type remains in the residual flat queue
   (the common case once harmonics is pre-passed away) — cost only varies by `eff_idx`,
   so with one residual type every item shares the same key and the sort is a
   guaranteed no-op that still pays full comparison cost. Now skipped via the existing
   `_count_flat_queue_only_effectors` helper.

5. **Env-var re-parsing in the per-satellite density hot path** —
   `src/core/types/runtime_types.jl` (new `SharedBuffers.density_static_config` field),
   `src/simulation/engine/setup.jl` (`_initialize_density_static_config!`),
   `src/simulation/engine/effector_sampling.jl`,
   `src/simulation/callbacks/density_callbacks/runtime.jl`.
   `_sample_atmosphere_from_planet_frame` was calling `_gram_track_cache_config()`
   (~10 env-var parses), `_gram_runtime_stats_enabled()`, and
   `_gram_track_cache_target_use_j2()`/`_uses_j2_gravity_effector()` fresh on every
   call — once per satellite per RHS call, millions of times over a large multi-sat
   mission. Same for the vacuum-cache accessors inside `_density_state_from_kinematics!`.
   All now computed once per solve and cached; both hot-path functions fall back to the
   raw (slow but correct) accessors if the cache wasn't populated, for safety with
   non-standard test harnesses.

Combined effect of (3)+(4)+(5): serial 568.2s → 434.6s (−23.5%), parallel 490.4s →
395.2s (−19.4%). The relative serial-vs-parallel gap didn't widen (1.16x → 1.10x)
because (5) benefits both modes about equally (it's on the shared aero/density path
regardless of RHS route), while (3)'s structural win competes against other
now-proportionally-larger serial-only costs (ODE solver bookkeeping, event/callback
overhead for 2048 satellites) that routing can't touch.

## GRAMSuite Native Library — Unresolved

**Bug:** the vendored native GRAM atmosphere model
(`GRAMAtmosphereModel`, `data/GRAMSuite.jl`) hangs indefinitely when 2+ satellites in
one constellation solve need the vacuum-predicted GRAM cache
(`SPACEAGORA_VACUUM_GRAM_CACHE`,
`src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`) to *rebuild*
(mission duration exceeds the cache horizon, default 600s).

- 1 satellite, any mission length: fine (tested up to 3600s / ~6 rebuilds).
- 2+ satellites, mission short enough that no rebuild is needed: fine.
- 2+ satellites, mission long enough to need a rebuild: hangs indefinitely, reproduces
  identically regardless of exact count (2/4/8 all hang the same way).
- Reproduces under fully single-threaded serial execution (zero Julia-level
  concurrency), so it's not a `Threads` race — it's sequential-call state corruption
  in the native library, most likely triggered by the vacuum cache's forward
  look-ahead queries (jumping ahead by up to the cache horizon for one satellite's
  rebuild, then a different satellite's rebuild querying an earlier time on the same
  underlying native state) violating an assumption of monotonically-increasing
  elapsed-time queries somewhere in the native atmosphere/perturbation model.

**Two targeted fixes tried and ruled out:**
- `SPACEAGORA_GRAM_PER_SAT_INSTANCES=1` (deepcopy'd model per satellite) — still
  hangs identically. The shared state isn't at the Julia object level (likely Fortran
  COMMON-block/module-global state `deepcopy` can't touch).
- `gram_perturbation_scales=0.0` (zeroing perturbation output) — still hangs. Whatever
  breaks happens inside `update!` itself, before output scaling is applied.

**Current workaround (in `leo_2048_constellation_gram_scaling.jl`, not a fix):** set
`SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S` and `_DEVIATION_M` both well past what the
mission can reach, so only the (proven-safe, satellite-count-independent) initial
build ever runs and no rebuild triggers.

**To pick this back up:** the fastest path is probably attaching to a hung repro
process (small N, e.g. 2 satellites, mission length just past the default 600s
horizon) for a real native stack trace, or reading the compiled
`data/GRAMSuite.jl/GRAM Suite 2.0/Build/lib/libGRAM.so` / underlying Fortran source
directly — both go beyond source-reading-plus-instrumented-Julia-repros and were
explicitly out of scope for this pass.

## Bottleneck Profile (re-profiled 2026-07-07, after the routing/sort/cache fixes)

Re-profiled the real 2048-sat/8-thread/3600s run with
`benchmarks/studies/leo_2048_constellation_gram_profile.jl` (`Profile.@profile` +
`Profile.print(groupby=:thread, C=true)`). The aggregate number was still
**"51% utilization across all threads and tasks"** — identical to the pre-fix
profile — which initially looked like the fixes hadn't moved the needle at all.

**That aggregate number is a red herring.** Breaking it down `groupby=:thread` shows
why: `Threads.nthreads()` (=8, the `--threads=8` "default" pool) is a *separate* set
of OS threads from thread 1 ("interactive"), and Julia auto-spawns additional parallel
GC helper threads ("foreign: gc") that Profile always reports at 100% utilization
regardless of real work, since it has no notion of "idle" for foreign threads. The
16-thread breakdown decomposed exactly as:

| Thread(s) | Role | Utilization |
|---|---|---|
| Thread 1 (interactive) | main/driver thread | 99% |
| Threads 2-9 (default ×8) | the actual `--threads=8` worker pool | **3% each** |
| Threads 10-16 (foreign: gc ×7) | Julia's auto-spawned GC helpers | 100% each (artifact — always reads 100%) |

`(0.99 + 8×0.03 + 7×1.00) / 16 = 51.4%` — matches exactly. **Both the pre- and
post-routing-fix "51%" numbers were dominated by this artifact**, not by real
application behavior; they were never a meaningful before/after comparison.

**Root cause of the real finding (8 of 9 compute threads idle):**
`auto_thread_min_budget(:density_callback)` in `src/parallel/policy/env_config.jl`
(line ~140) hardcodes `max(default_budget, 16)` — a **16-thread minimum** before
`SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto` will ever enable threading. In
`thread_policy_decision` (`src/parallel/policy/adaptive_decision.jl:12`):
```julia
auto_budget_allowed = mode != :auto || budget >= auto_thread_min_budget(source)
```
With `mode=:auto` and `budget=8 < 16`, this is `false`, forcing `use_threads=false`
*unconditionally, regardless of num_sats=2048*. Traced into
`_prefill_environment_samples!` (`dynamics_rhs.jl:1260-1262`): `worker_allotment=1`,
so the entire per-satellite GRAM atmosphere sample loop
(`_sample_atmosphere_from_planet_frame` — the expensive native MSIS/HWM/Ephemeris/
SPICE calls) runs as a plain serial loop on the calling thread for all 2048
satellites, every RHS call. Harmonics gravity and the aero-wrench residual queue have
no such gate and thread fine — only the single most expensive piece of the RHS never
parallelized in this configuration. This is the real explanation for the 1.10x (not
~8x) speedup: only the cheap parts of the workload were ever multithreaded.

**Fix validated:** set `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on` instead of `auto`
when the thread budget is below 16 (`mode==:on` short-circuits `auto_budget_allowed`
to `true` unconditionally). Re-ran the same scenario/script with `on`:

| Metric | `auto` (broken) | `on` (fixed) |
|---|---|---|
| Wall clock (8-thread parallel) | 395.0s | **361.7s** (−8.4%) |
| Threads 2-9 utilization (each) | 3% | **13%** (>4x) |
| Thread 1 utilization | 99% | 92% |
| Aggregate (misleading) | 51% | 56% |

Confirms the gate was real and the fix helps, but the improvement is capped well
short of proportional. The new dominant cost in the `on` profile is threads queuing
on native GRAM's own required global lock (`gram_use_global_lock()`/`_with_gram_lock`
in `GRAMSuite.jl`, showing up as `Base/lock.jl:335`/`wait_no_relock`, tens of
thousands of samples) — the vendored library isn't reentrant, so only one thread can
execute inside it at a time no matter the Julia thread budget. This is an inherent
ceiling for non-surrogate GRAM, not a further policy bug. Surrogate density models
(e.g. `ExponentialAtmosphereModel`, which don't share a native lock) would not hit
this ceiling and are the next thing to compare if more scaling is wanted.

Secondary, smaller-impact findings from the profile, not yet acted on:
- SPICE frame/ephemeris calls (`sxform_`, `frmchg_`, `tisbod_`, DAF file access) are a
  real but secondary cost (~single-digit % of samples); worth checking whether planet
  rotation gets recomputed per Runge-Kutta stage rather than cached per accepted step.
- Some GC/work-stealing-queue activity beyond what's explained by the env-var-parsing
  fix above; not isolated further.

## Suggested Next Checks

1. ~~Re-profile the fixed scenario~~ — done (see above); root cause found (density
   callback's 16-thread auto-mode gate) and empirically validated fix in hand
   (`SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on`).
2. Decide whether to land the `on`-instead-of-`auto` density-callback change as the
   new default for thread budgets below 16 (or just document it as a required env
   override for this class of benchmark), since `auto` silently serializes the most
   expensive part of the RHS below that threshold.
3. Compare against a surrogate density model (e.g. `ExponentialAtmosphereModel`) at
   the same satellite count/thread budget to quantify how much of the remaining gap
   is GRAM's global-lock ceiling specifically vs. other overhead — surrogate models
   aren't subject to that lock.
4. Audit/rerun the existing paper benchmark cases
   (`multi_16_aero_surrogate_cached`, `multi_64_high_fidelity`) now that the aero-cache
   race is fixed — their prior numbers may have been silently exposed to it.
5. Run the thread-scaling sweep (`leo_2048_constellation_gram_thread_scaling.jl`,
   1-64 threads) when there's a multi-hour window available — now also worth running
   with `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on` given the finding above, since
   `auto`'s 16-thread gate would otherwise mean the 16/32/64-thread points look
   artificially different in kind (gate open) from the 1/2/4/8-thread points (gate
   closed) rather than a clean scaling curve.
6. Decide whether to invest in the GRAMSuite native rebuild-hang root cause (see
   above) or continue accepting the horizon-widening workaround.
