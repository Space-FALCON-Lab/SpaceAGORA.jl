# Thread Allocation and GRAM Concurrency — Handoff

## Status

Follow-up to the constellation-size scaling work in this folder (`README.md`),
done on a lower-power machine than these scripts were originally sized for (11
physical cores). Covers three new studies (thread-count scaling by
density-model mode, a Monte Carlo thread-allocation grid search, and an
outer-parallel-route comparison) plus a real, previously-mischaracterized bug
found and root-caused along the way: **GRAM model construction is not
thread-safe**, independent of and in addition to the already-known
isolated-CSPICE kernel issue (`project_gram_mars_isolated_cspice_bug`
memory).

## Entry points

- Thread-count scaling by mode, monolithic route (one coupled state vector):
  `julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_thread_scaling_by_mode.jl`
- Same script, ensemble route (`run_constellation_ensemble`, one worker per
  satellite):
  `SPACEAGORA_TS_ROUTE=ensemble julia --project=. .../leo_thread_scaling_by_mode.jl`
  Override `SPACEAGORA_TS_N_SATS` / `SPACEAGORA_TS_THREADS` (comma list) as needed.
- Monte Carlo thread-allocation grid search (8 samples, 4 sats/sample):
  `julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/mc_multisat_thread_allocation.jl`
  Results already saved to `mc_multisat_thread_allocation_summary.csv` /
  `_samples.csv` in this folder — do not need to be rerun to analyze.
- Constellation-size scaling (`leo_constellation_size_scaling.jl`) unchanged in
  behavior from `README.md`; `leo_constellation_size_scaling_worker.jl` gained
  an optional third `route` arg (`monolithic|ensemble`, default `monolithic`)
  and a `SPACEAGORA_TS_SURROGATE_DENSITY_CALLBACK_OVERRIDE` env hook (see
  below) for the newer scripts to reuse.

## Finding 1 — the density-callback auto-gate silently caps the surrogate too

`auto_thread_min_budget(:density_callback)` (`src/parallel/policy/env_config.jl:140-154`)
defaults to a **16-thread minimum** before `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto`
ever threads density sampling — regardless of density-model type. This was
already known to hurt native GRAM (see
`benchmarks/studies/LEO_GRAM_CONSTELLATION_HANDOFF.md`), but it turns out to
gate the **lock-free surrogate** too, even though nothing about the surrogate
requires a 16-thread floor. Monolithic route, N=512 sats, 600s mission:

| Threads | gravity-only | surrogate (`auto`) | surrogate (forced `on`) | native GRAM |
|---|---|---|---|---|
| 1 | 1.00x | 1.00x | 1.00x | 1.00x |
| 2 | 2.55x | 1.88x | **2.40x** | 1.09x |
| 4 | **3.50x** (peak) | 1.88x | **2.84x** (peak) | 0.96x |
| 8 | 2.61x | 1.49x | 1.98x | 0.91x |

Forcing `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on` (bypassing the gate) raises
the surrogate's peak from 1.88x to 2.84x. Native GRAM never clears ~1.1x at any
thread count — its ceiling is the vendored library's own required global lock
(structural, not a policy bug). Plot: `leo_thread_scaling_by_mode.png`.

## Finding 2 — ensemble route: no-GRAM/surrogate scale well, native GRAM crashes

Ensemble route (`run_constellation_ensemble`, one worker per satellite), N=64
sats, 600s mission:

| Threads | gravity-only | surrogate | native GRAM |
|---|---|---|---|
| 1 | 1.00x (15.5s) | 1.00x | 1.00x (15.5s) |
| 2 | 1.80x | 1.82x | **crashed** |
| 4 | **2.92x** (peak) | **2.88x** (peak) | **crashed** |
| 8 | 1.58x | 2.67x | **crashed** |

The surrogate does *better* here than under the monolithic route (2.88x vs.
2.84x best-case) because single-satellite ensemble members never hit the
density-callback gate in the first place. Native GRAM crashes at every
threads≥2 — see Finding 4 for the actual root cause (not what was originally
suspected). Plot: `leo_thread_scaling_by_mode_ensemble.png`.

## Finding 3 — Monte Carlo thread allocation: outer-only wins, mixed is catastrophic

8-sample MC campaign, 4 satellites/sample, gravity-only (no_gram) vs GRAM
surrogate, sweeping every `(outer_workers, inner_threads)` split of each total
thread budget T ∈ {1,2,4,8}. Full raw data:
`mc_multisat_thread_allocation_summary.csv` / `_samples.csv`.

**Gravity-only — best allocation is always all-outer, near-linear through T=4:**

| T | best split | wall time | speedup |
|---|---|---|---|
| 1 | 1×1 | 0.132s | 1.00x |
| 2 | 2×1 | 0.069s | 1.92x |
| 4 | 4×1 | 0.035s | 3.78x |
| 8 | 8×1 | 0.030s | 4.42x |

T=8's shortfall from linear is a floor effect (30ms is close to fixed
per-campaign task-spawn/GC/scheduler overhead at this workload size), not a
real bottleneck — expect the crossover to push higher with more
satellites/longer missions.

**Any inner-thread allocation is bad; outer+inner combined is catastrophic**
(not just "no benefit" — genuine pathological contention):

| Config | Wall time | vs. serial |
|---|---|---|
| outer=1, inner=4 (T=4) | 0.75s | 5.7x slower |
| outer=1, inner=8 (T=8) | 0.72s | 5.5x slower |
| **outer=2, inner=2 (T=4)** | **20.4s** | **155x slower** |
| **outer=2, inner=4 (T=8)** | **185.8s** | **1400x slower** |

Likely mechanism: `run_monte_carlo`'s outer workers already occupy threads via
`Threads.@spawn`, blocked in `Threads.@sync` waiting on children; when each
*also* spawns nested inner-parallel tasks onto the same fixed-size default
thread pool, parent-waiting-on-child scheduling contention compounds instead
of just adding overhead. Gets worse, not just flat, as both axes grow —
consistent across two independent (T, split) data points. **Recommendation for
the paper: give the full thread budget to outer parallelism; never split it
with inner parallelism for this workload shape.**

**GRAM surrogate**: every combo with `outer_workers≥2` crashed — including the
pure-outer split that was gravity-only's best allocation. Only `outer=1`
survived (serial dispatch, any inner_threads), and inner-only there is worse
than serial (0.99s serial vs. 2.14–3.71s for inner=2/4/8) — same
"not-enough-per-sample-work" story as gravity-only's inner-only numbers. See
Finding 4 for why this crashes even though the surrogate has no native lock.

## Finding 4 — root cause of the crashes: unlocked GRAM model *construction*, not sampling

Full detail in memory: `project_gram_concurrent_construction_crash` (corrects
an earlier same-session hypothesis that blamed "native GRAM's own required
global lock" / the isolated-CSPICE atmosphere-sampling path — that path is
already correctly locked everywhere).

**Actual bug**: `GRAMSuite.GRAMAtmosphereModel`'s keyword constructor
(`data/GRAMSuite.jl/src/GRAMSuite.jl:408`) calls `gram.initialize!(spice_root)`,
`gram.create_atmosphere(...)`, `gram.set_start_time!(...)` etc. with **zero
locking**, unlike every density-query call site (`getDensity`/`density_state`/
`point_density_state`), which are all correctly wrapped in
`GRAM_LOCK`/`_with_gram_lock`. `gram.initialize!` furnishes SPICE kernels into
GRAM's private, process-global CSPICE instance; CSPICE kernel-furnishing is not
reentrant.

Confirmed with a minimal, simulation-free repro (not kept in-tree, reproducible
inline below): two threads each constructing a fresh
`GRAMAtmosphereModelSurrogate(planet_name="earth")` concurrently crash
identically (`SPICE(BADSUBSCRIPT)`/`trcpkg` or `libc++abi` abort), even after a
serial warmup construction first (rules out the unrelated Julia 1.12
world-age/dynamic-module-loading warning as a confound):

```julia
include("examples/common.jl"); ensure_gramsuite_loaded!()
_ = SimulationModel.GRAMAtmosphereModelSurrogate(planet_name="earth")  # serial warmup
Threads.@sync for i in 1:2
    Threads.@spawn SimulationModel.GRAMAtmosphereModelSurrogate(planet_name="earth")
end
# crashes with --threads=2
```

`Base.deepcopy_internal` for these model types **is** properly
`lock(GRAM_LOCK)`-wrapped and does not re-run `gram.initialize!` — so
duplicating an *already-constructed* model across threads (what
`run_constellation_ensemble` does for member isolation) is safe. The hazard is
specifically building *new* models concurrently, which is exactly what
per-sample/per-member config construction does in a naively-written Monte
Carlo/ensemble worker (including `mc_multisat_thread_allocation_worker.jl`'s
`build_sample_config`, deliberately reproducing the hazard for this study).

**Why this never showed up in earlier separate-process runs**: process
isolation means each worker process loads its own independent copy of
`libGRAM.dylib`/CSPICE — no shared memory, so "concurrent" construction across
processes never actually races on anything. `Threads.@spawn`-based outer
parallelism (`run_monte_carlo`, `run_constellation_ensemble`) shares one
process's address space, hence one loaded copy of that native global state,
which is what exposes the bug. The paper benchmark suite's `outer_process`
mode (`benchmarks/studies/parallelization_performance/execution.jl:74-96`,
`Distributed.addprocs`) is genuinely safe for this reason and is the existing
safe fallback for outer-parallel GRAM/surrogate campaigns until the
constructor is fixed.

**Not yet fixed** — `data/GRAMSuite.jl` is a git submodule; the fix would wrap
the constructor's `gram.initialize!`/`create_atmosphere`/`set_start_time!`
calls the same way `_with_gram_lock` already protects
`_gram_density_state_native`.

## Finding 5 — follow-up session: constructor lock fix, process backend, and what changed

Direct follow-up to items 1–2 above, done in a later session. Two library
changes, both now in place (submodule + parent repo, see
`project_parallelization_tier1_status` memory for exact commits):

- **GRAM constructor lock** (item 1): `_with_gram_lock`'s `lock_obj===nothing`
  branch (the one `GRAMAtmosphereModel`'s constructor now runs under) wraps the
  whole constructor body in the same `GRAM_LOCK`/`SPICE_LOCK` used everywhere
  else, closing the concurrent-construction hazard from Finding 4. A second,
  related fix in the same submodule: several `gram.xxx(...)` calls inside the
  constructor and `_gram_apply_user_ephemeris_state!` were missing their own
  `Base.invokelatest` wrap around the *call* (not just the `getfield` lookup) —
  needed specifically for the new process-backend path below, where the first
  construction on a fresh worker can happen inside an already-compiled
  Distributed closure. See `feedback_julia_distributed_gotchas` memory, point 5.
- **`SpaceAGORA.ParallelProcess`** (item 2, formalized): a shared,
  auto-bootstrapping `Distributed` worker pool (`src/parallel/process/`),
  wired into `run_monte_carlo`/`run_constellation_ensemble`'s adaptive
  `threads=:auto` routing as a real `:process` route — no more silent
  downgrade to `:threads`/`:none` for native-GRAM workloads.

**Re-verification of Finding 2's ensemble crash and Finding 3's MC surrogate
crash**: both confirmed fixed. `mc_multisat_thread_allocation.jl` now includes
a `"standard"` (native GRAM) mode in its full outer×inner threads grid
(previously excluded entirely because it crashed) — every combo up to
outer=8×inner=1 now completes 8/8 successfully; see
`mc_multisat_thread_allocation_summary.csv` (`mode=standard,outer_backend=threads`
rows). Native GRAM stays lock-limited (flat ~2.1–6.2s) except one point
(outer=4×inner=2: 146s) that clearly still hits Finding 3's nested-contention
pathology, just far less severely than gravity-only's worst case (see below) —
it degrades badly but does not hang past the 300s harness timeout or crash.

**Finding 3's pathology, revisited**: still present and unfixed (as expected —
out of scope for this follow-up), but its severity turns out to be
**workload-dependent, not universal**, refining the original characterization:

| Mode | outer=2×inner=2 | outer=2×inner=4 | outer=4×inner=2 |
|---|---|---|---|
| no_gram (gravity-only) | timed out (>300s) on 2 of 3 runs, 2.7–92.6s on the third | 26.6–57.1s (vs. ~0.03s outer=4×inner=1 baseline) | timed out (>300s) every run observed, once left running 16+ min before being force-killed |
| surrogate (lock-free) | 1.67s (mild, ~3x over best split) | 4.75s (mild) | 41.6s (moderate) |
| standard (native GRAM, lock-limited) | 3.0s (negligible) | 6.2s (mild) | 146.3s (moderate-severe, but bounded) |

Gravity-only is consistently the *worst* affected, not the mildest as might be
expected from having no shared lock at all — and its timing is highly
non-deterministic run-to-run (the same outer=2×inner=2 point ranged from 2.7s
to a 300s timeout across repeated runs in the same session). This looks more
like a scheduler/contention livelock than a fixed-cost slowdown; worth flagging
for whoever picks up Finding 3 as a real fix rather than a "give outer
parallelism the full budget" workaround; the workaround itself still holds
(never split the thread budget between outer and inner for any mode tested).

**New benchmark: `leo_ensemble_outer_backend_scaling.jl`/`_worker.jl`** — fixed
N=64-satellite constellation, native GRAM, ensemble route, sweeping
`outer_workers ∈ {1,2,4,8}` × `{threads, process}`. Clean run (see
`leo_ensemble_outer_backend_scaling_summary.csv`):

| outer_workers | threads | process |
|---|---|---|
| 1 | 1.00x (14.5s) | 1.00x (15.2s) |
| 2 | 0.98x | **1.20x** |
| 4 | 0.98x | **1.44x** |
| 8 | 0.93x | 0.75x (20.2s — regresses) |

`threads` stays flat at ~1x regardless of worker count (GRAM's process-wide
lock caps it, as in Finding 1/2) — `process` clearly breaks past that ceiling
up to 4 workers, then degrades at 8. The 8-worker regression is plausibly
this machine's own limit (11 physical cores; 8 process workers + 1 coordinator
= 9 concurrent OS processes, each doing real work, leaves little headroom) —
not retested on a larger machine.

**New benchmark: `mc_multisat_thread_allocation.jl`'s process backend** — less
clean. The full sweep's own timing for `outer_backend=process` at
`outer_workers≥2` is contaminated by a benchmark-harness issue (a
precompilation-cache invalidation this session accidentally triggered while
debugging an unrelated subprocess-hang bug in the sweep controller itself —
every fresh process-backend point ended up paying a ~190–380s one-time
recompile cost inside the harness's 300s per-point timeout, so most of those
rows are recorded as `TIMEOUT` in the CSV, not real measurements). Direct,
uncapped single-point runs (bypassing the sweep controller) give the real
numbers instead:

| outer_workers | wall_time_s | n_ok |
|---|---|---|
| 1 (serial baseline) | 2.23s | 8/8 |
| 2 | 2.45s | 8/8 |
| 4 | 2.43s | 8/8 |
| 8 | 2.45s | 8/8 |

Unlike the ensemble study above, this **does not show scaling** with worker
count — every value clusters around the serial baseline. `Distributed`'s own
concurrency was verified working correctly in this environment (an 8-way
`@sync`/`@async remotecall_fetch(sleep, w, 1.0)` probe across the same 8
warm workers took ~2s, not ~8s), so the flat timing isn't a routing or
dispatch bug — `select_outer_route!` also confirmed picks `:process`
unconditionally for native GRAM regardless of sample count (see
`default_outer_route`'s `_is_native_gram_point_density` short-circuit). Root
cause not confirmed; leading hypothesis is that each sample's
`build_sample_config` constructs a fresh `GRAMAtmosphereModel` (reads the
MERRA2 data file) on every dispatch, and that per-process file read may
serialize across the OS the same way the in-process `GRAM_LOCK` used to —
consistent with this 8-samples/~0.3s-each shape being small enough that
per-sample construction cost dominates, unlike the ensemble study's N=64,
no-per-member-reconstruction shape. **Flagged for follow-up, not resolved.**

## Finding 6 — Finding 3 root-caused and fixed; Finding 5's process-scaling hypothesis refuted, real cause found

Direct follow-up session to Finding 5's open items 6 and 7.

**Finding 3 (nested outer+inner contention): root-caused and fixed.** It was not
a fundamental Julia scheduler limitation with no available fix, as the
"livelock-like, out of scope" framing in Finding 5 assumed. `outer_parallel_active()`
turns out **not** to be a uniformly-enforced kill-switch across this codebase: it
correctly gates the density/control/thermal callback paths and the dynamic-effector
thread decision, but three multi-worker auto/forced-routing branches inside
`_rhs_execution_plan` (`src/simulation/engine/setup.jl`) never checked it at all —
the harmonics-batch/flat-constellation route (the exact mechanism
`SPACEAGORA_HARMONICS_BATCH_ENABLED=1` triggers, used by default in this study),
plus two structurally identical siblings (`single_invsq_flat`, and the
"many heavy effectors" flat-queue heuristic). Since this route fires on *every*
RHS/ODE step rather than once per sample, an outer worker already blocked in
`Base.@sync` was repeatedly spawning its own nested batch of `Threads.@spawn`-blocked
children throughout the whole integration — exactly the observed pattern. Fixed by
adding the same `outer_active && !allow_with_outer → force serial` guard already
used elsewhere in the same file to all three branches (new env var
`SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER`, default `false`). Directly
re-verified: the two worst repro points (`no_gram` outer=2×inner=2, previously
2.7s–300s+ non-deterministic; outer=4×inner=2, previously once left running 16+
minutes) now run consistently at ~0.04–0.07s across repeated trials — matching the
outer-only baseline, no residual slowdown. Full test suite re-run clean (4059
pass, 1 expected-broken, same as the pre-fix baseline — no regressions from
touching `setup.jl`).

**Finding 5 item 7's hypothesis (per-sample GRAM/MERRA2 reconstruction
serializing across processes) is refuted, and the real cause is bigger and more
interesting.** A direct probe (construct `GRAMAtmosphereModel` N=1 serially vs.
N=8 concurrently across warm `SpaceAGORA.ParallelProcess` pool workers) showed
construction itself is fast either way (sub-millisecond, no cross-process
contention) — ruling this out entirely. A second, corrected probe running the
*full* `build_sample_config` + `run_simulation` pipeline (properly wrapped with
the same `withenv(env_pairs(...))` the official worker script uses — an earlier,
uncorrected version of this probe hit a known, documented hang from skipping that
wrapper: native GRAM's vacuum-predicted cache hangs on rebuild with 2+ satellites
without it) found the actual mechanism: a worker's **first-ever** `run_simulation`
call pays a massive, one-time JIT/specialization cost — one already-warm worker
(reused from earlier serial calls) ran a sample in 0.44s, while seven genuinely
cold workers, dispatched to concurrently in the same probe, each took ~58s for
their first call, all succeeding but ~130x slower. `SpaceAGORA.ParallelProcess`'s
existing per-worker warmup (`_warm_gram_wrapper!`) only exercises GRAM's native
wrapper module (construction + one density query) — it does not touch the much
larger JIT surface of a full `run_simulation` call (ODE solver dispatch, RHS
assembly, effector/harmonics dispatch, etc.), which is specialized per exact
`SimulationConfiguration` type-parameter combination. This is a real, previously
unknown-to-this-project gap, but **not a small/safe fix**: a generic pool-level
warmup can't usefully pre-specialize for every type combination a future caller
might use (Julia specialization requires an exact type match, not "close
enough"), so a correct fix needs to be shaped by the specific caller (e.g. the
study's own worker-warmup step running one full representative call per new
worker, matching that study's exact types) — left for follow-up per the
"investigate and report, don't chase into an unclear fix" call for this item.
This also means `mc_multisat_thread_allocation`'s own two-pass
warmup-then-timed design may not be fully absorbing this cost across all workers
in the timed pass either (unconfirmed) — worth checking directly before trusting
its process-backend numbers at face value even after a clean sweep re-run.

**Full sweep re-run after the Finding 6 fixes landed** (see item 8 below,
now done) confirms the threads-backend fix comprehensively across every mode,
not just the two directly-repro'd points:

| mode | point (outer×inner) | before Finding 6 | after Finding 6 |
|---|---|---|---|
| no_gram | 2×2 | 2.7s–300s+ (non-deterministic) | 0.068s |
| no_gram | 2×4 | 26.6–57.1s | 0.069s |
| no_gram | 4×2 | timeout every run (once 16+ min) | 0.050s |
| surrogate | 2×2 | 1.67s | 0.623s |
| surrogate | 4×2 | 41.6s | 0.528s |
| standard | 4×2 | 146.3s | 2.286s |

Every nested `threads`-backend point across all three modes now lands within
a normal, bounded range matching (or beating) its own mode's outer-only
baseline — no timeouts, no multi-second-to-multi-minute outliers anywhere in
the 30 `threads`-backend rows. `process`-backend rows at `outer_workers>=2`
are still all `TIMEOUT`/`FAILED` in the raw CSV, consistent with the
cold-worker JIT cost described above (not a regression, and not something
Finding 6's fixes were meant to address) — `outer_workers=1` succeeds in all
three modes. Raw data: `mc_multisat_thread_allocation_summary.csv` /
`_samples.csv` (overwritten with this re-run's numbers).

## Finding 7 — Finding 6 item 7's cold-worker JIT cost: fixed, with a real diagnosis correction along the way

Direct follow-up to Finding 6's item 7. Implemented an optional `warmup_fn`
parameter on `ensure_process_workers!` (`src/parallel/process/worker_pool.jl`):
a zero-argument closure, `remotecall_fetch`'d best-effort on each *newly added*
worker only, after the existing generic bootstrap. `_run_campaign_with_route_env`'s
`:process` branch (`src/simulation/campaigns/adaptive_routing.jl`) now passes
`warmup_fn=() -> f(first(spec.seeds))` — reusing the campaign's own real
per-sample closure, so a newly grown worker's first-ever JIT/specialization
cost is paid once during pool growth instead of landing inside the caller's
first real (measured) dispatch. Fully backward compatible (default `nothing`);
no existing call site (library or benchmark script) needed changes.

**Validating this took much longer than the fix itself, and surfaced its own
lesson worth recording.** Every early verification attempt showed a confusing
"exactly one sample stays slow, the rest are fast" pattern that looked like a
bug in the new warmup mechanism. It wasn't — it was three stacked test-harness
mistakes, each masking the last:

1. A test closure referencing a benchmark-script-local helper
   (`build_sample_config`) as `warmup_fn` fails with `UndefVarError` on a bare
   worker (the closure-resolvability rule already documented in
   `feedback_julia_distributed_gotchas` memory point 1 applies equally to a
   warmup closure, not just real-dispatch ones) — but `ensure_process_workers!`'s
   own `try`/`catch` swallowed this correctly and silently, so the campaign
   kept running *without* the intended warmup ever having happened.
2. With that closure instead written using bare `SpaceAGORA` type names
   (`Earth(...)`), those *also* don't resolve on a worker bootstrapped only via
   `using SpaceAGORA` — `Earth` etc. live in the `SimulationModel` submodule
   and need either full qualification (`SpaceAGORA.SimulationModel.Earth`) or
   the worker to separately `using SpaceAGORA.SimulationModel`.
3. Most consequentially: a test passed `density_family="standard"` to
   `campaign_route_features` (matching the benchmark study's own `mode`
   naming), but the routing logic's `_is_native_gram_point_density` check
   actually matches the string `"gram_point"` (what `_campaign_density_family`
   produces for a real `GRAMAtmosphereModel`). With the wrong string, routing
   silently fell through to `:none` on this `:small`-classified machine —
   meaning **the entire campaign, across every early test, ran serially on the
   coordinator process**, never touching the process pool or the new
   `warmup_fn` code at all. The "one slow sample" was just the coordinator's
   own first-ever JIT cost; the "fast" ones were the same, now-warm,
   coordinator reused for the rest.

Adding a temporary trace print directly in `_campaign_route_plan` (confirming
`route=:none` when it should have been `:process`) is what surfaced #3
directly rather than continuing to guess. Once corrected (self-contained
closure, fully-qualified names, correct `density_family="gram_point"`), the
fix validated cleanly: all 4 newly-added workers' `warmup_fn` calls succeeded
sequentially (~57s each, matching the previously-measured cold cost) during
pool growth, and all 4 real, timed dispatch samples then came back uniformly
fast (~3.8s each, no asymmetry) — proving the JIT cost was fully absorbed
before the caller's measured work began. Full test suite re-run clean (no
regressions from touching `adaptive_routing.jl`/`worker_pool.jl`).

**Scope note, unchanged from Finding 6**: this fixes the cold-start case (a
caller dispatching to a fresh pool with no warmup discipline of their own). It
does not and cannot address the separate CPU-oversubscription ceiling on this
11-core dev machine (Finding 6's pass-2 instrumented data) — that remains a
hardware-capacity fact, not a code bug.

## Suggested next checks

1. ~~Implement the constructor lock fix in `data/GRAMSuite.jl` (submodule) and
   re-verify the MC surrogate outer≥2 combos and the ensemble native-GRAM
   crash both go away.~~ **Done — see Finding 5.**
2. ~~Until fixed, route any outer-parallel GRAM/surrogate campaign through
   `outer_process` (separate OS processes) instead of thread-based
   `run_monte_carlo`/`run_constellation_ensemble`.~~ **Formalized as
   `SpaceAGORA.ParallelProcess` — see Finding 5.**
3. Re-run the thread-count-by-mode and MC allocation studies at full scale
   (1024+ sats, full 1–64 thread ladder) now that the scripts exist and are
   debugged — tonight's numbers used cut-down sizes (512/64/4 sats) to fit
   this machine's lower core count and available time.
4. ~~Consider lowering `SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET`'s
   16-thread default...~~ **Done — see Finding 6: the gate is now
   density-model-aware (`:density_callback` vs. `:density_callback_lockfree`
   sources), native GRAM's floor is unchanged, the lock-free surrogate now
   gets the general 4-thread default.**
5. If repeating the 16-thread oversubscription point from the earlier
   (superseded) version of the thread-scaling plot on this specific machine,
   use 8 or 11 (physical core count) instead of 16 — 16 threads on an
   11-core machine is pure oversubscription and collapses all three modes
   uniformly, which isn't a meaningful data point.
6. ~~Root-cause Finding 3's nested outer+inner contention properly...~~
   **Done — see Finding 6.**
7. ~~Root-cause why `mc_multisat_thread_allocation`'s process backend doesn't
   scale with outer_workers...~~ **Done — see Finding 7. Original hypothesis
   (per-sample GRAM reconstruction) refuted; real cause (per-worker first-call
   JIT cost, ~57-70s) identified and fixed via `ensure_process_workers!`'s new
   `warmup_fn`.**
8. ~~Re-run `mc_multisat_thread_allocation.jl`'s full sweep now that Finding
   6's fixes are in place...~~ **Done — see Finding 6's before/after table.
   All 30 `threads`-backend rows now bounded; `process`-backend rows at
   `outer_workers>=2` remain `TIMEOUT`/`FAILED`, consistent with the
   cold-worker JIT cause (item 7 as it stood then), not a regression.**
9. ~~Investigate whether `mc_multisat_thread_allocation_worker.jl`'s existing
   warmup-then-timed two-pass design ... actually JIT-warms every worker...~~
   **Done — see Finding 6's pass-1/pass-2 instrumented data. Job distribution
   during warmup is confirmed fair (all 8 workers get exactly one job each);
   the flat timed-pass numbers are CPU oversubscription on this 11-core
   machine (each of 8 concurrent real-compute samples taking ~2.4-2.8s instead
   of ~0.3s), not unfair warmup distribution or residual JIT cost.**
10. Consider whether `mc_multisat_thread_allocation_worker.jl` and
    `mc_multisat_process_backend.jl`'s own hand-rolled two-pass
    warmup-then-timed design (and `ensure_process_workers_with_sample_logic!`)
    are now partially redundant with `ensure_process_workers!`'s built-in
    `warmup_fn` (Finding 7) and could be simplified — not required, flagged as
    optional cleanup only.
