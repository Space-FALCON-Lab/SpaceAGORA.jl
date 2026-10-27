# Native GRAM under threads: the global lock and the isolated pool

Native GRAM calls normally serialize on `RuntimeServices.GRAM_LOCK`, which is
the same object as `SPICE_LOCK`; the shared lock protects the native GRAM and
CSPICE state. The isolated pool instead evaluates an eligible density batch
across workers, each with its own model copy and its own lock.

This page covers the current controls and how to test them. The original
measurements, and their limits, are kept in the lab's research records (see
[Record a result](../src/user/studies_benchmarks.md#record-a-result)). A
historical result does not establish performance or scientific equivalence at a
later revision.

## Configuration and eligibility

The settings are read in
`src/orchestration/simulation/callbacks/density_callbacks/config.jl`.

| Setting | Default | Meaning |
|---|---|---|
| `SPACEAGORA_GRAM_ISOLATED_POOL` | `auto` | `auto` uses the pool only when the batch size, pool width and wind-history checks allow it and Julia has more than one thread. `off` keeps the locked route. `on` accepts independent instance histories. |
| `SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD` | `1024` | Minimum number of items in one callback batch that actually reach native GRAM before `auto` pools them. |
| `SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS` | `min(4, Threads.nthreads())` | Maximum number of worker instances; the pool needs at least two. |

These defaults came from a nominal-wind study and do not promise a speedup
elsewhere; each extra instance also holds its own native memory. The threshold
counts native calls, not configured spacecraft: items that return vacuum or use
the above-interface density approximation do not justify building a native
pool. The pooled native call applies the same `-30.0` m altitude floor as the
locked scalar path.

The density callback asks for the pool's width with `lock_free=true`. That
skips the 16-thread minimum budget the shared-lock path uses for native GRAM
(`SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET`), so the pool is
reachable with fewer than 16 threads. It is still bounded by the available
thread budget, the worker cap and its other eligibility checks.

## Wind history and identity

Under `auto`, the pool declines a batch that requests winds when the native core
reports history-dependent winds, before it builds or warms any worker. Unknown
cores are treated as history-dependent. The pinned GRAM wrapper's default wind
mode is perturbed, so wind queries reach the automatic pool only when nominal
winds are selected explicitly. Calls with `wind=false` stay eligible. The check
applies even when winds are requested only for diagnostics and the force model
does not use them.

`SPACEAGORA_GRAM_ISOLATED_POOL=on` bypasses the wind-history check. Separate
instances then follow separate stochastic histories, and winds and trajectories
can depend on pool width or thread count. It is not a speed switch that keeps
perturbed-wind runs reproducible. Routing never changes an atmosphere mode or
seed.

The recorded bit-identity comparisons cover specific nominal-wind cases. They
do not qualify perturbed winds, other model versions or later code without a
fresh comparison. Two tools in `benchmarks/studies/gram_thread_scaling/` support
one: `density_grid_parity.jl` compares exact floating-point bits on a density
grid, and `dump_states.jl` writes saved times and full spacecraft states for a
byte comparison.

## Which paths the pool can reach

| Path | Relationship to the pool |
|---|---|
| Density callback batch (`density_callbacks/runtime.jl`) | Evaluates eligible batches through the isolated pool, falling back to the locked route. |
| Vacuum-predicted look-ahead cache (`SPACEAGORA_VACUUM_GRAM_CACHE`) | Builds its own knots with scalar calls on the shared model. The callback's staged density batch is a separate operation. |
| GRAM track cache (`SPACEAGORA_GRAM_TRACK_CACHE`) | Its fallback can use the pool when the GRAM driver lacks a native `generate_trajectory`; the native trajectory path is a locked call. |
| Constellation RHS aerodynamic sampling | Uses the model chosen by `_density_model_for_sat`, not this pool. |

Two other instance-isolation controls serve different callers:

- `SPACEAGORA_GRAM_PER_SAT_INSTANCES` gives each spacecraft its own model for
  RHS sampling. Setup warms the fresh native instances one at a time.
- `SPACEAGORA_GRAM_LOCK_SCOPE=model` makes scalar calls take the wrapper
  instance's lock instead of the global one. A single shared model still
  serializes on its one instance lock. The isolated pool supplies its own locks
  and ignores this setting.

Model-scoped scalar locks are not counted in the shared native-lock counters,
so a small shared-lock total does not by itself show a speedup.

## Benchmark configuration and interpretation

Before timing a GRAM case, check how many members actually call native GRAM,
whether they are inside the entry interface, and whether the intended cache is
built. Otherwise a large constellation can mostly measure vacuum or a density
approximation while keeping a GRAM case name.

The density traces of the parallelization benchmark (traces 4, 5 and 6) fly
`ppc_p6_gram_constellation` from
`benchmarks/studies/parallelization_performance/cases.jl`: members in a
300–480 km band, below the 600 km entry interface `PPC_P6_GRAM_EI_KM`, so every
member is inside the atmosphere from the first step. The three traces share that
constellation and differ only in the density model and its access path; trace
6 flies trace 5's members as one-spacecraft samples. Rows measured before
traces 4, 5 and 6 moved to this shared constellation describe different work and
must not be pooled with current rows; a case name alone does not identify
equivalent work.

For a new comparison, record the revision, wrapper and data versions, wind mode,
seed, density path (including `SPACEAGORA_DENSITY_FREEZE_PER_STEP`), the actual
number of native items and the pool width. Alternate the locked and pooled arms,
keep compilation and setup out of the timed region, and compare repeated
timings together with exact state dumps. Include a one-thread control, where
the pool must decline. Separate serial native work from lock wait: fewer lock
acquisitions do not prove that contention was the limiting cost. Keep the raw
rows and machine details with the study's record.

## Refresh de-phasing: not available

Cache horizons set where the interpolation knots fall. Jittering them changes
the positions and times queried from GRAM, and so the spline and the densities
the RHS reads; it is a scientific change, not a scheduling-only optimization.
Staggered refreshes are neither implemented nor qualified. Batching or
prefetching would need its own validation that the sampled points and state
ownership stay the same.

## Pool-worker memory growth

A Julia heap-size hint does not bound memory owned by a native atmosphere.
Changing the epoch of a GRAM model, or copying one, can build a new native image,
and dropping the small Julia wrapper does not free that image until its
finalizer runs. A sequence of discarded models can therefore grow native memory
with no comparable growth in the Julia heap.

The parallelization benchmark caches its GRAM model at `ppc_initial_time()`,
the epoch its configuration builder also uses, so the run's epoch alignment
returns the same model instead of rebuilding it. The GRAM extension
(`ext/SpaceAGORAGRAMSuiteExt.jl`) counts the native atmospheres it creates and
finalizes. `run_simulation` calls `collect_unreferenced_gram_atmospheres!()`
after every run, successful or not, and the extension requests a full
collection once the number of additional live atmospheres reaches
`SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT` (default `8`; `0` or `off` disables it).
This triggers a collection; it is not a hard limit on resident memory. Live
models stay owned by their callers, and the allocator can keep resident memory
high after objects are freed.

`worker_growth.jl` separates Julia allocations, native growth, model
construction and finalization, and its state dumps support forward and reverse
sample-order comparisons. `worker_rss.jl` runs explicit collections, which you
must account for when reading its numbers. The tests in
`test/unit/environment/gram_native_collect_tests.jl` and
`test/unit/environment/gram_isolated_pool_tests.jl` cover these contracts;
claims about native memory or trajectories still need the actual wrapper,
inputs and scenarios. See
[Collector settings](../../benchmarks/studies/heap_contention/FINDINGS.md#collector-settings)
for the separate Julia heap controls.
