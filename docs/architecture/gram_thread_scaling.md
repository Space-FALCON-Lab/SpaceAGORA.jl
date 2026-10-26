# Native GRAM under threads: the global lock and the isolated pool

The original investigation, including its measurements and limitations, is
preserved at this same relative path in the private archive described by the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).
This page keeps current methods and configuration; historical results do not
qualify performance or scientific equivalence at a later revision.

Native GRAM calls normally serialize on `RuntimeServices.GRAM_LOCK`, the same
object as `SPICE_LOCK`. The shared lock protects the native GRAM/CSPICE state.
The isolated pool instead uses independent model copies, each with its own
lock, to evaluate an eligible density batch across workers.

## Configuration and eligibility

The settings live in `src/orchestration/simulation/callbacks/density_callbacks/config.jl`.

| Setting | Default | Meaning |
|---|---|---|
| `SPACEAGORA_GRAM_ISOLATED_POOL` | `auto` | Use the pool only when its batch, width and wind-history guards allow it. `off` retains the locked route; explicit `on` opts into independent instance histories. |
| `SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD` | `1024` | Minimum number of items that actually reach native GRAM in an automatic batch. |
| `SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS` | `min(4, Threads.nthreads())` | Maximum number of worker instances; at least two are required to use the pool. |

These defaults came from the archived nominal-wind study, not a general speedup
guarantee. The threshold counts native calls rather than configured spacecraft:
items that return vacuum or use the above-interface density approximation do
not justify constructing a native pool. The pooled native call applies the same
`-30.0` m altitude floor as the locked scalar path.

The density callback requests its pool width with `lock_free=true`. This avoids
the shared-lock path's default 16-thread budget floor
(`SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET`); the pool remains subject
to the available thread budget, worker cap and its other eligibility checks.
The earlier claim that the isolated pool is unreachable below 16 threads
applies to the historical implementation, not this path.

## Wind history and identity

Automatic pooling declines a wind-requesting batch when the native core reports
history-dependent winds, before constructing or warming workers. Unknown cores
are treated conservatively. The pinned wrapper's `auto` wind mode resolves to
`perturbed`, so nominal winds must be selected explicitly for automatic pooling
of wind queries. Calls with `wind=false` remain eligible. The guard also applies
when winds are requested for diagnostics but not used by the force model.

Explicit `SPACEAGORA_GRAM_ISOLATED_POOL=on` bypasses this wind-history guard.
Separate instances then follow separate stochastic histories; winds and
trajectories can depend on pool width or thread count. It is not a
reproducibility-preserving speed switch for perturbed winds. Routing does not
change an atmosphere mode or seed.

The archived bit-identity comparisons cover specific nominal-wind cases only.
They do not qualify perturbed winds, other model versions or the current code
without a fresh comparison. `density_grid_parity.jl` compares exact floating-point
bits on a grid; `dump_states.jl` writes saved times and full spacecraft states
for byte comparison. Both are in `benchmarks/studies/gram_thread_scaling/`.

## Which paths the pool can reach

| Path | Relationship to the pool |
|---|---|
| Density callback batch (`density_callbacks/runtime.jl`) | Evaluates eligible batches through the isolated pool, with a locked fallback. |
| Vacuum-predicted look-ahead cache (`SPACEAGORA_VACUUM_GRAM_CACHE`) | Its own knot construction uses scalar calls on the shared model. The callback's staged density batch is a separate operation. |
| GRAM track cache (`SPACEAGORA_GRAM_TRACK_CACHE`) | Its fallback can use the pool when the driver lacks native `generate_trajectory`; the native trajectory path is a locked call. |
| Constellation RHS aero sampling | Uses the model selected by `_density_model_for_sat`, not this callback pool. |

The other instance-isolation controls serve different callers:

- `SPACEAGORA_GRAM_PER_SAT_INSTANCES` supplies per-spacecraft models for RHS
  sampling. Initialization includes serial warm-up of fresh native instances.
- `SPACEAGORA_GRAM_LOCK_SCOPE=model` selects the wrapper's instance lock for
  scalar calls. One shared model still serializes on one instance lock. The
  isolated pool supplies its own locks and does not use this selector.

Model-scoped scalar locks are not recorded in the shared native-lock counters.
A small shared-lock total therefore does not, by itself, demonstrate a speedup.

## Benchmark configuration and interpretation

Before timing a GRAM case, check how many members actually call native GRAM,
whether they are inside the entry interface, and whether the intended cache is
built. A large constellation can otherwise mostly measure vacuum or a density
approximation while retaining a GRAM case name.

Current P6 density traces use `ppc_p6_gram_constellation` in
`benchmarks/studies/parallelization_performance/cases.jl`: a 300–480 km band
inside the 600 km entry interface. Traces 4, 5 and 6 share that constellation;
trace 6's individual samples correspond to members of trace 5. The
`p6_case_audit.jl` tool checks the case definitions. Historical rows predating
the trace 5/6 redefinition at `917ab6aa7`, or the trace 4 redefinition at
`be774d5b6`, must not be pooled with measurements of the current definitions.
Their case names alone do not identify equivalent work.

For a new comparison, record the revision, wrapper and data versions, wind
mode, seed, density path, actual native-item count and pool width, including `SPACEAGORA_DENSITY_FREEZE_PER_STEP`. Alternate the
locked and pooled arms, separate compilation and setup from the timed region,
and compare repeated timings alongside exact state dumps. Run a one-thread
control, where the pool must decline. Distinguish serial native work from lock
wait: reducing lock acquisitions does not prove that contention was the limiting
cost. Retain the raw rows and machine details with the resulting study.

## Refresh de-phasing: not available

Cache horizons determine the locations of interpolation knots. Jittering those
horizons changes the positions and times queried from GRAM, hence the spline and
the densities read by the RHS. It is a scientific change, not a scheduling-only
optimization. Staggered refreshes are not implemented or qualified here.
Batching or prefetching would need separate validation that the sampled points
and state ownership remain unchanged.

## Pool-worker memory growth

A Julia heap-size hint does not bound memory owned by a native atmosphere.
Changing the epoch or copying a GRAM model can construct a native image; releasing
the small Julia wrapper does not release that image until its finalizer runs.
A sequence of discarded models can therefore grow native memory without a
comparable increase in the Julia heap.

The benchmark harness caches its model at `ppc_initial_time()`, shared with its
configuration builder, to avoid rebuilding it solely for epoch alignment.
The GRAM extension also counts native atmospheres created and finalized.
`collect_unreferenced_gram_atmospheres!()` runs after successful or failed solves
and requests a full collection when the count of additional live atmospheres
reaches `SPACEAGORA_GRAM_NATIVE_COLLECT_LIMIT` (default `8`; `0` or `off`
disables it). This is a collection trigger, not a hard resident-memory limit.
Live models remain owned by their callers, and allocator behavior can keep
resident memory high after objects are freed.

Use `worker_growth.jl` to distinguish Julia allocations, native growth, model
construction and finalization. Its state dumps also support forward/reverse
sample-order comparisons. `worker_rss.jl` includes explicit collections, which
must be accounted for when interpreting its measurements. The tests in
`test/unit/environment/gram_native_collect_tests.jl` and
`gram_isolated_pool_tests.jl` cover the corresponding contracts; native memory
and trajectory claims require the actual wrapper, inputs and scenarios under
test. See [Collector settings](../../benchmarks/studies/heap_contention/FINDINGS.md#collector-settings) for the
separate Julia heap controls.
