---
id: simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang
label: _gram_isolated_pool_batch_eval!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _gram_isolated_pool_batch_eval!
  lines:
  - 155
  - 216
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: batch_request
  type: Tuple{AbstractVector,AbstractVector,AbstractVector,GRAMAtmosphereModel,AbstractVector,AbstractVector,AbstractVector}
  units: m, rad, rad
  required: true
  description: Output density, temperature and wind buffers together with the template
    GRAM model and the per-spacecraft altitude, latitude and longitude arrays to evaluate.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: pooled
  type: Bool
  units: n/a
  description: True when the isolated worker pool serviced the batch and the output
    buffers were filled; false when the caller must fall back to `getDensityBatch!`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# _gram_isolated_pool_batch_eval!

## Purpose
`_gram_isolated_pool_batch_eval!` evaluates an entire batch of GRAM atmosphere queries in parallel by giving each worker thread its own private copy of the GRAM model. GRAM's native state is not reentrant, so the ordinary path serialises behind `GRAM_LOCK`; the isolated pool trades memory for concurrency by removing the shared state entirely.

## Model & Assumptions
Two methods are defined. The generic `@inline` method accepts any density model and immediately returns `false`, so non-GRAM atmospheres pay nothing and the caller falls through to the standard batch path. The specialised method for `GRAMAtmosphereModel` first checks that the pool is enabled for this batch size, then validates that every output and input array has the same length as the altitude vector, refusing the batch on any mismatch instead of writing partial results. Worker count comes from `ParallelPolicy.thread_worker_count` capped by the configured maximum workers, and a count of one declines the batch because a single worker offers no advantage over the locked path.

## Design & Implementation
`_ensure_gram_isolated_pool!` lazily materialises the pool on `p.shared_buffers`, deep-copying the template model once per worker alongside a matching `ReentrantLock` vector, and rebuilds only when the requested worker count differs from the cached pool size. Work is then distributed with `threaded_foreach_worker_persistent(:rhs_gram_batch, n, workers)`, whose closure receives both a worker id and an element index: the worker id selects the private model and lock, and results are written with `@inbounds` into the caller's buffers at that index. Elapsed time is not recorded here; the calling density callback owns policy observation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `batch_request` | Tuple{AbstractVector,AbstractVector,AbstractVector,GRAMAtmosphereModel,AbstractVector,AbstractVector,AbstractVector} | m, rad, rad | yes | Output density, temperature and wind buffers together with the template GRAM model and the per-spacecraft altitude, latitude and longitude arrays to evaluate. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `pooled` | Bool | n/a | — | True when the isolated worker pool serviced the batch and the output buffers were filled; false when the caller must fall back to `getDensityBatch!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:307-307`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:243-243`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:307-307`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:204-204`
- `callees` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:197-197`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:201-201`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:185-185`
- `callees` → [[simulation.config__gram_isolated_pool_enabled|_gram_isolated_pool_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:186-186`
- `callees` → [[simulation.model_selection__ensure_gram_isolated_pool_bang|_ensure_gram_isolated_pool!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:200-200`
- `callees` → [[simulation.model_selection__gram_batch_elapsed_time|_gram_batch_elapsed_time]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:207-207`
- `callees` → [[simulation.model_selection__gram_isolated_pool_density_state|_gram_isolated_pool_density_state]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:208-208`
<!-- vulcan:connections:end -->

## Limitations
Deep-copying GRAM models multiplies the memory and native-handle footprint by the worker count, which is why the maximum worker count is a separate configurable ceiling. The per-worker lock still serialises any residual shared native state a model copy retains. The boolean return is the only failure signal, so callers must honour it and run the serial batch when it is false, as the density callback does.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl:170-216`.
