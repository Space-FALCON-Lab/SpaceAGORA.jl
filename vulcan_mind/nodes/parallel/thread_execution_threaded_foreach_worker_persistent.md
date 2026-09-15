---
id: parallel.thread_execution_threaded_foreach_worker_persistent
label: threaded_foreach_worker_persistent
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_foreach_worker_persistent
  lines:
  - 127
  - 127
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
- id: f
  type: F
  units: n/a
  required: true
  description: Positional argument `f`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Any
  units: n/a
  description: 'Return value of `threaded_foreach_worker_persistent`. Returns `threaded_foreach_worker(num_items,
    allotment, f)` or `_threaded_foreach_persistent!(pool, num_items, workers, f)`.
    Type parameters: `{F <: Function}`.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# threaded_foreach_worker_persistent

## Purpose
Worker-aware persistent-pool loop: invokes `f(worker_id, idx)` over all items on long-lived tasks associated with `source`, used for batches that need per-worker resources, such as the GRAM isolated-model pool.

## Design & Implementation
Returns for `num_items <= 0`. Computes `workers = _thread_worker_count(num_items, allotment)` and falls back to `threaded_foreach_worker(num_items, allotment, f)` when `workers <= 1` or `callback_persistent_workers_enabled()` is false. Otherwise obtains the pool through `_persistent_worker_pool_for(source)` and runs `_threaded_foreach_persistent!(pool, num_items, workers, f)`. The worker loop inside the pool supplies its own 1-based `worker_id` to `f`. A `do`-block ordering method is provided.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `threaded_foreach_worker_persistent`. Returns `threaded_foreach_worker(num_items, allotment, f)` or `_threaded_foreach_persistent!(pool, num_items, workers, f)`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:860-860`
- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:261-261`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1100-1100`
- [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1272-1272`
- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:201-201`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/thread_execution.jl:156-156`
- `callees` → [[parallel.thread_execution__persistent_worker_pool_for|_persistent_worker_pool_for]] · `callers` · call · `src/parallel/policy/thread_execution.jl:138-138`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:134-134`
- `callees` → [[parallel.thread_execution__threaded_foreach_persistent_bang|_threaded_foreach_persistent!]] · `callers` · call · `src/parallel/policy/thread_execution.jl:139-139`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callers` · call · `src/parallel/policy/thread_execution.jl:136-136`
<!-- vulcan:connections:end -->

## Limitations
The `worker_id` passed ranges over `1:workers` of the pool, which may exceed the number of resources a caller prepared if it sized them using a different allotment; `_ensure_gram_isolated_pool!` avoids this by sizing from the same `thread_worker_count`. Batches on one pool are serialised by `run_lock`. Disabled persistent workers silently change the scheduling path, which can alter timing but not results.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 127.
