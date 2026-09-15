---
id: parallel.thread_execution_threaded_foreach_persistent
label: threaded_foreach_persistent
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_foreach_persistent
  lines:
  - 94
  - 94
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
  description: 'Return value of `threaded_foreach_persistent`. Returns `threaded_foreach(num_items,
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

# threaded_foreach_persistent

## Purpose
Runs `f(idx)` over `1:num_items` using a persistent index-only worker pool identified by `source`, falling back to the spawn-per-call `threaded_foreach` when parallelism is not warranted or persistent workers are disabled.

## Design & Implementation
Returns immediately for `num_items <= 0`. Computes `workers = _thread_worker_count(num_items, allotment)`; if `workers <= 1` or `!callback_persistent_workers_enabled()` it delegates to `threaded_foreach(num_items, allotment, f)`. Otherwise it fetches the pool with `_persistent_pool_for(source)` and executes `_threaded_foreach_persistent!(pool, num_items, workers, f)`. A second method places `f` first to support `do`-block syntax. The function returns `nothing`.

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
| out | `result` | Any | n/a | — | Return value of `threaded_foreach_persistent`. Returns `threaded_foreach(num_items, allotment, f)` or `_threaded_foreach_persistent!(pool, num_items, workers, f)`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:291-291`
- [[simulation.thermal_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:75-75`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:110-110`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:291-291`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:75-75`

**Downstream**

- `callees` → [[parallel.thread_execution__persistent_pool_for|_persistent_pool_for]] · `callers` · call · `src/parallel/policy/thread_execution.jl:105-105`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:101-101`
- `callees` → [[parallel.thread_execution__threaded_foreach_persistent_bang|_threaded_foreach_persistent!]] · `callers` · call · `src/parallel/policy/thread_execution.jl:106-106`
- `callees` → [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callers` · call · `src/parallel/policy/thread_execution.jl:103-103`
<!-- vulcan:connections:end -->

## Limitations
The fallback path re-evaluates `_thread_worker_count` inside `threaded_foreach`, a minor duplicated cost. Persistent worker enablement is read from the environment on every call. Exceptions from workers surface only after every worker completes the batch. The pool lookup is keyed on the current policy scope, so calls from a task outside the intended scope get a different pool.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 94.
