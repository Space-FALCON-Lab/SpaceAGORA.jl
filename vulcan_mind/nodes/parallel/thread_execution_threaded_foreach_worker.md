---
id: parallel.thread_execution_threaded_foreach_worker
label: threaded_foreach_worker
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_foreach_worker
  lines:
  - 199
  - 199
inputs:
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
  type: Nothing
  units: n/a
  description: 'Return value of `threaded_foreach_worker`. Returns `nothing`. Type
    parameters: `{F <: Function}`.'
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

# threaded_foreach_worker

## Purpose
Spawn-per-call parallel loop that invokes `f(worker_id, idx)` for every item, giving each callback a stable worker index so it can use per-worker scratch buffers or model instances without locking.

## Design & Implementation
Returns for `num_items <= 0`; computes `workers = _thread_worker_count(num_items, allotment)` and, when that is 1, runs `f(1, idx)` serially. Otherwise it reads `inner_scheduler_mode()`. In `:dynamic` mode a shared `Threads.Atomic{Int}` cursor hands out chunks of `inner_dynamic_chunk_size()` indices via `atomic_add!` to `workers` spawned tasks under `Threads.@sync`. In static mode each task handles the strided range `worker_id:workers:num_items`. Both paths use `@inbounds`. A `do`-block ordering method exists.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `threaded_foreach_worker`. Returns `nothing`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:253-253`
- [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:136-136`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:204-204`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:204-204`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:204-204`
- `callees` → [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:210-210`
- `callees` → [[parallel.env_config_inner_scheduler_mode|inner_scheduler_mode]] · `callers` · call · `src/parallel/policy/thread_execution.jl:208-208`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:201-201`
- `callees` → [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callers` · feedback · `src/parallel/policy/thread_execution.jl:241-241`
- `callees` → [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `callers` · feedback · `src/parallel/policy/thread_execution.jl:242-242`
<!-- vulcan:connections:end -->

## Limitations
Tasks are spawned fresh each call, costing several microseconds per worker; hot loops should use the persistent or spin-barrier variants. Static striding gives poor balance when per-item cost is heterogeneous. `worker_id` is a logical index, not `Threads.threadid()`, so it is stable across task migration but cannot be used to index thread-local storage. Exceptions inside any task propagate through `@sync` as a `CompositeException`.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 199.
