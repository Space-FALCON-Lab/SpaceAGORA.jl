---
id: parallel.thread_execution_threaded_reduce
label: threaded_reduce
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_reduce
  lines:
  - 277
  - 277
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
- id: init
  type: I
  units: n/a
  required: true
  description: Positional argument `init`.
- id: body_bang
  type: B
  units: n/a
  required: true
  description: Positional argument `body!`.
- id: combine_bang
  type: C
  units: n/a
  required: true
  description: Positional argument `combine!`.
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
  description: 'Return value of `threaded_reduce`. Returns `acc0` or `result`. Type
    parameters: `{I <: Function, B <: Function, C <: Function}`.'
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

# threaded_reduce

## Purpose
Parallel map-reduce over `1:num_items` with mutable per-worker accumulators: each worker folds its items into a private accumulator via `body!`, and the results are merged in worker order with `combine!`, avoiding contention on a shared total.

## Design & Implementation
`workers = _thread_worker_count(num_items, allotment)`; `acc0 = init()` is created up front and returned directly when `num_items <= 0`, or filled serially when `workers <= 1`. In parallel mode a `Vector{typeof(acc0)}(undef, workers)` holds partials, with worker 1 reusing `acc0`. Under `Threads.@sync`, each spawned task builds `local_acc` (via `init()` for workers 2..n), processes either dynamic chunks handed out by an `Atomic{Int}` cursor or the static stride `worker_id:workers:num_items`, and stores `partials[worker_id]`. After the sync, `combine!(result, partials[k])` is applied for `k = 2:workers` and `result` is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `init` | I | n/a | yes | Positional argument `init`. |
| in | `body_bang` | B | n/a | yes | Positional argument `body!`. |
| in | `combine_bang` | C | n/a | yes | Positional argument `combine!`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `threaded_reduce`. Returns `acc0` or `result`. Type parameters: `{I <: Function, B <: Function, C <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`

**Downstream**

- `callees` → [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:300-300`
- `callees` → [[parallel.env_config_inner_scheduler_mode|inner_scheduler_mode]] · `callers` · call · `src/parallel/policy/thread_execution.jl:298-298`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:284-284`
<!-- vulcan:connections:end -->

## Limitations
The final result depends on which worker processed which item unless `body!` and `combine!` are exactly associative and commutative; for floating-point sums this makes the outcome non-deterministic across runs and worker counts, which is why `threaded_collect!` is recommended for integrator-facing reductions. `init()` is called once per extra worker per call, allocating accumulators every batch. `partials` slots are assigned from tasks; the type must be stable or the vector write throws.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 277.
