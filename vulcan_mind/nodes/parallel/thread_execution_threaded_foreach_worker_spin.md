---
id: parallel.thread_execution_threaded_foreach_worker_spin
label: threaded_foreach_worker_spin
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_foreach_worker_spin
  lines:
  - 158
  - 158
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
  type: Nothing
  units: n/a
  description: 'Return value of `threaded_foreach_worker_spin`. Returns `nothing`
    or `_spin_barrier_dispatch!(pool, num_items, workers, f)`. Type parameters: `{F
    <: Function}`.'
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

# threaded_foreach_worker_spin

## Purpose
Lowest-latency variant of the worker-aware loop, dispatching `f(worker_id, idx)` to a spin-barrier pool whose workers busy-poll an atomic generation counter, reducing dispatch overhead from microseconds to tens of nanoseconds for the harmonics SIMD batch at high thread counts.

## Design & Implementation
Returns for `num_items <= 0`; computes `workers = _thread_worker_count(num_items, allotment)` and runs `f(1, idx)` serially when that is 1. Otherwise it obtains `_spin_barrier_pool_for(source)` and calls `_spin_barrier_dispatch!(pool, num_items, workers, f)`. Unlike the channel-based persistent variants there is no `callback_persistent_workers_enabled()` check; opt-in is instead governed by callers consulting `harmonics_batch_spin_barrier_enabled()` (`SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER`, default false). A `do`-block ordering method exists.

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
| out | `result` | Nothing | n/a | — | Return value of `threaded_foreach_worker_spin`. Returns `nothing` or `_spin_barrier_dispatch!(pool, num_items, workers, f)`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:168-168`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:168-168`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:168-168`
- `callees` → [[parallel.context__spin_barrier_dispatch_bang|_spin_barrier_dispatch!]] · `callers` · call · `src/parallel/policy/thread_execution.jl:173-173`
- `callees` → [[parallel.context__spin_barrier_pool_for|_spin_barrier_pool_for]] · `callers` · call · `src/parallel/policy/thread_execution.jl:172-172`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:165-165`
<!-- vulcan:connections:end -->

## Limitations
Spinning workers consume full CPU cores while idle between dispatches, so this mode is inappropriate on shared machines or when the pool is idle for long periods. There is no fallback to the channel pool if the spin pool cannot be created. Worker count is bounded by the pool created at first use, and the caller must not request more workers than exist. Exception propagation semantics depend entirely on `_spin_barrier_dispatch!`.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 158.
