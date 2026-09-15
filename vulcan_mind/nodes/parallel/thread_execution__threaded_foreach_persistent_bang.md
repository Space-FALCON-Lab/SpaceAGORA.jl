---
id: parallel.thread_execution__threaded_foreach_persistent_bang
label: _threaded_foreach_persistent!
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: _threaded_foreach_persistent!
  lines:
  - 57
  - 57
inputs:
- id: pool
  type: _PersistentForeachPool
  units: n/a
  required: true
  description: Positional argument `pool`.
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  description: 'Return value of `_threaded_foreach_persistent!`; mutates `pool` in
    place. Returns `nothing`. Type parameters: `{F <: Function}`.'
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

# _threaded_foreach_persistent!

## Purpose
Dispatches one batch of `num_items` work items across the first `workers` long-lived tasks of a persistent pool, blocks until all have finished, and rethrows the first worker exception, providing the execution core for both persistent foreach variants.

## Design & Implementation
Reads `inner_scheduler_mode()` and `inner_dynamic_chunk_size()`, allocates a shared `Threads.Atomic{Int}(1)` cursor, then under `lock(pool.run_lock)` puts a `NamedTuple` request `(num_items, active_workers, scheduler, chunk, next_index, f)` into `pool.request_channels[worker_id]` for each of the `workers` workers. It then `take!`s exactly `workers` completions from `pool.done_channel`; each completion is either `nothing` or a captured exception, and the first non-`nothing` value is remembered and thrown after all workers report. Returns `nothing`. The `run_lock` guarantees one batch in flight per pool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool` | _PersistentForeachPool | n/a | yes | Positional argument `pool`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_threaded_foreach_persistent!`; mutates `pool` in place. Returns `nothing`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:106-106`
- [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:139-139`

**Downstream**

- `callees` → [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:64-64`
- `callees` → [[parallel.env_config_inner_scheduler_mode|inner_scheduler_mode]] · `callers` · call · `src/parallel/policy/thread_execution.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
Only the first error is surfaced; subsequent worker errors are discarded. If a worker task has died, `take!` blocks forever because fewer than `workers` completions will ever arrive. Callers must ensure `workers <= length(pool.request_channels)`, which is not checked here. The static scheduler ignores `chunk`; dynamic mode may leave later workers idle when `num_items` is small relative to `chunk`. Each call allocates the request tuples and the atomic.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 57.
