---
id: parallel.types__persistentforeachpool
label: _PersistentForeachPool
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: _PersistentForeachPool
  lines:
  - 112
  - 112
inputs:
- id: workers
  type: Int
  units: n/a
  required: true
  description: Field `workers`.
- id: request_channels
  type: Vector{Channel{Any}}
  units: n/a
  required: true
  description: Field `request_channels`.
- id: done_channel
  type: Channel{Any}
  units: n/a
  required: true
  description: Field `done_channel`.
- id: run_lock
  type: ReentrantLock
  units: n/a
  required: false
  description: Field `run_lock` (default `ReentrantLock()`).
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
  type: _PersistentForeachPool
  units: n/a
  description: Constructed `_PersistentForeachPool` (keyword constructor via @kwdef).
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

# _PersistentForeachPool

## Purpose
Long-lived worker pool for the persistent `foreach` dispatch path. Instead of spawning tasks per call, `workers` Julia tasks block on per-worker request channels, execute one shared work request, and report completion on a common done channel, amortising task-creation cost across many small parallel loops.

## Design & Implementation
`Base.@kwdef mutable struct` with `workers::Int`, `request_channels::Vector{Channel{Any}}` (one `Channel{Any}(1)` per worker), `done_channel::Channel{Any}` of capacity `workers`, and a `run_lock::ReentrantLock` that serialises dispatches. Pools are cached in `_persistent_foreach_pools` and `_persistent_foreach_worker_pools`, both `Dict{Tuple{UInt, Symbol}, _PersistentForeachPool}` guarded by `_persistent_foreach_lock`; the second variant passes `worker_id` to `f(worker_id, idx)` so callers can use per-worker scratch buffers. Workers receive a named tuple request (`num_items`, `active_workers`, `scheduler`, `chunk`, `f`, `next_index`) and run either a `:dynamic` chunked loop driven by `Threads.atomic_add!(next_index, chunk)` or a static stride `worker_id:active_workers:num_items`. Exceptions are caught as `Base.CapturedException` and put on `done_channel`; `_shutdown_persistent_foreach_pool!` sends `:stop` to every request channel.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workers` | Int | n/a | yes | Field `workers`. |
| in | `request_channels` | Vector{Channel{Any}} | n/a | yes | Field `request_channels`. |
| in | `done_channel` | Channel{Any} | n/a | yes | Field `done_channel`. |
| in | `run_lock` | ReentrantLock | n/a | no | Field `run_lock` (default `ReentrantLock()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentForeachPool | n/a | — | Constructed `_PersistentForeachPool` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.context__create_persistent_foreach_pool|_create_persistent_foreach_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:124-124`
- [[parallel.context__create_persistent_foreach_worker_pool|_create_persistent_foreach_worker_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:105-105`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Channel{Any}` boxes every request and result, so dispatch latency is in the microsecond range (the file comment cites 1 to 5 microseconds for a condvar wake), which is why the spin-barrier pool exists. Pools are never shrunk; a pool created for `workers` threads persists until shutdown. If a worker task dies outside the `try` block (for example an `InterruptException` while blocked on `take!`), the coordinator waits forever on `done_channel`.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 112.
