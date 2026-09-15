---
id: parallel.context__create_persistent_foreach_pool
label: _create_persistent_foreach_pool
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _create_persistent_foreach_pool
  lines:
  - 120
  - 120
inputs:
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  description: Return value of `_create_persistent_foreach_pool`.
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

# _create_persistent_foreach_pool

## Purpose
Constructs a `_PersistentForeachPool` whose long-lived tasks run `_persistent_foreach_worker_loop`, the index-only variant that calls `f(idx)`, and starts those tasks immediately so later dispatches pay only channel-handoff cost.

## Design & Implementation
Clamps `workers = max(2, workers)` so a pool always has at least two request lanes. Allocates one `Channel{Any}(1)` request channel per worker and a single `Channel{Any}(workers)` done channel sized so every worker can post a completion without blocking. Builds the pool struct via keyword constructor `(workers, request_channels, done_channel)` and then `Threads.@spawn`s `_persistent_foreach_worker_loop(worker_id, request_channels[worker_id], done_channel)` for each `worker_id`. Returns the pool; the spawned tasks are not stored, so they are kept alive only by the channels they block on.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentForeachPool | n/a | — | Return value of `_create_persistent_foreach_pool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.thread_execution__persistent_pool_for|_persistent_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:52-52`

**Downstream**

- `callees` → [[parallel.context__persistent_foreach_worker_loop|_persistent_foreach_worker_loop]] · `callers` · call · `src/parallel/policy/context.jl:130-130`
- `callees` → [[parallel.types__persistentforeachpool|_PersistentForeachPool]] · `callers` · call · `src/parallel/policy/context.jl:124-124`
<!-- vulcan:connections:end -->

## Limitations
Spawned tasks run until they receive `:stop`; if the pool is dropped without `_shutdown_persistent_foreach_pool!` they leak forever. Requests and completions are typed `Any`, so every `take!` involves dynamic typing and boxing. The `max(2, ...)` floor means even a 1-thread budget creates two worker tasks that will share one OS thread. The tasks are spawned on the default threadpool with no affinity.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 120.
