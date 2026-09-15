---
id: parallel.context__create_persistent_foreach_worker_pool
label: _create_persistent_foreach_worker_pool
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _create_persistent_foreach_worker_pool
  lines:
  - 101
  - 101
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
  description: Return value of `_create_persistent_foreach_worker_pool`.
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

# _create_persistent_foreach_worker_pool

## Purpose
Constructs a `_PersistentForeachPool` whose tasks run `_persistent_foreach_worker_loop_w`, the worker-aware variant that calls `f(worker_id, idx)`, enabling per-worker resources such as isolated GRAM model instances.

## Design & Implementation
Structurally identical to `_create_persistent_foreach_pool`: `workers = max(2, workers)`, a `Channel{Any}(1)` per worker for requests, a `Channel{Any}(workers)` done channel, keyword construction of `_PersistentForeachPool`, and one `Threads.@spawn` per worker invoking `_persistent_foreach_worker_loop_w(worker_id, request_channels[worker_id], done_channel)`. The only difference is which loop function the tasks execute, which determines whether `worker_id` is forwarded to the user closure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentForeachPool | n/a | — | Return value of `_create_persistent_foreach_worker_pool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.thread_execution__persistent_worker_pool_for|_persistent_worker_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:122-122`

**Downstream**

- `callees` → [[parallel.context__persistent_foreach_worker_loop_w|_persistent_foreach_worker_loop_w]] · `callers` · call · `src/parallel/policy/context.jl:111-111`
- `callees` → [[parallel.types__persistentforeachpool|_PersistentForeachPool]] · `callers` · call · `src/parallel/policy/context.jl:105-105`
<!-- vulcan:connections:end -->

## Limitations
The same pool struct type is used for both variants, so nothing at the type level prevents dispatching an index-only closure to a worker-aware pool or vice versa; a mismatch produces a `MethodError` inside the worker that is captured and rethrown after the batch. Tasks leak if the pool is never shut down. Because `workers` is fixed at creation, a later larger thread budget cannot be exploited without destroying the scope.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 101.
