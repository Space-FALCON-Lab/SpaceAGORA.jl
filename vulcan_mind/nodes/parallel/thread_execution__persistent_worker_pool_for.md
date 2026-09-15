---
id: parallel.thread_execution__persistent_worker_pool_for
label: _persistent_worker_pool_for
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: _persistent_worker_pool_for
  lines:
  - 118
  - 118
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
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
  description: Return value of `_persistent_worker_pool_for`.
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

# _persistent_worker_pool_for

## Purpose
Returns the lazily created `_PersistentForeachPool` whose workers invoke `f(worker_id, idx)` (worker-aware variant) for a given `source` in the active policy scope, backing `threaded_foreach_worker_persistent`.

## Design & Implementation
Identical in structure to `_persistent_pool_for` but keyed into the separate `_persistent_foreach_worker_pools` dictionary and constructed with `_create_persistent_foreach_worker_pool(_default_thread_pool_size())`. The lookup is guarded by `_persistent_foreach_lock`, and `get!` ensures only one pool is created per `(scope_id, source)` key. The worker-aware pool differs in that its worker loop passes its own 1-based `worker_id` to `f`, allowing per-worker scratch buffers or model instances (as in the GRAM isolated pool).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentForeachPool | n/a | — | Return value of `_persistent_worker_pool_for`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:138-138`

**Downstream**

- `callees` → [[parallel.context__create_persistent_foreach_worker_pool|_create_persistent_foreach_worker_pool]] · `callers` · call · `src/parallel/policy/thread_execution.jl:122-122`
- `callees` → [[parallel.env_config__default_thread_pool_size|_default_thread_pool_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:122-122`
- `callees` → [[parallel.thread_execution__persistent_pool_key|_persistent_pool_key]] · `callers` · call · `src/parallel/policy/thread_execution.jl:119-119`
<!-- vulcan:connections:end -->

## Limitations
Maintaining two parallel dictionaries (index-only and worker-aware) doubles resident worker tasks when the same `source` is used with both APIs. Pool size is frozen at first creation. No mechanism verifies that the `f` closure type is stable across calls, so each distinct closure type incurs a new dynamic dispatch inside the worker loop. Pools persist until the owning scope is destroyed.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 118.
