---
id: parallel.thread_execution__persistent_pool_for
label: _persistent_pool_for
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: _persistent_pool_for
  lines:
  - 48
  - 48
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
  description: Return value of `_persistent_pool_for`.
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

# _persistent_pool_for

## Purpose
Returns the lazily created `_PersistentForeachPool` (index-only worker variant) associated with a named `source` within the active policy scope, so repeated `threaded_foreach_persistent` calls from the same call site reuse long-lived worker tasks instead of spawning fresh ones.

## Design & Implementation
Builds the lookup key with `_persistent_pool_key(source)`, which pairs `_active_policy_scope_id()` with the `source` symbol. Under `lock(_persistent_foreach_lock)` it calls `get!(_persistent_foreach_pools, key) do ... end`, constructing a pool of `_default_thread_pool_size()` workers via `_create_persistent_foreach_pool` only on first access. The `return` inside the `do` block returns the pool from the closure and thus from the function. Return type is annotated `_PersistentForeachPool`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _PersistentForeachPool | n/a | — | Return value of `_persistent_pool_for`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:105-105`

**Downstream**

- `callees` → [[parallel.context__create_persistent_foreach_pool|_create_persistent_foreach_pool]] · `callers` · call · `src/parallel/policy/thread_execution.jl:52-52`
- `callees` → [[parallel.env_config__default_thread_pool_size|_default_thread_pool_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:52-52`
- `callees` → [[parallel.thread_execution__persistent_pool_key|_persistent_pool_key]] · `callers` · call · `src/parallel/policy/thread_execution.jl:49-49`
<!-- vulcan:connections:end -->

## Limitations
Pools are never evicted by this function; every distinct `(scope, source)` pair keeps its worker tasks alive until `_destroy_persistent_foreach_scope!` runs, so unbounded `source` symbols leak tasks. The pool size is fixed at creation from `_default_thread_pool_size()`, so a later change to the thread budget does not resize an existing pool. The global lock serialises all pool lookups, which matters only when many sources are first-touched concurrently.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 48.
