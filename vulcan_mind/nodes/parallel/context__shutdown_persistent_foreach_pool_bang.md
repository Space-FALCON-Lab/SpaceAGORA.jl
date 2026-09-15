---
id: parallel.context__shutdown_persistent_foreach_pool_bang
label: _shutdown_persistent_foreach_pool!
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _shutdown_persistent_foreach_pool!
  lines:
  - 139
  - 139
inputs:
- id: pool
  type: _PersistentForeachPool
  units: n/a
  required: true
  description: Positional argument `pool`.
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
  description: Return value of `_shutdown_persistent_foreach_pool!`; mutates `pool`
    in place.
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

# _shutdown_persistent_foreach_pool!

## Purpose
Signals every worker task of a channel-based persistent pool to exit by posting the `:stop` sentinel on each request channel, allowing the tasks to terminate and be garbage-collected.

## Design & Implementation
Acquires `pool.run_lock` so the stop signals cannot interleave with an in-flight batch's requests, then iterates `pool.request_channels` with `@inbounds` and `put!(channel, :stop)` on each. Because request channels have capacity 1 and the workers block on `take!`, each `put!` completes as soon as the corresponding worker consumes any prior request. Returns `nothing` after releasing the lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool` | _PersistentForeachPool | n/a | yes | Positional argument `pool`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_shutdown_persistent_foreach_pool!`; mutates `pool` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__destroy_persistent_foreach_scope_bang|_destroy_persistent_foreach_scope!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:324-324`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
If a worker task has already died (for example from an uncaught error outside the `try`), its channel's single slot fills and `put!` blocks forever, wedging the shutdown and the scope exit. The function does not wait for workers to actually exit, so callers cannot assume the threads are free on return. It is not idempotent: a second call would block on the full channels of already-stopped workers.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 139.
