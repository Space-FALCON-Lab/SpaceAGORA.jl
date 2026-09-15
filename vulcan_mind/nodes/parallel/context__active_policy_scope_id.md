---
id: parallel.context__active_policy_scope_id
label: _active_policy_scope_id
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _active_policy_scope_id
  lines:
  - 17
  - 17
inputs:
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
  type: UInt
  units: n/a
  description: Return value of `_active_policy_scope_id`.
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

# _active_policy_scope_id

## Purpose
Returns the `UInt` identifier of the currently active policy scope, the first element of every persistent-pool dictionary key, so that pools created inside one `with_policy_context` block are invisible to other scopes.

## Design & Implementation
An `@inline` composition: `_policy_scope_id(_active_policy_context())`. `_active_policy_context` returns the task-local or global `PolicyContext`, and `_policy_scope_id` maps it to `UInt(objectid(ctx))`. No allocation is performed and the return type is `UInt`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | UInt | n/a | — | Return value of `_active_policy_scope_id`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.thread_execution__persistent_pool_key|_persistent_pool_key]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:45-45`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/context.jl:18-18`
- `callees` → [[parallel.context__policy_scope_id|_policy_scope_id]] · `callers` · call · `src/parallel/policy/context.jl:18-18`
<!-- vulcan:connections:end -->

## Limitations
The identifier inherits the caveats of `objectid`: it is unique only while the context object is alive, and a garbage-collected context's id may be reused by a later allocation. Since `_destroy_persistent_foreach_scope!` removes all keys for a scope at exit, reuse does not normally leak pools, but a scope that exits without running the `finally` block (process abort) leaves stale entries. Cost is dominated by the `try`/`catch` in `_active_policy_context`.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 17.
