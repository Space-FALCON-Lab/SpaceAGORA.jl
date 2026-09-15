---
id: parallel.context__active_policy_context
label: _active_policy_context
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _active_policy_context
  lines:
  - 1
  - 1
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
  type: PolicyContext
  units: n/a
  description: Return value of `_active_policy_context`.
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

# _active_policy_context

## Purpose
Resolves the `PolicyContext` in effect for the calling task: the one installed by `with_policy_context` in task-local storage if present, otherwise the process-global default. Every persistent pool lookup derives its scope id from this object.

## Design & Implementation
Attempts `Base.task_local_storage(_policy_context_tls_key)` inside a `try` block, since that call throws `KeyError` when the key is absent; any exception is mapped to `nothing`. If the fetched value `isa PolicyContext` it is returned, otherwise `_global_policy_context[]` (a `Ref`) is dereferenced and returned. The function is `@inline` with a `::PolicyContext` return annotation, so the union is resolved before returning.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PolicyContext | n/a | — | Return value of `_active_policy_context`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.context__active_policy_scope_id|_active_policy_scope_id]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:18-18`
- [[parallel.observation_tracking_record_route_discard_bang|record_route_discard!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:3-3`
- [[parallel.policy_telemetry__adaptive_state_for|_adaptive_state_for]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:2-2`
- [[parallel.policy_telemetry__record_policy_decision_bang|_record_policy_decision!]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:30-30`
- [[parallel.policy_telemetry_reset_policy_telemetry_bang|reset_policy_telemetry!]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:78-78`
- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:126-126`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:28-28`
- [[parcore.policy_telemetry_policy_telemetry_snapshot|policy_telemetry_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:89-89`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Using `try`/`catch` for the missing-key case is comparatively slow because exception setup runs on every call from a task without a context; `haskey`-style probing would be cheaper. Task-local storage is inherited by child tasks spawned with `Threads.@spawn` only if Julia copies TLS (it does not by default), so worker tasks created inside a scope observe the global context instead, which is why pools are looked up by the dispatching task, not the workers. Any non-`PolicyContext` value stored under the key is silently ignored.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 1.
