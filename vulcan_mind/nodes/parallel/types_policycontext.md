---
id: parallel.types_policycontext
label: PolicyContext
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: PolicyContext
  lines:
  - 89
  - 89
inputs:
- id: telemetry
  type: PolicyTelemetry
  units: n/a
  required: false
  description: Field `telemetry` (default `PolicyTelemetry()`).
- id: adaptive_state
  type: Dict{Symbol, AdaptiveControllerState}
  units: n/a
  required: false
  description: Field `adaptive_state` (default `Dict{Symbol, AdaptiveControllerState}()`).
- id: decision_signature
  type: Dict{Symbol, String}
  units: n/a
  required: false
  description: Field `decision_signature` (default `Dict{Symbol, String}()`).
- id: decision_allotment
  type: Dict{Symbol, Int64}
  units: n/a
  required: false
  description: Field `decision_allotment` (default `Dict{Symbol, Int64}()`).
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
  description: Constructed `PolicyContext` (keyword constructor via @kwdef).
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

# PolicyContext

## Purpose
Bundle of mutable policy state that the parallel dispatch layer consults on every threading decision. It groups the `PolicyTelemetry` counters with the per-source `AdaptiveControllerState` map and the last decision signature and allotment per source, so a decision and its later timing observation can be matched.

## Design & Implementation
`Base.@kwdef mutable struct` whose fields default to fresh empty containers: `telemetry::PolicyTelemetry`, `adaptive_state::Dict{Symbol, AdaptiveControllerState}`, `decision_signature::Dict{Symbol, String}` and `decision_allotment::Dict{Symbol, Int64}`. A global default lives in `_global_policy_context::Ref{PolicyContext}`; a task-scoped override is installed with `Base.task_local_storage(_policy_context_tls_key, ctx)` in `context.jl` so nested or concurrent simulations can isolate their statistics. `adaptive_decision.jl` writes `decision_signature[source]` and `decision_allotment[source]` after each decision; `observation_tracking.jl` reads them back with `get(..., "")` defaults; `policy_telemetry.jl` empties both dicts on reset.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `telemetry` | PolicyTelemetry | n/a | no | Field `telemetry` (default `PolicyTelemetry()`). |
| in | `adaptive_state` | Dict{Symbol, AdaptiveControllerState} | n/a | no | Field `adaptive_state` (default `Dict{Symbol, AdaptiveControllerState}()`). |
| in | `decision_signature` | Dict{Symbol, String} | n/a | no | Field `decision_signature` (default `Dict{Symbol, String}()`). |
| in | `decision_allotment` | Dict{Symbol, Int64} | n/a | no | Field `decision_allotment` (default `Dict{Symbol, Int64}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PolicyContext | n/a | — | Constructed `PolicyContext` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.context_with_policy_context|with_policy_context]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:333-333`

**Downstream**

- `callees` → [[parcore.types_policytelemetry|PolicyTelemetry]] · `callers` · call · `src/parallel/policy/types.jl:90-90`
<!-- vulcan:connections:end -->

## Limitations
The dicts are unsynchronised; callers must hold `_policy_telemetry_lock` for every read and write, and nothing in the type checks that. Because the context is looked up from task-local storage, a `Threads.@spawn`ed task does not inherit the parent's override and silently falls back to the global context. Only one outstanding decision per `source` is remembered, so overlapping dispatches of the same source overwrite each other's signature.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 89.
