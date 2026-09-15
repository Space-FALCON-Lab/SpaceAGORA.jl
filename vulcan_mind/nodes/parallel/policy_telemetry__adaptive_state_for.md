---
id: parallel.policy_telemetry__adaptive_state_for
label: _adaptive_state_for
kind: function
source:
  file: src/parallel/policy/policy_telemetry.jl
  symbol: _adaptive_state_for
  lines:
  - 1
  - 1
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
  type: AdaptiveControllerState
  units: n/a
  description: Return value of `_adaptive_state_for`.
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

# _adaptive_state_for

## Purpose
Returns the adaptive controller state belonging to one threading decision source, creating it on first use so call sites never have to initialise it.

## Design & Implementation
Resolves the live policy context through `_active_policy_context()`, then uses `get!` with a zero-argument closure so a fresh `AdaptiveControllerState()` is constructed only when `source` has no entry yet. Keying by source symbol is what lets the density, control and multibody loops adapt their allotments independently instead of sharing one controller. Declared `@inline` with a concrete return type because it sits on the decision path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AdaptiveControllerState | n/a | — | Return value of `_adaptive_state_for`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:48-48`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:48-48`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:2-2`
- `callees` → [[parallel.types_adaptivecontrollerstate|AdaptiveControllerState]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:4-4`
<!-- vulcan:connections:end -->

## Limitations
The lookup is unsynchronised, so two threads first touching the same source concurrently can race on the dictionary; callers reach it through the telemetry lock in practice, but nothing in this function enforces that.

## Provenance
Mapped from `src/parallel/policy/policy_telemetry.jl` line 1.
