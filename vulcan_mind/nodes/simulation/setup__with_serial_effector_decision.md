---
id: simulation.setup__with_serial_effector_decision
label: _with_serial_effector_decision
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _with_serial_effector_decision
  lines:
  - 749
  - 749
inputs:
- id: effector_decision
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector_decision`.
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
  type: Any
  units: n/a
  description: Return value of `_with_serial_effector_decision`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _with_serial_effector_decision

## Purpose
Rewrites an effector threading decision to serial while preserving its `mode` and `policy_applied` fields, so telemetry still reports what the policy would have chosen even though nested threading is structurally disabled under satellite batching.

## Design & Implementation
Takes any `effector_decision` named tuple and returns `(use_threads=false, allotment=1, mode=effector_decision.mode, policy_applied=effector_decision.policy_applied)`. `@inline`, no allocation beyond the tuple, never throws for inputs carrying the two fields. Applied by `_rhs_execution_plan_uncached` when `_satellite_batch_saturates_pool` is true.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector_decision` | Any | n/a | yes | Positional argument `effector_decision`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_with_serial_effector_decision`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1077-1077`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `policy_applied` may remain `true` while `use_threads` is `false`, telemetry consumers must not infer that the policy itself chose serial. Inputs lacking `mode` or `policy_applied` raise a `FieldError`.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 749.
