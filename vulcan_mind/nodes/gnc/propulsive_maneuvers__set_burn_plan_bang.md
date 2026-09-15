---
id: gnc.propulsive_maneuvers__set_burn_plan_bang
label: _set_burn_plan!
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _set_burn_plan!
  lines:
  - 82
  - 82
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: plan
  type: PropulsiveBurnPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
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
  description: Return value of `_set_burn_plan!`; mutates `p` in place. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _set_burn_plan!

## Purpose
Stores a computed `PropulsiveBurnPlan` at index `i` of the shared burn-plan buffer.

## Design & Implementation
Resolves the buffer through `_burn_plan_buffer(p)`, returns `nothing` when the buffer is absent or `i` lies outside `1:length(plans)`, and otherwise performs `plans[i] = plan`, mutating the shared array in place. It always returns `nothing`, so the caller cannot tell a successful store from a skipped one. `calcControlEffect!` calls it after computing the apoapsis-centred burn window.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `plan` | PropulsiveBurnPlan | n/a | yes | Positional argument `plan`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_set_burn_plan!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:558-558`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers__burn_plan_buffer|_burn_plan_buffer]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
Silent failure is the main hazard: when the parameter object has no `maneuver_burn_plans` field the plan is discarded and scheduling appears to succeed, since the thruster model's `start_burn_time` and `stop_burn_time` were already written just before the call. The store is an unlocked write to process-shared state, safe only because each thread owns a distinct `i`.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 82.
