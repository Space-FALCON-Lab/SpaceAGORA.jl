---
id: gnc.propulsive_maneuvers__effective_burn_window
label: _effective_burn_window
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _effective_burn_window
  lines:
  - 149
  - 149
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
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
  description: Return value of `_effective_burn_window`. Returns `start_time, stop_time`
    or `plan.start_burn_s, plan.stop_burn_s`.
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

# _effective_burn_window

## Purpose
Determines the burn window actually in force for spacecraft `i`, preferring the thruster model's own schedule and falling back to the active burn plan.

## Design & Implementation
Calls `_model_burn_window`; if that yields a window where both endpoints are finite and `stop_time > start_time`, it is returned unchanged. Otherwise it consults `_active_burn_plan(p, i)` and returns `(plan.start_burn_s, plan.stop_burn_s)` when a valid plan exists. With neither available it returns the model window as-is, which may be `(NaN, NaN)` or the cleared sentinel `(-1.0, -1.0)`. Both `calcControlForceTorque` and `calcControlMassFlowRate` use it to gate thrust on `t >= start_time && t <= stop_time`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_effective_burn_window`. Returns `start_time, stop_time` or `plan.start_burn_s, plan.stop_burn_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:378-378`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers__active_burn_plan|_active_burn_plan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:154-154`
- `callees` → [[gnc.propulsive_maneuvers__model_burn_window|_model_burn_window]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:150-150`
<!-- vulcan:connections:end -->

## Limitations
Precedence favours the model over the plan, so a stale model schedule shadows a freshly planned burn until the model arrays are cleared. The returned window is used in a closed-interval comparison with no tolerance, whereas `calcControlEffect!` applies a 1e-9 second pad, so the two paths can disagree about membership exactly at an endpoint. A degenerate window with `stop_time == start_time` fires for a single instant.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 149.
