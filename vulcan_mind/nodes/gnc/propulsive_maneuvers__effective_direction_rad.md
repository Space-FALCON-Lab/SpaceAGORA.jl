---
id: gnc.propulsive_maneuvers__effective_direction_rad
label: _effective_direction_rad
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _effective_direction_rad
  lines:
  - 161
  - 161
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
  type: Float64
  units: n/a
  description: Return value of `_effective_direction_rad`.
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

# _effective_direction_rad

## Purpose
Resolves the thrust direction angle in radians for spacecraft `i`, preferring the active burn plan over the thruster model.

## Design & Implementation
Consults `_active_burn_plan(p, i)` first and returns `plan.direction_rad` when a valid plan exists. Otherwise it bounds-checks `i` against `length(controlModel.direction)`, returning `NaN` if out of range, and returns `Float64(controlModel.direction[i])`. The consumer `calcControlForceTorque` reduces this angle to a sign through `cos(direction_rad) >= 0.0`, so only the hemisphere matters: near zero means prograde, near `π` means retrograde.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effective_direction_rad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:383-383`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:169-169`
- `callees` → [[gnc.propulsive_maneuvers__active_burn_plan|_active_burn_plan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:162-162`
<!-- vulcan:connections:end -->

## Limitations
Note the precedence is the opposite of `_effective_burn_window`: the plan wins here but loses there, so a plan and a model schedule that disagree yield a window from one source and a direction from the other. Because the downstream consumer only takes the sign of the cosine, an angle of 1.2 radians thrusts fully prograde rather than off-axis, silently discarding the intended obliquity.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 161.
