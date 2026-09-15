---
id: gnc.propulsive_maneuvers__effective_thrust_isp
label: _effective_thrust_isp
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _effective_thrust_isp
  lines:
  - 172
  - 172
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
  description: Return value of `_effective_thrust_isp`. Returns `plan.thrust_n, plan.isp_s`
    or `NaN, NaN` or `Float64(controlModel.thrust[i]), Float64(controlModel.Isp[i])`.
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

# _effective_thrust_isp

## Purpose
Resolves the thrust magnitude in newtons and the specific impulse in seconds for spacecraft `i`, preferring the active burn plan over the thruster model arrays.

## Design & Implementation
Returns `(plan.thrust_n, plan.isp_s)` when `_active_burn_plan(p, i)` yields a valid plan. Otherwise it requires `i` to be in range for both `controlModel.thrust` and `controlModel.Isp`, returning `(NaN, NaN)` if either check fails, and otherwise `(Float64(controlModel.thrust[i]), Float64(controlModel.Isp[i]))`. Callers include `_validated_burn_plan`, `calcControlForceTorque` and `calcControlMassFlowRate`, all of which reject non-finite or non-positive values themselves.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effective_thrust_isp`. Returns `plan.thrust_n, plan.isp_s` or `NaN, NaN` or `Float64(controlModel.thrust[i]), Float64(controlModel.Isp[i])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__validated_burn_plan|_validated_burn_plan]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:205-205`
- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:380-380`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:180-180`
- `callees` → [[gnc.propulsive_maneuvers__active_burn_plan|_active_burn_plan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:173-173`
<!-- vulcan:connections:end -->

## Limitations
Thrust and specific impulse are constants per spacecraft, so throttling, tank blowdown and the dependence of specific impulse on back-pressure or thruster duty cycle are not represented. The `(NaN, NaN)` sentinel is returned for an out-of-range index rather than throwing, so a mis-sized model array degrades into zero thrust instead of an error.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 172.
