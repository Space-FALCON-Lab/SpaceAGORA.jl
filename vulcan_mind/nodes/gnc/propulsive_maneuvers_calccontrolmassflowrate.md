---
id: gnc.propulsive_maneuvers_calccontrolmassflowrate
label: calcControlMassFlowRate
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: calcControlMassFlowRate
  lines:
  - 336
  - 336
inputs:
- id: controlModel
  type: AbstractControlEffectorModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
- id: u
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `calcControlMassFlowRate`.
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

# calcControlMassFlowRate

## Purpose
Returns the propellant mass flow rate in kilograms per second contributed by a control effector, negative because mass is leaving the vehicle.

## Design & Implementation
Three methods exist. Two fallbacks — one on `AbstractControlEffectorModel` and one untyped — return `0.0`, so any effector that does not consume propellant needs no override. The `BaseThrusterModel` method bounds-checks `i` against `controlModel.thrust`, then gates on the effective burn window from `_effective_burn_window`, requiring `t >= start_time && t <= stop_time`. It rejects a non-finite or non-positive specific impulse, then calls `calcControlForceTorque` and takes `norm(force)` as the applied thrust, returning `0.0` if that is not finite and positive. The result is `-applied_thrust / (Isp * _STANDARD_GRAVITY_MPS2)`, the rocket-equation flow rate with sign convention that mass decreases.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | AbstractControlEffectorModel | n/a | yes | Positional argument `controlModel`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `calcControlMassFlowRate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:400-400`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:395-395`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:345-345`
- `callees` → [[gnc.rpo_mpc_control_model_calcreactionwheeltorque|calcReactionWheelTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:345-345`
<!-- vulcan:connections:end -->

## Limitations
It recomputes the full force vector through `calcControlForceTorque` purely to recover its magnitude, duplicating the burn-window test and the velocity normalisation on every right-hand-side evaluation. The flow rate is not limited by remaining propellant, so a burn scheduled without a propellant check can drive the integrated mass below dry mass. Using the norm of the force discards direction, which is correct only while thrust magnitude equals the model's `thrust` entry.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 336.
