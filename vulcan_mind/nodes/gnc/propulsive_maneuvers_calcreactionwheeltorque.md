---
id: gnc.propulsive_maneuvers_calcreactionwheeltorque
label: calcReactionWheelTorque
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: calcReactionWheelTorque
  lines:
  - 353
  - 353
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
  type: Nothing
  units: n/a
  description: Return value of `calcReactionWheelTorque`. Returns `nothing`.
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

# calcReactionWheelTorque

## Purpose
Optional hook letting an effector declare which part of its commanded body torque is produced through a reaction wheel, so that wheel momentum can be accumulated rather than treated as an ideal external torque.

## Design & Implementation
Two methods are defined here, one on `AbstractControlEffectorModel` and one fully untyped, and both return `nothing`. Returning `nothing` is the documented signal that the effector does not drive a reaction wheel, so every effector type that does not override the method is handled by the default. Effectors that do react against stored wheel angular momentum are expected to add a method returning the wheel-borne torque component.

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
| out | `result` | Nothing | n/a | — | Return value of `calcReactionWheelTorque`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:345-345`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:197-197`

**Downstream**

- `callees` → [[gnc.momentum_manager_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callers` · feedback · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- `callees` → [[gnc.robot_arm_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- `callees` → [[gnc.targeting_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
<!-- vulcan:connections:end -->

## Limitations
The default is `nothing` rather than a zero vector, so every consumer must branch on the type of the returned value instead of accumulating unconditionally. There is no counterpart hook reporting wheel momentum or a saturation limit, so an overriding effector cannot signal that its wheels have run out of authority, and no method here validates that the returned torque is a subset of the torque reported by `calcControlForceTorque`.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 353.
