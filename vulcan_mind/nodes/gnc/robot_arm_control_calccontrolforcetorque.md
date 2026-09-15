---
id: gnc.robot_arm_control_calccontrolforcetorque
label: calcControlForceTorque
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: calcControlForceTorque
  lines:
  - 205
  - 205
inputs:
- id: model
  type: RobotArmControlEffector
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: Any
  units: n/a
  description: Return value of `calcControlForceTorque`. Returns `model.held.base_force_ii,
    model.held.base_torque_body`.
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

# calcControlForceTorque

## Purpose
Standard control-effector hook that reports the force and torque this effector applies to the spacecraft base at integrator stage time `t`, returning the values cached in `model.held` by the last `calcControlEffect!` so that the dynamics sees a piecewise-constant command between control ticks.

## Design & Implementation
Signature `calcControlForceTorque(model::RobotArmControlEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)`. If `i != model.spacecraft_idx` it returns a tuple of two zero `SVector{3,Float64}`s; otherwise it returns `(model.held.base_force_ii, model.held.base_torque_body)`, an inertial-frame force in newtons and a body-frame torque in newton-metres. The arguments `u`, `p` and `t` are accepted for interface conformity and unused.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RobotArmControlEffector | n/a | yes | Positional argument `model`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `calcControlForceTorque`. Returns `model.held.base_force_ii, model.held.base_torque_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:386-386`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `calcControlEffect!` always stores zero base force and torque, this hook currently contributes nothing to the base dynamics; the arm's reaction wrench reaches the spacecraft only through `RobotArmReactionEffector`. The hook cannot signal staleness of the held command. The zero-tuple return for other spacecraft allocates nothing but is constructed on every call.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 205.
