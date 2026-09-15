---
id: dynamics.robot_arm_reaction_effector_robotarmreactioneffectors
label: RobotArmReactionEffectors
kind: module
source:
  file: src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl
  symbol: RobotArmReactionEffectors
  lines:
  - 2
  - 2
inputs:
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
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# RobotArmReactionEffectors

## Purpose

`RobotArmReactionEffectors` is the module that couples robot-arm joint loads back into the spacecraft force and torque stack. It defines and exports `RobotArmReactionEffector`, a mutable `AbstractForceTorqueModel` holding the arm plan, the compliant joint actuators and the impedance gains, and it adds a `calcForceTorque` method so the arm participates in the same effector dispatch as gravity, aerodynamic and thruster models.

## Design & Implementation

The module imports `SVector`/`SMatrix` from StaticArrays, `AbstractForceTorqueModel` from `AbstractTypes`, `CompliantJointActuator` from `ClothMultibody`, and `RobotArmPlan` plus `robot_arm_plan_sample` from `RobotArmPlanning`; it imports `calcForceTorque` from the parent `DynamicEffectors` so its method extends the existing generic function rather than shadowing it. `RobotArmReactionEffector` is declared with `Base.@kwdef mutable struct` so every field has a default and the plan can be swapped in at runtime: `spacecraft_idx = 1` selects the owning satellite, `plan` is `Union{Nothing, RobotArmPlan}` and starts `nothing`, `updated_at_s`, `force_scale` and `torque_scale` all start at `0.0`, the impedance gains default to `k_translation_n_m = 5.0e3` N/m, `c_translation_n_s_m = 30.0` N·s/m, `k_rotation_n_m_rad = 15.0` N·m/rad and `c_rotation_n_m_s_rad = 0.5` N·m·s/rad, and `joint_actuators` is an empty `Vector{CompliantJointActuator}`. Only `RobotArmReactionEffector` is exported.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The module is a coupling seam that is not yet load-bearing: its `calcForceTorque` method returns zero wrenches unconditionally, so the gains, `joint_actuators`, `force_scale`, `torque_scale` and `updated_at_s` fields are stored but never read by any code in this file, and `robot_arm_plan_sample` is imported without being called. The four impedance gains are typed `Any` rather than `Float64`, which defeats concrete field inference for the struct and forces boxed access once they are used. Only one spacecraft index is served per effector instance, so a multi-arm vehicle needs one instance per arm.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl` line 2.
