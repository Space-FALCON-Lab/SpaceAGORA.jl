---
id: gnc.robot_arm_control_robotarmcontroleffector
label: RobotArmControlEffector
kind: struct
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: RobotArmControlEffector
  lines:
  - 21
  - 21
inputs:
- id: spacecraft_idx
  type: Int
  units: n/a
  required: false
  description: Field `spacecraft_idx` (default `1`).
- id: plan
  type: Union{Nothing, RobotArmPlan}
  units: n/a
  required: false
  description: Field `plan` (default `nothing`).
- id: controller
  type: Any
  units: n/a
  required: false
  description: Field `controller` (default `nothing`).
- id: updated_at_s
  type: Float64
  units: n/a
  required: false
  description: Field `updated_at_s` (default `0.0`).
- id: joint_kp
  type: Float64
  units: n/a
  required: false
  description: Field `joint_kp` (default `0.0`).
- id: joint_kd
  type: Float64
  units: n/a
  required: false
  description: Field `joint_kd` (default `0.0`).
- id: control_dt_s
  type: Float64
  units: n/a
  required: false
  description: Field `control_dt_s` (default `0.1`).
- id: k_translation_n_m
  type: Float64
  units: n/a
  required: false
  description: Field `k_translation_n_m` (default `5.0e3`).
- id: c_translation_n_s_m
  type: Float64
  units: n/a
  required: false
  description: Field `c_translation_n_s_m` (default `30.0`).
- id: k_rotation_n_m_rad
  type: Float64
  units: n/a
  required: false
  description: Field `k_rotation_n_m_rad` (default `15.0`).
- id: c_rotation_n_m_s_rad
  type: Float64
  units: n/a
  required: false
  description: Field `c_rotation_n_m_s_rad` (default `0.5`).
- id: held
  type: RobotArmHeldActuation
  units: n/a
  required: false
  description: Field `held` (default `RobotArmHeldActuation()`).
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
  type: RobotArmControlEffector
  units: n/a
  description: Constructed `RobotArmControlEffector` (keyword constructor via @kwdef).
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

# RobotArmControlEffector

## Purpose
Control effector that tracks a `RobotArmPlan` for a single spacecraft and publishes joint torque commands through a held actuation record. It is the control-tier counterpart of `RobotArmReactionEffector`, which applies the resulting reactions in the dynamics.

## Design & Implementation
A `Base.@kwdef mutable struct RobotArmControlEffector <: AbstractControlEffectorModel` with fields `spacecraft_idx::Int = 1`, `plan::Union{Nothing,RobotArmPlan} = nothing`, `controller::Any = nothing` (expected to be a `RobotArmJointMPCController` from `init_robot_arm_joint_mpc`), `updated_at_s::Float64 = 0.0` (simulation time at which the plan clock starts), `joint_kp`, `joint_kd` (both default 0.0), `control_dt_s = 0.1` (MPC sample time), base compliance gains `k_translation_n_m = 5.0e3`, `c_translation_n_s_m = 30.0`, `k_rotation_n_m_rad = 15.0`, `c_rotation_n_m_s_rad = 0.5`, and `held::RobotArmHeldActuation`. `calcControlEffect!` fills `held`; `calcControlForceTorque` and `calcControlMassFlowRate` implement the standard control-effector hooks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft_idx` | Int | n/a | no | Field `spacecraft_idx` (default `1`). |
| in | `plan` | Union{Nothing, RobotArmPlan} | n/a | no | Field `plan` (default `nothing`). |
| in | `controller` | Any | n/a | no | Field `controller` (default `nothing`). |
| in | `updated_at_s` | Float64 | n/a | no | Field `updated_at_s` (default `0.0`). |
| in | `joint_kp` | Float64 | n/a | no | Field `joint_kp` (default `0.0`). |
| in | `joint_kd` | Float64 | n/a | no | Field `joint_kd` (default `0.0`). |
| in | `control_dt_s` | Float64 | n/a | no | Field `control_dt_s` (default `0.1`). |
| in | `k_translation_n_m` | Float64 | n/a | no | Field `k_translation_n_m` (default `5.0e3`). |
| in | `c_translation_n_s_m` | Float64 | n/a | no | Field `c_translation_n_s_m` (default `30.0`). |
| in | `k_rotation_n_m_rad` | Float64 | n/a | no | Field `k_rotation_n_m_rad` (default `15.0`). |
| in | `c_rotation_n_m_s_rad` | Float64 | n/a | no | Field `c_rotation_n_m_s_rad` (default `0.5`). |
| in | `held` | RobotArmHeldActuation | n/a | no | Field `held` (default `RobotArmHeldActuation()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmControlEffector | n/a | — | Constructed `RobotArmControlEffector` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[gnc.robot_arm_control_robotarmheldactuation|RobotArmHeldActuation]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
`joint_kp`, `joint_kd` and the four base compliance gains are stored but not referenced anywhere in this file, so setting them has no effect on the torques produced here. With `controller === nothing` the effector silently commands zero torque rather than warning that no controller is attached. `controller::Any` means a wrong object type is only detected by the `isa` check falling through to the zero-torque path. The effector is mutable and per-spacecraft, so it must be deep-copied for concurrent runs.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 21.
