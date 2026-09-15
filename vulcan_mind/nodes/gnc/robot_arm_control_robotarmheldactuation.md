---
id: gnc.robot_arm_control_robotarmheldactuation
label: RobotArmHeldActuation
kind: struct
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: RobotArmHeldActuation
  lines:
  - 2
  - 2
inputs:
- id: joint_torque_nm
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `joint_torque_nm` (default `Float64[]`).
- id: base_force_ii
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `base_force_ii` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`).
- id: base_torque_body
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `base_torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`).
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
  type: RobotArmHeldActuation
  units: n/a
  description: Constructed `RobotArmHeldActuation` (keyword constructor via @kwdef).
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

# RobotArmHeldActuation

## Purpose
Zero-order-hold container for the most recent robot-arm command, so that the control effector can compute joint torques once per controller update in `calcControlEffect!` and the dynamics RHS can read the same command at every intermediate integrator stage without recomputation.

## Design & Implementation
A `Base.@kwdef mutable struct` with three fields: `joint_torque_nm::Vector{Float64}` (one entry per arm joint, default empty), `base_force_ii::SVector{3,Float64}` (inertial-frame force on the spacecraft base, N) and `base_torque_body::SVector{3,Float64}` (body-frame torque, N m), both defaulting to zero vectors. `RobotArmControlEffector.held` stores one instance; `calcControlEffect!` replaces it wholesale with a new object rather than mutating fields, and `calcControlForceTorque` returns the base force and torque pair.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `joint_torque_nm` | Vector{Float64} | n/a | no | Field `joint_torque_nm` (default `Float64[]`). |
| in | `base_force_ii` | SVector{3, Float64} | n/a | no | Field `base_force_ii` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `base_torque_body` | SVector{3, Float64} | n/a | no | Field `base_torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmHeldActuation | n/a | — | Constructed `RobotArmHeldActuation` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:196-196`
- [[gnc.robot_arm_control_robotarmcontroleffector|RobotArmControlEffector]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:33-33`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only `joint_torque_nm` is ever populated by the controller in this file; `base_force_ii` and `base_torque_body` are always written as zero, so the reaction of the arm on the base must come from the separate `RobotArmReactionEffector` dynamic effector. Being mutable and replaced by reference, a reader that captured the previous `held` object will keep seeing stale values. There is no timestamp, so consumers cannot tell how old the held command is.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 2.
