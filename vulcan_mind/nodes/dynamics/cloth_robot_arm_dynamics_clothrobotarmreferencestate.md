---
id: dynamics.cloth_robot_arm_dynamics_clothrobotarmreferencestate
label: ClothRobotArmReferenceState
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: ClothRobotArmReferenceState
  lines:
  - 22
  - 22
inputs:
- id: t_s
  type: Float64
  units: n/a
  required: true
  description: Field `t_s`.
- id: state
  type: ClothArmState
  units: n/a
  required: true
  description: Field `state`.
- id: end_effector_position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `end_effector_position`.
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
  type: ClothRobotArmReferenceState
  units: n/a
  description: Constructed `ClothRobotArmReferenceState`.
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

# ClothRobotArmReferenceState

## Purpose
Immutable sample of a robot-arm plan at one instant: time, full-kinematic `ClothArmState`, and end-effector position, used as the tracking reference for cloth coupling.

## Design & Implementation
Three fields: `t_s::Float64` (seconds), `state::ClothArmState` (joint angles, rates, and forward-kinematic link poses produced by `cloth_fk_state`), and `end_effector_position::SVector{3,Float64}` (metres, world frame). Constructed only by `cloth_reference_state`, which converts the caller's `t_s::Real` to `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_s` | Float64 | n/a | yes | Field `t_s`. |
| in | `state` | ClothArmState | n/a | yes | Field `state`. |
| in | `end_effector_position` | SVector{3, Float64} | n/a | yes | Field `end_effector_position`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothRobotArmReferenceState | n/a | — | Constructed `ClothRobotArmReferenceState`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_reference_state|cloth_reference_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:32-32`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct carries no joint acceleration even though its docstring mentions one; only what `cloth_fk_state` returns is stored. There is no validation that `end_effector_position` matches the state's kinematics. Being immutable, per-step reference updates allocate a fresh instance.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 22.
