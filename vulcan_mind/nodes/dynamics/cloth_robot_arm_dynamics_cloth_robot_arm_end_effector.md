---
id: dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_end_effector
label: cloth_robot_arm_end_effector
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_robot_arm_end_effector
  lines:
  - 427
  - 427
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `cloth_robot_arm_end_effector`. Returns `state.r +
    _rot(state.q) * tip_body`.
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

# cloth_robot_arm_end_effector

## Purpose
Computes the world-frame tip position of the last link from a packed compliant state vector `x`, for tracking-error evaluation.

## Design & Implementation
With `n = length(plan.model.links)`, returns `plan.base_pose.position` when `n == 0`. Otherwise unpacks `state = compliant_state_parts(x, n)`, forms `tip_body = link.vector_parent - link.com_offset_parent` for the last link, and returns `state.r + _rot(state.q) * tip_body`. The docstring mentions orientation but only position is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_robot_arm_end_effector`. Returns `state.r + _rot(state.q) * tip_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:501-501`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:433-433`
- `callees` → [[dynamics.cloth_multibody_compliant_state_parts|compliant_state_parts]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:430-430`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:433-433`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:433-433`
<!-- vulcan:connections:end -->

## Limitations
Relies on `compliant_state_parts(x, n)` returning the last body's `r` and `q`; if it returns all bodies the indexing is implicit in that helper. Uses the plan's base position for the empty case rather than a live base. No orientation output despite the docstring.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 427.
