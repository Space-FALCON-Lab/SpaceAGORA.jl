---
id: dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions
label: cloth_robot_arm_rest_quaternions
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_robot_arm_rest_quaternions
  lines:
  - 165
  - 165
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: t_s
  type: Real
  units: n/a
  required: true
  description: Positional argument `t_s`.
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
  description: Return value of `cloth_robot_arm_rest_quaternions`. Returns `rests`.
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

# cloth_robot_arm_rest_quaternions

## Purpose
Evaluates the plan's forward kinematics at `t_s` and returns the child-to-parent rest quaternion for every link, defining where each compliant joint's rotational spring is unloaded.

## Design & Implementation
Samples `robot_arm_plan_sample(plan, t_s)`, runs `cloth_fk(plan.model, plan.base_pose, sample.q)`, then walks links in order with `parent_q` starting at `plan.base_pose.quaternion` and updating to each `pose.link_quaternions[i]`. Each entry is `_rest_child_parent_quat(parent_q, child_q)`. Returns a `Vector{SVector{4,Float64}}` of length `length(plan.model.links)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_s` | Real | n/a | yes | Positional argument `t_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_robot_arm_rest_quaternions`. Returns `rests`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:242-242`
- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:476-476`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:357-357`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rest_child_parent_quat|_rest_child_parent_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:172-172`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat|_rest_child_parent_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:172-172`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:166-166`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:167-167`
<!-- vulcan:connections:end -->

## Limitations
Uses the plan's static `base_pose`, not the live spacecraft attitude, so in coupled mode the rest orientation is relative and correct only because the RHS composes it with the actual parent quaternion. Called once per RHS evaluation and per output sample with full FK each time; no memoisation. Assumes a serial chain where link i-1 is the parent of link i.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 165.
