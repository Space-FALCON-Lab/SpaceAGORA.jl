---
id: dynamics.cloth_robot_arm_dynamics_cloth_reference_state
label: cloth_reference_state
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_reference_state
  lines:
  - 29
  - 29
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
  type: ClothRobotArmReferenceState
  units: n/a
  description: Return value of `cloth_reference_state`.
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

# cloth_reference_state

## Purpose
Samples a `RobotArmPlan` at time `t_s` and evaluates forward kinematics to produce the `ClothRobotArmReferenceState` that the coupled dynamics track.

## Design & Implementation
Calls `robot_arm_plan_sample(plan, t_s)` to obtain joint position `q`, rate `dq`, and end-effector `ee`, then `cloth_fk_state(plan.model, plan.base_pose, sample.q; dq=sample.dq)` to compute link poses and velocities. Returns `ClothRobotArmReferenceState(Float64(t_s), state, sample.ee)`. The return type is annotated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_s` | Real | n/a | yes | Positional argument `t_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothRobotArmReferenceState | n/a | — | Return value of `cloth_reference_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:32-32`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_clothrobotarmreferencestate|ClothRobotArmReferenceState]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:32-32`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:30-30`
- `callees` → [[vehicle.robotics_cloth_fk_state|cloth_fk_state]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:31-31`
<!-- vulcan:connections:end -->

## Limitations
Sampling outside the plan's time range relies on whatever extrapolation `robot_arm_plan_sample` implements; nothing here clamps `t_s`. The forward kinematics run every call with no caching, so callers evaluating many times per step pay repeatedly. The base pose is the plan's static base, not a live spacecraft pose.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 29.
