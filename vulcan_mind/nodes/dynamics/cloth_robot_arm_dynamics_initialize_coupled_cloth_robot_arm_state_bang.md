---
id: dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang
label: initialize_coupled_cloth_robot_arm_state!
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: initialize_coupled_cloth_robot_arm_state!
  lines:
  - 202
  - 202
inputs:
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: t_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `t_s` (default `0.0`).
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
  description: Return value of `initialize_coupled_cloth_robot_arm_state!`; mutates
    `sc_view` in place. Returns `nothing`.
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

# initialize_coupled_cloth_robot_arm_state!

## Purpose
Writes initial link positions, quaternions, velocities, and angular rates into a spacecraft state view so the arm starts rigidly attached to the moving base.

## Design & Implementation
Returns `nothing` immediately if `sc_view` lacks `arm_r`. Reads base attitude `q` (default `plan.base_pose.quaternion`), position `pos`, and optional `vel` and `ω` (defaulting to zero). Runs `cloth_fk` with `ClothArmBasePose(base_r, base_q)` at `t_s`, computes `ω_world = R_base * base_ω`, and for each link mutates `sc_view.arm_r[:,i]`, `arm_q[:,i]`, `arm_v[:,i] = base_v + ω_world × (r - base_r)`, and `arm_ω[:,i] = R_link' * ω_world`. The plan's joint rates at `t_s` are not included.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_s` | Real | n/a | no | Keyword argument `t_s` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `initialize_coupled_cloth_robot_arm_state!`; mutates `sc_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2346-2346`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:210-210`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:205-205`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:210-210`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:205-205`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:204-204`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:209-209`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:210-210`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:205-205`
- `callees` → [[vehicle.robotics_clotharmbasepose|ClothArmBasePose]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:209-209`
<!-- vulcan:connections:end -->

## Limitations
Links start with the base's rigid-body velocity only; any planned joint motion at `t_s` causes an initial spring transient. `sc_view.pos` is required (no `hasproperty` guard), so a view without `pos` throws. The base pose is taken from the live view while rest quaternions elsewhere use the plan base, which is consistent only because relative orientations are used.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 202.
