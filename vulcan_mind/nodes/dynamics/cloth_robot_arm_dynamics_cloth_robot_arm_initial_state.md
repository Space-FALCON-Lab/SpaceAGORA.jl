---
id: dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state
label: cloth_robot_arm_initial_state
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_robot_arm_initial_state
  lines:
  - 179
  - 179
inputs:
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
  type: Any
  units: n/a
  description: Return value of `cloth_robot_arm_initial_state`. Returns `compliant_state_vector(`.
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

# cloth_robot_arm_initial_state

## Purpose
Packs the compliant multibody state vector for a standalone arm simulation from the plan's pose at `t_s`, with all link velocities set to zero.

## Design & Implementation
Samples the plan, computes `cloth_fk` link centre-of-mass positions and quaternions, and calls `compliant_state_vector(positions, quaternions; velocities=zeros, angular_velocities=zeros)` where the zero vectors are `fill(SVector{3,Float64}(0,0,0), n)`. Default `t_s=0.0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_s` | Real | n/a | no | Keyword argument `t_s` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_robot_arm_initial_state`. Returns `compliant_state_vector(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:473-473`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody_compliant_state_vector|compliant_state_vector]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:182-182`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:180-180`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:181-181`
<!-- vulcan:connections:end -->

## Limitations
Velocities are always zero regardless of the plan's joint rates at `t_s`, so starting mid-plan introduces an initial tracking transient. The returned layout is whatever `compliant_state_vector` produces and must match `compliant_state_parts` used elsewhere.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 179.
