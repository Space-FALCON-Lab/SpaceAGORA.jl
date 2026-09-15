---
id: dynamics.cloth_robot_arm_dynamics_coupled_cloth_robot_arm_state_shape
label: coupled_cloth_robot_arm_state_shape
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: coupled_cloth_robot_arm_state_shape
  lines:
  - 191
  - 191
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
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
  description: Return value of `coupled_cloth_robot_arm_state_shape`. Returns `(`.
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

# coupled_cloth_robot_arm_state_shape

## Purpose
Returns a NamedTuple of zero-filled arrays describing the per-link state blocks a spacecraft state vector must reserve for a coupled arm.

## Design & Implementation
With `n = length(plan.model.links)`, returns `(arm_r=zeros(3,n), arm_q=zeros(4,n), arm_v=zeros(3,n), arm_ω=zeros(3,n))`: link centre-of-mass positions (m, world), scalar-last attitude quaternions, world velocities (m/s), and body angular rates (rad/s). The state-layout machinery uses these shapes to allocate views named `arm_r`, `arm_q`, `arm_v`, `arm_ω`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `coupled_cloth_robot_arm_state_shape`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2317-2317`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The quaternion block is initialised to all zeros, which is not a valid attitude; `initialize_coupled_cloth_robot_arm_state!` must run before integration or `_unit_quat` will substitute identity silently. A plan with zero links yields empty 3 x 0 arrays, which downstream code must tolerate.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 191.
