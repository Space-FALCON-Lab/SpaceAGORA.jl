---
id: dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_actuators
label: cloth_robot_arm_actuators
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_robot_arm_actuators
  lines:
  - 274
  - 274
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: torque_limit_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `torque_limit_n_m` (default `Inf`).
- id: kp_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `kp_n_m_rad` (default `0.0`).
- id: kd_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `kd_n_m_s_rad` (default `0.0`).
- id: feedforward_torque_child_body
  type: Any
  units: n/a
  required: false
  description: Keyword argument `feedforward_torque_child_body` (default `(0.0, 0.0,
    0.0)`).
- id: efficiency
  type: Real
  units: n/a
  required: false
  description: Keyword argument `efficiency` (default `1.0`).
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
  description: Return value of `cloth_robot_arm_actuators`. Returns `CompliantJointActuator[`.
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

# cloth_robot_arm_actuators

## Purpose
Creates one `CompliantJointActuator` per arm joint with shared PD gains, torque limit, feed-forward torque, and efficiency, for tracking the planned rest orientation.

## Design & Implementation
Comprehension over `eachindex(plan.model.links)` calling `CompliantJointActuator(Symbol("cloth_arm_actuator_i"), i; torque_limit_n_m, kp_n_m_rad, kd_n_m_s_rad, feedforward_torque_child_body, efficiency)`. Defaults are `Inf` limit, zero gains, zero feed-forward `(0,0,0)`, and `efficiency=1.0`, which makes the default actuators inert.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `torque_limit_n_m` | Any | n/a | no | Keyword argument `torque_limit_n_m` (default `Inf`). |
| in | `kp_n_m_rad` | Any | n/a | no | Keyword argument `kp_n_m_rad` (default `0.0`). |
| in | `kd_n_m_s_rad` | Any | n/a | no | Keyword argument `kd_n_m_s_rad` (default `0.0`). |
| in | `feedforward_torque_child_body` | Any | n/a | no | Keyword argument `feedforward_torque_child_body` (default `(0.0, 0.0, 0.0)`). |
| in | `efficiency` | Real | n/a | no | Keyword argument `efficiency` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_robot_arm_actuators`. Returns `CompliantJointActuator[`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:459-459`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody_compliantjointactuator|CompliantJointActuator]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:283-283`
<!-- vulcan:connections:end -->

## Limitations
All joints receive identical parameters; per-joint tuning requires constructing actuators manually. Gains are not validated for sign, so negative `kp` produces a destabilising actuator. The feed-forward torque is constant in time.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 274.
