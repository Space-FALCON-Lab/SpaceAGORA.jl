---
id: dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan
label: simulate_cloth_robot_arm_plan
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: simulate_cloth_robot_arm_plan
  lines:
  - 437
  - 437
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: dt_s
  type: Real
  units: n/a
  required: false
  description: 'Keyword argument `dt_s` (default `(length(plan.t_ref_s) >= 2 ? plan.t_ref_s[2]
    - plan.t_ref_s[1] : 0.1)`).'
- id: duration_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `duration_s` (default `plan.t_ref_s[end]`).
- id: integrator
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `integrator` (default `:implicit_midpoint`).
- id: k_translation_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_translation_n_m` (default `5.0e3`).
- id: c_translation_n_s_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_translation_n_s_m` (default `30.0`).
- id: k_rotation_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_rotation_n_m_rad` (default `15.0`).
- id: c_rotation_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_rotation_n_m_s_rad` (default `0.5`).
- id: joint_actuators
  type: Union{Nothing, AbstractVector{CompliantJointActuator}}
  units: n/a
  required: false
  description: Keyword argument `joint_actuators` (default `nothing`).
- id: actuator_torque_limit_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `actuator_torque_limit_n_m` (default `Inf`).
- id: actuator_kp_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `actuator_kp_n_m_rad` (default `0.0`).
- id: actuator_kd_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `actuator_kd_n_m_s_rad` (default `0.0`).
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
  type: ClothRobotArmSimulation
  units: n/a
  description: Return value of `simulate_cloth_robot_arm_plan`. Returns `ClothRobotArmSimulation(`.
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

# simulate_cloth_robot_arm_plan

## Purpose
Runs a standalone compliant simulation of a robot-arm plan with a fixed base, stepping with implicit midpoint or RK4, and returns end-effector tracking and joint torque histories in a `ClothRobotArmSimulation`.

## Design & Implementation
Builds the model via `cloth_robot_arm_multibody` and actuators via `cloth_robot_arm_actuators` unless `joint_actuators` is supplied. Default `dt_s` is the plan's first time spacing (or 0.1 s) and `duration_s` is `plan.t_ref_s[end]`; `times = collect(0:dt:duration)` with the end appended if missing. For each step it recomputes rest quaternions at `t_prev` and calls `step_compliant_multibody_implicit_midpoint` or `step_compliant_multibody_rk4`, throwing `ArgumentError` for any other `integrator` symbol. Afterwards it evaluates `cloth_robot_arm_end_effector`, the reference `ee` from the plan, the error norm, and `compliant_joint_loads` per sample, filling 3 x n_joints x N_t torque arrays.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `dt_s` | Real | n/a | no | Keyword argument `dt_s` (default `(length(plan.t_ref_s) >= 2 ? plan.t_ref_s[2] - plan.t_ref_s[1] : 0.1)`). |
| in | `duration_s` | Real | n/a | no | Keyword argument `duration_s` (default `plan.t_ref_s[end]`). |
| in | `integrator` | Symbol | n/a | no | Keyword argument `integrator` (default `:implicit_midpoint`). |
| in | `k_translation_n_m` | Any | n/a | no | Keyword argument `k_translation_n_m` (default `5.0e3`). |
| in | `c_translation_n_s_m` | Any | n/a | no | Keyword argument `c_translation_n_s_m` (default `30.0`). |
| in | `k_rotation_n_m_rad` | Any | n/a | no | Keyword argument `k_rotation_n_m_rad` (default `15.0`). |
| in | `c_rotation_n_m_s_rad` | Any | n/a | no | Keyword argument `c_rotation_n_m_s_rad` (default `0.5`). |
| in | `joint_actuators` | Union{Nothing, AbstractVector{CompliantJointActuator}} | n/a | no | Keyword argument `joint_actuators` (default `nothing`). |
| in | `actuator_torque_limit_n_m` | Any | n/a | no | Keyword argument `actuator_torque_limit_n_m` (default `Inf`). |
| in | `actuator_kp_n_m_rad` | Any | n/a | no | Keyword argument `actuator_kp_n_m_rad` (default `0.0`). |
| in | `actuator_kd_n_m_s_rad` | Any | n/a | no | Keyword argument `actuator_kd_n_m_s_rad` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothRobotArmSimulation | n/a | — | Return value of `simulate_cloth_robot_arm_plan`. Returns `ClothRobotArmSimulation(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:466-466`
- `callees` → [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:504-504`
- `callees` → [[dynamics.cloth_multibody_compliantmultibodytrajectory|CompliantMultibodyTrajectory]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:518-518`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_actuators|cloth_robot_arm_actuators]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:459-459`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_end_effector|cloth_robot_arm_end_effector]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:501-501`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state|cloth_robot_arm_initial_state]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:473-473`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:451-451`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:476-476`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_clothrobotarmsimulation|ClothRobotArmSimulation]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:516-516`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:470-470`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:502-502`
<!-- vulcan:connections:end -->

## Limitations
The integrator symbol is validated inside the loop, so an invalid choice is only detected on the second time sample and after model construction. Rest quaternions are held at `t_prev` across each step (explicit in time) and recomputed twice per sample overall. The base is fixed at `plan.base_pose`; spacecraft coupling must go through `assign_coupled_cloth_robot_arm_rhs!` instead. All states are stored as `Vector{Float64}` copies, so memory grows linearly with `N_t`.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 437.
