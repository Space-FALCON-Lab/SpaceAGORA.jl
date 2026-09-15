---
id: gnc.robot_arm_planning_robot_arm_plan_sample
label: robot_arm_plan_sample
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: robot_arm_plan_sample
  lines:
  - 145
  - 145
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
  type: Tuple
  units: n/a
  description: Return value of `robot_arm_plan_sample`. Returns `(` or `(q=q, dq=dq,
    ddq=ddq, ee=SVector{3, Float64}(ee))`.
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

# robot_arm_plan_sample

## Purpose
Evaluates a stored `RobotArmPlan` at an arbitrary time, returning joint position, velocity, acceleration and end-effector position, so controllers can track the reference at their own rate.

## Design & Implementation
Takes `plan::RobotArmPlan` and `t_s::Real`, converted to `Float64`. Times at or before `plan.t_ref_s[1]` return copies of the first column of `q_ref`, `dq_ref`, `ddq_ref` and an `SVector{3}` from `ee_ref`; times at or after the last grid point return the last column. Otherwise `searchsortedfirst` finds the bracketing index `hi` (`lo = hi-1`), computes `α = (t - t_lo)/(t_hi - t_lo)` and linearly interpolates all four arrays. The result is a NamedTuple `(q, dq, ddq, ee)` with `ee::SVector{3,Float64}` in metres and joint quantities in rad, rad/s, rad/s^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_s` | Real | n/a | yes | Positional argument `t_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `robot_arm_plan_sample`. Returns `(` or `(q=q, dq=dq, ddq=ddq, ee=SVector{3, Float64}(ee))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_reference_state|cloth_reference_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:30-30`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state|cloth_robot_arm_initial_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:180-180`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:166-166`
- [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:204-204`
- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:502-502`
- [[gnc.robot_arm_control__robot_arm_control_reference_state|_robot_arm_control_reference_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:173-173`
- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:187-187`
- [[gnc.robot_arm_control_robot_arm_joint_mpc_reference_preview|robot_arm_joint_mpc_reference_preview]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:106-106`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_planning.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:146-146`
<!-- vulcan:connections:end -->

## Limitations
Linear interpolation between quintic samples means the returned `dq` and `ddq` are not the exact derivatives of the returned `q`; consistency depends on `dt_s` being small. The end-effector position is interpolated in Cartesian space rather than recomputed via forward kinematics, so it drifts from `cloth_fk(q)` between samples. A degenerate last interval from `_reference_times` gives `t_hi == t_lo` and a NaN `α`. Each call allocates new vectors.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 145.
