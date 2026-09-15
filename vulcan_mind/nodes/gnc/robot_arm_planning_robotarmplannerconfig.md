---
id: gnc.robot_arm_planning_robotarmplannerconfig
label: RobotArmPlannerConfig
kind: struct
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: RobotArmPlannerConfig
  lines:
  - 16
  - 16
inputs:
- id: dt_s
  type: Float64
  units: n/a
  required: false
  description: Field `dt_s` (default `0.1`).
- id: duration_s
  type: Float64
  units: n/a
  required: false
  description: Field `duration_s` (default `12.0`).
- id: ik_tol_m
  type: Float64
  units: n/a
  required: false
  description: Field `ik_tol_m` (default `1.0e-4`).
- id: ik_max_iters
  type: Int
  units: n/a
  required: false
  description: Field `ik_max_iters` (default `100`).
- id: ik_damping
  type: Float64
  units: n/a
  required: false
  description: Field `ik_damping` (default `1.0e-3`).
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
  type: RobotArmPlannerConfig
  units: n/a
  description: Constructed `RobotArmPlannerConfig` (keyword constructor via @kwdef).
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

# RobotArmPlannerConfig

## Purpose
Keyword-constructed configuration for the quintic joint-space robot-arm planner, grouping the reference time step, motion duration and the inverse-kinematics solver tolerances.

## Design & Implementation
`Base.@kwdef struct RobotArmPlannerConfig` with fields `dt_s::Float64 = 0.1` (reference sample spacing, s), `duration_s::Float64 = 12.0` (total motion time, s), `ik_tol_m::Float64 = 1.0e-4` (end-effector position tolerance for `cloth_ik`, m), `ik_max_iters::Int = 100` (damped least-squares iteration cap) and `ik_damping::Float64 = 1.0e-3` (Levenberg-Marquardt style damping). `plan_robot_arm_motion` passes `dt_s` and `duration_s` to `_reference_times` and forwards the three `ik_*` values to `cloth_ik`; the HYPR planner receives the whole config as `planner_config`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dt_s` | Float64 | n/a | no | Field `dt_s` (default `0.1`). |
| in | `duration_s` | Float64 | n/a | no | Field `duration_s` (default `12.0`). |
| in | `ik_tol_m` | Float64 | n/a | no | Field `ik_tol_m` (default `1.0e-4`). |
| in | `ik_max_iters` | Int | n/a | no | Field `ik_max_iters` (default `100`). |
| in | `ik_damping` | Float64 | n/a | no | Field `ik_damping` (default `1.0e-3`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmPlannerConfig | n/a | — | Constructed `RobotArmPlannerConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:179-179`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:77-77`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No field is validated at construction; non-positive `dt_s` or `duration_s` are rejected only later by `_reference_times`, and a zero `ik_max_iters` or negative damping is passed straight through to `cloth_ik`. There is no joint-rate or acceleration limit, so `duration_s` must be chosen by the user to keep the quintic peak rate `1.875*|Δq|/duration_s` feasible.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 16.
