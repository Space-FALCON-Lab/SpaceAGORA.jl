---
id: gncz.robot_arm_planning_plan_robot_arm_motion
label: plan_robot_arm_motion
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: plan_robot_arm_motion
  lines:
  - 72
  - 142
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning module namespace exporting the planner configuration,
    the plan record, and both planner backends.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: plan
  type: RobotArmPlan
  units: s, rad, m
  description: Reference time grid with joint position, rate, and acceleration matrices,
    end-effector track, goal configuration, and final position error.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# plan_robot_arm_motion

## Purpose
`plan_robot_arm_motion` is the public planning entry point of the robot-arm module. It resolves a Cartesian target into a joint-space reference trajectory and selects between the direct quintic blend and the obstacle-aware hybrid planner, so callers depend on one function rather than on a backend choice.

## Theory & Math
The default backend applies a quintic scalar blend $a(s) = 10s^3 - 15s^4 + 6s^5$ over normalised time $s = t/T$, whose first and second derivatives vanish at both ends. Joint references are then $q(t) = q_0 + a(s)\,\Delta q$ with $\dot q = \dot a \Delta q$ and $\ddot q = \ddot a \Delta q$, giving a motion that starts and ends at rest with zero acceleration and needs no separate smoothing.

## Model & Assumptions
The goal configuration comes from damped least-squares inverse kinematics seeded at the start configuration, with tolerance, iteration cap, and damping taken from the planner configuration. The quintic path is a straight line in joint space, so it is only safe in an uncluttered workspace; anything else must use the hybrid backend, which is dispatched by symbol and returns the plan field of its richer result. Any other planner symbol raises an argument error naming the two supported values.

## Design & Implementation
The reference grid is built once from the step and duration, then joint position, rate, acceleration, and end-effector matrices are filled column by column, with forward kinematics evaluated at each column to record the Cartesian track. The final position error is the norm of the difference between the last end-effector column and the target, and it is stored in the plan so a caller can detect an inverse kinematics shortfall.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning module namespace exporting the planner configuration, the plan record, and both planner backends. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `plan` | RobotArmPlan | s, rad, m | — | Reference time grid with joint position, rate, and acceleration matrices, end-effector track, goal configuration, and final position error. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.robot_arm_planning__quintic_scalar|_quintic_scalar]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:118-118`
- `callees` → [[gnc.robot_arm_planning__quintic_scalar_ddot|_quintic_scalar_ddot]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:120-120`
- `callees` → [[gnc.robot_arm_planning__quintic_scalar_dot|_quintic_scalar_dot]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:119-119`
- `callees` → [[gnc.robot_arm_planning__reference_times|_reference_times]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:108-108`
- `callees` → [[gnc.robot_arm_planning_robotarmplan|RobotArmPlan]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:128-128`
- `callees` → [[gnc.robot_arm_planning_robotarmplannerconfig|RobotArmPlannerConfig]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:77-77`
- `callees` → [[gncz.config_robotarmhyprconfig|RobotArmHYPRConfig]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:90-90`
- `callees` → [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:84-84`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:124-124`
- `callees` → [[vehicle.robotics_cloth_ik|cloth_ik]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:99-99`
<!-- vulcan:connections:end -->

## Limitations
The quintic backend performs no collision checking and no rate limiting, so its plan can violate joint limits reached by interpolation. Duration is prescribed rather than derived, and the final error is reported rather than corrected.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl:1-174`.
