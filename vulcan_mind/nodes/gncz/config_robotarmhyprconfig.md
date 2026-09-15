---
id: gncz.config_robotarmhyprconfig
label: RobotArmHYPRConfig
kind: struct
source:
  file: src/gnc/robotics/robot_arm_hypr/config.jl
  symbol: RobotArmHYPRConfig
  lines:
  - 11
  - 76
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning namespace where the planner configuration, obstacle
    type, and result record are declared and exported.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: config
  type: RobotArmHYPRConfig
  units: mixed
  description: Immutable configuration governing swarm size, cost weights, scheduling,
    culling, warmstart, refinement, early stopping, and retiming.
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

# RobotArmHYPRConfig

## Purpose
`RobotArmHYPRConfig` is the single tuning surface of the robot-arm hybrid planner. Every stage of the planner reads its behaviour from this record, which keeps the algorithm free of embedded constants and makes a planning run reproducible from one value.

## Model & Assumptions
The fields group into seven concerns. Problem size fixes the number of interior waypoints, particles, iterations, path samples, and the curve type. Cost weighting balances path length, smoothness, and obstacle penalty against a safe distance margin. Swarm dynamics set inertia and the two acceleration coefficients, the initial spread, and a velocity cap expressed as a fraction of the joint range. A schedule anneals inertia and both coefficients toward end fractions over a transition fraction of the run, clamped to configured bounds. Culling replaces the worst particles after a start iteration with noise-perturbed copies. Warmstart controls an optional rapidly exploring random tree search with its own iteration, step, goal bias, collision sampling, connection, shortcut, and runtime limits. Refinement, early stopping, and retiming close the list, the last covering joint velocity and acceleration caps, reaction scaling, base force and torque limits, and a cloth physics model of the mounting compliance.

## Design & Implementation
The struct is keyword-constructed and immutable, so every field has a working default and a caller overrides only what it cares about. A separate validation routine in the same file checks the ranges before the planner allocates anything.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning namespace where the planner configuration, obstacle type, and result record are declared and exported. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config` | RobotArmHYPRConfig | mixed | — | Immutable configuration governing swarm size, cost weights, scheduling, culling, warmstart, refinement, early stopping, and retiming. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:180-180`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:90-90`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The flat design means several fields are inert unless the switch that governs them is enabled, and nothing in the type expresses that coupling. Defaults were tuned for one arm scale and are not dimensionless, so they do not transfer unchanged to a much larger or smaller manipulator.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/config.jl:1-158`.
