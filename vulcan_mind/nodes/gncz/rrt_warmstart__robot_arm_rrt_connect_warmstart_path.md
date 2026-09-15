---
id: gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path
label: _robot_arm_rrt_connect_warmstart_path
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_connect_warmstart_path
  lines:
  - 212
  - 290
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning namespace providing the tree record, the extend and
    connect primitives, and the segment safety test.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: warmstart
  type: Tuple
  units: mixed
  description: Shortcut warmstart path or nothing, paired with diagnostics recording
    attempt, success, iterations, cost, raw cost, and point count.
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

# _robot_arm_rrt_connect_warmstart_path

## Purpose
`_robot_arm_rrt_connect_warmstart_path` produces a feasible seed for the particle swarm by running a bidirectional rapidly exploring random tree between the start and goal joint configurations. Seeding the swarm inside the free space is what lets the planner solve cluttered problems that random initialisation would fail.

## Theory & Math
Two trees are grown, one from the start and one from the goal. On each iteration the parity of the counter selects which tree extends toward a sample, drawn either uniformly within the joint limits or, with the configured goal bias probability, at the opposite root. After the chosen tree extends by one step, the other tree attempts a greedy connect toward the new node, and a reached status means the two trees meet and the joined branch forms a path.

## Model & Assumptions
The routine returns immediately with empty diagnostics when warmstart is disabled. Before growing anything it tests the straight segment from start to goal, and if that is safe it returns the direct two-point path, which makes uncluttered problems free. Segment safety is judged by sampling the segment at the configured angular spacing and calling the capsule clearance test, so safety is discrete rather than certified.

## Design & Implementation
A wall-clock budget is checked at the top of every iteration after the first, so an unsolvable problem degrades to the direct path instead of hanging. A successful path is passed through a randomised shortcut pass, and both the shortcut and raw path are scored so the diagnostics show how much the shortcut gained. When no path is found the routine returns nothing while still reporting the attempt.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning namespace providing the tree record, the extend and connect primitives, and the segment safety test. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `warmstart` | Tuple | mixed | — | Shortcut warmstart path or nothing, paired with diagnostics recording attempt, success, iterations, cost, raw cost, and point count. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:239-239`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:227-227`
- `callees` → [[gnc.rrt_warmstart__robot_arm_empty_rrt_warmstart_diagnostics|_robot_arm_empty_rrt_warmstart_diagnostics]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:221-221`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_connect_bang|_robot_arm_rrt_connect!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:265-265`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:261-261`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_join_paths|_robot_arm_rrt_join_paths]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:269-269`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_path_score|_robot_arm_rrt_path_score]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:230-230`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_random_state|_robot_arm_rrt_random_state]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:257-257`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe|_robot_arm_rrt_segment_is_safe]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:229-229`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_shortcut_path|_robot_arm_rrt_shortcut_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:277-277`
- `callees` → [[gnc.rrt_warmstart_robotarmrrtconnecttree|RobotArmRRTConnectTree]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:244-244`
<!-- vulcan:connections:end -->

## Limitations
A trapped extension only skips the iteration, so a narrow passage can consume the whole budget. Sampled safety checking can miss a thin obstacle, and the trees are discarded after one call rather than being reused across planning invocations.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:1-290`.
