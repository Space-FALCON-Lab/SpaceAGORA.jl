---
id: gncz.robot_arm_hypr_robot_arm_hypr
label: robot_arm_hypr
kind: struct
source:
  file: src/gnc/robotics/robot_arm_hypr.jl
  symbol: robot_arm_hypr
  lines:
  - 1
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning module namespace into which this file splices the
    five robot-arm HYPR planner source files.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: hypr_planner_scope
  type: Module
  units: n/a
  description: The configuration types, swarm and retiming routines, clearance statistics,
    warmstart search, and planner entry point brought into the enclosing module.
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

# robot_arm_hypr

## Purpose
`robot_arm_hypr.jl` is the aggregation point of the robot-arm hybrid planner. It contains no declarations of its own and exists to splice the five implementation files into the enclosing planning module in an order that satisfies their parse-time dependencies.

## Model & Assumptions
The order is load bearing rather than cosmetic. Configuration comes first because the sphere obstacle type, the planner configuration, and the result record are referenced in the signatures of everything that follows. Swarm scheduling and retiming come next, since the planner core calls the retiming routine and the control point helpers. Clearance follows because it depends on both the obstacle type and the segment distance helper defined in the retiming file. The warmstart search comes fourth, needing clearance to test segment safety, and the planner core comes last as the only consumer of all four.

## Design & Implementation
Each include joins the directory of the current file with the subdirectory name, so the split survives being loaded from any working directory. Because these are includes rather than a nested module, every symbol lands directly in the planning module namespace and the private names beginning with an underscore stay internal to it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning module namespace into which this file splices the five robot-arm HYPR planner source files. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `hypr_planner_scope` | Module | n/a | — | The configuration types, swarm and retiming routines, clearance statistics, warmstart search, and planner entry point brought into the enclosing module. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Splitting by include gives no encapsulation. A name collision between two of the five files would be silently resolved by load order, and reordering the includes to satisfy a new dependency can break an existing one without an obvious error. Because the file declares nothing, its behaviour can only be understood by reading the files it pulls in.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr.jl:1-5`.
