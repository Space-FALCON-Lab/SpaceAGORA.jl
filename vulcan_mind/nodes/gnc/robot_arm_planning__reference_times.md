---
id: gnc.robot_arm_planning__reference_times
label: _reference_times
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: _reference_times
  lines:
  - 59
  - 59
inputs:
- id: dt_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_s`.
- id: duration_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `duration_s`.
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
  description: Return value of `_reference_times`. Returns `t`.
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

# _reference_times

## Purpose
Builds the monotone time grid `t_ref_s` on which a `RobotArmPlan` stores its joint references, ensuring the grid starts at 0 and ends exactly at `duration_s`.

## Design & Implementation
Takes `dt_s::Float64` and `duration_s::Float64` and throws `ArgumentError("dt_s must be positive.")` or `ArgumentError("duration_s must be positive.")` when either is non-positive. It materialises `collect(0.0:dt_s:duration_s)`; because a `StepRangeLen` may stop short of `duration_s` when the duration is not an integer multiple of `dt_s`, it appends `duration_s` if the last element is smaller (or the vector is empty). Returns a `Vector{Float64}` in seconds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `duration_s` | Float64 | n/a | yes | Positional argument `duration_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_reference_times`. Returns `t`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:336-336`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:201-201`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:108-108`
- [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:355-355`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:64-64`
<!-- vulcan:connections:end -->

## Limitations
The appended final sample can be arbitrarily close to the previous one (down to floating-point epsilon), producing a nearly degenerate last interval that `robot_arm_plan_sample` later divides by. No maximum sample count is enforced, so a tiny `dt_s` with a long duration allocates a very large grid. The grid is always uniform except possibly the last step.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 59.
