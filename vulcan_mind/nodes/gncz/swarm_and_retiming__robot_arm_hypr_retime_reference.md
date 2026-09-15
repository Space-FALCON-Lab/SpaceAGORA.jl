---
id: gncz.swarm_and_retiming__robot_arm_hypr_retime_reference
label: _robot_arm_hypr_retime_reference
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_retime_reference
  lines:
  - 347
  - 430
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning namespace providing the base wrench models, the reaction
    scale rule, and the nominal reference time grid.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: t_ref
  type: Vector{Float64}
  units: s
  description: Monotonic reference time grid whose per-segment durations satisfy the
    configured joint rate, acceleration, and base wrench limits.
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

# _robot_arm_hypr_retime_reference

## Purpose
`_robot_arm_hypr_retime_reference` stretches the time grid of a joint-space reference until the motion respects joint rate limits, joint acceleration limits, and the force and torque the arm may impose on its mounting base. It is what turns a geometrically valid path into one a free-flying or compliantly mounted robot can actually execute.

## Theory & Math
Each segment carries a scale factor applied to the nominal step. A velocity limit requires $\max_j |\Delta q_j| / (v_{max}\Delta t) \le 1$, so the scale is raised to that ratio. An acceleration limit uses the central second difference $\ddot q \approx (q_{k+1} - 2q_k + q_{k-1})/\Delta t^2$, and because acceleration scales with the square of duration the node scale is $\sqrt{\ddot q / a_{max}}$, applied to both adjoining segments. Base wrench limits are enforced by iteration, scaling by the margin times the square root of the worst demand ratio until every node ratio drops to unity or nothing changes.

## Model & Assumptions
Retiming is skipped entirely when disabled or when the reference has fewer than two columns, returning the nominal grid. When no base wrench limit is finite, an additional reaction scale rule shapes the profile from the residual demand. All scales are clamped between the configured minimum and maximum.

## Design & Implementation
The wrench ratios come from either a rigid base model or a cloth compliance model, selected by configuration, and the loop mixes a global scale with per-node corrections so a single hot node does not have to be fixed one segment at a time. The iteration terminates on a change threshold as well as on satisfaction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning namespace providing the base wrench models, the reaction scale rule, and the nominal reference time grid. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `t_ref` | Vector{Float64} | s | — | Monotonic reference time grid whose per-segment durations satisfy the configured joint rate, acceleration, and base wrench limits. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:338-338`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:203-203`

**Downstream**

- `callees` → [[gnc.robot_arm_planning__reference_times|_reference_times]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:355-355`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios|_robot_arm_hypr_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:399-399`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_reaction_scale|_robot_arm_hypr_reaction_scale]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:387-387`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_reference_times_from_scales|_robot_arm_hypr_reference_times_from_scales]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:398-398`
<!-- vulcan:connections:end -->

## Limitations
The finite-difference rate estimates are as coarse as the reference spacing, and uniform stretching cannot exploit slack in one joint while another saturates. The wrench iteration is capped and may exit with a residual violation.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:1-439`.
