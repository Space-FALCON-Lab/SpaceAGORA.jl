---
id: gncz.clearance_robot_arm_clearance_stats_from_samples
label: robot_arm_clearance_stats_from_samples
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/clearance.jl
  symbol: robot_arm_clearance_stats_from_samples
  lines:
  - 2
  - 44
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning namespace providing the cloth arm model, the forward
    kinematics routine, and the spherical obstacle type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: clearance_stats
  type: NamedTuple
  units: m
  description: Minimum clearance, violation count, violated fraction of checks, and
    the squared-deficit penalty summed over all link and obstacle pairs.
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

# robot_arm_clearance_stats_from_samples

## Purpose
`robot_arm_clearance_stats_from_samples` scores a sampled joint-space path against the spherical obstacles in the workspace. It produces both the diagnostic minimum clearance a reviewer wants to see and the smooth penalty term that the particle swarm optimiser needs in order to be pushed away from collisions.

## Theory & Math
Each link is treated as a capsule between its joint origin and its tip. For a link with radius $r_l$ and an obstacle of radius $r_o$ centred at $c$, the clearance is $c_{lo} = \mathrm{dist}(c, \overline{ab}) - r_o - r_l$, where the segment distance is the standard clamped projection of the centre onto the link axis. Violations are counted whenever the clearance falls below the safe distance, and the penalty accumulates the squared deficit $\sum \max(0, s - c_{lo})^2$, which is zero and continuously differentiable at the boundary.

## Model & Assumptions
Forward kinematics is recomputed from the cloth arm model and its base pose at every sampled configuration, so the statistics reflect the actual link placement rather than an interpolation of endpoint poses. With no obstacles the routine returns immediately with infinite clearance and zero penalty, which lets the planner run unobstructed cases without a special code path.

## Design & Implementation
The triple loop over samples, links, and obstacles runs inbounds and counts every pair it evaluates, so the violated fraction is normalised by the number of checks rather than by the number of samples.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RobotArmPlanning namespace providing the cloth arm model, the forward kinematics routine, and the spherical obstacle type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `clearance_stats` | NamedTuple | m | — | Minimum clearance, violation count, violated fraction of checks, and the squared-deficit penalty summed over all link and obstacle pairs. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:15-15`
- [[gnc.rrt_warmstart__robot_arm_rrt_path_score|_robot_arm_rrt_path_score]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:178-178`
- [[gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe|_robot_arm_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:55-55`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:15-15`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_segment_distance|_robot_arm_segment_distance]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:27-27`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
Only spheres are supported as obstacles, and the capsule model ignores link cross-section shape and any payload carried by the end effector. Collision is judged only at sampled configurations, so a thin obstacle can be tunnelled between two consecutive samples.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/clearance.jl:1-44`.
