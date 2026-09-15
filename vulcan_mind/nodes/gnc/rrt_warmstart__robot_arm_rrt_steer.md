---
id: gnc.rrt_warmstart__robot_arm_rrt_steer
label: _robot_arm_rrt_steer
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_steer
  lines:
  - 26
  - 26
inputs:
- id: q_near
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `q_near`.
- id: q_target
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `q_target`.
- id: step_size_rad
  type: Real
  units: n/a
  required: true
  description: Positional argument `step_size_rad`.
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
  description: Return value of `_robot_arm_rrt_steer`. Returns `hypr_rrt_steer(Float64.(collect(q_near)),
    Float64.(collect(q_target)), step_size`.
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

# _robot_arm_rrt_steer

## Purpose
`_robot_arm_rrt_steer` computes the new tree node obtained by moving from a nearest node toward a target configuration by at most one RRT step, and reports whether the target was reached, the tree advanced, or the move was degenerate. It is the local planner of the warm start.

## Design & Implementation
Signature `(q_near::AbstractVector{<:Real}, q_target::AbstractVector{<:Real}, step_size_rad::Real)`. Both configurations are converted to `Vector{Float64}` and passed to `hypr_rrt_steer(q_near, q_target, step_size)`, which forms `direction = q_target - q_near` and `distance = norm(direction)`; if `distance <= 1e-10` it returns `(q_near, :trapped)`, if `distance <= max(step_size, 1e-9)` it returns `(q_target, :reached)`, and otherwise `(q_near + (step/distance) * direction, :advanced)`. The step is a straight line in joint space of length `step_size_rad`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_near` | AbstractVector{<:Real} | n/a | yes | Positional argument `q_near`. |
| in | `q_target` | AbstractVector{<:Real} | n/a | yes | Positional argument `q_target`. |
| in | `step_size_rad` | Real | n/a | yes | Positional argument `step_size_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_steer`. Returns `hypr_rrt_steer(Float64.(collect(q_near)), Float64.(collect(q_target)), step_size`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:69-69`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_steer|hypr_rrt_steer]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations
The step is Euclidean in joint space with no per-joint scaling, so a single large-range joint dominates the direction. The default trap tolerance `1e-10` and step floor `1e-9` are hard-coded in the shared helper and not exposed through `RobotArmHYPRConfig`. Joint limits are not enforced on the steered point; only the random samples respect limits, so a step toward a biased goal can leave the feasible box if the goal itself does. Each call allocates two converted vectors and a direction vector.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 26.
