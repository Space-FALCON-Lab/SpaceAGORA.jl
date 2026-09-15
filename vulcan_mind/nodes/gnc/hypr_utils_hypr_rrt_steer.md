---
id: gnc.hypr_utils_hypr_rrt_steer
label: hypr_rrt_steer
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_steer
  lines:
  - 166
  - 166
inputs:
- id: q_near
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_near`.
- id: q_target
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_target`.
- id: step_size
  type: Any
  units: n/a
  required: true
  description: Positional argument `step_size`.
- id: trap_tol
  type: Real
  units: n/a
  required: false
  description: Keyword argument `trap_tol` (default `1.0e-10`).
- id: step_floor
  type: Real
  units: n/a
  required: false
  description: Keyword argument `step_floor` (default `1.0e-9`).
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
  description: Return value of `hypr_rrt_steer`. Returns `q_near + (step / distance)
    * direction, :advanced`.
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

# hypr_rrt_steer

## Purpose
Advances from a tree node toward a target by at most one step, reporting whether the move was trapped, reached the target, or merely advanced.

## Design & Implementation
Computes the direction and its norm. If the distance is within `trap_tol` it returns the start and `:trapped`; if within the step (floored at `step_floor`) it returns the target and `:reached`; otherwise it returns the point one step along the direction and `:advanced`. The status symbol drives the RRT-Connect extend loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_near` | Any | n/a | yes | Positional argument `q_near`. |
| in | `q_target` | Any | n/a | yes | Positional argument `q_target`. |
| in | `step_size` | Any | n/a | yes | Positional argument `step_size`. |
| in | `trap_tol` | Real | n/a | no | Keyword argument `trap_tol` (default `1.0e-10`). |
| in | `step_floor` | Real | n/a | no | Keyword argument `step_floor` (default `1.0e-9`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_steer`. Returns `q_near + (step / distance) * direction, :advanced`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_steer|rpo_rrt_steer]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:64-64`
- [[gnc.rrt_warmstart__robot_arm_rrt_steer|_robot_arm_rrt_steer]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:27-27`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:169-169`
<!-- vulcan:connections:end -->

## Limitations
The step is purely geometric with no collision check, which the caller performs afterwards; the `step_floor` guard against a zero step means a caller asking for zero motion silently gets 1e-9 instead.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 166.
