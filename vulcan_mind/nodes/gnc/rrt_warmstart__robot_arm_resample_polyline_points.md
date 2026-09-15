---
id: gnc.rrt_warmstart__robot_arm_resample_polyline_points
label: _robot_arm_resample_polyline_points
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_resample_polyline_points
  lines:
  - 134
  - 134
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
- id: n_points
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_points`.
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
  description: Return value of `_robot_arm_resample_polyline_points`. Returns `result`.
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

# _robot_arm_resample_polyline_points

## Purpose
`_robot_arm_resample_polyline_points` re-parameterises a joint-space polyline to a fixed number of waypoints spaced uniformly by arc length, so an RRT path with an arbitrary number of nodes can seed a HYPR optimiser that expects a fixed decision-variable count.

## Design & Implementation
Signature `(path, n_points::Int)`; throws `ArgumentError` if `n_points < 2`. The input is copied to `pts = Matrix{Float64}(path)`; if it already has `n_points` columns a copy is returned. The result matrix gets the first and last columns copied verbatim; a single-column input is padded with zeros beyond the endpoints and returned. Otherwise it builds `cumulative` arc lengths via `norm(pts[:, j] - pts[:, j-1])`, and if `total <= 1e-12` fills all interior columns with the first point. For each interior index `j` it computes `target_s = total (j-1)/(n_points-1)`, advances a segment pointer `seg` while `cumulative[seg+1] < target_s`, and interpolates linearly within that segment with `α = clamp((target_s - cumulative[seg]) / max(seg_len, 1e-12), 0, 1)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `n_points` | Int | n/a | yes | Positional argument `n_points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_resample_polyline_points`. Returns `result`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:242-242`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For a single-column input the interior columns are left as zeros rather than copies of the point, an inconsistency with the degenerate-length branch. Uniform arc-length spacing in joint space does not equal uniform end-effector spacing. The segment pointer only moves forward, which is correct for monotone cumulative lengths but silently mis-assigns if `pts` contains `NaN`. Interpolation is linear, so curvature of the original polyline at waypoints is not preserved.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 134.
