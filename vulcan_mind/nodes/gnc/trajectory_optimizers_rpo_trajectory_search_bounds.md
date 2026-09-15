---
id: gnc.trajectory_optimizers_rpo_trajectory_search_bounds
label: rpo_trajectory_search_bounds
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_search_bounds
  lines:
  - 111
  - 111
inputs:
- id: start
  type: Any
  units: n/a
  required: true
  description: Positional argument `start`.
- id: goal
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal`.
- id: search_margin
  type: Any
  units: n/a
  required: true
  description: Positional argument `search_margin`.
- id: n_waypoints
  type: Integer
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
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
  description: Return value of `rpo_trajectory_search_bounds`. Returns `repeat(lo,
    1, n_waypoints), repeat(hi, 1, n_waypoints)`.
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

# rpo_trajectory_search_bounds

## Purpose
Produces per-waypoint lower and upper bound matrices defining the axis-aligned box the optimizers may move internal waypoints within: the bounding box of `start` and `goal` expanded by `search_margin` on every side.

## Design & Implementation
Clamps `margin = max(0.0, Float64(search_margin))`, computes `lo = min.(start, goal) .- margin` and `hi = max.(start, goal) .+ margin`, and returns `(repeat(lo, 1, n_waypoints), repeat(hi, 1, n_waypoints))` so each column of the `3 x n_waypoints` bound matrices is identical. The planners subsequently widen these with `lo .= min.(lo, theta)` and `hi .= max.(hi, theta)` so a seed path outside the box is never clipped.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `search_margin` | Any | n/a | yes | Positional argument `search_margin`. |
| in | `n_waypoints` | Integer | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_search_bounds`. Returns `repeat(lo, 1, n_waypoints), repeat(hi, 1, n_waypoints)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:373-373`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:241-241`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:112-112`
<!-- vulcan:connections:end -->

## Limitations
The box is axis-aligned in the RTN frame; a route that must loop around a station well outside the start-goal corridor is only reachable if `search_margin` is large. A negative `search_margin` is silently treated as zero, which collapses the box to the corridor itself when start and goal share a coordinate.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 111.
