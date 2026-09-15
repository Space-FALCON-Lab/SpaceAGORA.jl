---
id: gnc.trajectory_optimizers_rpo_clamp_internal_waypoints
label: rpo_clamp_internal_waypoints
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_clamp_internal_waypoints
  lines:
  - 119
  - 119
inputs:
- id: theta
  type: Any
  units: n/a
  required: true
  description: Positional argument `theta`.
- id: lo
  type: Any
  units: n/a
  required: true
  description: Positional argument `lo`.
- id: hi
  type: Any
  units: n/a
  required: true
  description: Positional argument `hi`.
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
  description: 'Return value of `rpo_clamp_internal_waypoints`. Returns `size(theta,
    2) == 0 ? copy(theta) : clamp.(theta, lo, hi)`.'
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

# rpo_clamp_internal_waypoints

## Purpose
One-line projection of the internal waypoint matrix `theta` onto the box `[lo, hi]` after each optimizer update, keeping CHOMP descent steps and STOMP rollouts within the search corridor.

## Design & Implementation
Defined as `rpo_clamp_internal_waypoints(theta, lo, hi) = size(theta, 2) == 0 ? copy(theta) : clamp.(theta, lo, hi)`. The zero-column branch avoids broadcasting a `3 x 0` matrix against bound matrices and still returns a fresh copy so callers can safely mutate the result. `lo` and `hi` must be the same `3 x n` shape as `theta` (as produced by `rpo_trajectory_search_bounds`).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Any | n/a | yes | Positional argument `theta`. |
| in | `lo` | Any | n/a | yes | Positional argument `lo`. |
| in | `hi` | Any | n/a | yes | Positional argument `hi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_clamp_internal_waypoints`. Returns `size(theta, 2) == 0 ? copy(theta) : clamp.(theta, lo, hi)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:418-418`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:283-283`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Elementwise clamping is a projection onto a box, not onto the feasible obstacle-free set, so a clamped point can still lie inside the station geometry. Always allocates a new matrix; there is no in-place variant. Mismatched bound shapes raise `DimensionMismatch` from the broadcast.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 119.
