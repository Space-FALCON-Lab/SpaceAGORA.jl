---
id: gnc.trajectory_optimizers_rpo_trajectory_internal_points
label: rpo_trajectory_internal_points
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_internal_points
  lines:
  - 62
  - 62
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
  description: Return value of `rpo_trajectory_internal_points`. Returns `theta`.
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

# rpo_trajectory_internal_points

## Purpose
Builds the initial decision variable for the trajectory optimizers: a `3 x n_waypoints` matrix of internal waypoints placed at equal parameter spacing on the straight segment between `start` and `goal`, excluding the endpoints themselves.

## Design & Implementation
Allocates `theta = zeros(3, n_waypoints)` and for each column `i` sets `alpha = i / (n_waypoints + 1)` and `theta[:, i] = (1 - alpha) * start + alpha * goal` with broadcast arithmetic. `start` and `goal` may be any 3-element indexable; the result is always `Float64`. Called directly by `rpo_trajectory_internal_points_from_seed` as the fallback when no usable seed path is supplied.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `n_waypoints` | Integer | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_internal_points`. Returns `theta`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_trajectory_internal_points_from_seed|rpo_trajectory_internal_points_from_seed]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:74-74`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing checks that `start` and `goal` have length 3; a mismatch raises a `DimensionMismatch` from the broadcast rather than a descriptive error. A straight-line seed that passes through the station geometry starts the optimizer inside the obstacle potential, which is why callers prefer an RRT-derived seed where available.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 62.
