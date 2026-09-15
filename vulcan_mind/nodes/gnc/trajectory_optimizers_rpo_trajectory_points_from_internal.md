---
id: gnc.trajectory_optimizers_rpo_trajectory_points_from_internal
label: rpo_trajectory_points_from_internal
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_points_from_internal
  lines:
  - 101
  - 101
inputs:
- id: theta
  type: Any
  units: n/a
  required: true
  description: Positional argument `theta`.
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
  description: Return value of `rpo_trajectory_points_from_internal`. Returns `points`.
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

# rpo_trajectory_points_from_internal

## Purpose
Assembles the full waypoint trajectory matrix consumed by the cost functions by prepending `start` and appending `goal` to the `3 x n` internal waypoint matrix `theta`.

## Design & Implementation
Reads `n_waypoints = size(theta, 2)`, allocates `points = zeros(3, n_waypoints + 2)`, writes `start` into column 1 and `goal` into the last column, and copies `theta` into columns `2:end-1` only when `n_waypoints > 0`. It is called many times per iteration inside `rpo_chomp_numeric_gradient` (twice per perturbed coordinate) and once per rollout in STOMP.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Any | n/a | yes | Positional argument `theta`. |
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_points_from_internal`. Returns `points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_chomp_numeric_gradient|rpo_chomp_numeric_gradient]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:214-214`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:390-390`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:256-256`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every call allocates a fresh matrix, so the finite-difference gradient allocates `6 * n_waypoints` matrices per iteration. The result is always `Float64` even if `theta` were a different element type. No shape check on `start`/`goal`; a wrong length fails at the broadcast assignment.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 101.
