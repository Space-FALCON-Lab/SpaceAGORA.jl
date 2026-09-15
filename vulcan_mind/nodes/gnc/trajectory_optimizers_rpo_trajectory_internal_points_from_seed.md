---
id: gnc.trajectory_optimizers_rpo_trajectory_internal_points_from_seed
label: rpo_trajectory_internal_points_from_seed
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_internal_points_from_seed
  lines:
  - 72
  - 72
inputs:
- id: seed_path
  type: Any
  units: n/a
  required: true
  description: Positional argument `seed_path`.
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
  description: Return value of `rpo_trajectory_internal_points_from_seed`. Returns
    `theta`.
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

# rpo_trajectory_internal_points_from_seed

## Purpose
Resamples an arbitrary seed path (for example an RRT-Connect result) into exactly `n_waypoints` internal optimizer waypoints spaced uniformly in arc length, so CHOMP and STOMP can be warm-started from a collision-free route.

## Design & Implementation
Early exits: `n_waypoints <= 0` returns `zeros(3, 0)`; `seed_path === nothing`, a seed without 3 rows, fewer than 2 columns, or total length `<= 1.0e-12` all fall back to `rpo_trajectory_internal_points`. Otherwise the seed is copied to `Matrix{Float64}`, its first and last columns are overwritten with `start` and `goal`, and cumulative Euclidean arc length is accumulated. For each `i` the target arc length is `s = total_len * i / (n_waypoints + 1)`; `searchsortedfirst` locates the segment (clamped to `[2, end]`), and the point is linearly interpolated within it with `alpha` set to 0 when the segment span is degenerate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `seed_path` | Any | n/a | yes | Positional argument `seed_path`. |
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `n_waypoints` | Integer | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_internal_points_from_seed`. Returns `theta`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:372-372`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:240-240`

**Downstream**

- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_internal_points|rpo_trajectory_internal_points]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:74-74`
<!-- vulcan:connections:end -->

## Limitations
Overwriting the seed endpoints silently distorts a seed whose ends differ from `start`/`goal` instead of raising. `searchsortedfirst` assumes `cumulative` is non-decreasing, which holds because norms are non-negative. Converting with `Matrix{Float64}(seed_path)` throws if the seed contains non-numeric entries.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 72.
