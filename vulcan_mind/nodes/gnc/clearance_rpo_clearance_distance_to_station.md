---
id: gnc.clearance_rpo_clearance_distance_to_station
label: rpo_clearance_distance_to_station
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/clearance.jl
  symbol: rpo_clearance_distance_to_station
  lines:
  - 9
  - 9
inputs:
- id: p_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `p_body`.
- id: geometry
  type: RPOReferenceGeometry
  units: n/a
  required: true
  description: Positional argument `geometry`.
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
  description: Return value of `rpo_clearance_distance_to_station`. Returns `distance
    - geometry.station.keepout_radius_m - maximum(geometry.chaser.half_exte`.
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

# rpo_clearance_distance_to_station

## Purpose
Reports how much clear space a chaser body point has from the target station, accounting for both the station keep-out sphere and the chaser's own size.

## Design & Implementation
Takes the squared nearest-station distance from `nearest_station_distance_sq`, square-roots it once, then subtracts the station `keepout_radius_m` and the largest of the chaser's body half-extents. Subtracting the maximum half-extent treats the chaser conservatively as a bounding sphere, so a positive result is safe for any attitude. `@inline` because it runs per path sample.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `geometry` | RPOReferenceGeometry | n/a | yes | Positional argument `geometry`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_clearance_distance_to_station`. Returns `distance - geometry.station.keepout_radius_m - maximum(geometry.chaser.half_exte`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.clearance_rpo_path_clearance_stats|rpo_path_clearance_stats]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:23-23`
- [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:48-48`
- [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:134-134`
- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:516-516`
- [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:131-131`
- [[gnc.trajectory_optimizers_rpo_soft_obstacle_cost_from_samples|rpo_soft_obstacle_cost_from_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:175-175`
- [[gnc.trajectory_optimizers_rpo_stomp_waypoint_state_cost|rpo_stomp_waypoint_state_cost]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:348-348`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:218-218`

**Downstream**

- `callees` → [[gnc.mesh_distance_nearest_station_distance_sq|nearest_station_distance_sq]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
The bounding-sphere treatment is conservative and will report a violation for attitudes where the real geometry would clear; it cannot exploit a slender chaser's actual orientation.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/clearance.jl` line 9.
