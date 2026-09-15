---
id: gnc.path_sampling_rpo_inflated_obstacle_radius_m
label: rpo_inflated_obstacle_radius_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_inflated_obstacle_radius_m
  lines:
  - 66
  - 66
inputs:
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `safe_distance_m`.
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
  description: Return value of `rpo_inflated_obstacle_radius_m`. Returns `geometry.station.keepout_radius_m
    +`.
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

# rpo_inflated_obstacle_radius_m

## Purpose
Gives the effective radius of the station obstacle once the chaser's own size and the requested safety margin are folded in, so clearance tests can treat the chaser as a point.

## Design & Implementation
Adds three terms: `geometry.station.keepout_radius_m`, the largest chaser body half-extent from `geometry.chaser.half_extents_body`, and `safe_distance_m` clamped at zero. Using the maximum half-extent bounds the chaser by a sphere, which is conservative for any attitude.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_inflated_obstacle_radius_m`. Returns `geometry.station.keepout_radius_m +`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_adaptive_sampling_min_ds_m|rpo_adaptive_sampling_min_ds_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:81-81`
- [[gnc.rrt_connect_rpo_rrt_collision_min_ds_m|rpo_rrt_collision_min_ds_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:79-79`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
The bounding-sphere inflation over-approximates a slender chaser, so paths that would clear the station edge-on are rejected; a negative safety margin is silently clamped to zero rather than rejected.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 66.
