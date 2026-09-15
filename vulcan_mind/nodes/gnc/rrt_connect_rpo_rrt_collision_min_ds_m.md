---
id: gnc.rrt_connect_rpo_rrt_collision_min_ds_m
label: rpo_rrt_collision_min_ds_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_collision_min_ds_m
  lines:
  - 68
  - 68
inputs:
- id: sample_ds_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `sample_ds_m`.
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
- id: safe_distance_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `safe_distance_fraction`.
- id: guard_fraction
  type: Real
  units: n/a
  required: true
  description: Positional argument `guard_fraction`.
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
  description: Return value of `rpo_rrt_collision_min_ds_m`. Returns `max(min(min_ds,
    Float64(guard_fraction) * inflated_radius), 1.0e-9)`.
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

# rpo_rrt_collision_min_ds_m

## Purpose
Computes the minimum collision-check sample spacing (metres) for a segment so that the adaptive sampler cannot step over the safety shell of the station geometry.

## Design & Implementation
Starts from `min_ds = max(sample_ds_m, 1e-9)`. If `safe_distance_m > 0` it tightens to `min(min_ds, safe_distance_fraction * safe_distance_m)`. It then queries `rpo_inflated_obstacle_radius_m(geometry, safe_distance_m)`; if that radius is positive the spacing is further tightened to `min(min_ds, guard_fraction * inflated_radius)` and floored at `1e-9`. With the default fractions of 0.5 this guarantees at least two samples across any safety margin or obstacle radius.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sample_ds_m` | Real | n/a | yes | Positional argument `sample_ds_m`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `safe_distance_fraction` | Real | n/a | yes | Positional argument `safe_distance_fraction`. |
| in | `guard_fraction` | Real | n/a | yes | Positional argument `guard_fraction`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_collision_min_ds_m`. Returns `max(min(min_ds, Float64(guard_fraction) * inflated_radius), 1.0e-9)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:102-102`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:75-75`
- `callees` → [[gnc.path_sampling_rpo_inflated_obstacle_radius_m|rpo_inflated_obstacle_radius_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
The `1e-9` m floor can produce enormous sample counts when a caller passes a tiny `sample_ds_m`, making segment checks very slow rather than failing fast. The obstacle radius is a single scalar summary of the geometry, so thin protrusions smaller than the inflated radius are not specifically protected. Fractions greater than 1 defeat the guard but are not rejected.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 68.
