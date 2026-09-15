---
id: gnc.rrt_connect_rpo_rrt_segment_is_safe
label: rpo_rrt_segment_is_safe
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_segment_is_safe
  lines:
  - 85
  - 85
inputs:
- id: q_from
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_from`.
- id: q_to
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_to`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
- id: sample_ds_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `sample_ds_m`.
- id: adaptive_enable
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `adaptive_enable` (default `true`).
- id: max_sample_ds_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `max_sample_ds_m` (default `0.50`).
- id: far_clearance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `far_clearance_m` (default `1.0`).
- id: sampling_power
  type: Real
  units: n/a
  required: false
  description: Keyword argument `sampling_power` (default `1.0`).
- id: safe_distance_fraction
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_fraction` (default `0.5`).
- id: obstacle_guard_fraction
  type: Real
  units: n/a
  required: false
  description: Keyword argument `obstacle_guard_fraction` (default `0.5`).
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
  type: Bool
  units: n/a
  description: Return value of `rpo_rrt_segment_is_safe`. Returns `true`.
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

# rpo_rrt_segment_is_safe

## Purpose
Tests whether the straight segment from `q_from` to `q_to` keeps at least `safe_distance_m` clearance from the station geometry at every sample point, either with adaptive or uniform sampling.

## Design & Implementation
Two methods. The keyword method converts endpoints to `SVector{3,Float64}` and, when `adaptive_enable` is true, computes `min_ds` via `rpo_rrt_collision_min_ds_m` and calls `rpo_adaptive_segment_samples` with `max_ds_m = max(max_sample_ds_m, min_ds)`, `far_clearance_m`, and `power`. Otherwise it builds a uniform `3 x (n+1)` matrix with `n = max(1, ceil(dist / max(sample_ds_m, 1e-6)))` by linear interpolation. Each column is checked with `rpo_clearance_distance_to_station(q, geometry) + 1e-9 >= safe`; the first violation returns `false`. The second method `(q_from, q_to, geometry, settings; safe_distance_m)` unpacks the `collision_*` fields from either settings struct.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_from` | Any | n/a | yes | Positional argument `q_from`. |
| in | `q_to` | Any | n/a | yes | Positional argument `q_to`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `sample_ds_m` | Real | n/a | yes | Keyword argument `sample_ds_m`. |
| in | `adaptive_enable` | Bool | n/a | no | Keyword argument `adaptive_enable` (default `true`). |
| in | `max_sample_ds_m` | Real | n/a | no | Keyword argument `max_sample_ds_m` (default `0.50`). |
| in | `far_clearance_m` | Real | n/a | no | Keyword argument `far_clearance_m` (default `1.0`). |
| in | `sampling_power` | Real | n/a | no | Keyword argument `sampling_power` (default `1.0`). |
| in | `safe_distance_fraction` | Real | n/a | no | Keyword argument `safe_distance_fraction` (default `0.5`). |
| in | `obstacle_guard_fraction` | Real | n/a | no | Keyword argument `obstacle_guard_fraction` (default `0.5`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `rpo_rrt_segment_is_safe`. Returns `true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:164-164`
- [[gnc.rrt_connect_rpo_rrt_shortcut_path|rpo_rrt_shortcut_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:227-227`
- [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:260-260`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:507-507`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:327-327`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:100-100`
- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:131-131`
- `callees` → [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:109-109`
- `callees` → [[gnc.rrt_connect_rpo_rrt_collision_min_ds_m|rpo_rrt_collision_min_ds_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:102-102`
<!-- vulcan:connections:end -->

## Limitations
Safety is only verified at discrete samples, so a segment can graze an obstacle between samples if the spacing exceeds the obstacle feature size; the `1e-9` tolerance also admits points marginally inside the margin. The uniform path allocates a dense matrix per call, and the adaptive path allocates similarly, which dominates planner runtime. Endpoints are re-checked on every call even when the tree already validated them.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 85.
