---
id: gnc.rrt_connect_rporrtstarsettings
label: RPORRTStarSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: RPORRTStarSettings
  lines:
  - 18
  - 18
inputs:
- id: n_iters
  type: Int
  units: n/a
  required: false
  description: Field `n_iters` (default `1000`).
- id: step_size_m
  type: Float64
  units: n/a
  required: false
  description: Field `step_size_m` (default `0.75`).
- id: goal_sample_rate
  type: Float64
  units: n/a
  required: false
  description: Field `goal_sample_rate` (default `0.05`).
- id: collision_sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `collision_sample_ds_m` (default `0.10`).
- id: adaptive_collision_sampling_enable
  type: Bool
  units: n/a
  required: false
  description: Field `adaptive_collision_sampling_enable` (default `true`).
- id: collision_max_sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `collision_max_sample_ds_m` (default `0.50`).
- id: collision_far_clearance_m
  type: Float64
  units: n/a
  required: false
  description: Field `collision_far_clearance_m` (default `1.0`).
- id: collision_sampling_power
  type: Float64
  units: n/a
  required: false
  description: Field `collision_sampling_power` (default `1.0`).
- id: collision_safe_distance_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `collision_safe_distance_fraction` (default `0.5`).
- id: collision_obstacle_guard_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `collision_obstacle_guard_fraction` (default `0.5`).
- id: neighbor_radius_m
  type: Float64
  units: n/a
  required: false
  description: Field `neighbor_radius_m` (default `2.0`).
- id: shortcut_iters
  type: Int
  units: n/a
  required: false
  description: Field `shortcut_iters` (default `80`).
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
  type: RPORRTStarSettings
  units: n/a
  description: Constructed `RPORRTStarSettings` (keyword constructor via @kwdef).
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

# RPORRTStarSettings

## Purpose
Keyword-constructible settings struct for the RPO RRT* comparison planner, mirroring the RRT-Connect settings but replacing the connect step budget with a rewiring `neighbor_radius_m`.

## Design & Implementation
Defined with `Base.@kwdef`; defaults are `n_iters=1000`, `step_size_m=0.75`, `goal_sample_rate=0.05`, `collision_sample_ds_m=0.10`, `adaptive_collision_sampling_enable=true`, `collision_max_sample_ds_m=0.50`, `collision_far_clearance_m=1.0`, `collision_sampling_power=1.0`, `collision_safe_distance_fraction=0.5`, `collision_obstacle_guard_fraction=0.5`, `neighbor_radius_m::Float64=2.0`, `shortcut_iters=80`. `neighbor_radius_m` bounds the near-set query in `rpo_rrt_star_add_node!`. The second `rpo_rrt_segment_is_safe` method accepts either settings type because it only reads the shared `collision_*` fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_iters` | Int | n/a | no | Field `n_iters` (default `1000`). |
| in | `step_size_m` | Float64 | n/a | no | Field `step_size_m` (default `0.75`). |
| in | `goal_sample_rate` | Float64 | n/a | no | Field `goal_sample_rate` (default `0.05`). |
| in | `collision_sample_ds_m` | Float64 | n/a | no | Field `collision_sample_ds_m` (default `0.10`). |
| in | `adaptive_collision_sampling_enable` | Bool | n/a | no | Field `adaptive_collision_sampling_enable` (default `true`). |
| in | `collision_max_sample_ds_m` | Float64 | n/a | no | Field `collision_max_sample_ds_m` (default `0.50`). |
| in | `collision_far_clearance_m` | Float64 | n/a | no | Field `collision_far_clearance_m` (default `1.0`). |
| in | `collision_sampling_power` | Float64 | n/a | no | Field `collision_sampling_power` (default `1.0`). |
| in | `collision_safe_distance_fraction` | Float64 | n/a | no | Field `collision_safe_distance_fraction` (default `0.5`). |
| in | `collision_obstacle_guard_fraction` | Float64 | n/a | no | Field `collision_obstacle_guard_fraction` (default `0.5`). |
| in | `neighbor_radius_m` | Float64 | n/a | no | Field `neighbor_radius_m` (default `2.0`). |
| in | `shortcut_iters` | Int | n/a | no | Field `shortcut_iters` (default `80`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPORRTStarSettings | n/a | — | Constructed `RPORRTStarSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_rrt_star_settings|_rpo_comparison_rrt_star_settings]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:216-216`
- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:34-34`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:494-494`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `neighbor_radius_m` is a fixed constant rather than the shrinking radius of asymptotically optimal RRT*, so optimality guarantees do not hold. There is no validation that `neighbor_radius_m >= step_size_m`; a smaller radius makes the near set frequently empty and degrades to plain RRT. Duplicated fields with `RPORRTConnectSettings` must be kept in sync manually.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 18.
