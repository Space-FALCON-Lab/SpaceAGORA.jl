---
id: gnc.rrt_connect_rporrtconnectsettings
label: RPORRTConnectSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: RPORRTConnectSettings
  lines:
  - 2
  - 2
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
- id: connect_max_steps
  type: Int
  units: n/a
  required: false
  description: Field `connect_max_steps` (default `10_000`).
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
  type: RPORRTConnectSettings
  units: n/a
  description: Constructed `RPORRTConnectSettings` (keyword constructor via @kwdef).
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

# RPORRTConnectSettings

## Purpose
Keyword-constructible settings struct for the RPO bidirectional RRT-Connect planner, holding iteration budgets, step size, goal bias, adaptive collision-sampling parameters, and shortcut smoothing count.

## Design & Implementation
Defined with `Base.@kwdef`. Fields and defaults: `n_iters::Int=1000`, `step_size_m::Float64=0.75`, `goal_sample_rate::Float64=0.05`, `collision_sample_ds_m::Float64=0.10`, `adaptive_collision_sampling_enable::Bool=true`, `collision_max_sample_ds_m::Float64=0.50`, `collision_far_clearance_m::Float64=1.0`, `collision_sampling_power::Float64=1.0`, `collision_safe_distance_fraction::Float64=0.5`, `collision_obstacle_guard_fraction::Float64=0.5`, `connect_max_steps::Int=10_000`, `shortcut_iters::Int=80`. All distances are metres in the RTN relative frame. The struct is immutable and is read by `rpo_rrt_segment_is_safe`, `rpo_rrt_extend!`, `rpo_rrt_connect!`, and `rpo_rrt_shortcut_path`.

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
| in | `connect_max_steps` | Int | n/a | no | Field `connect_max_steps` (default `10_000`). |
| in | `shortcut_iters` | Int | n/a | no | Field `shortcut_iters` (default `80`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPORRTConnectSettings | n/a | — | Constructed `RPORRTConnectSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_rrt_connect_settings|_rpo_comparison_rrt_connect_settings]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:197-197`
- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:33-33`
- [[gnc.pso_path_planning_rpo_pso_rrt_warmstart_path|rpo_pso_rrt_warmstart_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:126-126`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:427-427`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:581-581`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:314-314`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No constructor validation: a non-positive `step_size_m` or `collision_sample_ds_m` is accepted and only partially guarded by `max(..., 1e-9)` inside consumers. `goal_sample_rate` outside [0,1] is clamped at use sites rather than rejected. The field set duplicates `RPORRTStarSettings` except for `connect_max_steps` versus `neighbor_radius_m`, and `rpo_rrt_star_plan_path` must copy fields one by one to build a shortcut settings object.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 2.
