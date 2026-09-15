---
id: gnc.rrt_connect_rpo_rrt_star_plan_path
label: rpo_rrt_star_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_star_plan_path
  lines:
  - 488
  - 488
inputs:
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
- id: settings
  type: RPORRTStarSettings
  units: n/a
  required: false
  description: Keyword argument `settings` (default `RPORRTStarSettings()`).
- id: max_runtime_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `max_runtime_s` (default `Inf`).
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
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
  description: Return value of `rpo_rrt_star_plan_path`. Returns `(`.
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

# rpo_rrt_star_plan_path

## Purpose
Plans an RPO relative-motion path from `start_rtn` to `goal_rtn` around station geometry using single-tree RRT*, then shortcuts and post-refines the best goal-connected path and reports normalized cost components.

## Design & Implementation
Builds `local_cfg = rpo_pso_config(cfg; curve_type=:polyline)` and bounds from `rpo_pso_bounds`. If the direct segment is safe it returns immediately with a two-column path. Otherwise it iterates up to `settings.n_iters` (breaking when `max_runtime_s` elapses, measured with `time_ns()`), sampling the goal with probability `goal_sample_rate` or a random state, steering from the nearest node, and calling `rpo_rrt_star_add_node!`. When the new node is within `step_size_m` of the goal, cheaper than `best_goal_cost`, and the goal edge is safe, it records the candidate path and its cost in `history`/`cost_history`. After the loop it converts the settings into an `RPORRTConnectSettings` for `rpo_rrt_shortcut_path`, runs `rpo_post_refine_path`, and returns a NamedTuple with `path`, `raw_path`, `cost`, `components`, `config`, `iterations`, `path_found`, and other fields matching the PSO planner interface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `settings` | RPORRTStarSettings | n/a | no | Keyword argument `settings` (default `RPORRTStarSettings()`). |
| in | `max_runtime_s` | Real | n/a | no | Keyword argument `max_runtime_s` (default `Inf`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_star_plan_path`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:344-344`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:505-505`
- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:514-514`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:574-574`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:498-498`
- `callees` → [[gnc.rrt_connect_rpo_rrt_nearest_index|rpo_rrt_nearest_index]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:547-547`
- `callees` → [[gnc.rrt_connect_rpo_rrt_random_state|rpo_rrt_random_state]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:546-546`
- `callees` → [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:507-507`
- `callees` → [[gnc.rrt_connect_rpo_rrt_shortcut_path|rpo_rrt_shortcut_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:596-596`
- `callees` → [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:552-552`
- `callees` → [[gnc.rrt_connect_rpo_rrt_steer|rpo_rrt_steer]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:548-548`
- `callees` → [[gnc.rrt_connect_rpo_rrt_tree_path|rpo_rrt_tree_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:573-573`
- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:581-581`
- `callees` → [[gnc.rrt_connect_rporrtconnecttree|RPORRTConnectTree]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:533-533`
- `callees` → [[gnc.rrt_connect_rporrtstarsettings|RPORRTStarSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:494-494`
- `callees` → [[gncy.pso_helpers_rpo_pso_bounds|rpo_pso_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:501-501`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:597-597`
<!-- vulcan:connections:end -->

## Limitations
When no goal connection is found the function returns the unsafe `direct_path` with `path_found=false`; callers must check that flag. The runtime check is skipped on iteration 1, and each iteration performs O(N) nearest and near queries plus multiple segment checks, so wall time grows super-linearly. `cost_history` reflects raw candidate costs while the final `cost` is post-refined, so the two are not directly comparable.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 488.
