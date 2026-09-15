---
id: gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path
label: rpo_rrt_connect_bezier_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_connect_bezier_plan_path
  lines:
  - 421
  - 421
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
  type: RPORRTConnectSettings
  units: n/a
  required: false
  description: Keyword argument `settings` (default `RPORRTConnectSettings()`).
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
  description: Return value of `rpo_rrt_connect_bezier_plan_path`. Returns `(`.
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

# rpo_rrt_connect_bezier_plan_path

## Purpose
Runs the RRT-Connect polyline planner and then refits the result as a Bezier curve with fixed endpoints, selecting the control-point count that minimizes obstacle penalty and cost.

## Design & Implementation
Calls `rpo_rrt_connect_plan_path` with the same keywords, then builds `bezier_cfg = rpo_pso_config(cfg; curve_type=:bezier)` and resamples the polyline with `rpo_sample_path(..., base_ds_m=bezier_cfg.sample_ds_m, curve_type=:polyline)`. It tries control-point counts from `base_controls = max(2, n_waypoints + 2)` up to `max_controls = max(base_controls, min(size(samples,2), max(base_controls+6, size(path,2))))`, fitting each with `rpo_fit_bezier_fixed_endpoints` and scoring with `rpo_normalized_path_cost_components`. A candidate replaces the best when `J_obs` drops by more than `1e-9`, or ties on `J_obs` with a lower `total`; the loop stops early once `violation_count == 0`. The winner is passed through `rpo_post_refine_path` and returned alongside the raw polyline fields from the base plan.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `settings` | RPORRTConnectSettings | n/a | no | Keyword argument `settings` (default `RPORRTConnectSettings()`). |
| in | `max_runtime_s` | Real | n/a | no | Keyword argument `max_runtime_s` (default `Inf`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_connect_bezier_plan_path`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:323-323`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:454-454`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:442-442`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:441-441`
- `callees` → [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:453-453`
- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:427-427`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:466-466`
- `callees` → [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:431-431`
<!-- vulcan:connections:end -->

## Limitations
The search is greedy over increasing control counts and stops at the first violation-free fit, which is not necessarily the lowest-cost one. Fitting is least-squares on samples, so the Bezier curve can cut corners into the safety margin and rely on `rpo_post_refine_path` to recover. Runtime is roughly `max_controls - base_controls + 1` fits plus cost evaluations on top of the RRT-Connect cost.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 421.
