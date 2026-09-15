---
id: gnc.planner_comparison_rpoplannercomparisonconfig
label: RPOPlannerComparisonConfig
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: RPOPlannerComparisonConfig
  lines:
  - 30
  - 30
inputs:
- id: planners
  type: Vector{Symbol}
  units: n/a
  required: false
  description: Field `planners` (default `[:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier,
    :rrt_star, :chomp, :stomp]`).
- id: pso_config
  type: RPOPSOConfig
  units: n/a
  required: false
  description: Field `pso_config` (default `rpo_740_mpc_final_pso_config(safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M)`).
- id: rrt_connect
  type: RPORRTConnectSettings
  units: n/a
  required: false
  description: Field `rrt_connect` (default `RPORRTConnectSettings()`).
- id: rrt_star
  type: RPORRTStarSettings
  units: n/a
  required: false
  description: Field `rrt_star` (default `RPORRTStarSettings()`).
- id: chomp
  type: RPOCHOMPSettings
  units: n/a
  required: false
  description: Field `chomp` (default `RPOCHOMPSettings(n_iters=pso_config.n_iters)`).
- id: stomp
  type: RPOSTOMPSettings
  units: n/a
  required: false
  description: Field `stomp` (default `RPOSTOMPSettings(n_iters=pso_config.n_iters)`).
- id: optimizer
  type: RPOTrajectoryOptimizerSettings
  units: n/a
  required: false
  description: Field `optimizer` (default `RPOTrajectoryOptimizerSettings()`).
- id: tracking
  type: RPOLQMPCTrackingSettings
  units: n/a
  required: false
  description: Field `tracking` (default `RPOLQMPCTrackingSettings()`).
- id: safe_distance_m
  type: Float64
  units: n/a
  required: false
  description: 'Field `safe_distance_m` (default `pso_config.safe_distance_m > 0.0
    ? pso_config.safe_distance_m : RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M`).'
- id: output_dir
  type: String
  units: n/a
  required: false
  description: Field `output_dir` (default `joinpath(pwd(), "rpo_planner_comparison")`).
- id: write_plotly_outputs
  type: Bool
  units: n/a
  required: false
  description: Field `write_plotly_outputs` (default `true`).
- id: write_failed_path_outputs
  type: Bool
  units: n/a
  required: false
  description: Field `write_failed_path_outputs` (default `true`).
- id: rng_seed
  type: Int
  units: n/a
  required: false
  description: Field `rng_seed` (default `740`).
- id: show_progress
  type: Bool
  units: n/a
  required: false
  description: Field `show_progress` (default `true`).
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
  type: RPOPlannerComparisonConfig
  units: n/a
  description: Constructed `RPOPlannerComparisonConfig` (keyword constructor via @kwdef).
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

# RPOPlannerComparisonConfig

## Purpose
Top-level configuration for an RPO planner comparison batch: which planners to run, the shared HYPR/PSO base config, per-planner settings for RRT-Connect, RRT*, CHOMP, and STOMP, the LQ-MPC tracking settings, keep-out distance, output location, and reproducibility seed.

## Design & Implementation
`Base.@kwdef struct` with `planners` defaulting to all seven canonical symbols; `pso_config = rpo_740_mpc_final_pso_config(safe_distance_m = 0.5)` (the constant `RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M`); `rrt_connect`, `rrt_star`, `optimizer`, `tracking` at their defaults; `chomp` and `stomp` constructed with `n_iters = pso_config.n_iters` so iteration budgets match HYPR; `safe_distance_m` derived from `pso_config.safe_distance_m` when positive, else 0.5; `output_dir = joinpath(pwd(), "rpo_planner_comparison")`; `write_plotly_outputs = true`; `write_failed_path_outputs = true`; `rng_seed = 740`; `show_progress = true`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planners` | Vector{Symbol} | n/a | no | Field `planners` (default `[:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, :stomp]`). |
| in | `pso_config` | RPOPSOConfig | n/a | no | Field `pso_config` (default `rpo_740_mpc_final_pso_config(safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M)`). |
| in | `rrt_connect` | RPORRTConnectSettings | n/a | no | Field `rrt_connect` (default `RPORRTConnectSettings()`). |
| in | `rrt_star` | RPORRTStarSettings | n/a | no | Field `rrt_star` (default `RPORRTStarSettings()`). |
| in | `chomp` | RPOCHOMPSettings | n/a | no | Field `chomp` (default `RPOCHOMPSettings(n_iters=pso_config.n_iters)`). |
| in | `stomp` | RPOSTOMPSettings | n/a | no | Field `stomp` (default `RPOSTOMPSettings(n_iters=pso_config.n_iters)`). |
| in | `optimizer` | RPOTrajectoryOptimizerSettings | n/a | no | Field `optimizer` (default `RPOTrajectoryOptimizerSettings()`). |
| in | `tracking` | RPOLQMPCTrackingSettings | n/a | no | Field `tracking` (default `RPOLQMPCTrackingSettings()`). |
| in | `safe_distance_m` | Float64 | n/a | no | Field `safe_distance_m` (default `pso_config.safe_distance_m > 0.0 ? pso_config.safe_distance_m : RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M`). |
| in | `output_dir` | String | n/a | no | Field `output_dir` (default `joinpath(pwd(), "rpo_planner_comparison")`). |
| in | `write_plotly_outputs` | Bool | n/a | no | Field `write_plotly_outputs` (default `true`). |
| in | `write_failed_path_outputs` | Bool | n/a | no | Field `write_failed_path_outputs` (default `true`). |
| in | `rng_seed` | Int | n/a | no | Field `rng_seed` (default `740`). |
| in | `show_progress` | Bool | n/a | no | Field `show_progress` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlannerComparisonConfig | n/a | — | Constructed `RPOPlannerComparisonConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_config_with_fixed_safe_distance|_rpo_comparison_config_with_fixed_safe_distance]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:49-49`
- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:544-544`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_740_mpc_final_pso_config|rpo_740_mpc_final_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:32-32`
- `callees` → [[gnc.planner_comparison_rpolqmpctrackingsettings|RPOLQMPCTrackingSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:38-38`
- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:33-33`
- `callees` → [[gnc.rrt_connect_rporrtstarsettings|RPORRTStarSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:34-34`
- `callees` → [[gnc.trajectory_optimizers_rpochompsettings|RPOCHOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:35-35`
- `callees` → [[gnc.trajectory_optimizers_rpostompsettings|RPOSTOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:36-36`
- `callees` → [[gnc.trajectory_optimizers_rpotrajectoryoptimizersettings|RPOTrajectoryOptimizerSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:37-37`
<!-- vulcan:connections:end -->

## Limitations
Because `output_dir` is evaluated from `pwd()` at construction, two configs built in different directories write to different places. `_rpo_comparison_config_with_fixed_safe_distance` overrides `safe_distance_m` back to 0.5 m in both the batch runner and single-path runner, so a user-set value is silently discarded. `planners` is not validated until `normalize_rpo_comparison_planner_type` throws at run time.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 30.
