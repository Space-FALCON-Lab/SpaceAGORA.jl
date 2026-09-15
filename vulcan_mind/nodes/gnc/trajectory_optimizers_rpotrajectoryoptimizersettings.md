---
id: gnc.trajectory_optimizers_rpotrajectoryoptimizersettings
label: RPOTrajectoryOptimizerSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: RPOTrajectoryOptimizerSettings
  lines:
  - 22
  - 22
inputs:
- id: w_obs_scale
  type: Float64
  units: n/a
  required: false
  description: Field `w_obs_scale` (default `10.0`).
- id: obstacle_margin_m
  type: Float64
  units: n/a
  required: false
  description: Field `obstacle_margin_m` (default `0.5`).
- id: no_change_iters
  type: Int
  units: n/a
  required: false
  description: Field `no_change_iters` (default `8`).
- id: no_change_tol
  type: Float64
  units: n/a
  required: false
  description: Field `no_change_tol` (default `1.0e-7`).
- id: match_hypr_iters
  type: Bool
  units: n/a
  required: false
  description: Field `match_hypr_iters` (default `false`).
- id: match_hypr_runtime
  type: Bool
  units: n/a
  required: false
  description: Field `match_hypr_runtime` (default `false`).
- id: runtime_limit_s
  type: Float64
  units: n/a
  required: false
  description: Field `runtime_limit_s` (default `30.0`).
- id: runtime_max_iters
  type: Int
  units: n/a
  required: false
  description: Field `runtime_max_iters` (default `100_000`).
- id: force_full_iters
  type: Bool
  units: n/a
  required: false
  description: Field `force_full_iters` (default `true`).
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
  type: RPOTrajectoryOptimizerSettings
  units: n/a
  description: Constructed `RPOTrajectoryOptimizerSettings` (keyword constructor via
    @kwdef).
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

# RPOTrajectoryOptimizerSettings

## Purpose
Shared `Base.@kwdef` settings struct governing obstacle weighting, stall detection, and runtime/iteration matching that both `rpo_chomp_plan_path` and `rpo_stomp_plan_path` consume. It lets the comparison harness run CHOMP and STOMP under the same stopping rules as HYPR.

## Design & Implementation
Fields: `w_obs_scale::Float64 = 10.0` multiplies `cfg.w_obs` when the planners derive `local_cfg` via `rpo_pso_config`; `obstacle_margin_m::Float64 = 0.5` is the soft-potential width passed to `rpo_chomp_obstacle_potential` (the planners take `max(obstacle_margin_m, safe_distance_m)`); `no_change_iters::Int = 8` and `no_change_tol::Float64 = 1.0e-7` define stall termination, where an iteration is unchanged if `abs(previous_obj - best_obj) <= no_change_tol * max(1, abs(previous_obj))`; `match_hypr_iters::Bool = false`, `match_hypr_runtime::Bool = false`, `runtime_limit_s::Float64 = 30.0`, and `runtime_max_iters::Int = 100_000` are read by the planner-comparison driver (`runtime_limited_iters`) rather than by the optimizers themselves; `force_full_iters::Bool = true` disables the stall early-exit so every configured iteration is run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `w_obs_scale` | Float64 | n/a | no | Field `w_obs_scale` (default `10.0`). |
| in | `obstacle_margin_m` | Float64 | n/a | no | Field `obstacle_margin_m` (default `0.5`). |
| in | `no_change_iters` | Int | n/a | no | Field `no_change_iters` (default `8`). |
| in | `no_change_tol` | Float64 | n/a | no | Field `no_change_tol` (default `1.0e-7`). |
| in | `match_hypr_iters` | Bool | n/a | no | Field `match_hypr_iters` (default `false`). |
| in | `match_hypr_runtime` | Bool | n/a | no | Field `match_hypr_runtime` (default `false`). |
| in | `runtime_limit_s` | Float64 | n/a | no | Field `runtime_limit_s` (default `30.0`). |
| in | `runtime_max_iters` | Int | n/a | no | Field `runtime_max_iters` (default `100_000`). |
| in | `force_full_iters` | Bool | n/a | no | Field `force_full_iters` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOTrajectoryOptimizerSettings | n/a | — | Constructed `RPOTrajectoryOptimizerSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:37-37`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:363-363`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:232-232`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
With the default `force_full_iters = true`, `no_change_iters` and `no_change_tol` have no effect; callers must flip the flag explicitly. `runtime_limit_s` is not enforced by the optimizers unless the caller forwards it as `max_runtime_s`. No field range checks exist, so a negative `w_obs_scale` would invert obstacle avoidance.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 22.
