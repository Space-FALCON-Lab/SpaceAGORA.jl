---
id: gnc.planner_comparison_rpo_plan_comparison_path
label: rpo_plan_comparison_path
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_plan_comparison_path
  lines:
  - 258
  - 258
inputs:
- id: planner_type
  type: Any
  units: n/a
  required: true
  description: Positional argument `planner_type`.
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
  type: RPOPlannerComparisonConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
- id: runtime_limit_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `runtime_limit_s` (default `Inf`).
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
  description: Return value of `rpo_plan_comparison_path`. Returns `merge(plan, (`.
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

# rpo_plan_comparison_path

## Purpose
Runs a single planner on one start-goal case with the comparison's shared configuration, timing the run and returning the planner's native result tuple merged with normalised bookkeeping fields (`planner_type`, `planner_label`, `planner_compute_time`, `planner_plan_time`, `planner_refinement_time`, `planner_iteration_count`).

## Design & Implementation
Signature `rpo_plan_comparison_path(planner_type, start_rtn, goal_rtn, geometry, cfg::RPOPlannerComparisonConfig; rng = Random.default_rng(), runtime_limit_s = Inf)`. It pins the keep-out distance, normalises the planner symbol, derives `base_cfg = rpo_pso_config(cfg.pso_config; safe_distance_m = cfg.safe_distance_m)`, defines the `runtime_limited_iters` closure, and dispatches: `:hypr` -> `rpo_pso_plan_path`; `:pso_unrefined` -> same with `refinement_enable = false`; `:rrt_connect` / `:rrt_connect_bezier` / `:rrt_star` -> the respective planners with iteration-capped settings and `max_runtime_s`; `:chomp` and `:stomp` -> an RRT-Connect seed followed by `rpo_chomp_plan_path` or `rpo_stomp_plan_path` with `cfg.optimizer`. Elapsed time is `(time_ns() - t0) / 1e9`; iteration count is `length(plan.cost_history)` for PSO variants and `plan.iterations` otherwise.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planner_type` | Any | n/a | yes | Positional argument `planner_type`. |
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPlannerComparisonConfig | n/a | yes | Positional argument `cfg`. |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `runtime_limit_s` | Real | n/a | no | Keyword argument `runtime_limit_s` (default `Inf`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_plan_comparison_path`. Returns `merge(plan, (`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:575-575`

**Downstream**

- `callees` → [[gnc.planner_comparison__rpo_comparison_config_with_fixed_safe_distance|_rpo_comparison_config_with_fixed_safe_distance]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:267-267`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:269-269`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:268-268`
<!-- vulcan:connections:end -->

## Limitations
`planner_refinement_time` is always 0.0 and `planner_plan_time == planner_compute_time`, so refinement cost is never separated. For CHOMP/STOMP the seed RRT time is folded into the optimizer's time. `runtime_limit_s` is not passed to the HYPR branches, so HYPR always runs its own budget. Unsupported planners throw `ArgumentError` from normalisation.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 258.
