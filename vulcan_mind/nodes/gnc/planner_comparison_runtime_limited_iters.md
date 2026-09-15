---
id: gnc.planner_comparison_runtime_limited_iters
label: runtime_limited_iters
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: runtime_limited_iters
  lines:
  - 271
  - 271
inputs:
- id: default_iters
  type: Any
  units: n/a
  required: true
  description: Positional argument `default_iters`.
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
  description: Return value of `runtime_limited_iters`. Returns `isfinite(Float64(runtime_limit_s))
    && cfg.optimizer.runtime_max_iters > 0 ?`.
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

# runtime_limited_iters

## Purpose
Local closure inside `rpo_plan_comparison_path` that decides how many iterations a non-HYPR planner receives: the runtime-limited cap, the HYPR-matched count, or the planner's own default.

## Design & Implementation
Defined at line 271 as `runtime_limited_iters(default_iters)`. It returns `cfg.optimizer.runtime_max_iters` when `isfinite(Float64(runtime_limit_s)) && cfg.optimizer.runtime_max_iters > 0` (a finite wall-clock limit is active, so the iteration cap is set high and time governs termination); otherwise `base_cfg.n_iters` when `cfg.optimizer.match_hypr_iters` is true; otherwise `default_iters`. It captures `runtime_limit_s`, `cfg`, and `base_cfg` from the enclosing scope and is applied to RRT-Connect, RRT*, CHOMP, STOMP, and the RRT seed iteration counts.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `default_iters` | Any | n/a | yes | Positional argument `default_iters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `runtime_limited_iters`. Returns `isfinite(Float64(runtime_limit_s)) && cfg.optimizer.runtime_max_iters > 0 ?`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:271-271`
- `callees` → [[gnc.planner_comparison__rpo_comparison_rrt_connect_seed_path|_rpo_comparison_rrt_connect_seed_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:370-370`
- `callees` → [[gnc.planner_comparison__rpo_comparison_rrt_connect_settings|_rpo_comparison_rrt_connect_settings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:301-301`
- `callees` → [[gnc.planner_comparison__rpo_comparison_rrt_star_settings|_rpo_comparison_rrt_star_settings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:343-343`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:288-288`
- `callees` → [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:323-323`
- `callees` → [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:344-344`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:281-281`
- `callees` → [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:419-419`
- `callees` → [[gnc.trajectory_optimizers_rpochompsettings|RPOCHOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:364-364`
- `callees` → [[gnc.trajectory_optimizers_rpostompsettings|RPOSTOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:401-401`
- `callees` → [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:276-276`
- `callees` → [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:302-302`
- `callees` → [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:380-380`
<!-- vulcan:connections:end -->

## Limitations
When a finite runtime limit is present, the cap is `runtime_max_iters` (default 100_000) regardless of `match_hypr_iters`, so iteration matching is overridden by time matching. The closure is recreated on every call to the enclosing function. `default_iters` is returned untyped, so callers convert with `Int` downstream.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 271.
