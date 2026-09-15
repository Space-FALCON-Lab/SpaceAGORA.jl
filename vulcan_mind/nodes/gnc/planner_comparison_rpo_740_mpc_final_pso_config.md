---
id: gnc.planner_comparison_rpo_740_mpc_final_pso_config
label: rpo_740_mpc_final_pso_config
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_740_mpc_final_pso_config
  lines:
  - 93
  - 93
inputs:
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
- id: goal_collision_margin_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `goal_collision_margin_m` (default `0.0`).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  description: Return value of `rpo_740_mpc_final_pso_config`. Returns `rpo_pso_config(cfg;
    kwargs...)`.
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

# rpo_740_mpc_final_pso_config

## Purpose
Builds the tuned `RPOPSOConfig` used for the published 740 MPC comparison case: a 100-particle, 60-iteration HYPR configuration with 5 kg spacecraft mass, 0.00625 m/s^2 retiming acceleration, adaptive downscaling, early stopping, culling, scheduling, probing, re-exploration, RRT warm-start, and refinement all enabled.

## Design & Implementation
Signature `rpo_740_mpc_final_pso_config(; safe_distance_m = 0.0, goal_collision_margin_m = 0.0, kwargs...)`. Constructs `RPOPSOConfig` with roughly 90 explicit keyword values, notably `n_waypoints = 5`, `n_particles = 100`, `n_iters = 60`, `w_obs = 1e6`, `retime_dt_s = 0.1`, `retime_max_speed_mps = 0.25`, `retime_max_steps = 20_000`, `adaptive_n_particles_min/max = 60/160`, `adaptive_n_iters_min/max = 10/60`, `early_stopping_patience = 10`, `cull_start_iter = 50`, `refinement_start_iter = 60`, and `rrt_warmstart_enable = true`. The two keep-out arguments are converted to `Float64`, then `rpo_pso_config(cfg; kwargs...)` applies further overrides and validation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `goal_collision_margin_m` | Real | n/a | no | Keyword argument `goal_collision_margin_m` (default `0.0`). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_740_mpc_final_pso_config`. Returns `rpo_pso_config(cfg; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:32-32`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:120-120`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:191-191`
- `callees` → [[gncy.pso_parameters_rpopsoconfig|RPOPSOConfig]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:94-94`
<!-- vulcan:connections:end -->

## Limitations
Any keyword in `kwargs` overrides the tuned values without warning, and unknown keywords raise `MethodError` from the config constructor. With `refinement_start_iter = 60` equal to `n_iters`, in-loop refinement never triggers; only the final refinement pass runs. The constants are specific to a 5 kg vehicle and 0.0125 m/s^2 class thrusters.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 93.
