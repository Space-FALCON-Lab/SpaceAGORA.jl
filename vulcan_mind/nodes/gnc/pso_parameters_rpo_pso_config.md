---
id: gnc.pso_parameters_rpo_pso_config
label: rpo_pso_config
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: rpo_pso_config
  lines:
  - 583
  - 583
inputs:
- id: base
  type: RPOPSOConfig
  units: n/a
  required: false
  description: Positional argument `base` (default `RPOPSOConfig()`).
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
  description: Return value of `rpo_pso_config`. Returns `validate_rpo_pso_config(RPOPSOConfig(;
    values...))`.
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

# rpo_pso_config

## Purpose
Primary factory for `RPOPSOConfig`: starts from an existing flattened config (default-constructed if omitted) or an `RPOPSOConfigurator`, applies keyword overrides with alias support, syncs sampling density with the keep-out distance, and validates the result.

## Design & Implementation
Method 1: `rpo_pso_config(base::RPOPSOConfig = RPOPSOConfig(); kwargs...)` computes `values = merge(_rpo_pso_config_tuple(base), _rpo_pso_normalize_kwargs(kwargs))`, then `values = _rpo_pso_sync_sample_ds_with_safe_distance(values)`, and returns `validate_rpo_pso_config(RPOPSOConfig(; values...))`. Method 2: `rpo_pso_config(configurator::RPOPSOConfigurator; kwargs...)` delegates to `RPOPSOConfig(configurator; kwargs...)`, which flattens the groups and then calls method 1. Both return a new immutable config; nothing is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `base` | RPOPSOConfig | n/a | no | Positional argument `base` (default `RPOPSOConfig()`). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_config`. Returns `validate_rpo_pso_config(RPOPSOConfig(; values...))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_config_with_fixed_safe_distance|_rpo_comparison_config_with_fixed_safe_distance]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:51-51`
- [[gnc.planner_comparison_rpo_740_mpc_final_pso_config|rpo_740_mpc_final_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:191-191`
- [[gnc.planner_comparison_rpo_plan_comparison_path|rpo_plan_comparison_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:269-269`
- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:465-465`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:288-288`
- [[gnc.pso_parameters_push_bang|push!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:579-579`
- [[gnc.pso_path_planning_cfg_for_current|cfg_for_current]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:286-286`
- [[gnc.pso_refinement_rpo_refinement_config|rpo_refinement_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:4-4`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:19-19`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:441-441`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:498-498`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:368-368`
- [[gncy.pso_adaptive_policy_rpo_adaptive_pso_config|rpo_adaptive_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:63-63`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:189-189`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:318-318`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:236-236`

**Downstream**

- `callees` → [[gnc.pso_parameters__rpo_pso_config_tuple|_rpo_pso_config_tuple]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:584-584`
- `callees` → [[gnc.pso_parameters__rpo_pso_normalize_kwargs|_rpo_pso_normalize_kwargs]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:584-584`
- `callees` → [[gnc.pso_parameters__rpo_pso_sync_sample_ds_with_safe_distance|_rpo_pso_sync_sample_ds_with_safe_distance]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:585-585`
- `callees` → [[gnc.pso_parameters_validate_rpo_pso_config|validate_rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:586-586`
- `callees` → [[gncy.pso_parameters_rpopsoconfig|RPOPSOConfig]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:583-583`
<!-- vulcan:connections:end -->

## Limitations
Throws `ArgumentError` from validation for out-of-range fields and `MethodError` for unknown keyword names. Every call allocates a ~115-element NamedTuple twice, so it should be called at configuration time rather than per planning iteration. The sample-spacing sync silently overrides an explicit `sample_ds_m` whenever `safe_distance_m > 0`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 583.
