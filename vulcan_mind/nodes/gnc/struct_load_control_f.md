---
id: gnc.struct_load_control_f
label: f
kind: function
source:
  file: src/gnc/control/struct_load_control.jl
  symbol: f
  lines:
  - 80
  - 80
inputs:
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  description: Return value of `f`. Returns `q * _energy_depletion_struct_drag_area(         spacecraft,         temperature,`.
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

# f

## Purpose
`f(alpha)` is the residual closure that `_energy_depletion_struct_load_root_alpha` hands to `Roots.find_zero` so that a panel angle can be found at which the aerodynamic drag load on the spacecraft exactly equals the configured structural load budget. It exists only inside that function's scope and is the bridge between the drag-area model and the bisection solver.

## Design & Implementation
For a candidate panel angle `alpha` (rad) it evaluates `q * _energy_depletion_struct_drag_area(spacecraft, temperature, speed_ratio, controlled_panel_links, alpha, config) - drag_limit`, where `q = env.dynamic_pressure` (Pa), the drag area is the sum over `spacecraft.links` of `max(0, CD) * ref_area` (m^2) with the controlled links held at `alpha` and the remaining links at their own clamped `link.α`, and `drag_limit = config.structural_load_limit_pa * reference_drag_area` with the reference computed at `config.max_alpha_rad`. The sign of the result is therefore positive when the load exceeds the limit. The caller has already checked that `f(min_alpha) <= 0 < f(max_alpha)`, so the bracket `(min_alpha, max_alpha)` is valid for `Roots.Bisection()`; any exception from the solver is swallowed and `min_alpha` returned as the protective fallback. Nothing is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha` | Any | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `f`. Returns `q * _energy_depletion_struct_drag_area(         spacecraft,         temperature,`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_step_compliant_multibody_rk4|step_compliant_multibody_rk4]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:637-637`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:112-112`
- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `callers` · call · `src/parallel/policy/context.jl:46-46`
- [[grp.src_parallel_process|parallel/process/]] · `members_out` → `callers` · call · `src/parallel/process/worker_pool.jl:186-186`
- [[grp.src_parallel_routing|parallel/routing/]] · `members_out` → `callers` · call · `src/parallel/routing/env_mapping.jl:187-187`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- [[parallel.context__persistent_foreach_worker_loop|_persistent_foreach_worker_loop]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:46-46`
- [[parallel.context__persistent_foreach_worker_loop_w|_persistent_foreach_worker_loop_w]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:86-86`
- [[parallel.context__spin_barrier_dispatch_bang|_spin_barrier_dispatch!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:243-243`
- [[parallel.context__spin_barrier_worker_loop_w|_spin_barrier_worker_loop_w]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:178-178`
- [[parallel.env_mapping_with_parallel_profile|with_parallel_profile]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:187-187`
- [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:254-254`
- [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:262-262`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:204-204`
- [[parallel.thread_execution_threaded_foreach_worker_spin|threaded_foreach_worker_spin]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:168-168`
- [[parallel.worker_pool__warm_gram_wrapper_bang|_warm_gram_wrapper!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:186-186`
- [[parcore.context_with_policy_context|with_policy_context]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:337-337`
- [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:6-6`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- [[simulation.from_env__parse_or_default|_parse_or_default]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- [[simulation.from_env__with_engine_env_overrides|_with_engine_env_overrides]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:280-280`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:214-214`
- [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:79-79`
- [[simulation.registry__gram_runtime_stats_update_bang|_gram_runtime_stats_update!]] · `callees` → `callers` · call · `src/simulation/callbacks/registry.jl:78-78`

**Downstream**

- `callees` → [[gnc.struct_load_control__energy_depletion_struct_drag_area|_energy_depletion_struct_drag_area]] · `callers` · call · `src/gnc/control/struct_load_control.jl:80-80`
<!-- vulcan:connections:end -->

## Limitations
Each evaluation recomputes `aerodynamic_coefficient_fM` for every link, so the bisection cost scales with link count times solver iterations. Monotonicity of drag area in `alpha` is assumed but not verified; a non-monotone CD(alpha) could put multiple roots in the bracket and bisection will return an arbitrary one. The `catch` clause discards the exception type, hiding genuine errors (for example NaN coefficients) behind the min-drag fallback. `speed_ratio` is floored at `eps(Float64)` by the caller, so an unphysical near-zero ratio still produces a numeric, if meaningless, result.

## Provenance
Mapped from `src/gnc/control/struct_load_control.jl` line 80.
