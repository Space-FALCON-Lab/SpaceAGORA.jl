---
id: gnc.tracking_executor_f
label: f
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: f
  lines:
  - 44
  - 44
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `f`. Returns `q * aerodynamic_coefficient_fM(x, m.body,
    T_p, S, m.aerodynamics, MonteCarlo)[2]`.
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
Residual closure defined inside `control_struct_load` (line 44) expressing the difference between the panel drag force at candidate angle `x` and the structural drag limit, so a root of `f` is the largest angle that keeps drag at the limit.

## Design & Implementation
`f(x) = q * aerodynamic_coefficient_fM(x, m.body, T_p, S, m.aerodynamics, MonteCarlo)[2] * area_tot - drag_limit`. It captures the dynamic pressure `q` (Pa), the body model, freestream temperature `T_p`, speed ratio `S`, the aerodynamics configuration, the Monte Carlo flag, the reference area `area_tot` (m^2) and `drag_limit = max_dyn_press * CD90 * area_tot` (N). Index `[2]` picks the drag coefficient from the `(CL, CD)` tuple returned by `aerodynamic_coefficient_fM`. It is handed to `find_zero(f, (0, pi/2), Roots.Bisection())` and only evaluated when the limit lies between `drag_min` and `drag_max`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `f`. Returns `q * aerodynamic_coefficient_fM(x, m.body, T_p, S, m.aerodynamics, MonteCarlo)[2]`. |
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

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:60-60`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:44-44`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:59-59`
- `callees` → [[gnc.tracking_executor__control_exception_fallback|_control_exception_fallback]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:56-56`
<!-- vulcan:connections:end -->

## Limitations
Each evaluation calls the full free-molecular coefficient routine, so bisection to default tolerance costs tens of aerodynamic evaluations per control step. The closure assumes drag is monotone in angle over `(0, pi/2)`; if the coefficient model is non-monotone, bisection may return an unintended root. The bracket ignores `max_α`, so the root may exceed the physical panel limit and later be clamped to 0.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 44.
