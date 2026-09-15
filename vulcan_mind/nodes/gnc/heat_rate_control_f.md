---
id: gnc.heat_rate_control_f
label: f
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: f
  lines:
  - 57
  - 57
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
  description: Return value of `f`. Returns `begin`.
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
The residual whose root is the largest angle of attack that keeps free-molecular heat rate at the thermal limit, defined locally inside `_energy_depletion_heatrate_root_alpha` for the Newton solve.

## Theory & Math
With $s = S \sin\alpha$, the residual is

$$
f(\alpha) = L\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(e^{-s^2} + \sqrt{\pi}\, s\,(1+\operatorname{erf} s)\right) - \tfrac{1}{2} e^{-s^2}\right] - \dot{q}_{\text{lim}}
$$

where $L = 10^{-4}\,\text{taf}\,\rho R T_p \sqrt{R T_p / 2\pi}$, $S$ is the molecular speed ratio, $\gamma$ the ratio of specific heats, $T_w$ and $T_p$ the wall and freestream temperatures in kelvin, and $\dot{q}_{\text{lim}}$ the heat-rate limit in W/cm².

## Design & Implementation
A closure over `S`, `gamma`, `T_w`, `T_p`, the precomputed leading factor `L` and `thermal_limit`. It evaluates the Maxwellian flat-plate heat-rate expression in `alpha` — the Gaussian term in `S sin(alpha)`, the error-function term, and the enthalpy bracket in `gamma` and the wall-to-freestream temperature ratio — multiplies by `L` and subtracts the limit, so a zero crossing is exactly the constrained angle. It is paired with an analytic derivative `df` and handed to `Roots.find_zero` with `Roots.Newton()`, falling back through three further initial guesses and finally bisection on the whole bracket.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha` | Any | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `f`. Returns `begin`. |
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

- *none*
<!-- vulcan:connections:end -->

## Limitations
Being a closure it is reconstructed on every guidance call rather than compiled once, and `L` is captured at definition so a change in density mid-solve is not seen; the bracketing bisection fallback assumes the residual changes sign between `min_alpha` and `max_alpha`, which the earlier limit checks guarantee only when both endpoints straddle the limit.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl` line 57.
