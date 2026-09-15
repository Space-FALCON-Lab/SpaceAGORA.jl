---
id: simulation.config__policy_env_config
label: _policy_env_config
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _policy_env_config
  lines:
  - 222
  - 222
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_policy_env_config`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _policy_env_config

## Purpose
Run-scoped accessor for the `PolicyDecisionEnvConfig` snapshot, returning `nothing` when the run has not installed one so that the policy layer falls back to live environment reads.

## Design & Implementation
Returns `Union{Nothing, PolicyDecisionEnvConfig}`. It mirrors `_callback_env_config`: guard on `p !== nothing`, `hasproperty(p, :shared_buffers)` and `hasproperty(sb, :policy_env_config)`, then dereference `sb.policy_env_config[]`. Unlike its sibling it returns whatever the `Ref` holds — including `nothing` — rather than rebuilding a snapshot, because the downstream `thread_policy_decision` accepts `env=nothing` and reads the environment itself.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_policy_env_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
- [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:258-258`
- [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:334-334`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:633-633`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1051-1051`

**Downstream**

- `callees` → [[simulation.config__density_batch_enabled|_density_batch_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:232-232`
- `callees` → [[simulation.config__gram_isolated_pool_enabled|_gram_isolated_pool_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:242-242`
<!-- vulcan:connections:end -->

## Limitations
Returning `nothing` transfers the cost of environment parsing into `thread_policy_decision` on every call, so the absence of a snapshot degrades performance rather than raising an error. Callers must also fall back to `_callback_outer_parallel_hint()` for the outer-parallel signal, duplicating that logic at three call sites in this file.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 222.
