---
id: simulation.model_selection__gram_isolated_pool_batch_model_for_callback
label: _gram_isolated_pool_batch_model_for_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _gram_isolated_pool_batch_model_for_callback
  lines:
  - 46
  - 46
inputs:
- id: density_models
  type: AbstractVector{<:AbstractDensityModel}
  units: n/a
  required: true
  description: Positional argument `density_models`.
- id: fallback_model
  type: AbstractDensityModel
  units: n/a
  required: true
  description: Positional argument `fallback_model`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  description: 'Return value of `_gram_isolated_pool_batch_model_for_callback`. Returns
    `model isa EnvironmentModels.GRAMAtmosphereModel ? model : nothing`.'
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

# _gram_isolated_pool_batch_model_for_callback

## Purpose
Decides whether the GRAM isolated worker-pool batch path is applicable for a density callback over `num_sats` satellites, returning the `GRAMAtmosphereModel` template to clone or `nothing` if the pool should not be used.

## Design & Implementation
Three methods share the same gating logic. The pool is used only when `_gram_isolated_pool_enabled` (either the env-var-driven form or the `CallbackEnvConfig`-driven form with `env.gram_isolated_pool_mode` and `env.gram_isolated_pool_threshold`) returns true, when `density_models` is empty (no per-satellite models already exist), and when `fallback_model isa EnvironmentModels.GRAMAtmosphereModel`. Each failed condition short-circuits with `return nothing`. The `(p, num_sats)` method obtains the env config via `_callback_env_config(p)` and reads `p.shared_buffers.density_models` and `p.args.environment_model.density_model`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_models` | AbstractVector{<:AbstractDensityModel} | n/a | yes | Positional argument `density_models`. |
| in | `fallback_model` | AbstractDensityModel | n/a | yes | Positional argument `fallback_model`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_isolated_pool_batch_model_for_callback`. Returns `model isa EnvironmentModels.GRAMAtmosphereModel ? model : nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:267-267`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:267-267`

**Downstream**

- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:71-71`
- `callees` → [[simulation.config__gram_isolated_pool_enabled|_gram_isolated_pool_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
The pool is refused whenever any per-satellite density models are configured, even if they are all GRAM instances. The non-env method re-reads environment variables on every call via `_gram_isolated_pool_enabled(num_sats)`, which is slower than the cached `CallbackEnvConfig` route. `num_sats <= 0` returns `nothing` in `:on` mode. Nothing here validates that `Threads.nthreads() > 1`, which is left to the enabled check in `:auto` mode.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 46.
