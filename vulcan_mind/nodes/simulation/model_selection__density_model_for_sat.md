---
id: simulation.model_selection__density_model_for_sat
label: _density_model_for_sat
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _density_model_for_sat
  lines:
  - 1
  - 1
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
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  description: Return value of `_density_model_for_sat`. Returns `density_models[sat_idx]`
    or `fallback_model`.
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

# _density_model_for_sat

## Purpose
Resolves which `AbstractDensityModel` should evaluate atmospheric density for a single satellite in a multi-satellite simulation, allowing per-satellite model instances (for example isolated GRAM copies) while falling back to the environment-wide model.

## Design & Implementation
Two `@inline` methods. The three-argument form takes `density_models::AbstractVector{<:AbstractDensityModel}`, `fallback_model::AbstractDensityModel`, and `sat_idx::Int`; it returns `density_models[sat_idx]` when `sat_idx <= length(density_models)` and `fallback_model` otherwise. The two-argument form `(p, sat_idx)` unpacks the integrator parameter object, reading `p.shared_buffers.density_models` and `p.args.environment_model.density_model` and delegating. No state is mutated and no allocation occurs; the return type is a union of the vector eltype and the fallback type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_models` | AbstractVector{<:AbstractDensityModel} | n/a | yes | Positional argument `density_models`. |
| in | `fallback_model` | AbstractDensityModel | n/a | yes | Positional argument `fallback_model`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_density_model_for_sat`. Returns `density_models[sat_idx]` or `fallback_model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query|_aero_link_atmosphere_query]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:508-508`
- [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1948-1948`
- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:224-224`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:224-224`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A `sat_idx` of zero or negative is not rejected; because the bounds test only checks the upper end, `density_models[0]` would raise a `BoundsError`. Satellites beyond the populated vector silently share the fallback model, which for a stateful GRAM model reintroduces the cross-thread contention that per-satellite instances were meant to avoid. The function assumes `p.shared_buffers` and `p.args.environment_model` exist with those exact field names.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 1.
