---
id: dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query
label: _aero_link_atmosphere_query
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _aero_link_atmosphere_query
  lines:
  - 506
  - 506
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: pos_pp_link
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_pp_link`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: SimulationModel.getDensity
  units: n/a
  description: Return value of `_aero_link_atmosphere_query`. Returns `SimulationModel.getDensity(density_model,
    alt, lat, lon, t, true, p)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _aero_link_atmosphere_query

## Purpose
Direct, uncached density/temperature/wind query at a single link's planet-fixed position for per-link atmosphere sampling.

## Design & Implementation
Converts `pos_pp_link` to `(alt, lat, lon)` via `rtolatlong(pos_pp_link, planet)`, resolves the satellite's density model with `SimulationCallbacks._density_model_for_sat(p, sat_idx)`, and returns `SimulationModel.getDensity(density_model, alt, lat, lon, t, true, p)`. It intentionally bypasses the satellite-level GRAM tracking cache because that cache assumes trajectory continuity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `pos_pp_link` | SVector{3, Float64} | n/a | yes | Positional argument `pos_pp_link`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.getDensity | n/a | — | Return value of `_aero_link_atmosphere_query`. Returns `SimulationModel.getDensity(density_model, alt, lat, lon, t, true, p)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:551-551`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:507-507`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:509-509`
- `callees` → [[simulation.model_selection__density_model_for_sat|_density_model_for_sat]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:508-508`
<!-- vulcan:connections:end -->

## Limitations
Costs one raw density-model call per non-root link per wrench evaluation, which for GRAM-backed models can dominate the RHS. The returned tuple is destructured as `(rho, T, wind)` by the caller, so the `getDensity` return contract must stay in that order.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 506.
