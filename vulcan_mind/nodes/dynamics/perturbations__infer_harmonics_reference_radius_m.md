---
id: dynamics.perturbations__infer_harmonics_reference_radius_m
label: _infer_harmonics_reference_radius_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _infer_harmonics_reference_radius_m
  lines:
  - 223
  - 223
inputs:
- id: coefficients_file
  type: String
  units: n/a
  required: true
  description: Positional argument `coefficients_file`.
- id: planet
  type: AbstractPlanet
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
  type: Float64
  units: n/a
  description: Return value of `_infer_harmonics_reference_radius_m`.
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

# _infer_harmonics_reference_radius_m

## Purpose
Chooses the reference radius a coefficient file was defined against when the caller did not specify one.

## Design & Implementation
Returns 1,738.0 km for the lunar `lp165p.csv`, 3,396.0 km for Mars, and the planet's `Rp_e` otherwise. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `coefficients_file` | String | n/a | yes | Positional argument `coefficients_file`. |
| in | `planet` | AbstractPlanet | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_infer_harmonics_reference_radius_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:906-906`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Filename matching is a heuristic; a renamed lunar file gets the planet's equatorial radius and the field is evaluated at the wrong scale.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 223.
