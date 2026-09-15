---
id: dynamics.perturbations__harmonics_model_cache_key
label: _harmonics_model_cache_key
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_model_cache_key
  lines:
  - 675
  - 675
inputs:
- id: L
  type: Int64
  units: n/a
  required: true
  description: Positional argument `L`.
- id: M
  type: Int64
  units: n/a
  required: true
  description: Positional argument `M`.
- id: coefficients_file
  type: String
  units: n/a
  required: true
  description: Positional argument `coefficients_file`.
- id: coefficient_normalization
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `coefficient_normalization`.
- id: coefficients_normalized
  type: Bool
  units: n/a
  required: true
  description: Positional argument `coefficients_normalized`.
- id: j2_source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `j2_source`.
- id: reference_radius_m
  type: Union{Nothing, Float64}
  units: n/a
  required: true
  description: Positional argument `reference_radius_m`.
- id: include_central
  type: Bool
  units: n/a
  required: true
  description: Positional argument `include_central`.
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
  type: Any
  units: n/a
  description: Return value of `_harmonics_model_cache_key`. Returns `(L, M, coefficients_file,
    coefficient_normalization, coefficients_normalized,`.
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

# _harmonics_model_cache_key

## Purpose
Forms the key under which constructed harmonics models are memoised so repeated construction with the same parameters reuses the loaded coefficient matrices.

## Design & Implementation
Returns a tuple of degree, order, file path, normalisation symbol and flag, J2 source, reference radius, `include_central` and `objectid(planet)`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `L` | Int64 | n/a | yes | Positional argument `L`. |
| in | `M` | Int64 | n/a | yes | Positional argument `M`. |
| in | `coefficients_file` | String | n/a | yes | Positional argument `coefficients_file`. |
| in | `coefficient_normalization` | Symbol | n/a | yes | Positional argument `coefficient_normalization`. |
| in | `coefficients_normalized` | Bool | n/a | yes | Positional argument `coefficients_normalized`. |
| in | `j2_source` | Symbol | n/a | yes | Positional argument `j2_source`. |
| in | `reference_radius_m` | Union{Nothing, Float64} | n/a | yes | Positional argument `reference_radius_m`. |
| in | `include_central` | Bool | n/a | yes | Positional argument `include_central`. |
| in | `planet` | AbstractPlanet | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_harmonics_model_cache_key`. Returns `(L, M, coefficients_file, coefficient_normalization, coefficients_normalized,`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- `callees` → [[dynamics.perturbations__canonical_harmonics_normalization|_canonical_harmonics_normalization]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:725-725`
- `callees` → [[dynamics.perturbations__convert_harmonics_coefficients_to_full_bang|_convert_harmonics_coefficients_to_full!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:835-835`
- `callees` → [[dynamics.perturbations__infer_harmonics_reference_radius_m|_infer_harmonics_reference_radius_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:906-906`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- `callees` → [[dynamics.perturbations_gravitationalharmonicsmodel|GravitationalHarmonicsModel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:685-685`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:886-886`
<!-- vulcan:connections:end -->

## Limitations
Using `objectid(planet)` means two equal but distinct planet objects do not share a cached model.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 675.
