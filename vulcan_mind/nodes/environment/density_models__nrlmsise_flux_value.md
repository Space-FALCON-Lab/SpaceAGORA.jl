---
id: environment.density_models__nrlmsise_flux_value
label: _nrlmsise_flux_value
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_flux_value
  lines:
  - 334
  - 334
inputs:
- id: value
  type: Real
  units: n/a
  required: true
  description: Positional argument `value`.
- id: name
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `name`.
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
  description: Return value of `_nrlmsise_flux_value`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _nrlmsise_flux_value

## Purpose
Validates a solar flux input as finite and non-negative, naming the field in the error.

## Design & Implementation
Converts to `Float64` and throws `ArgumentError` mentioning `name` if the value is non-finite or negative. `@inline`. Used for both `f107` and `f107a` at construction and at provider resolution.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Real | n/a | yes | Positional argument `value`. |
| in | `name` | AbstractString | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nrlmsise_flux_value`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:570-570`
- [[environment.density_models__nrlmsise_space_indices_f107|_nrlmsise_space_indices_f107]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:519-519`
- [[environment.density_models__nrlmsise_space_indices_f107a|_nrlmsise_space_indices_f107a]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:523-523`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:335-335`
- `callees` → [[environment.density_models__gram_lock_scope|_gram_lock_scope]] · `callers` · call · `src/environment/atmosphere/density_models.jl:404-404`
- `callees` → [[environment.density_models__nrlmsise_ap_value|_nrlmsise_ap_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:371-371`
- `callees` → [[environment.density_models_nrlmsise00atmospheremodel|NRLMSISE00AtmosphereModel]] · `callers` · call · `src/environment/atmosphere/density_models.jl:342-342`
- `callees` → [[environment.density_models_nrlmsise00spaceindicesprovider|NRLMSISE00SpaceIndicesProvider]] · `callers` · call · `src/environment/atmosphere/density_models.jl:372-372`
<!-- vulcan:connections:end -->

## Limitations
No upper bound, so an obviously wrong flux such as 10,000 solar flux units passes.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 334.
