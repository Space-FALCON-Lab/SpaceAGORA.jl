---
id: environment.density_models__nrlmsise_resolved_indices
label: _nrlmsise_resolved_indices
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_resolved_indices
  lines:
  - 615
  - 615
inputs:
- id: model
  type: NRLMSISE00AtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: instant
  type: DateTime
  units: n/a
  required: true
  description: Positional argument `instant`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat`.
- id: lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon`.
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
  type: Tuple
  units: n/a
  description: Return value of `_nrlmsise_resolved_indices`. Returns `(f107a=model.f107a,
    f107=model.f107, ap=model.ap)` or `_nrlmsise_provider_indices(result)`.
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

# _nrlmsise_resolved_indices

## Purpose
Determines the indices for one evaluation: the model's fixed values, or the provider's result under whichever call signature it supports.

## Design & Implementation
Returns the fixed triple if `index_provider` is `nothing`. Otherwise it uses `applicable` to try the four-argument `(instant, h, lat, lon)` form first, then the single-argument form, raising if neither applies, and normalises the result through `_nrlmsise_provider_indices`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | NRLMSISE00AtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `lat` | Float64 | n/a | yes | Positional argument `lat`. |
| in | `lon` | Float64 | n/a | yes | Positional argument `lon`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_nrlmsise_resolved_indices`. Returns `(f107a=model.f107a, f107=model.f107, ap=model.ap)` or `_nrlmsise_provider_indices(result)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_density_state|_nrlmsise_density_state]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:645-645`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `callers` · call · `src/environment/atmosphere/density_models.jl:635-635`
<!-- vulcan:connections:end -->

## Limitations
`applicable` is evaluated on every call, which is a runtime method-table query; and a provider that accepts both signatures always gets the four-argument one.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 615.
