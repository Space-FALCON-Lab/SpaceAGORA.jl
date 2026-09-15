---
id: environment.density_models__nrlmsise_density_state
label: _nrlmsise_density_state
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_density_state
  lines:
  - 638
  - 638
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_nrlmsise_density_state`.
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

# _nrlmsise_density_state

## Purpose
Evaluates NRLMSISE-00 at one point and returns density, temperature and a zero wind vector in the common density-state format.

## Design & Implementation
Resolves indices, calls `AtmosphericModels.nrlmsise00` with the instant, altitude, latitude, longitude, fluxes, Ap and the anomalous-oxygen flag, and returns `total_density`, `temperature` and zero wind. `@inline`.

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
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_nrlmsise_density_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__gram_point_density|_gram_point_density]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1078-1078`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_resolved_indices|_nrlmsise_resolved_indices]] · `callers` · call · `src/environment/atmosphere/density_models.jl:645-645`
<!-- vulcan:connections:end -->

## Limitations
The underlying call allocates its output struct; and no horizontal wind model accompanies it, so runs with `wind=true` still see zero wind.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 638.
