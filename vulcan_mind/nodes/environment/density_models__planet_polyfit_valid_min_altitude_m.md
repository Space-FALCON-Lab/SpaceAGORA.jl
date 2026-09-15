---
id: environment.density_models__planet_polyfit_valid_min_altitude_m
label: _planet_polyfit_valid_min_altitude_m
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _planet_polyfit_valid_min_altitude_m
  lines:
  - 218
  - 218
inputs:
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
  description: Return value of `_planet_polyfit_valid_min_altitude_m`.
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

# _planet_polyfit_valid_min_altitude_m

## Purpose
Reads a planet's lower polyfit validity altitude if it declares one, otherwise the 50 km default.

## Design & Implementation
Tests `hasproperty(planet, :polyfit_valid_min_altitude_m)` and converts the value to `Float64`, else returns `_POLYFIT_DEFAULT_VALID_MIN_ALTITUDE_M`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | AbstractPlanet | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_planet_polyfit_valid_min_altitude_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__planet_polyfit_valid_max_altitude_m|_planet_polyfit_valid_max_altitude_m]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:235-235`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:220-220`
<!-- vulcan:connections:end -->

## Limitations
A planet declaring the property with a non-numeric value fails at the `Float64` conversion with an unhelpful message.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 218.
