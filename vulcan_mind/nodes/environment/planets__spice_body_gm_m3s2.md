---
id: environment.planets__spice_body_gm_m3s2
label: _spice_body_gm_m3s2
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _spice_body_gm_m3s2
  lines:
  - 245
  - 245
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  description: Return value of `_spice_body_gm_m3s2`.
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

# _spice_body_gm_m3s2

## Purpose
Reads a body's gravitational parameter from the SPICE kernel pool and converts it from km^3/s^2 to SI m^3/s^2.

## Design & Implementation
`@inline` function taking `planet_name::String`. Under `lock(SPICE_LOCK)` it calls `bodvrd(_spice_body_pool_name(planet_name), "GM")`, which returns a vector; the first element is multiplied by `1e9` and returned as `Float64`. Requires that a GM text kernel such as `gm_de440.tpc` was previously furnished.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_spice_body_gm_m3s2`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_backed_planet_kwargs|_spice_backed_planet_kwargs]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:281-281`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__spice_body_pool_name|_spice_body_pool_name]] · `callers` · call · `src/environment/ephemerides/planets.jl:247-247`
<!-- vulcan:connections:end -->

## Limitations
Throws a SPICE error if the `GM` variable is absent from the pool; `_spice_backed_planet_kwargs` relies on catching that. The `1e9` factor is a hard-coded unit conversion with no assertion that the kernel actually stores km^3/s^2.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 245.
