---
id: environment.planets__spice_body_radii_m
label: _spice_body_radii_m
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _spice_body_radii_m
  lines:
  - 238
  - 238
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
  type: NTuple{3,
  units: n/a
  description: Return value of `_spice_body_radii_m`.
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

# _spice_body_radii_m

## Purpose
Fetches a body's triaxial ellipsoid radii from the loaded PCK and converts them to metres for use in the planet struct's `Rp_e`, `Rp_p`, and `Rp_m` fields.

## Design & Implementation
`@inline` function returning `NTuple{3, Float64}`. It takes `SPICE_LOCK`, calls `bodvrd(_spice_body_pool_name(planet_name), "RADII")` to get the three radii in kilometres, and returns each element multiplied by `1e3`. Order follows the PCK convention: two equatorial axes then the polar axis.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NTuple{3, | n/a | — | Return value of `_spice_body_radii_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_backed_planet_kwargs|_spice_backed_planet_kwargs]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:271-271`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__spice_body_pool_name|_spice_body_pool_name]] · `callers` · call · `src/environment/ephemerides/planets.jl:240-240`
<!-- vulcan:connections:end -->

## Limitations
Requires a text PCK (pck00010/pck00011) to be furnished first; otherwise `bodvrd` raises a SPICE error that propagates. No check that exactly three values were returned; a malformed kernel with fewer entries would throw `BoundsError` on indexing.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 238.
