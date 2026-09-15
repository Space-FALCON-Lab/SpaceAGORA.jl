---
id: environment.planets__spice_body_pool_name
label: _spice_body_pool_name
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _spice_body_pool_name
  lines:
  - 234
  - 234
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
  type: String
  units: n/a
  description: Return value of `_spice_body_pool_name`.
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

# _spice_body_pool_name

## Purpose
Maps a planet display name to the identifier used for kernel-pool variable lookups with `bodvrd`.

## Design & Implementation
`@inline` function returning `"MOON"` when `planet_name == "Moon"` and `uppercase(planet_name)` otherwise, so `"Mars"` becomes `"MARS"` and `"Titan"` becomes `"TITAN"`. Used by `_spice_body_radii_m` and `_spice_body_gm_m3s2`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_spice_body_pool_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_gm_m3s2|_spice_body_gm_m3s2]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:247-247`
- [[environment.planets__spice_body_radii_m|_spice_body_radii_m]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:240-240`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The Moon special case is redundant with `uppercase` and exists only for clarity. Names with spaces or non-ASCII characters are not sanitised, and no validation that the resulting name is a body SPICE knows is performed here.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 234.
