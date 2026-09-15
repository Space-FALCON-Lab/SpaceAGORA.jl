---
id: environment.simple_ephemerides__spice_body_fixed_frame
label: _spice_body_fixed_frame
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _spice_body_fixed_frame
  lines:
  - 75
  - 75
inputs:
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
  type: String
  units: n/a
  description: Return value of `_spice_body_fixed_frame`.
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

# _spice_body_fixed_frame

## Purpose
Returns the NAIF body-fixed reference frame name to use for a given planet when querying SPICE, so that `pxform` calls target the highest-fidelity frame available for that body.

## Design & Implementation
A nested ternary on `planet.name`: `"Moon"` yields `"MOON_PA_DE421"` (the DE421 principal-axis lunar frame), `"Earth"` yields the constant `_EARTH_HIGH_PREC_BODY_FIXED_FRAME = "ITRF93"`, and any other planet yields `"IAU_" * uppercase(planet.name)`, following the NAIF convention for IAU rotational-element frames such as `IAU_MARS` or `IAU_VENUS`. The result is a `String` and the function is `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_spice_body_fixed_frame`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__body_fixed_to_j2000_state|_body_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:84-84`
- [[core.reference_system__j2000_to_body_fixed_state|_j2000_to_body_fixed_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:67-67`
- [[environment.simple_ephemerides__spice_planet_frame_lpi|_spice_planet_frame_lpi]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:88-88`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `IAU_` fallback assumes the planet name matches NAIF's spelling once upper-cased; bodies with multi-word or non-standard names produce a frame that `pxform` will reject with a SPICE error. `ITRF93` and `MOON_PA_DE421` require their respective binary PCK kernels to be furnished; this function does not verify kernel availability, leaving `_spice_planet_frame_lpi` to catch the Earth failure. Comparison is case-sensitive on `planet.name`.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 75.
