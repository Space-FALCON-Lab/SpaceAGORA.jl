---
id: core.reference_system__spice_body_fixed_frame
label: _spice_body_fixed_frame
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _spice_body_fixed_frame
  lines:
  - 43
  - 43
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
- core
charts:
- core
origin: agent
---

# _spice_body_fixed_frame

## Purpose
Names the SPICE body-fixed frame appropriate for a planet, choosing higher-precision frames for the Moon and Earth.

## Design & Implementation
Returns `MOON_PA_DE421` for the Moon, `ITRF93` for Earth, and `IAU_` followed by the upper-cased planet name otherwise. `@inline` with a `::String` return.

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
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The Moon and Earth frames require their respective binary PCK kernels to be furnished; if they are not, `sxform` throws and only the Earth callers have a fallback.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 43.
