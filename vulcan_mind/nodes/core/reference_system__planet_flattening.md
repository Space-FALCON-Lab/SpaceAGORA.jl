---
id: core.reference_system__planet_flattening
label: _planet_flattening
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _planet_flattening
  lines:
  - 87
  - 87
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
  type: Float64
  units: n/a
  description: Return value of `_planet_flattening`.
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

# _planet_flattening

## Purpose
Computes the planet's geometric flattening from its equatorial and polar radii, used by the geodetic conversions.

## Theory & Math
$$
f = \frac{R_e - R_p}{R_e}
$$

## Design & Implementation
Returns `(Rp_e - Rp_p) / Rp_e`. `@inline` with a `::Float64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_planet_flattening`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_rtolatlong|rtolatlong]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:357-357`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system__body_fixed_to_j2000_state|_body_fixed_to_j2000_state]] · `callers` · call · `src/core/interfaces/reference_system.jl:117-117`
- `callees` → [[core.reference_system__j2000_to_body_fixed_state|_j2000_to_body_fixed_state]] · `callers` · call · `src/core/interfaces/reference_system.jl:92-92`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/core/interfaces/reference_system.jl:91-91`
- `callees` → [[core.reference_system_r_pintor_i|r_pintor_i]] · `callers` · call · `src/core/interfaces/reference_system.jl:116-116`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/core/interfaces/reference_system.jl:108-108`
<!-- vulcan:connections:end -->

## Limitations
Assumes an oblate spheroid with `Rp_e ≥ Rp_p`; a prolate or triaxial body model would produce a negative or meaningless flattening.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 87.
