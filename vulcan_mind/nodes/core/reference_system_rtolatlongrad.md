---
id: core.reference_system_rtolatlongrad
label: rtolatlongrad
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtolatlongrad
  lines:
  - 392
  - 392
inputs:
- id: r_p
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_p`.
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
  type: SVector
  units: n/a
  description: Return value of `rtolatlongrad`. Returns `SVector{3, Float64}([r, lat,
    lon])`.
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

# rtolatlongrad

## Purpose
Converts a planet-fixed position to planetocentric radius, latitude and longitude, the spherical companion to the geodetic `rtolatlong`.

## Design & Implementation
Latitude is `asin(z / |r|)`, longitude `atan(y, x)`, and the first element is the full radial distance. Returns an `SVector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_p` | Any | n/a | yes | Positional argument `r_p`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `rtolatlongrad`. Returns `SVector{3, Float64}([r, lat, lon])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The variable `r` is assigned twice with different meanings — first the equatorial distance, then the norm — and only the second is used; planetocentric latitude differs from geodetic by up to the flattening angle, so mixing the two conventions between callers is a real source of altitude error.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 392.
