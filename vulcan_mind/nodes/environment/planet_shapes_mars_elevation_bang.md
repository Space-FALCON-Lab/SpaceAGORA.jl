---
id: environment.planet_shapes_mars_elevation_bang
label: Mars_elevation!
kind: function
source:
  file: src/environment/ephemerides/planet_shapes.jl
  symbol: Mars_elevation!
  lines:
  - 97
  - 97
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: Clm
  type: Any
  units: n/a
  required: true
  description: Positional argument `Clm`.
- id: Slm
  type: Any
  units: n/a
  required: true
  description: Positional argument `Slm`.
- id: latitude
  type: Any
  units: n/a
  required: true
  description: Positional argument `latitude`.
- id: longitude
  type: Any
  units: n/a
  required: true
  description: Positional argument `longitude`.
- id: A
  type: Any
  units: n/a
  required: true
  description: Positional argument `A`.
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
  type: Any
  units: n/a
  description: Return value of `Mars_elevation!`; mutates `args` in place. Returns
    `harmonic`.
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

# Mars_elevation!

## Purpose

Evaluates Mars topography at a given geodetic latitude and longitude. Because the MOLA-derived Mars coefficient set already expresses the shape as a distance from the centre of mass, the harmonic sum is itself the planetary radius and is returned unmodified.

## Design & Implementation

Reads `topo_degree` and `topo_order` from the `args` dictionary under the `:topo_degree` and `:topo_order` keys, then delegates to `calculate_topography_harmonics!(Clm, Slm, latitude, longitude, A, topo_degree, topo_order)`, which fills the caller's Legendre workspace `A` in place and returns the summed harmonic. The bang in the name refers to that workspace mutation, not to any change of `Clm` or `Slm`, which are read only. Latitude and longitude are radians; the result is metres of radius from the Mars centre of mass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `Clm` | Any | n/a | yes | Positional argument `Clm`. |
| in | `Slm` | Any | n/a | yes | Positional argument `Slm`. |
| in | `latitude` | Any | n/a | yes | Positional argument `latitude`. |
| in | `longitude` | Any | n/a | yes | Positional argument `longitude`. |
| in | `A` | Any | n/a | yes | Positional argument `A`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `Mars_elevation!`; mutates `args` in place. Returns `harmonic`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`

**Downstream**

- `callees` → [[envana.env_planet_shapes_calculate_topography_harmonics_bang|calculate_topography_harmonics!]] · `callers` · call · `src/environment/ephemerides/planet_shapes.jl:124-124`
<!-- vulcan:connections:end -->

## Limitations

The scratch array `A` is written on every call, so a single workspace cannot be shared between threads sampling different ground points. Unlike `Venus_elevation!` no mean radius is added, which is correct only for a coefficient set already referenced to the centre of mass — feeding Venus-style normalised coefficients here silently returns a number near zero instead of a radius. `args` is untyped and the two keys are assumed present; degree and order are not checked against the dimensions of `Clm`, `Slm` or `A`.

## Provenance
Mapped from `src/environment/ephemerides/planet_shapes.jl` line 97.
