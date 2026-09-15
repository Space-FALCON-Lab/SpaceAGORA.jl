---
id: environment.planet_shapes_earth_elevation_bang
label: Earth_elevation!
kind: function
source:
  file: src/environment/ephemerides/planet_shapes.jl
  symbol: Earth_elevation!
  lines:
  - 161
  - 161
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: Clm
  type: AbstractArray{Float64}
  units: n/a
  required: true
  description: Positional argument `Clm`.
- id: Slm
  type: AbstractArray{Float64}
  units: n/a
  required: true
  description: Positional argument `Slm`.
- id: latitude
  type: Float64
  units: n/a
  required: true
  description: Positional argument `latitude`.
- id: longitude
  type: Float64
  units: n/a
  required: true
  description: Positional argument `longitude`.
- id: A
  type: AbstractArray{Float64}
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
  description: Return value of `Earth_elevation!`; mutates `args` in place. Returns
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

# Earth_elevation!

## Purpose

Evaluates Earth topography at a given latitude and longitude, returning the harmonic sum directly because the Earth coefficient set — like the Mars one — already gives distance from the centre of mass rather than a height above a reference sphere.

## Design & Implementation

Reads `args[:topo_degree]` and `args[:topo_order]`, forwards `Clm`, `Slm`, `latitude`, `longitude` and the Legendre workspace `A` to `calculate_topography_harmonics!`, and returns the resulting scalar. That callee evaluates `fnALF_IDR!` into `A` at `sin(latitude)` and accumulates the doubly indexed sum with `@turbo`, using a Chebyshev-style recursion for `cos mλ` and `sin mλ` instead of repeated trigonometric calls. Signature and behaviour deliberately mirror `Mars_elevation!` so the three planetary shape models are interchangeable at their call sites.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `Clm` | AbstractArray{Float64} | n/a | yes | Positional argument `Clm`. |
| in | `Slm` | AbstractArray{Float64} | n/a | yes | Positional argument `Slm`. |
| in | `latitude` | Float64 | n/a | yes | Positional argument `latitude`. |
| in | `longitude` | Float64 | n/a | yes | Positional argument `longitude`. |
| in | `A` | AbstractArray{Float64} | n/a | yes | Positional argument `A`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `Earth_elevation!`; mutates `args` in place. Returns `harmonic`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`

**Downstream**

- `callees` → [[envana.env_planet_shapes_calculate_topography_harmonics_bang|calculate_topography_harmonics!]] · `callers` · call · `src/environment/ephemerides/planet_shapes.jl:188-188`
<!-- vulcan:connections:end -->

## Limitations

No mean radius is added, so the interpretation of the return value depends entirely on the convention of the loaded coefficient file; there is no runtime check that the right file was loaded. The `A` workspace is mutated on each call and therefore not safe to share across threads. Requested degree and order are trusted, and the `@turbo` accumulation in the callee reassociates the floating-point sum, so results are not bit-reproducible against a plain scalar loop.

## Provenance
Mapped from `src/environment/ephemerides/planet_shapes.jl` line 161.
