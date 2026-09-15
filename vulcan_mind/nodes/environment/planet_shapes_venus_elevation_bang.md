---
id: environment.planet_shapes_venus_elevation_bang
label: Venus_elevation!
kind: function
source:
  file: src/environment/ephemerides/planet_shapes.jl
  symbol: Venus_elevation!
  lines:
  - 128
  - 128
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
  description: Return value of `Venus_elevation!`; mutates `args` in place. Returns
    `elevation`.
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

# Venus_elevation!

## Purpose

Evaluates Venus topography at a given latitude and longitude and converts it to an absolute planetary radius. The Venus coefficient set expresses shape relative to the mean radius, so this variant adds that reference sphere back in.

## Design & Implementation

Pulls `topo_degree` and `topo_order` from `args`, calls `calculate_topography_harmonics!(Clm, Slm, latitude, longitude, A, topo_degree, topo_order)` — which mutates the preallocated Legendre workspace `A` — and returns `6051.8e3 + harmonic` in metres. The mean radius is a hard-coded literal rather than a lookup, which is the only structural difference from `Mars_elevation!` and `Earth_elevation!`. Arguments are typed as `AbstractArray{Float64}` and `Float64`, so the method will not accept a mixed-precision coefficient table.

## Theory & Math

The elevation is the reference sphere plus the truncated surface harmonic expansion in fully normalized associated Legendre functions,

$$R(\varphi, \lambda) = R_0 + \sum_{n=0}^{N} \sum_{m=0}^{M} \bar{P}_{nm}(\sin\varphi)\left[\bar{C}_{nm}\cos m\lambda + \bar{S}_{nm}\sin m\lambda\right]$$

with $R_0 = 6\,051.8$ km, $\varphi$ the latitude and $\lambda$ the longitude, both in radians.

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
| out | `result` | Any | n/a | — | Return value of `Venus_elevation!`; mutates `args` in place. Returns `elevation`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`

**Downstream**

- `callees` → [[envana.env_planet_shapes_calculate_topography_harmonics_bang|calculate_topography_harmonics!]] · `callers` · call · `src/environment/ephemerides/planet_shapes.jl:155-155`
<!-- vulcan:connections:end -->

## Limitations

The 6051.8 km mean radius is embedded in the source, so it cannot be reconciled with whatever reference radius the loaded coefficient file actually used; a mismatched coefficient set produces a radius biased by the difference. Truncating at `topo_degree`/`topo_order` gives a smoothed shape whose error grows in steep terrain, and no bounds are enforced between the requested degree and the size of `Clm`, `Slm` or `A`. The mutated workspace `A` makes the call non-reentrant on a shared buffer.

## Provenance
Mapped from `src/environment/ephemerides/planet_shapes.jl` line 128.
