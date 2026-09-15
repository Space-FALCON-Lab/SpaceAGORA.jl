---
id: dynamics.perturbations__interp_vec3_linear
label: _interp_vec3_linear
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _interp_vec3_linear
  lines:
  - 1474
  - 1474
inputs:
- id: p0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p0`.
- id: p1
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p1`.
- id: tau
  type: Float64
  units: n/a
  required: true
  description: Positional argument `tau`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_interp_vec3_linear`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _interp_vec3_linear

## Purpose
Linear interpolation between two ephemeris positions, used in the first and last intervals of a cache table where Catmull-Rom interpolation lacks the outer neighbour it needs.

## Design & Implementation
Returns `p0 + tau (p1 - p0)` for `tau` in the unit interval, on static three-vectors. Declared `@inline` so the cache lookup's branch between linear and cubic interpolation costs nothing extra.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p0` | SVector{3, Float64} | n/a | yes | Positional argument `p0`. |
| in | `p1` | SVector{3, Float64} | n/a | yes | Positional argument `p1`. |
| in | `tau` | Float64 | n/a | yes | Positional argument `tau`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_interp_vec3_linear`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1104-1104`
- [[dynamics.perturbations__srp_sun_position_from_cache_j2000_m|_srp_sun_position_from_cache_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1467-1467`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Velocity is discontinuous at the boundary between a linear end segment and the cubic interior, which a very fine-stepped integrator can see as a small force jump twice per mission.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1474.
