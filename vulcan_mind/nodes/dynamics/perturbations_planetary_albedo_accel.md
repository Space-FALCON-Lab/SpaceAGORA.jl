---
id: dynamics.perturbations_planetary_albedo_accel
label: planetary_albedo_accel
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: planetary_albedo_accel
  lines:
  - 1183
  - 1183
inputs:
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: pos_primary_sun
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_primary_sun`.
- id: planet_radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `planet_radius_m`.
- id: p_srp_unscaled
  type: Float64
  units: n/a
  required: true
  description: Positional argument `p_srp_unscaled`.
- id: reflection_coefficient
  type: Float64
  units: n/a
  required: true
  description: Positional argument `reflection_coefficient`.
- id: reference_area_m2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `reference_area_m2`.
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass_kg`.
- id: planet_albedo
  type: Float64
  units: n/a
  required: true
  description: Positional argument `planet_albedo`.
- id: AU_m
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `AU_m` (default `149_597_870_700.0`).
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
  description: Return value of `planetary_albedo_accel`.
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

# planetary_albedo_accel

## Purpose
Acceleration from sunlight reflected by the planet, modelled as a Lambertian sphere with a single albedo.

## Design & Implementation
Returns zero for a non-positive radius or albedo, or a spacecraft below the surface. Scales the solar pressure by `(R/r)²`, the albedo and the Lambert phase function at the Sun-planet-spacecraft angle, and applies the cannonball kernel radially outward. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `pos_primary_sun` | SVector{3, Float64} | n/a | yes | Positional argument `pos_primary_sun`. |
| in | `planet_radius_m` | Float64 | n/a | yes | Positional argument `planet_radius_m`. |
| in | `p_srp_unscaled` | Float64 | n/a | yes | Positional argument `p_srp_unscaled`. |
| in | `reflection_coefficient` | Float64 | n/a | yes | Positional argument `reflection_coefficient`. |
| in | `reference_area_m2` | Float64 | n/a | yes | Positional argument `reference_area_m2`. |
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `planet_albedo` | Float64 | n/a | yes | Positional argument `planet_albedo`. |
| in | `AU_m` | Float64 | n/a | no | Keyword argument `AU_m` (default `149_597_870_700.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `planetary_albedo_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1340-1340`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__cannonball_radiation_accel|_cannonball_radiation_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1218-1218`
- `callees` → [[dynamics.perturbations__lambert_phase_function|_lambert_phase_function]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1211-1211`
<!-- vulcan:connections:end -->

## Limitations
Uniform albedo and a phase function that ignores the terminator's actual position on a rotating planet.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1183.
