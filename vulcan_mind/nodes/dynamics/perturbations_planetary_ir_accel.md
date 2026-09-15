---
id: dynamics.perturbations_planetary_ir_accel
label: planetary_ir_accel
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: planetary_ir_accel
  lines:
  - 1227
  - 1227
inputs:
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: planet_radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `planet_radius_m`.
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
- id: planet_ir_flux_w_m2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `planet_ir_flux_w_m2`.
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
  description: Return value of `planetary_ir_accel`.
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

# planetary_ir_accel

## Purpose
Acceleration from the planet's thermal infrared emission, treated as a uniform isotropic radiator and applied radially outward through the cannonball kernel.

## Theory & Math
$$
P_{IR} = \frac{\Phi_{IR}}{c}\left(\frac{R}{r}\right)^2
$$

## Design & Implementation
Returns zero for a non-positive radius or flux, or a spacecraft range not exceeding the radius. Otherwise the pressure is `flux / c` scaled by `(R / r)²` for the solid-angle falloff, and `_cannonball_radiation_accel` is applied along the spacecraft position vector. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `planet_radius_m` | Float64 | n/a | yes | Positional argument `planet_radius_m`. |
| in | `reflection_coefficient` | Float64 | n/a | yes | Positional argument `reflection_coefficient`. |
| in | `reference_area_m2` | Float64 | n/a | yes | Positional argument `reference_area_m2`. |
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `planet_ir_flux_w_m2` | Float64 | n/a | yes | Positional argument `planet_ir_flux_w_m2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `planetary_ir_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1353-1353`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__cannonball_radiation_accel|_cannonball_radiation_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1248-1248`
<!-- vulcan:connections:end -->

## Limitations
No day-night asymmetry and no latitude dependence; a real planet's infrared emission is stronger on the day side and varies with surface temperature.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1227.
