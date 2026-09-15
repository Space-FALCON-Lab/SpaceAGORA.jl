---
id: dynamics.perturbations_srp_cannonball_accel
label: srp_cannonball_accel
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: srp_cannonball_accel
  lines:
  - 1111
  - 1111
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
  description: Return value of `srp_cannonball_accel`.
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

# srp_cannonball_accel

## Purpose
Direct solar radiation pressure acceleration with inverse-square distance scaling and eclipse shadowing.

## Theory & Math
$$
P = P_{1\,\text{AU}} \left(\frac{\text{AU}}{r_{\odot}}\right)^2 \nu,\qquad \vec{a} = \frac{C_r A P}{m}\hat{r}_{\odot \to sc}
$$

## Design & Implementation
Validates mass, area and reflectivity, forms the Sun-to-spacecraft vector, computes the eclipse fraction, scales the one-AU pressure by `(AU/r)²`, and applies the cannonball kernel along the Sun-to-spacecraft direction. `@inline`.

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
| in | `AU_m` | Float64 | n/a | no | Keyword argument `AU_m` (default `149_597_870_700.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `srp_cannonball_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1328-1328`
- [[dynx.coupled_perturbations_srp|srp]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1271-1271`

**Downstream**

- `callees` → [[dynamics.perturbations__cannonball_radiation_accel|_cannonball_radiation_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1139-1139`
- `callees` → [[dynamics.perturbations_eclipse_area_calc|eclipse_area_calc]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1137-1137`
<!-- vulcan:connections:end -->

## Limitations
Eclipse is a scalar fraction, so penumbral pressure is reduced uniformly rather than directionally.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1111.
