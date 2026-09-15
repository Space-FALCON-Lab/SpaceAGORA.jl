---
id: dynamics.perturbations__cannonball_radiation_accel
label: _cannonball_radiation_accel
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _cannonball_radiation_accel
  lines:
  - 1148
  - 1148
inputs:
- id: direction_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `direction_ii`.
- id: pressure_n_m2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `pressure_n_m2`.
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
  description: Return value of `_cannonball_radiation_accel`.
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

# _cannonball_radiation_accel

## Purpose
The shared cannonball radiation-pressure kernel: converts a pressure, reflectivity, area and mass into an acceleration along a given direction, used by direct SRP, albedo and infrared.

## Theory & Math
$$
\vec{a} = \frac{C_r A P}{m}\,\hat{d}
$$

## Design & Implementation
Returns zero unless mass, area, pressure and the direction magnitude are finite and positive and the reflection coefficient is finite and non-negative. Otherwise scales the direction by `Cr A P / (m |d|)` so the result is the unit direction times the cannonball acceleration. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `direction_ii` | SVector{3, Float64} | n/a | yes | Positional argument `direction_ii`. |
| in | `pressure_n_m2` | Float64 | n/a | yes | Positional argument `pressure_n_m2`. |
| in | `reflection_coefficient` | Float64 | n/a | yes | Positional argument `reflection_coefficient`. |
| in | `reference_area_m2` | Float64 | n/a | yes | Positional argument `reference_area_m2`. |
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_cannonball_radiation_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_planetary_albedo_accel|planetary_albedo_accel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1218-1218`
- [[dynamics.perturbations_planetary_ir_accel|planetary_ir_accel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1248-1248`
- [[dynamics.perturbations_srp_cannonball_accel|srp_cannonball_accel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1139-1139`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Cannonball means the area is attitude-independent; a flat panel's true SRP varies with incidence and this kernel cannot express that.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1148.
