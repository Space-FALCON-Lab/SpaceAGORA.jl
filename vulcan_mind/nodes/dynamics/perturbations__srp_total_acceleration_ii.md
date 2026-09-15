---
id: dynamics.perturbations__srp_total_acceleration_ii
label: _srp_total_acceleration_ii
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _srp_total_acceleration_ii
  lines:
  - 1315
  - 1315
inputs:
- id: model
  type: SolarRadiationPressureModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  description: Return value of `_srp_total_acceleration_ii`.
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

# _srp_total_acceleration_ii

## Purpose
Sums the enabled radiation-pressure components — direct solar, planetary albedo and planetary infrared — into one inertial acceleration.

## Design & Implementation
Returns zero if the area is zero. Adds `srp_cannonball_accel` with the solar constant pressure 4.56e-6 Pa at one AU when `direct`, `planetary_albedo_accel` when `albedo`, and `planetary_ir_accel` when `ir`, each with the model's `Cr`, `A` and the planet's radius. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SolarRadiationPressureModel | n/a | yes | Positional argument `model`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `pos_primary_sun` | SVector{3, Float64} | n/a | yes | Positional argument `pos_primary_sun`. |
| in | `mass_kg` | Float64 | n/a | yes | Positional argument `mass_kg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_srp_total_acceleration_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1365-1365`
- `callees` → [[dynamics.perturbations_gravity_backbone_kick_acceleration_ii|gravity_backbone_kick_acceleration_ii]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1387-1387`
- `callees` → [[dynamics.perturbations_planetary_albedo_accel|planetary_albedo_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1340-1340`
- `callees` → [[dynamics.perturbations_planetary_ir_accel|planetary_ir_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1353-1353`
- `callees` → [[dynamics.perturbations_srp_cannonball_accel|srp_cannonball_accel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1328-1328`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1365-1365`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1365-1365`
<!-- vulcan:connections:end -->

## Limitations
One reflectivity coefficient is used for all three sources although the spectral response to solar and infrared differs.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1315.
