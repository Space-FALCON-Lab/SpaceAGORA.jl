---
id: dynamics.perturbations__magnetic_field_inertial
label: _magnetic_field_inertial
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _magnetic_field_inertial
  lines:
  - 2004
  - 2004
inputs:
- id: model
  type: Union{MagneticTorqueRodModel, EddyCurrentDampingModel}
  units: n/a
  required: true
  description: Positional argument `model`.
- id: l_pi
  type: SMatrix{3, 3, Float64, 9}
  units: n/a
  required: true
  description: Positional argument `l_pi`.
- id: pos_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_pp`.
- id: lat_rad
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat_rad`.
- id: lon_rad
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon_rad`.
- id: alt_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt_m`.
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
  description: Return value of `_magnetic_field_inertial`.
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

# _magnetic_field_inertial

## Purpose
Returns the local magnetic field in the inertial frame in tesla, from either the tilted dipole or IGRF, so torque-rod and eddy-damping models share one unit contract.

## Design & Implementation
For `:igrf` it calls `igrf` at the model's epoch year and geodetic coordinates, converts the NED result to planet-fixed with `ned_to_ecef`, rotates to inertial through `l_pi'` and scales from nanotesla. Otherwise it calls `get_magnetic_field_dipole`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | Union{MagneticTorqueRodModel, EddyCurrentDampingModel} | n/a | yes | Positional argument `model`. |
| in | `l_pi` | SMatrix{3, 3, Float64, 9} | n/a | yes | Positional argument `l_pi`. |
| in | `pos_pp` | SVector{3, Float64} | n/a | yes | Positional argument `pos_pp`. |
| in | `lat_rad` | Float64 | n/a | yes | Positional argument `lat_rad`. |
| in | `lon_rad` | Float64 | n/a | yes | Positional argument `lon_rad`. |
| in | `alt_m` | Float64 | n/a | yes | Positional argument `alt_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_magnetic_field_inertial`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2056-2056`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1911-1911`
- [[dynamics.perturbations_magnetictorquerodmodel|MagneticTorqueRodModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1998-1998`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2020-2020`
<!-- vulcan:connections:end -->

## Limitations
IGRF is Earth-only and the epoch is a single fixed year per model; the dipole branch uses Earth constants regardless of planet.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2004.
