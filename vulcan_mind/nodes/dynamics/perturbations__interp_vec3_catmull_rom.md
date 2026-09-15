---
id: dynamics.perturbations__interp_vec3_catmull_rom
label: _interp_vec3_catmull_rom
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _interp_vec3_catmull_rom
  lines:
  - 1482
  - 1482
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
- id: p2
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p2`.
- id: p3
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p3`.
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
  description: Return value of `_interp_vec3_catmull_rom`.
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

# _interp_vec3_catmull_rom

## Purpose
Interpolates a position between two ephemeris samples with a Catmull-Rom cubic using the two neighbours on either side, giving smooth velocity across the table.

## Theory & Math
$$
p(\tau) = \tfrac{1}{2}\left[ 2p_1 + (-p_0 + p_2)\tau + (2p_0 - 5p_1 + 4p_2 - p_3)\tau^2 + (-p_0 + 3p_1 - 3p_2 + p_3)\tau^3 \right]
$$

## Design & Implementation
Evaluates the standard Catmull-Rom blend of `p0` through `p3` at parameter `tau` in the unit interval. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p0` | SVector{3, Float64} | n/a | yes | Positional argument `p0`. |
| in | `p1` | SVector{3, Float64} | n/a | yes | Positional argument `p1`. |
| in | `p2` | SVector{3, Float64} | n/a | yes | Positional argument `p2`. |
| in | `p3` | SVector{3, Float64} | n/a | yes | Positional argument `p3`. |
| in | `tau` | Float64 | n/a | yes | Positional argument `tau`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_interp_vec3_catmull_rom`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1108-1108`
- [[dynamics.perturbations__srp_sun_position_from_cache_j2000_m|_srp_sun_position_from_cache_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1471-1471`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · feedback · `src/dynamics/coupled/perturbations.jl:1500-1500`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
<!-- vulcan:connections:end -->

## Limitations
Assumes uniform sample spacing; the clamped last interval of the ephemeris tables is slightly non-uniform, introducing a small error there.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1482.
