---
id: dynamics.perturbations__nbody_body_position_from_cache_j2000_m
label: _nbody_body_position_from_cache_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_body_position_from_cache_j2000_m
  lines:
  - 1071
  - 1071
inputs:
- id: cache
  type: NBodyEphemerisCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: body_name_spice
  type: String
  units: n/a
  required: true
  description: Positional argument `body_name_spice`.
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_nbody_body_position_from_cache_j2000_m`.
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

# _nbody_body_position_from_cache_j2000_m

## Purpose
Looks up a third body's position in the ephemeris table at a time, interpolating between samples and returning `nothing` when the cache does not cover the request.

## Design & Implementation
Returns `nothing` on primary mismatch, unknown body, fewer than two samples, or a time outside the table. Finds the bracketing index, returns the last sample at the end, and otherwise interpolates with Catmull-Rom when two neighbours exist on each side, falling back to linear near the ends. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | NBodyEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `body_name_spice` | String | n/a | yes | Positional argument `body_name_spice`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_nbody_body_position_from_cache_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:936-936`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__accumulate_nbody_flat_batch_bang|_accumulate_nbody_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:731-731`
- [[simulation.effector_sampling_sample_third_body_ephemerides|sample_third_body_ephemerides]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:203-203`

**Downstream**

- `callees` → [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1108-1108`
- `callees` → [[dynamics.perturbations__interp_vec3_linear|_interp_vec3_linear]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1104-1104`
<!-- vulcan:connections:end -->

## Limitations
A `nothing` return triggers a live SPICE call in the caller, so a table that ends just before the mission end degrades silently for its last seconds.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1071.
