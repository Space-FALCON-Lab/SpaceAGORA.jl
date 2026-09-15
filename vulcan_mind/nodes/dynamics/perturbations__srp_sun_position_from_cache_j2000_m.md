---
id: dynamics.perturbations__srp_sun_position_from_cache_j2000_m
label: _srp_sun_position_from_cache_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _srp_sun_position_from_cache_j2000_m
  lines:
  - 1443
  - 1443
inputs:
- id: cache
  type: SRPSunEphemerisCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  description: Return value of `_srp_sun_position_from_cache_j2000_m`.
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

# _srp_sun_position_from_cache_j2000_m

## Purpose
Looks up the Sun's position in the SRP ephemeris table with the same interpolation and coverage rules as the N-body cache.

## Design & Implementation
Returns `nothing` for fewer than two samples or a time outside the table, the last sample at the end, and Catmull-Rom or linear interpolation inside depending on neighbour availability. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | SRPSunEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_srp_sun_position_from_cache_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__accumulate_srp_flat_batch_bang|_accumulate_srp_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:790-790`
- [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:171-171`

**Downstream**

- `callees` → [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1471-1471`
- `callees` → [[dynamics.perturbations__interp_vec3_linear|_interp_vec3_linear]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1467-1467`
<!-- vulcan:connections:end -->

## Limitations
Shares the silent-fallback behaviour of the N-body lookup.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1443.
