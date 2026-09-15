---
id: dynamics.perturbations__spice_query_name
label: _spice_query_name
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _spice_query_name
  lines:
  - 58
  - 58
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  type: Any
  units: n/a
  description: 'Return value of `_spice_query_name`. Returns `key` or `key in _SPICE_FORCE_BARYCENTER_BODIES
    ? key * "_barycenter" : key`.'
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

# _spice_query_name

## Purpose
Converts a body name to the SPICE target name, appending `_barycenter` for outer planets whose planet-centre kernels are typically unavailable.

## Design & Implementation
Canonicalises the name, returns it if it already ends in `_barycenter`, and appends the suffix if the name is in `_SPICE_FORCE_BARYCENTER_BODIES`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_spice_query_name`. Returns `key` or `key in _SPICE_FORCE_BARYCENTER_BODIES ? key * "_barycenter" : key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:921-921`
- [[dynx.coupled_perturbations_srp|srp]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1268-1268`
- [[simulation.dynamics_rhs__accumulate_nbody_flat_batch_bang|_accumulate_nbody_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:720-720`
- [[simulation.dynamics_rhs__accumulate_srp_flat_batch_bang|_accumulate_srp_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:781-781`
- [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:165-165`
- [[simulation.effector_sampling_sample_third_body_ephemerides|sample_third_body_ephemerides]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:195-195`
- [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1481-1481`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1819-1819`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1766-1766`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1713-1713`

**Downstream**

- `callees` → [[dynamics.perturbations__canonical_spice_name|_canonical_spice_name]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:59-59`
<!-- vulcan:connections:end -->

## Limitations
The forced-barycentre set is a fixed constant; a user with full planet-centre kernels cannot opt out.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 58.
