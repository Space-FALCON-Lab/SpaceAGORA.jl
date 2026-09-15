---
id: simulation.setup__prewarmed_nbody_ephemeris_lookup
label: _prewarmed_nbody_ephemeris_lookup
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _prewarmed_nbody_ephemeris_lookup
  lines:
  - 1575
  - 1575
inputs:
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
- id: body_query_names
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `body_query_names`.
- id: et_start
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et_start`.
- id: mission_end_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_end_s`.
- id: dt_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_s`.
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
  description: Return value of `_prewarmed_nbody_ephemeris_lookup`. Returns `_ephemeris_reuse_lookup(_NBODY_EPHEMERIS_PREWARMED_CACHE,
    reuse_key)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _prewarmed_nbody_ephemeris_lookup

## Purpose
Looks up a prewarmed N-body table by its build parameters.

## Design & Implementation
Forms the reuse key from primary, bodies, epoch, end and step and queries `_NBODY_EPHEMERIS_PREWARMED_CACHE` via `_ephemeris_reuse_lookup`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `body_query_names` | Vector{String} | n/a | yes | Positional argument `body_query_names`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_prewarmed_nbody_ephemeris_lookup`. Returns `_ephemeris_reuse_lookup(_NBODY_EPHEMERIS_PREWARMED_CACHE, reuse_key)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1820-1820`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1715-1715`

**Downstream**

- `callees` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1577-1577`
- `callees` → [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1576-1576`
<!-- vulcan:connections:end -->

## Limitations
Prewarmed and reuse registries are separate dictionaries with the same key type, so a table built by a run is not found by a later prewarm lookup and vice versa without explicit registration.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1575.
