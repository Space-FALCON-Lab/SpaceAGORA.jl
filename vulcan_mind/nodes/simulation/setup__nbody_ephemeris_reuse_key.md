---
id: simulation.setup__nbody_ephemeris_reuse_key
label: _nbody_ephemeris_reuse_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_reuse_key
  lines:
  - 279
  - 279
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
  type: NBodyEphemerisReuseKey
  units: n/a
  description: Return value of `_nbody_ephemeris_reuse_key`.
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

# _nbody_ephemeris_reuse_key

## Purpose
Constructs the reuse-dictionary key for an `NBodyEphemerisCache`, combining primary body, ordered third-body names, and the quantised time grid.

## Design & Implementation
Signature `(primary_body_name::String, body_query_names::Vector{String}, et_start::Float64, mission_end_s::Float64, dt_s::Float64)::NBodyEphemerisReuseKey`. Returns `(primary_body_name, _body_query_names_reuse_key(body_query_names), _cache_time_key(et_start), _cache_time_key(mission_end_s), _cache_time_key(dt_s))`, matching `NBodyEphemerisReuseKey = Tuple{String, String, Int64, Int64, Int64}`. The same key type indexes both `_NBODY_EPHEMERIS_REUSE_CACHE` and `_NBODY_EPHEMERIS_PREWARMED_CACHE`.

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
| out | `result` | NBodyEphemerisReuseKey | n/a | — | Return value of `_nbody_ephemeris_reuse_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1826-1826`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1719-1719`
- [[simulation.setup__prewarmed_nbody_ephemeris_lookup|_prewarmed_nbody_ephemeris_lookup]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1576-1576`
- [[simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang|_register_prewarmed_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1589-1589`

**Downstream**

- `callees` → [[simulation.setup__body_query_names_reuse_key|_body_query_names_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:282-282`
- `callees` → [[simulation.setup__cache_time_key|_cache_time_key]] · `callers` · call · `src/simulation/engine/setup.jl:283-283`
<!-- vulcan:connections:end -->

## Limitations
Body-name order and case are significant, so equivalent but differently ordered configurations do not share caches. The ephemerides model is not part of the key, which is safe only because `_validate_ephemerides_support!` forbids N-body effectors with the simple model. The key ignores `max_samples`, so a cache built under a tighter cap could be reused by a run expecting finer sampling.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 279.
