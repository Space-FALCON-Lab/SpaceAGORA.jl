---
id: simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang
label: _register_prewarmed_nbody_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _register_prewarmed_nbody_ephemeris_cache!
  lines:
  - 1580
  - 1580
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
- id: cache
  type: SimulationModel.NBodyEphemerisCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: replace
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `replace` (default `false`).
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
  type: SimulationModel.NBodyEphemerisCache
  units: n/a
  description: Return value of `_register_prewarmed_nbody_ephemeris_cache!`; mutates
    `primary_body_name` in place.
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

# _register_prewarmed_nbody_ephemeris_cache!

## Purpose
Stores an N-body table in the prewarmed registry under its parameter key, optionally replacing an existing entry.

## Design & Implementation
Forms the key and calls `_ephemeris_explicit_cache_store!` with the `replace` flag. Returns the stored cache. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `body_query_names` | Vector{String} | n/a | yes | Positional argument `body_query_names`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `cache` | SimulationModel.NBodyEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `replace` | Bool | n/a | no | Keyword argument `replace` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.NBodyEphemerisCache | n/a | — | Return value of `_register_prewarmed_nbody_ephemeris_cache!`; mutates `primary_body_name` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__load_nbody_ephemeris_cache_bang|_load_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1669-1669`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1722-1722`

**Downstream**

- `callees` → [[simulation.setup__ephemeris_explicit_cache_store_bang|_ephemeris_explicit_cache_store!]] · `callers` · call · `src/simulation/engine/setup.jl:1590-1590`
- `callees` → [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1589-1589`
<!-- vulcan:connections:end -->

## Limitations
With `replace=false` an existing entry wins and the caller's table is discarded, which is surprising if the caller intended to refresh.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1580.
