---
id: simulation.setup__write_nbody_ephemeris_cache_file_bang
label: _write_nbody_ephemeris_cache_file!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _write_nbody_ephemeris_cache_file!
  lines:
  - 1612
  - 1612
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
- id: cache
  type: SimulationModel.NBodyEphemerisCache
  units: n/a
  required: true
  description: Positional argument `cache`.
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
  type: String
  units: n/a
  description: Return value of `_write_nbody_ephemeris_cache_file!`; mutates `path`
    in place.
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

# _write_nbody_ephemeris_cache_file!

## Purpose
Serializes an N-body ephemeris cache and its parameters to a file atomically for later prewarm loading.

## Design & Implementation
Builds the payload with `_nbody_ephemeris_cache_payload` and writes it through `_atomic_write_file` with `serialize`, so a crash mid-write leaves no partial file. Returns the path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `cache` | SimulationModel.NBodyEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_write_nbody_ephemeris_cache_file!`; mutates `path` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1746-1746`

**Downstream**

- `callees` → [[io.io_serialization__atomic_write_file|_atomic_write_file]] · `callers` · call · `src/simulation/engine/setup.jl:1620-1620`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_payload|_nbody_ephemeris_cache_payload]] · `callers` · call · `src/simulation/engine/setup.jl:1619-1619`
<!-- vulcan:connections:end -->

## Limitations
Julia's `Serialization` format is version-specific, so a file written by one Julia release may not load in another; no digest is recorded alongside.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1612.
