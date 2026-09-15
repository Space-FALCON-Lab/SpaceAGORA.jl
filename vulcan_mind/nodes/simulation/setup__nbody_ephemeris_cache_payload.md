---
id: simulation.setup__nbody_ephemeris_cache_payload
label: _nbody_ephemeris_cache_payload
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_payload
  lines:
  - 1593
  - 1593
inputs:
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
  type: Any
  units: n/a
  description: Return value of `_nbody_ephemeris_cache_payload`. Returns `(`.
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

# _nbody_ephemeris_cache_payload

## Purpose
Packages an N-body cache and its build parameters into the named tuple that is serialized to disk.

## Design & Implementation
Returns a tuple with `schema_version`, a `created_utc` stamp, the primary name, copies of body names, times and positions, and `et_start`, `mission_end_s` and `dt_s` as `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | SimulationModel.NBodyEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_nbody_ephemeris_cache_payload`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__write_nbody_ephemeris_cache_file_bang|_write_nbody_ephemeris_cache_file!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1619-1619`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:1604-1604`
<!-- vulcan:connections:end -->

## Limitations
The `created_utc` stamp makes two payloads of identical content differ byte-for-byte, so files cannot be compared by hash.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1593.
