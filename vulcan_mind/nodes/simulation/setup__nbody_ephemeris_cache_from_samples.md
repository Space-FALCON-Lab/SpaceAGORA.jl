---
id: simulation.setup__nbody_ephemeris_cache_from_samples
label: _nbody_ephemeris_cache_from_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_from_samples
  lines:
  - 1525
  - 1525
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
- id: ets
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `ets`.
- id: positions_j2000_m
  type: Matrix{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `positions_j2000_m`.
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
  description: Return value of `_nbody_ephemeris_cache_from_samples`.
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

# _nbody_ephemeris_cache_from_samples

## Purpose
Assembles an `NBodyEphemerisCache` from its primary name, body names, time vector and position matrix, adding the derived index map so both the builder and the file loader construct caches identically.

## Design & Implementation
Calls `_nbody_ephemeris_body_index_by_name` on the body names and passes all five fields to the `NBodyEphemerisCache` constructor. Keeping the index-map derivation here means neither the SPICE-sampling builder nor the deserialiser can forget it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `body_query_names` | Vector{String} | n/a | yes | Positional argument `body_query_names`. |
| in | `ets` | Vector{Float64} | n/a | yes | Positional argument `ets`. |
| in | `positions_j2000_m` | Matrix{SVector{3, Float64}} | n/a | yes | Positional argument `positions_j2000_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.NBodyEphemerisCache | n/a | — | Return value of `_nbody_ephemeris_cache_from_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__cache_from_nbody_ephemeris_payload|_cache_from_nbody_ephemeris_payload]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1655-1655`
- [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1572-1572`

**Downstream**

- `callees` → [[core.runtime_types_nbodyephemeriscache|NBodyEphemerisCache]] · `callers` · call · `src/simulation/engine/setup.jl:1532-1532`
- `callees` → [[simulation.setup__nbody_ephemeris_body_index_by_name|_nbody_ephemeris_body_index_by_name]] · `callers` · call · `src/simulation/engine/setup.jl:1531-1531`
<!-- vulcan:connections:end -->

## Limitations
It performs no validation that the matrix has one row per time and one column per body; the file loader validates before calling and the builder constructs consistent inputs by construction.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1525.
