---
id: core.runtime_types_nbodyephemeriscache
label: NBodyEphemerisCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: NBodyEphemerisCache
  lines:
  - 473
  - 473
inputs:
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Field `primary_body_name`.
- id: body_query_names
  type: Vector{String}
  units: n/a
  required: true
  description: Field `body_query_names`.
- id: body_index_by_name
  type: Dict{String, Int}
  units: n/a
  required: true
  description: Field `body_index_by_name`.
- id: ets
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ets`.
- id: positions_j2000_m
  type: Matrix{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `positions_j2000_m`.
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
  type: NBodyEphemerisCache
  units: n/a
  description: Constructed `NBodyEphemerisCache`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# NBodyEphemerisCache

## Purpose
A precomputed table of third-body positions relative to the primary over the mission, serving the N-body gravity effector without per-call SPICE queries.

## Design & Implementation
An immutable struct with the primary body name, the list of body query names, a `Dict` mapping name to column index, a sorted vector of ephemeris times and a matrix of J2000 positions indexed by time and body. Installed in `SharedBuffers.nbody_ephemeris_cache` at setup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `primary_body_name` | String | n/a | yes | Field `primary_body_name`. |
| in | `body_query_names` | Vector{String} | n/a | yes | Field `body_query_names`. |
| in | `body_index_by_name` | Dict{String, Int} | n/a | yes | Field `body_index_by_name`. |
| in | `ets` | Vector{Float64} | n/a | yes | Field `ets`. |
| in | `positions_j2000_m` | Matrix{SVector{3, Float64}} | n/a | yes | Field `positions_j2000_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyEphemerisCache | n/a | — | Constructed `NBodyEphemerisCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.setup__nbody_ephemeris_cache_from_samples|_nbody_ephemeris_cache_from_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1532-1532`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every configured body is tabulated on the same time grid, so a fast-moving moon and a slow outer planet share one resolution; the name lookup goes through a `Dict` on each query.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 473.
