---
id: core.runtime_types_srpsunephemeriscache
label: SRPSunEphemerisCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: SRPSunEphemerisCache
  lines:
  - 468
  - 468
inputs:
- id: ets
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ets`.
- id: positions_j2000_m
  type: Vector{SVector{3, Float64}}
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
  type: SRPSunEphemerisCache
  units: n/a
  description: Constructed `SRPSunEphemerisCache`.
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

# SRPSunEphemerisCache

## Purpose
A precomputed table of Sun position relative to the primary body over the mission, so solar radiation pressure need not query SPICE on every RHS evaluation.

## Design & Implementation
An immutable struct holding a sorted vector of ephemeris times and a parallel vector of J2000 positions in metres. Built once at setup and installed in `SharedBuffers.srp_sun_ephemeris_cache`; queries interpolate between neighbouring entries and fall back to a live lookup outside the covered span.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ets` | Vector{Float64} | n/a | yes | Field `ets`. |
| in | `positions_j2000_m` | Vector{SVector{3, Float64}} | n/a | yes | Field `positions_j2000_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SRPSunEphemerisCache | n/a | — | Constructed `SRPSunEphemerisCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1787-1787`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Coverage and resolution are fixed at build time, so a solve that runs past the planned end time or needs finer temporal resolution silently falls back to per-call SPICE queries with the associated lock contention.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 468.
