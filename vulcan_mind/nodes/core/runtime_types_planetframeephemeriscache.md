---
id: core.runtime_types_planetframeephemeriscache
label: PlanetFrameEphemerisCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: PlanetFrameEphemerisCache
  lines:
  - 481
  - 481
inputs:
- id: ets
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ets`.
- id: quaternions
  type: Vector{SVector{4, Float64}}
  units: n/a
  required: true
  description: Field `quaternions`.
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
  type: PlanetFrameEphemerisCache
  units: n/a
  description: Constructed `PlanetFrameEphemerisCache`.
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

# PlanetFrameEphemerisCache

## Purpose
A precomputed table of planet orientation quaternions over the mission, so the inertial-to-planet-fixed rotation can be interpolated instead of fetched from SPICE per call.

## Design & Implementation
An immutable struct holding a sorted vector of ephemeris times and a parallel vector of scalar-last quaternions. Storing quaternions rather than matrices allows spherical interpolation between entries and halves the memory per sample.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ets` | Vector{Float64} | n/a | yes | Field `ets`. |
| in | `quaternions` | Vector{SVector{4, Float64}} | n/a | yes | Field `quaternions`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameEphemerisCache | n/a | — | Constructed `PlanetFrameEphemerisCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1892-1892`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Interpolating orientation is only accurate when the grid spacing is short relative to the planet's rotation period; the builder chooses the spacing, and this struct cannot express a variable step.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 481.
