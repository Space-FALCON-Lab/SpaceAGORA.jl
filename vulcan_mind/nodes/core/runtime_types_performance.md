---
id: core.runtime_types_performance
label: Performance
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Performance
  lines:
  - 436
  - 436
inputs:
- id: mass
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `mass` (default `[]`).
- id: heat_rate
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `heat_rate` (default `[]`).
- id: heat_load
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `heat_load` (default `[]`).
- id: T_r
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `T_r` (default `[]`).
- id: q
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `q` (default `[]`).
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
  type: Performance
  units: n/a
  description: Constructed `Performance` (keyword constructor via @kwdef).
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

# Performance

## Purpose
The time-series container for mass and thermal performance in the legacy `Solution` record.

## Design & Implementation
A `@kwdef mutable struct` with `mass`, per-link `heat_rate` and `heat_load` as vectors of vectors, the recovery temperature `T_r` and dynamic pressure `q`, all defaulting to empty vectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mass` | Vector{Float64} | n/a | no | Field `mass` (default `[]`). |
| in | `heat_rate` | Vector{Vector{Float64}} | n/a | no | Field `heat_rate` (default `[]`). |
| in | `heat_load` | Vector{Vector{Float64}} | n/a | no | Field `heat_load` (default `[]`). |
| in | `T_r` | Vector{Float64} | n/a | no | Field `T_r` (default `[]`). |
| in | `q` | Vector{Float64} | n/a | no | Field `q` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Performance | n/a | — | Constructed `Performance` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:505-505`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Heat rate and load are stored per link at every sample with no unit annotation, so the W/cm² and J/cm² conventions used by the guidance layer must be known from context.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 436.
