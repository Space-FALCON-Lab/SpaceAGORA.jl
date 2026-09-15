---
id: core.runtime_types_savecache
label: SaveCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: SaveCache
  lines:
  - 820
  - 820
inputs:
- id: rho_cache
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `ρ_cache` (default `[]`).
- id: heat_rate_cache
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_rate_cache` (default `[]`).
- id: drag_cache
  type: Vector{SVector{3,Float64}}
  units: n/a
  required: false
  description: Field `drag_cache` (default `[]`).
- id: lift_cache
  type: Vector{SVector{3,Float64}}
  units: n/a
  required: false
  description: Field `lift_cache` (default `[]`).
- id: cross_cache
  type: Vector{SVector{3,Float64}}
  units: n/a
  required: false
  description: Field `cross_cache` (default `[]`).
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
  type: SaveCache
  units: n/a
  description: Constructed `SaveCache` (keyword constructor via @kwdef).
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

# SaveCache

## Purpose
Per-step cache of expensive derived quantities so the save callback can record them without recomputing aerodynamics.

## Design & Implementation
A `@kwdef struct` of vectors — density, heat rate, and static drag, lift and cross-force vectors — one entry per satellite at the current step, all defaulting to empty. The RHS fills them during force evaluation and the save callback reads them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rho_cache` | Vector{Float64} | n/a | no | Field `ρ_cache` (default `[]`). |
| in | `heat_rate_cache` | Vector{Float64} | n/a | no | Field `heat_rate_cache` (default `[]`). |
| in | `drag_cache` | Vector{SVector{3,Float64}} | n/a | no | Field `drag_cache` (default `[]`). |
| in | `lift_cache` | Vector{SVector{3,Float64}} | n/a | no | Field `lift_cache` (default `[]`). |
| in | `cross_cache` | Vector{SVector{3,Float64}} | n/a | no | Field `cross_cache` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SaveCache | n/a | — | Constructed `SaveCache` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It is a single-step cache with no timestamp, so a save callback that fires at a time when the RHS has not yet run for the current step reads the previous step's values without any indication.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 820.
