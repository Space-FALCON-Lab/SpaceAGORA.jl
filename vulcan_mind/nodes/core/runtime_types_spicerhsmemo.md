---
id: core.runtime_types_spicerhsmemo
label: SpiceRhsMemo
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: SpiceRhsMemo
  lines:
  - 495
  - 495
inputs:
- id: lock
  type: ReentrantLock
  units: n/a
  required: false
  description: Field `lock` (default `ReentrantLock()`).
- id: et
  type: Float64
  units: n/a
  required: false
  description: Field `et` (default `NaN`).
- id: primary_body_name
  type: String
  units: n/a
  required: false
  description: Field `primary_body_name` (default `""`).
- id: body_positions_j2000_m
  type: Dict{String, SVector{3, Float64}}
  units: n/a
  required: false
  description: Field `body_positions_j2000_m` (default `Dict{String, SVector{3, Float64}}()`).
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
  type: SpiceRhsMemo
  units: n/a
  description: Constructed `SpiceRhsMemo` (keyword constructor via @kwdef).
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

# SpiceRhsMemo

## Purpose
A one-entry memo of body positions at the most recent ephemeris time, so multiple effectors evaluated at the same RHS time share one SPICE lookup per body.

## Design & Implementation
A `@kwdef mutable struct` with a `ReentrantLock`, the memoised `et` (initially `NaN`), the primary body name, and a `Dict` from body name to J2000 position. A query at the memoised time and primary reads the dictionary under the lock; a different time replaces the whole entry.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lock` | ReentrantLock | n/a | no | Field `lock` (default `ReentrantLock()`). |
| in | `et` | Float64 | n/a | no | Field `et` (default `NaN`). |
| in | `primary_body_name` | String | n/a | no | Field `primary_body_name` (default `""`). |
| in | `body_positions_j2000_m` | Dict{String, SVector{3, Float64}} | n/a | no | Field `body_positions_j2000_m` (default `Dict{String, SVector{3, Float64}}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpiceRhsMemo | n/a | — | Constructed `SpiceRhsMemo` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_sharedbuffers|SharedBuffers]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:739-739`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the single latest time is memoised, so the multi-stage RK evaluations of one step each miss when their stage times differ; the dictionary and lock make each hit heavier than a cache-table interpolation.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 495.
