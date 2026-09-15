---
id: core.runtime_types_gramtrackcache
label: GramTrackCache
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: GramTrackCache
  lines:
  - 511
  - 511
inputs:
- id: valid
  type: Bool
  units: n/a
  required: true
  description: Field `valid`.
- id: t0
  type: Float64
  units: n/a
  required: true
  description: Field `t0`.
- id: t1
  type: Float64
  units: n/a
  required: true
  description: Field `t1`.
- id: index_hint
  type: Int
  units: n/a
  required: true
  description: Field `index_hint`.
- id: times
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `times`.
- id: alts
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `alts`.
- id: lats
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `lats`.
- id: lons
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `lons`.
- id: rhos
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `rhos`.
- id: Ts
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `Ts`.
- id: winds
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `winds`.
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
  type: GramTrackCache
  units: n/a
  description: Constructed `GramTrackCache`.
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

# GramTrackCache

## Purpose
Per-satellite along-track cache of GRAM atmosphere samples, letting density lookups interpolate previously computed density, temperature, and wind over a short time horizon instead of calling the expensive GRAM model at every integrator stage.

## Design & Implementation
`mutable struct GramTrackCache` with `valid::Bool`, the validity window `t0`, `t1` (s), an `index_hint::Int` for monotone search reuse, and parallel knot vectors `times` (s), `alts` (m), `lats`, `lons` (rad), `rhos` (kg/m^3), `Ts` (K), and `winds::Vector{SVector{3,Float64}}` (m/s). Tolerances that decide whether a query is close enough to the cached track come from `GramTrackCacheConfig`. Slots are stored in `SharedBuffers.gram_density_cache[i]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `valid` | Bool | n/a | yes | Field `valid`. |
| in | `t0` | Float64 | n/a | yes | Field `t0`. |
| in | `t1` | Float64 | n/a | yes | Field `t1`. |
| in | `index_hint` | Int | n/a | yes | Field `index_hint`. |
| in | `times` | Vector{Float64} | n/a | yes | Field `times`. |
| in | `alts` | Vector{Float64} | n/a | yes | Field `alts`. |
| in | `lats` | Vector{Float64} | n/a | yes | Field `lats`. |
| in | `lons` | Vector{Float64} | n/a | yes | Field `lons`. |
| in | `rhos` | Vector{Float64} | n/a | yes | Field `rhos`. |
| in | `Ts` | Vector{Float64} | n/a | yes | Field `Ts`. |
| in | `winds` | Vector{SVector{3, Float64}} | n/a | yes | Field `winds`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GramTrackCache | n/a | — | Constructed `GramTrackCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.model_selection__gram_density_cache_for_sat_bang|_gram_density_cache_for_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:122-122`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validity depends on the actual trajectory staying within altitude and angular tolerances of the predicted track; the struct itself carries no measure of prediction error. `index_hint` assumes monotonically increasing query times and can degrade to a linear scan when the integrator rejects steps. Mutable and unlocked, so it must only be touched by the owning satellite's callback.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 511.
