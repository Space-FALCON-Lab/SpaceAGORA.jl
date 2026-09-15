---
id: core.runtime_types__typed_nothing_vector
label: _typed_nothing_vector
kind: function
source:
  file: src/core/types/runtime_types.jl
  symbol: _typed_nothing_vector
  lines:
  - 698
  - 698
inputs:
- id: _type
  type: Type{T}
  units: n/a
  required: true
  description: Positional argument `_type`.
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
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
  description: 'Return value of `_typed_nothing_vector`. Returns `out`. Type parameters:
    `{T}`.'
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

# _typed_nothing_vector

## Purpose
Allocates a concretely typed `Vector{Union{Nothing, T}}` pre-filled with `nothing`, used by `SharedBuffers` defaults so per-satellite cache slots are type-stable rather than `Vector{Any}`.

## Design & Implementation
`@inline function _typed_nothing_vector(::Type{T}, n::Int) where {T}` creates `Vector{Union{Nothing, T}}(undef, n)` and calls `fill!(out, nothing)` before returning it. Because the element type is a small `Union`, Julia stores it inline without boxing, and downstream code can branch on `=== nothing` with full inference. It is invoked for `GramTrackCache`, `VacuumPredictedGRAMCache`, `_HarmonicsWorkspaceMap`, `NBodyScratchWorkspace`, and `AeroScratchWorkspace` slots.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `_type` | Type{T} | n/a | yes | Positional argument `_type`. |
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_typed_nothing_vector`. Returns `out`. Type parameters: `{T}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_sharedbuffers|SharedBuffers]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:722-722`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Negative `n` throws `ArgumentError` from the vector constructor. The function allocates fresh storage each call and does not accept a preallocated destination. For non-isbits `T` the `fill!` is needed because `undef` slots would otherwise be uninitialised references.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 698.
