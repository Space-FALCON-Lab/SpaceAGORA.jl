---
id: core.runtime_types_nbodyscratchworkspace
label: NBodyScratchWorkspace
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: NBodyScratchWorkspace
  lines:
  - 555
  - 555
inputs:
- id: pos_primary_k_all
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `pos_primary_k_all`.
- id: body_force_ii
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `body_force_ii`.
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
  type: NBodyScratchWorkspace
  units: n/a
  description: Constructed `NBodyScratchWorkspace`.
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

# NBodyScratchWorkspace

## Purpose
Preallocated per-satellite scratch space for the N-body gravity effector so its per-body loop allocates nothing on the RHS path.

## Design & Implementation
An immutable struct of two vectors of static three-vectors: `pos_primary_k_all` for each body's position relative to the primary, and `body_force_ii` for each body's acceleration contribution, summed in body order so the result is deterministic across threads.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_primary_k_all` | Vector{SVector{3, Float64}} | n/a | yes | Field `pos_primary_k_all`. |
| in | `body_force_ii` | Vector{SVector{3, Float64}} | n/a | yes | Field `body_force_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyScratchWorkspace | n/a | — | Constructed `NBodyScratchWorkspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__make_nbody_scratch_workspace|_make_nbody_scratch_workspace]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:97-97`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Sized to the configured body count at construction; changing the body list mid-run would silently index out of range.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 555.
