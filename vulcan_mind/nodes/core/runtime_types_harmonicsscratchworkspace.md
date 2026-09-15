---
id: core.runtime_types_harmonicsscratchworkspace
label: HarmonicsScratchWorkspace
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: HarmonicsScratchWorkspace
  lines:
  - 560
  - 560
inputs:
- id: A
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `A`.
- id: R
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `R`.
- id: I
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `I`.
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
  type: HarmonicsScratchWorkspace
  units: n/a
  description: Constructed `HarmonicsScratchWorkspace`.
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

# HarmonicsScratchWorkspace

## Purpose
Preallocated storage for the Pines-style spherical-harmonic gravity evaluation so the recursion arrays are reused across RHS calls for one satellite and one model.

## Design & Implementation
Immutable `struct HarmonicsScratchWorkspace` with `A::Matrix{Float64}` (the derived Legendre function table sized degree+2 by order+2), and `R::Vector{Float64}` and `I::Vector{Float64}` (real and imaginary parts of the longitude recursion `(x + i y)^m`). Instances are stored per satellite in a `Dict{UInt, HarmonicsScratchWorkspace}` keyed by a model cache key inside `SharedBuffers.harmonics_workspaces[i]`, so different degree/order models do not share arrays.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `A` | Matrix{Float64} | n/a | yes | Field `A`. |
| in | `R` | Vector{Float64} | n/a | yes | Field `R`. |
| in | `I` | Vector{Float64} | n/a | yes | Field `I`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | HarmonicsScratchWorkspace | n/a | — | Constructed `HarmonicsScratchWorkspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__make_harmonics_scratch_workspace|_make_harmonics_scratch_workspace]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:248-248`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Array dimensions are fixed at construction; requesting a higher degree than allocated leads to `BoundsError` unless the caller allocates a new workspace. The struct carries no record of the degree it was sized for. No synchronisation is provided, so concurrent evaluation for the same satellite and model would race.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 560.
