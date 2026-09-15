---
id: core.runtime_types_forces
label: Forces
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Forces
  lines:
  - 444
  - 444
inputs:
- id: gravity_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `gravity_ii` (default `[[], [], []]`).
- id: drag_pp
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `drag_pp` (default `[[], [], []]`).
- id: drag_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `drag_ii` (default `[[], [], []]`).
- id: lift_pp
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `lift_pp` (default `[[], [], []]`).
- id: lift_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `lift_ii` (default `[[], [], []]`).
- id: force_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `force_ii` (default `[[], [], []]`).
- id: tau_body
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `τ_body` (default `[[], [], []]`).
- id: energy
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `energy` (default `[]`).
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
  type: Forces
  units: n/a
  description: Constructed `Forces` (keyword constructor via @kwdef).
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

# Forces

## Purpose
Legacy time-history container for force vectors and orbital energy, accumulating per-axis component lists across a simulation for post-processing.

## Design & Implementation
`@kwdef mutable struct Forces` whose vector fields `gravity_ii`, `drag_pp`, `drag_ii`, `lift_pp`, `lift_ii`, `force_ii`, and `τ_body` are each `Vector{Vector{Float64}}` initialised to three empty inner vectors (one per axis, N or N m), plus `energy::Vector{Float64}` (J/kg). Values are appended axis-by-axis from `IntermediateSolution` records. Embedded in `Solution.forces`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gravity_ii` | Vector{Vector{Float64}} | n/a | no | Field `gravity_ii` (default `[[], [], []]`). |
| in | `drag_pp` | Vector{Vector{Float64}} | n/a | no | Field `drag_pp` (default `[[], [], []]`). |
| in | `drag_ii` | Vector{Vector{Float64}} | n/a | no | Field `drag_ii` (default `[[], [], []]`). |
| in | `lift_pp` | Vector{Vector{Float64}} | n/a | no | Field `lift_pp` (default `[[], [], []]`). |
| in | `lift_ii` | Vector{Vector{Float64}} | n/a | no | Field `lift_ii` (default `[[], [], []]`). |
| in | `force_ii` | Vector{Vector{Float64}} | n/a | no | Field `force_ii` (default `[[], [], []]`). |
| in | `tau_body` | Vector{Vector{Float64}} | n/a | no | Field `τ_body` (default `[[], [], []]`). |
| in | `energy` | Vector{Float64} | n/a | no | Field `energy` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Forces | n/a | — | Constructed `Forces` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:506-506`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The axis-major layout (`[[x...], [y...], [z...]]`) is inefficient for per-sample access and easy to misuse by pushing whole vectors. No length invariant is enforced across fields. The inner default `[]` literals are typed `Vector{Float64}` only because of the outer annotation.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 444.
