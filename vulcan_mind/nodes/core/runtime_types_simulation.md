---
id: core.runtime_types_simulation
label: Simulation
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Simulation
  lines:
  - 455
  - 455
inputs:
- id: MC_seed
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `MC_seed` (default `[]`).
- id: drag_passage
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `drag_passage` (default `[]`).
- id: solution_states
  type: Int64
  units: n/a
  required: false
  description: Field `solution_states` (default `0`).
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
  type: Simulation
  units: n/a
  description: Constructed `Simulation` (keyword constructor via @kwdef).
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

# Simulation

## Purpose
Bookkeeping about the run itself within the legacy `Solution` record: the Monte Carlo seed history, drag-passage markers and the saved state count.

## Design & Implementation
A `@kwdef mutable struct` with `MC_seed` and `drag_passage` as `Int64` vectors and `solution_states` as a single counter, all defaulting to empty or zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `MC_seed` | Vector{Int64} | n/a | no | Field `MC_seed` (default `[]`). |
| in | `drag_passage` | Vector{Int64} | n/a | no | Field `drag_passage` (default `[]`). |
| in | `solution_states` | Int64 | n/a | no | Field `solution_states` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Simulation | n/a | — | Constructed `Simulation` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:507-507`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`solution_states` is a plain integer that must be kept in step with the vector lengths of the sibling structs by hand; nothing checks consistency.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 455.
