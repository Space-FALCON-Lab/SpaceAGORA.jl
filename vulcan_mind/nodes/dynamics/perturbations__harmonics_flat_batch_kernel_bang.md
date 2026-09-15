---
id: dynamics.perturbations__harmonics_flat_batch_kernel_bang
label: _harmonics_flat_batch_kernel!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_flat_batch_kernel!
  lines:
  - 376
  - 376
inputs:
- id: totals
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `totals`.
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
- id: work_items
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `work_items`.
- id: item_start
  type: Int
  units: n/a
  required: true
  description: Positional argument `item_start`.
- id: item_end
  type: Int
  units: n/a
  required: true
  description: Positional argument `item_end`.
- id: lpi
  type: SMatrix{3, 3, Float64, 9}
  units: n/a
  required: true
  description: Positional argument `lpi`.
- id: ws
  type: HarmonicsBatchWorkspace
  units: n/a
  required: true
  description: Positional argument `ws`.
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
  type: Nothing
  units: n/a
  description: Return value of `_harmonics_flat_batch_kernel!`; mutates `totals` in
    place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _harmonics_flat_batch_kernel!

## Purpose
The batched harmonics kernel: evaluates the same Pines recurrence for a block of satellites at once, with the inner loop across satellites so the coefficient reads are shared.

## Design & Implementation
For work items `item_start` through `item_end`, it stages each satellite's planet-frame position, direction cosines, inverse radius and mass into the batch workspace, then runs the recurrences with the satellite index as the innermost, contiguous dimension of the three-dimensional `A` array. Results are written into the `totals` matrix. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `totals` | Matrix{Float64} | n/a | yes | Positional argument `totals`. |
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `item_start` | Int | n/a | yes | Positional argument `item_start`. |
| in | `item_end` | Int | n/a | yes | Positional argument `item_end`. |
| in | `lpi` | SMatrix{3, 3, Float64, 9} | n/a | yes | Positional argument `lpi`. |
| in | `ws` | HarmonicsBatchWorkspace | n/a | yes | Positional argument `ws`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_harmonics_flat_batch_kernel!`; mutates `totals` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:948-948`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:508-508`
<!-- vulcan:connections:end -->

## Limitations
The batch size is fixed by the workspace allocation, so a work range larger than it would overflow; the caller partitions accordingly.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 376.
