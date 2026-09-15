---
id: simulation.persistence__append_saved_segment_bang
label: _append_saved_segment!
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _append_saved_segment!
  lines:
  - 10
  - 10
inputs:
- id: times_acc
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times_acc`.
- id: data_acc
  type: Vector
  units: n/a
  required: true
  description: Positional argument `data_acc`.
- id: saved_values
  type: Any
  units: n/a
  required: true
  description: Positional argument `saved_values`.
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
  description: Return value of `_append_saved_segment!`; mutates `times_acc` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _append_saved_segment!

## Purpose
Appends one solver segment's saved output onto the run-level time and data accumulators, stitching together the pieces produced when propagation is split across several `solve` calls.

## Design & Implementation
Forwards to `SimulationModel.IOOutputs._append_saved_segment!(times_acc, data_acc, saved_values)`, which mutates both `times_acc::Vector{Float64}` and `data_acc::Vector` in place. An empty segment returns `nothing` immediately. If `times_acc` is non-empty and its last entry equals `saved_values.t[1]` exactly, under `isapprox` with `atol=0.0` and `rtol=0.0`, the start index advances to 2 so the shared boundary sample is not duplicated. The remaining range is appended through `@view` slices of `saved_values.t` and `saved_values.saveval`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `times_acc` | Vector{Float64} | n/a | yes | Positional argument `times_acc`. |
| in | `data_acc` | Vector | n/a | yes | Positional argument `data_acc`. |
| in | `saved_values` | Any | n/a | yes | Positional argument `saved_values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_append_saved_segment!`; mutates `times_acc` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`
- [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:92-92`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Deduplication uses bit-exact equality, so a restart whose first sample differs by one unit in the last place leaves a duplicated timestamp in the results. Only the single leading sample is examined; a segment that re-saves several overlapping points is appended whole. Nothing enforces that times are monotonically increasing, and the two accumulators can be left at unequal length if the caller mutates only one of them elsewhere.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 10.
