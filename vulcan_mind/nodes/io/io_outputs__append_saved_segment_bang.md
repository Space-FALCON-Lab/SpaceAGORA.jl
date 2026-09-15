---
id: io.io_outputs__append_saved_segment_bang
label: _append_saved_segment!
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _append_saved_segment!
  lines:
  - 12
  - 12
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
  type: Nothing
  units: n/a
  description: Return value of `_append_saved_segment!`; mutates `times_acc` in place.
    Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# _append_saved_segment!

## Purpose
`_append_saved_segment!` concatenates one solver segment's `SavedValues` (from a `SavingCallback`) onto running accumulators of times and snapshots, de-duplicating the boundary sample that adjacent integration segments share. It lets the engine run a mission as several solver calls (for example across maneuvers or discontinuities) while producing one continuous time series.

## Design & Implementation
Signature `(times_acc::Vector{Float64}, data_acc::Vector, saved_values)`. It reads `seg_len = length(saved_values.t)` and returns immediately when the segment is empty. If `times_acc` already has entries and its last time is exactly equal (`isapprox` with `atol=0.0, rtol=0.0`, i.e. bitwise-equal floats) to `saved_values.t[1]`, the segment is appended from index 2 so the shared boundary sample is kept only once. Both `saved_values.t` and `saved_values.saveval` are appended via `@view` slices with `append!`, mutating the two accumulator vectors in place. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `times_acc` | Vector{Float64} | n/a | yes | Positional argument `times_acc`. |
| in | `data_acc` | Vector | n/a | yes | Positional argument `data_acc`. |
| in | `saved_values` | Any | n/a | yes | Positional argument `saved_values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_saved_segment!`; mutates `times_acc` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`
- [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:92-92`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Duplicate detection requires exact floating-point equality, so a boundary that differs by one ULP (common after a solver re-initialises at `t + eps`) produces a duplicated timestamp. Only the first sample of the new segment is checked; a segment that overlaps the accumulator by more than one sample is appended with all overlaps. No check enforces that `saved_values.t` is monotonically increasing or that `length(saveval) == length(t)`, so mismatched segments corrupt the accumulators silently. `data_acc` is an untyped `Vector`, so element types are not validated.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 12.
