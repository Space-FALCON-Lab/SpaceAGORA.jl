---
id: simulation.save_fields__save_drag
label: _save_drag
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_drag
  lines:
  - 33
  - 33
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `_save_drag`. Returns `drags`.
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

# _save_drag

## Purpose
Save-time getter for the aerodynamic drag force on each spacecraft, read out of the drag cache that the dynamics right-hand side populated during the step rather than recomputed.

## Design & Implementation
Marked `@inline`. It reaches `drag_cache = integrator.p.save_cache.drag_cache` and fills a `Vector{SVector{3, Float64}}` of length `num_sats`. Each entry is guarded by `i <= length(drag_cache)`, substituting the zero vector `SVector{3, Float64}(0.0, 0.0, 0.0)` when the cache is shorter than the constellation, which keeps the save path from throwing when a spacecraft never entered the atmospheric branch of the right-hand side.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_drag`. Returns `drags`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:180-180`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The values are stale by construction: they reflect the most recent right-hand-side evaluation, which for an adaptive solver may belong to a rejected stage rather than the accepted step now being saved. The zero-vector fallback is indistinguishable from genuinely zero drag in vacuum, so a truncated cache silently reads as a drag-free spacecraft. Nothing verifies that the cache was written at the current time.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 33.
