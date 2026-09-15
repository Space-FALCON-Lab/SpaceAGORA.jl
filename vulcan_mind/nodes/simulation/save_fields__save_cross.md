---
id: simulation.save_fields__save_cross
label: _save_cross
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_cross
  lines:
  - 53
  - 53
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
  description: Return value of `_save_cross`. Returns `crosses`.
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

# _save_cross

## Purpose
Save-time getter for the cross-flow, or side, component of the aerodynamic force on each spacecraft, completing the drag, lift and cross decomposition in the saved output.

## Design & Implementation
Marked `@inline`, reading `integrator.p.save_cache.cross_cache` and producing a `Vector{SVector{3, Float64}}` of length `num_sats`. As with the other force getters, the per-index access is guarded by `i <= length(cross_cache)` and falls back to the zero `SVector{3, Float64}`. Keeping this third component cached and saved separately lets downstream analysis verify that the three saved vectors sum to the total aerodynamic force the dynamics actually applied.

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
| out | `result` | Any | n/a | — | Return value of `_save_cross`. Returns `crosses`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:182-182`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Cross-flow force is exactly zero for symmetric-attitude cases, so the zero fallback used for a short cache is especially easy to mistake for a real result. The vector is a snapshot of the last right-hand-side evaluation, not an integral or average over the step. No consistency check ties the cross cache's length to the drag and lift caches, so the three can disagree on how many spacecraft they cover.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 53.
