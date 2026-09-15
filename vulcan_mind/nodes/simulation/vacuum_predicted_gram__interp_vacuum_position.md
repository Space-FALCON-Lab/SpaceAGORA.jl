---
id: simulation.vacuum_predicted_gram__interp_vacuum_position
label: _interp_vacuum_position
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _interp_vacuum_position
  lines:
  - 165
  - 165
inputs:
- id: cache
  type: VacuumPredictedGRAMCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_interp_vacuum_position`.
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

# _interp_vacuum_position

## Purpose
Linearly interpolates the vacuum-predicted inertial position between cache knots so `_query_vacuum_gram_cache!` can measure how far the actual spacecraft has strayed from the drag-free prediction and decide whether the cached splines are still trustworthy.

## Design & Implementation
Signature `_interp_vacuum_position(cache::VacuumPredictedGRAMCache, t::Float64)::SVector{3,Float64}`, `@inline`. Using `n = length(cache.vac_positions)`, it forms the segment index `idx = clamp(floor(Int, (t - cache.t0)/cache.h) + 1, 1, n - 1)` and offset `x = (t - (cache.t0 + (idx-1) cache.h))/cache.h`, then returns `(1 - x) vac_positions[idx] + x vac_positions[idx+1]` as a new static vector, reading with `@inbounds`. The knot positions were produced by successive `_vacuum_rk4_step` calls at build time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | VacuumPredictedGRAMCache | n/a | yes | Positional argument `cache`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_interp_vacuum_position`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:251-251`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Chord interpolation between RK4 knots introduces a sagitta error of order `v^2 h^2 / (8 r)` relative to the true arc, which at orbital speed and ~30 s spacing is tens of metres and eats into the 5 km default deviation budget. The clamp extrapolates outside the horizon instead of signalling. No `valid` check or length check is performed, so an empty cache leads to an out-of-bounds `@inbounds` access.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 165.
