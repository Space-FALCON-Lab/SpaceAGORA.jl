---
id: simulation.vacuum_predicted_gram_vacuumpredictedgramcache
label: VacuumPredictedGRAMCache
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: VacuumPredictedGRAMCache
  lines:
  - 38
  - 38
inputs:
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
  type: VacuumPredictedGRAMCache
  units: n/a
  description: Return value of `VacuumPredictedGRAMCache`. Returns `VacuumPredictedGRAMCache(`.
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

# VacuumPredictedGRAMCache

## Purpose
Per-satellite mutable cache that stores a drag-free prediction of the trajectory over a look-ahead horizon together with natural cubic spline fits of `log(ρ)` and temperature and linearly interpolated wind at uniformly spaced time knots. It lets the density callback answer repeated GRAM queries from the splines instead of calling the expensive density model at every RHS evaluation.

## Design & Implementation
The struct is declared in `src/core/types/runtime_types.jl` as `mutable struct VacuumPredictedGRAMCache` with fields `valid::Bool`, `t0`, `t1`, `h::Float64` (first and last knot times and knot spacing in seconds), `log_rhos`, `Ms_rho`, `Ts`, `Ms_T::Vector{Float64}` (knot values and spline second derivatives for log density and temperature), `winds::Vector{SVector{3,Float64}}`, `vac_alts::Vector{Float64}` and `vac_positions::Vector{SVector{3,Float64}}`. This file adds a zero-argument constructor returning an invalid cache with `t0 = t1 = 0.0`, `h = 1.0` and empty vectors. `_build_vacuum_gram_cache!` resizes and fills every vector in place and finally sets `valid = true`; `_query_vacuum_gram_cache!` reads them. Vectors are reused across rebuilds via `resize!`, so no allocation occurs after the first build for a fixed `n_pts`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VacuumPredictedGRAMCache | n/a | — | Return value of `VacuumPredictedGRAMCache`. Returns `VacuumPredictedGRAMCache(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.vacuum_predicted_gram__vacuum_gram_cache_for_sat_bang|_vacuum_gram_cache_for_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:53-53`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The cache is not thread-safe: a rebuild mutates `valid`, the knot vectors and `t0`/`h` non-atomically, so concurrent readers of the same satellite cache can observe a half-built state. `h = 1.0` in the empty constructor is a placeholder that would silently give wrong interpolation indices if a query bypassed the `valid` check. Nothing records which density model or planet the cache was built against, so swapping models mid-run reuses stale splines until the position-deviation check fails.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 38.
