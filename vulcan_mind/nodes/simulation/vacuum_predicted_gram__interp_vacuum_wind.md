---
id: simulation.vacuum_predicted_gram__interp_vacuum_wind
label: _interp_vacuum_wind
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _interp_vacuum_wind
  lines:
  - 172
  - 172
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
  description: Return value of `_interp_vacuum_wind`.
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

# _interp_vacuum_wind

## Purpose
Linearly interpolates the wind vector sampled from the density model at cache knots to the query time, completing the `(ρ, T, wind)` triple that the cached density lookup returns. Wind is interpolated linearly rather than with a spline because it is a three-component vector and its accuracy matters less than density.

## Design & Implementation
Signature `_interp_vacuum_wind(cache::VacuumPredictedGRAMCache, t::Float64)::SVector{3,Float64}`, `@inline`. It computes `idx` and `x` exactly as the altitude and position interpolators do, from `cache.t0`, `cache.h` and `n = length(cache.winds)`, loads `w0 = winds[idx]` and `w1 = winds[idx+1]` with `@inbounds`, and returns `(1 - x) w0 + x w1`. The winds were stored by `_build_vacuum_gram_cache!` from the third return value of `getDensity(density_model, alt, lat, lon, t, true, p)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | VacuumPredictedGRAMCache | n/a | yes | Positional argument `cache`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_interp_vacuum_wind`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:255-255`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The frame of the wind vector is whatever `getDensity` returns and is not transformed here; interpolating vectors expressed in a rotating planet-fixed frame across a time gap ignores frame rotation over that gap. Query times beyond the horizon are extrapolated. The function shares the unchecked `cache.valid` and empty-vector hazards of its siblings, relying entirely on `_query_vacuum_gram_cache!` to have validated the cache first.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 172.
