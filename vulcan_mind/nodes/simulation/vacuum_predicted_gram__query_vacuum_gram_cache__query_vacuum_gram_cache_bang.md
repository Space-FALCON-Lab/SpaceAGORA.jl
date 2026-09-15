---
id: simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang
label: _query_vacuum_gram_cache!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _query_vacuum_gram_cache!
  lines:
  - 238
  - 238
inputs:
- id: cache
  type: VacuumPredictedGRAMCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: density_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `density_model`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_ii`.
- id: alt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: n_pts
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_pts`.
- id: horizon_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `horizon_s`.
- id: deviation_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `deviation_m`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_query_vacuum_gram_cache!`; mutates `cache` in place.
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

# _query_vacuum_gram_cache!

## Purpose
Serves a density query `(ρ, T, wind)` for one satellite from the vacuum-predicted spline cache when the cache is valid, covers the query time, and the actual position is within `deviation_m` of the drag-free prediction; otherwise it rebuilds the cache from the current state. It is the entry point the density callback uses when `SPACEAGORA_VACUUM_GRAM_CACHE` is on.

## Design & Implementation
Signature `_query_vacuum_gram_cache!(cache, density_model, p, pos_ii, vel_ii::SVector{3,Float64}, alt, t::Float64, n_pts::Int, horizon_s, deviation_m::Float64)::Tuple{Float64, Float64, SVector{3,Float64}}`. The hit path requires `cache.valid && cache.t0 <= t <= cache.t1` and `norm(pos_ii - _interp_vacuum_position(cache, t)) <= deviation_m`; it then evaluates `log(ρ)` and `T` with `_eval_natural_cubic_spline` and wind with `_interp_vacuum_wind`, returning `exp(log_rho)`. On a miss it calls `_build_vacuum_gram_cache!(cache, density_model, p, pos_ii, vel_ii, t, n_pts, horizon_s)`, which propagates with `_vacuum_rk4_step`, samples `getDensity` at each knot, stores `log(max(rho, 1e-40))`, and fits both splines. Since the rebuild starts at `t`, the first knot is the current point and `(exp(log_rhos[1]), Ts[1], winds[1])` is returned directly. If the build left the cache invalid (`n_pts < 2`), it falls back to a direct `getDensity` call after converting `pos_ii` to planet-fixed coordinates with `_planet_lpi_at(p, t)` and `rtolatlong`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | VacuumPredictedGRAMCache | n/a | yes | Positional argument `cache`. |
| in | `density_model` | Any | n/a | yes | Positional argument `density_model`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Positional argument `vel_ii`. |
| in | `alt` | Float64 | n/a | yes | Positional argument `alt`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `n_pts` | Int | n/a | yes | Positional argument `n_pts`. |
| in | `horizon_s` | Float64 | n/a | yes | Positional argument `horizon_s`. |
| in | `deviation_m` | Float64 | n/a | yes | Positional argument `deviation_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_query_vacuum_gram_cache!`; mutates `cache` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:103-103`

**Downstream**

- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:272-272`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:273-273`
- `callees` → [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:270-270`
- `callees` → [[simulation.vacuum_predicted_gram__eval_natural_cubic_spline|_eval_natural_cubic_spline]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:253-253`
- `callees` → [[simulation.vacuum_predicted_gram__interp_vacuum_position|_interp_vacuum_position]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:251-251`
- `callees` → [[simulation.vacuum_predicted_gram__interp_vacuum_wind|_interp_vacuum_wind]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:255-255`
- `callees` → [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:261-261`
<!-- vulcan:connections:end -->

## Limitations
A rebuild costs `n_pts` synchronous density-model calls plus the RK4 propagation, so a trajectory that oscillates around the deviation threshold triggers repeated expensive rebuilds; there is no hysteresis or minimum time between rebuilds. The `alt` argument is accepted but never used. Backward time queries (`t < cache.t0`, as happen when an adaptive integrator rejects a step) always miss and rebuild from the earlier state. The floor `1e-40` on density keeps `log` finite but means a genuinely zero density (`NoAtmosphereModel`) is returned as `1e-40` rather than `0.0`. The mutation of the cache is not synchronised, so one cache must not be shared between threads.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 238.
