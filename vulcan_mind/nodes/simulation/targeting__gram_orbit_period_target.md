---
id: simulation.targeting__gram_orbit_period_target
label: _gram_orbit_period_target
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_orbit_period_target
  lines:
  - 144
  - 144
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: include_j2
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `include_j2` (default `true`).
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_gram_orbit_period_target`.
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

# _gram_orbit_period_target

## Purpose
Computes one full Keplerian orbital period from the current state and predicts the geodetic position at that time, giving the GRAM track cache a target for orbit-repeat prefetching.

## Theory & Math
$T = 2\pi\sqrt{a^{3}/\mu}$.

## Design & Implementation
Inside a `try`, extracts `a` and `e` from `rvtoorbitalelement(pos, vel, planet)` and returns `nothing` unless both are finite with `a > 0` and `0 <= e < 1`. Mean motion `n = sqrt(planet.μ / a^3)` yields `dt_orbit = 2π / n`, rejected if non-finite or `<= 1e-6`. It then delegates to `_gram_kepler_target(pos, vel, planet, dt_orbit; include_j2)` and returns `(dt_orbit, alt_end, lat_end, lon_end)` or `nothing` if that fails. Any exception returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `include_j2` | Bool | n/a | no | Keyword argument `include_j2` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_gram_orbit_period_target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:136-136`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:152-152`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:151-151`
- `callees` → [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:164-164`
<!-- vulcan:connections:end -->

## Limitations
The period is the osculating two-body value, so with J2 the true nodal or anomalistic period differs slightly, and with drag the next periapsis arrives earlier; the end point therefore drifts from the actual state after one revolution. The orbital elements are computed twice (here and inside `_gram_kepler_target`), doubling conversion cost. The bare `catch` suppresses diagnostic information.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 144.
