---
id: simulation.targeting__gram_linear_target
label: _gram_linear_target
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_linear_target
  lines:
  - 173
  - 173
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
- id: dt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt`.
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
  description: Return value of `_gram_linear_target`.
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

# _gram_linear_target

## Purpose
Simplest track-target predictor: extrapolates the inertial position linearly by `dt` along the current velocity and converts the result to geodetic coordinates, serving as the fallback when Keplerian prediction is refused.

## Design & Implementation
Computes `pos_target = pos + vel * dt` using `SVector{3,Float64}` arithmetic, rotates it into the planet-fixed frame via `r_intor_p!(pos_target, vel, planet)` (velocity is passed through but its rotated value is discarded), and calls `rtolatlong(rp_target, planet)` to obtain `(alt, lat, lon)`. The function is `@inline` and returns a `Tuple{Float64, Float64, Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `dt` | Float64 | n/a | yes | Positional argument `dt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_linear_target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- [[simulation.refresh__gram_kepler_or_linear_target|_gram_kepler_or_linear_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:54-54`

**Downstream**

- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:180-180`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:181-181`
<!-- vulcan:connections:end -->

## Limitations
Linear extrapolation ignores gravity, so the altitude error grows as `0.5 g dt^2` (about 4.9 km after 100 s at Earth surface gravity) and the prediction overshoots radially outward; it is only sensible for `dt` of a few seconds. `r_intor_p!` is evaluated at the current time, not `t + dt`, so planet rotation over `dt` is not accounted for in longitude. No validation of `dt` is performed.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 173.
