---
id: simulation.refresh__gram_kepler_or_linear_target
label: _gram_kepler_or_linear_target
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/refresh.jl
  symbol: _gram_kepler_or_linear_target
  lines:
  - 46
  - 46
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
- id: include_j2
  type: Bool
  units: n/a
  required: true
  description: Positional argument `include_j2`.
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
  description: Return value of `_gram_kepler_or_linear_target`.
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

# _gram_kepler_or_linear_target

## Purpose
Chooses the endpoint of a GRAM track-cache segment: it prefers a Keplerian propagation of the state `dt` seconds ahead and degrades to straight-line propagation when the Keplerian solve is not applicable. Returns `(altitude_m, latitude_rad, longitude_rad)` for the segment's far end.

## Design & Implementation
Marked `@inline` and fully type-annotated, taking `pos` and `vel` as `SVector{3, Float64}` in metres and metres per second, a `planet` model, `dt` in seconds and an `include_j2` flag. It calls `_gram_kepler_target(pos, vel, planet, dt; include_j2=include_j2)`; that routine signals inapplicability by returning `nothing`, in which case this routine falls through to `_gram_linear_target(pos, vel, planet, dt)`. The declared return type `Tuple{Float64, Float64, Float64}` forces both branches into the same concrete shape so callers inside the refresh loop stay allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `dt` | Float64 | n/a | yes | Positional argument `dt`. |
| in | `include_j2` | Bool | n/a | yes | Positional argument `include_j2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_kepler_or_linear_target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:89-89`

**Downstream**

- `callees` → [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:53-53`
- `callees` → [[simulation.targeting__gram_linear_target|_gram_linear_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:54-54`
<!-- vulcan:connections:end -->

## Limitations
The fallback is silent: a hyperbolic or otherwise degenerate state yields a linear extrapolation with no warning, which for long `dt` can place the target far from the true trajectory and mis-size the cached track. There is no validation that `dt` is finite or positive, and non-finite inputs propagate straight through into the returned tuple.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/refresh.jl` line 46.
