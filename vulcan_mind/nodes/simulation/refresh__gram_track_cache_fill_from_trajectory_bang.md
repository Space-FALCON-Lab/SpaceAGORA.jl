---
id: simulation.refresh__gram_track_cache_fill_from_trajectory_bang
label: _gram_track_cache_fill_from_trajectory!
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/refresh.jl
  symbol: _gram_track_cache_fill_from_trajectory!
  lines:
  - 1
  - 1
inputs:
- id: cache
  type: GramTrackCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: trajectory
  type: Any
  units: n/a
  required: true
  description: Positional argument `trajectory`.
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
  description: Return value of `_gram_track_cache_fill_from_trajectory!`; mutates
    `cache` in place. Returns `nothing`.
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

# _gram_track_cache_fill_from_trajectory!

## Purpose
Copies a GRAM-generated trajectory into a `GramTrackCache`, converting each sample into the internal SI and radian conventions used by the density callbacks. Mutates `cache` in place and returns `nothing`.

## Design & Implementation
It requires at least two points, throwing `ArgumentError` when `length(trajectory) < 2`. When the cached arrays are a different length it `resize!`s `times`, `alts`, `lats`, `lons`, `rhos`, `Ts` and `winds` together to `n`. The `@inbounds` loop then writes `elapsedTime` as seconds, `height` scaled by `1e3` from kilometres to metres, `latitude` and `longitude` through `deg2rad`, `dynamics.density` and `dynamics.temperature` verbatim, and packs `perturbedEWWind`, `perturbedNSWind` and `perturbedVerticalWind` into an `SVector{3, Float64}`. Longitude is wrapped into the half-open range by a single conditional add or subtract of `2π` instead of a trigonometric reduction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | GramTrackCache | n/a | yes | Positional argument `cache`. |
| in | `trajectory` | Any | n/a | yes | Positional argument `trajectory`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_gram_track_cache_fill_from_trajectory!`; mutates `cache` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:226-226`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
The longitude wrap is one-shot: an input outside roughly $(-3\pi, 3\pi]$ stays outside the intended range because only one correction is applied. The resize path assumes the seven cache vectors are always the same length, so an externally desynchronised cache is silently left inconsistent for the arrays that already matched. Field access is duck-typed on `pt.position`, `pt.dynamics` and `pt.winds`, so a GRAM driver with renamed fields fails only at runtime. The routine does not mark the cache valid or set its time bounds; the caller must do that.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/refresh.jl` line 1.
