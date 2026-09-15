---
id: simulation.targeting__gram_track_cache_target_spacing_m
label: _gram_track_cache_target_spacing_m
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_track_cache_target_spacing_m
  lines:
  - 37
  - 37
inputs:
- id: alt_tol_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt_tol_m`.
- id: ang_tol_rad
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ang_tol_rad`.
- id: radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `radius_m`.
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
  type: Float64
  units: n/a
  description: Return value of `_gram_track_cache_target_spacing_m`.
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

# _gram_track_cache_target_spacing_m

## Purpose
Chooses the along-track sample spacing (m) for a GRAM track cache so that interpolation between samples stays tighter than the altitude and angular tolerances at which cached values are accepted.

## Design & Implementation
Converts the angular tolerance to a linear distance `ang_tol_m = max(1.0, radius_m * max(ang_tol_rad, 1e-9))`, takes `tol_scale_m = min(max(1.0, alt_tol_m), ang_tol_m)` as the tighter of the two tolerances, and returns `max(1.0, 0.5 * tol_scale_m)`, i.e. half the tighter tolerance with a 1 m floor. All arguments are `Float64` and the function is `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alt_tol_m` | Float64 | n/a | yes | Positional argument `alt_tol_m`. |
| in | `ang_tol_rad` | Float64 | n/a | yes | Positional argument `ang_tol_rad`. |
| in | `radius_m` | Float64 | n/a | yes | Positional argument `radius_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gram_track_cache_target_spacing_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:171-171`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Half-tolerance spacing is a heuristic with no error analysis behind it; steep density gradients near periapsis may still be under-resolved when tolerances are loose. The 1 m floor and the `1e-9` rad floor are hard-coded. Because the spacing is derived from acceptance tolerances rather than density scale height, a large `alt_tol_m` yields coarse sampling irrespective of atmospheric structure.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 37.
