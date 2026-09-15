---
id: simulation.targeting__gram_expected_track_length_m
label: _gram_expected_track_length_m
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_expected_track_length_m
  lines:
  - 15
  - 15
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: alt0
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt0`.
- id: lat0
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat0`.
- id: lon0
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon0`.
- id: alt1
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt1`.
- id: lat1
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat1`.
- id: lon1
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon1`.
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
  description: Return value of `_gram_expected_track_length_m`.
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

# _gram_expected_track_length_m

## Purpose
Estimates the ground-track length (m) between two geodetic points at possibly different altitudes and the mean radius at which that track lies, so the GRAM track cache can decide how many interpolation samples to request.

## Theory & Math
$\cos\sigma = \sin\phi_0\sin\phi_1 + \cos\phi_0\cos\phi_1\cos\Delta\lambda$, $L = \sqrt{(R\,\sigma)^2 + (h_1-h_0)^2}$ with $R = R_p + \tfrac{1}{2}(h_0+h_1)$.

## Design & Implementation
Takes `planet`, start `(alt0, lat0, lon0)`, and end `(alt1, lat1, lon1)` in metres and radians. `radius_m = max(1.0, planet.Rp_m + 0.5 * (alt0 + alt1))` is the mean geocentric radius. The great-circle central angle is `acos` of the clamped spherical-law-of-cosines expression `sin(lat0)sin(lat1) + cos(lat0)cos(lat1)cos(Δlon)` with `Δlon` from `_angle_delta_rad`. Horizontal length is `radius_m * central_angle`, vertical is `alt1 - alt0`, and the total is `hypot(horizontal, vertical)`. Returns `(max(length_m, 1.0), radius_m)` as a `Tuple{Float64, Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `alt0` | Float64 | n/a | yes | Positional argument `alt0`. |
| in | `lat0` | Float64 | n/a | yes | Positional argument `lat0`. |
| in | `lon0` | Float64 | n/a | yes | Positional argument `lon0`. |
| in | `alt1` | Float64 | n/a | yes | Positional argument `alt1`. |
| in | `lat1` | Float64 | n/a | yes | Positional argument `lat1`. |
| in | `lon1` | Float64 | n/a | yes | Positional argument `lon1`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_expected_track_length_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:162-162`

**Downstream**

- `callees` → [[simulation.targeting__angle_delta_rad|_angle_delta_rad]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:26-26`
<!-- vulcan:connections:end -->

## Limitations
The spherical law of cosines loses precision for very short tracks (central angle below ~1e-4 rad) where the haversine form would be preferred; the 1 m floor masks the worst of this. A spherical planet of radius `Rp_m` is assumed, ignoring flattening. The vertical and horizontal components are combined as if the track were a straight line in a flat projection, which is inaccurate for long arcs with large altitude change.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 15.
