---
id: simulation.interpolation__angdiff_rad
label: _angdiff_rad
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _angdiff_rad
  lines:
  - 19
  - 19
inputs:
- id: a
  type: Float64
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: Float64
  units: n/a
  required: true
  description: Positional argument `b`.
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
  description: Return value of `_angdiff_rad`.
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

# _angdiff_rad

## Purpose
Computes the absolute short-arc angular separation in radians between two angles, used by the GRAM track cache to decide whether the current latitude and longitude are close enough to the interpolated cached ground track to reuse a cached atmosphere sample.

## Theory & Math
For angles $a, b$ the returned value is $|\,\mathrm{atan2}(\sin(b-a), \cos(b-a))\,| \in [0, \pi]$, the length of the shorter arc on the unit circle between the two directions.

## Design & Implementation
Takes `a` and `b` as `Float64` and returns `Float64`. It first forms `d = b - a` and takes the cheap path `abs(d)` when `d` already lies in the half-open interval `(-π, π]`, which the comment notes is the overwhelmingly common case for consecutive orbital cache points where `|b - a| << π`. Only when that test fails does it fall back to `abs(atan(sin(d), cos(d)))`, the full two-argument arctangent wrap, avoiding two transcendental calls per comparison on the hot path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Positional argument `a`. |
| in | `b` | Float64 | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_angdiff_rad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:140-140`
- [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:109-109`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The fallback path costs a `sin`, a `cos`, and an `atan` and is reached whenever the two angles differ by more than half a turn, so a cache whose stored longitudes are not pre-normalised pays that cost on every sample. The result is always non-negative, discarding the sign of the rotation, so it cannot be used where direction matters. NaN inputs fail the interval test and flow into the `atan` path, returning NaN, which then makes every tolerance comparison in `_gram_track_cache_ready` false and silently disables the cache.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 19.
