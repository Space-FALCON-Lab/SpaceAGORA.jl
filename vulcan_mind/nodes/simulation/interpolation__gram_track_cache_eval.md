---
id: simulation.interpolation__gram_track_cache_eval
label: _gram_track_cache_eval
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _gram_track_cache_eval
  lines:
  - 116
  - 116
inputs:
- id: cache
  type: GramTrackCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
- id: x
  type: Float64
  units: n/a
  required: true
  description: Positional argument `x`.
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
  type: Any
  units: n/a
  description: Return value of `_gram_track_cache_eval`. Returns `ρ, T, wind`.
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

# _gram_track_cache_eval

## Purpose
Evaluates the cached atmosphere at a located segment, returning the interpolated density, temperature, and three-component wind vector that substitute for a full GRAM model call.

## Design & Implementation
Takes the `cache`, the segment index `idx`, and the fraction `x` produced by `_gram_track_cache_segment`. Density and temperature come from `_lerp` on `cache.rhos` and `cache.Ts`. The wind is built component-wise as `SVector{3, Float64}` using `muladd(x, w1[k], mx * w0[k])` with `mx = 1.0 - x`, which the inline comment justifies: writing it in the `(1-x)*w0 + x*w1` form avoids constructing the intermediate `SVector` difference that `w0 + x*(w1 - w0)` would need. The endpoint fetches from `cache.winds` are wrapped in `@inbounds`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | GramTrackCache | n/a | yes | Positional argument `cache`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_track_cache_eval`. Returns `ρ, T, wind`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:171-171`

**Downstream**

- `callees` → [[simulation.interpolation__lerp|_lerp]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:117-117`
<!-- vulcan:connections:end -->

## Limitations
The `@inbounds` annotations remove the safety net, so an `idx` of `length(cache.winds)` — or any index not produced by `_gram_track_cache_segment` — reads past the end of the array and returns garbage or crashes rather than throwing a `BoundsError`. Nothing revalidates that `cache.valid` still holds or that `x` is within `[0, 1]`; the contract with the segment finder is enforced only by convention. Linear interpolation of density is a poor approximation across a scale height, so wide cache spacing biases the drag computation.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 116.
