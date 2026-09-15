---
id: simulation.interpolation__lerp
label: _lerp
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _lerp
  lines:
  - 14
  - 14
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
  description: Return value of `_lerp`. Returns `a + x * (b - a)`.
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

# _lerp

## Purpose
Scalar linear interpolation used throughout the GRAM track cache to blend cached density, temperature, and altitude samples between two bracketing time points.

## Design & Implementation
Defined `@inline _lerp(a::Float64, b::Float64, x::Float64) = a + x * (b - a)`. The form `a + x*(b - a)` is chosen over `(1-x)*a + x*b` because it is a single fused multiply-add on the hot path and returns exactly `a` when `x == 0.0`. All three arguments are concretely typed `Float64`, so the call site never triggers dynamic dispatch inside the interpolation loop of `_gram_track_cache_eval` and `_gram_track_cache_ready`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Positional argument `a`. |
| in | `b` | Float64 | n/a | yes | Positional argument `b`. |
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_lerp`. Returns `a + x * (b - a)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.interpolation__gram_track_cache_eval|_gram_track_cache_eval]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:117-117`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:136-136`
- [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:105-105`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no clamping of `x` to `[0, 1]`, so out-of-range fractions extrapolate silently; callers are responsible for computing `x` from a bracketing segment. The `a + x*(b - a)` form is not exact at `x == 1.0` in floating point, where it can return a value one ulp from `b`, and it loses precision when `a` and `b` are large and nearly equal. Passing an `Int` or a `Float32` fails to match the signature rather than promoting.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 14.
