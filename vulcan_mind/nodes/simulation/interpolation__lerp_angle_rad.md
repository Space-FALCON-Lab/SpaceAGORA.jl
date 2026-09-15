---
id: simulation.interpolation__lerp_angle_rad
label: _lerp_angle_rad
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _lerp_angle_rad
  lines:
  - 31
  - 31
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
  type: Float64
  units: n/a
  description: Return value of `_lerp_angle_rad`.
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

# _lerp_angle_rad

## Purpose
Interpolates between two angles along the shorter arc, so that cached latitude and longitude tracks blend correctly across the $\pm\pi$ branch cut instead of sweeping the long way round.

## Theory & Math
With $\delta = \mathrm{wrap}(b - a) \in (-\pi, \pi]$, the result is $\mathrm{wrap}(a + x\,\delta)$, where $\mathrm{wrap}$ maps into $(-\pi, \pi]$ by adding or subtracting $2\pi$ once.

## Design & Implementation
Computes `d = b - a`, then wraps it into the short arc with two conditional adjustments: subtract `2π` when `d > π`, add `2π` when `d < -π`. The interpolated value is `r = a + x * d`, and `r` is renormalised with the same conditional subtract-or-add so the result lands in `(-π, π]`. The comment explains the design intent explicitly: this avoids the second `atan2` that a naive wrap-then-atan2 implementation would need, keeping the routine trig-free. Called by `_gram_track_cache_ready` on the cached `lats` and `lons` arrays.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Positional argument `a`. |
| in | `b` | Float64 | n/a | yes | Positional argument `b`. |
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_lerp_angle_rad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:137-137`
- [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:106-106`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Both the input wrap and the output renormalisation are single-shot conditionals, so inputs more than one full turn apart, or an `x` far outside `[0, 1]`, leave the result outside `(-π, π]`; the function silently assumes its inputs are already normalised angles and that `x` is a proper segment fraction. The asymmetric renormalisation boundaries (`r > π` versus `r ≤ -π`) place exactly `-π` into the wrapped branch while `π` is kept. Near the poles, latitude interpolated this way is not the great-circle path.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 31.
