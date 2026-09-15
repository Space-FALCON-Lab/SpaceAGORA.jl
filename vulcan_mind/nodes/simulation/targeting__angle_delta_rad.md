---
id: simulation.targeting__angle_delta_rad
label: _angle_delta_rad
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _angle_delta_rad
  lines:
  - 13
  - 13
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
  type: Any
  units: n/a
  description: Return value of `_angle_delta_rad`. Returns `atan(sin(b - a), cos(b
    - a))`.
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

# _angle_delta_rad

## Purpose
Returns the signed shortest angular difference `b - a` wrapped into `(-π, π]`, used when comparing longitudes so that a track crossing the antimeridian is not measured as nearly a full revolution.

## Theory & Math
$\Delta = \operatorname{atan2}\big(\sin(b-a),\,\cos(b-a)\big) \in (-\pi, \pi]$, the principal value of $b - a$ modulo $2\pi$.

## Design & Implementation
An `@inline` one-liner `atan(sin(b - a), cos(b - a))`, which uses the two-argument arctangent to wrap the difference without explicit modular arithmetic or branching. Both arguments are `Float64` radians. Evaluating `sin` and `cos` of the raw difference is robust for any magnitude of input, including differences exceeding many revolutions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Positional argument `a`. |
| in | `b` | Float64 | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_angle_delta_rad`. Returns `atan(sin(b - a), cos(b - a))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.targeting__gram_expected_track_length_m|_gram_expected_track_length_m]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:26-26`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:200-200`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Costs three transcendental evaluations per call, which is heavier than a `rem2pi`-based wrap when the inputs are known to be already in range. Precision degrades for very large `|b - a|` (above ~1e8 rad) because argument reduction in `sin`/`cos` loses digits. The result for a difference of exactly π is `π` (not `-π`), and NaN inputs propagate to a NaN result rather than throwing.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 13.
