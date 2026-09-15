---
id: simulation.targeting__true_to_eccentric_anomaly
label: _true_to_eccentric_anomaly
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _true_to_eccentric_anomaly
  lines:
  - 62
  - 62
inputs:
- id: nu
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ν`.
- id: e
  type: Float64
  units: n/a
  required: true
  description: Positional argument `e`.
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
  description: Return value of `_true_to_eccentric_anomaly`.
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

# _true_to_eccentric_anomaly

## Purpose
Converts a true anomaly `ν` to the eccentric anomaly `E` for an elliptic orbit of eccentricity `e`, the first step in converting the current state to mean anomaly for Kepler propagation.

## Theory & Math
$E = \operatorname{atan2}\!\left(\sqrt{1-e^{2}}\,\sin\nu,\; e + \cos\nu\right) \bmod 2\pi$.

## Design & Implementation
Uses the quadrant-safe form `E = atan(sqrt(1 - e^2) * sin(ν), e + cos(ν))`, with `max(0.0, 1 - e^2)` protecting the square root from tiny negative rounding errors near `e = 1`. The result is wrapped into `[0, 2π)` via `mod(E, 2pi)`. The function is `@inline`, takes and returns `Float64`, and allocates nothing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nu` | Float64 | n/a | yes | Positional argument `ν`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_true_to_eccentric_anomaly`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:89-89`
- [[simulation.targeting__gram_periapsis_target|_gram_periapsis_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:129-129`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Valid only for `0 <= e < 1`; no check is made and hyperbolic eccentricities silently yield meaningless values. As `e` approaches 1 the numerator collapses and `E` becomes insensitive to `ν`. NaN propagates. The `mod` wrap can return `2π` exactly for tiny negative arguments because of floating-point rounding.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 62.
