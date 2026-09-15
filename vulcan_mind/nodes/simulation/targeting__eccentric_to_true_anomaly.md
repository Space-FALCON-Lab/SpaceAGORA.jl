---
id: simulation.targeting__eccentric_to_true_anomaly
label: _eccentric_to_true_anomaly
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _eccentric_to_true_anomaly
  lines:
  - 67
  - 67
inputs:
- id: E
  type: Float64
  units: n/a
  required: true
  description: Positional argument `E`.
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
  description: Return value of `_eccentric_to_true_anomaly`.
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

# _eccentric_to_true_anomaly

## Purpose
Converts an eccentric anomaly `E` to the true anomaly `ν` for an elliptic orbit of eccentricity `e`, the final step of Kepler propagation inside `_gram_kepler_target`.

## Theory & Math
$\nu = \operatorname{atan2}\!\left(\sqrt{1-e^{2}}\,\sin E,\; \cos E - e\right) \bmod 2\pi$.

## Design & Implementation
Uses the quadrant-safe form `ν = atan(sqrt(1 - e^2) * sin(E), cos(E) - e)`, with `max(0.0, 1 - e^2)` guarding against negative radicands from rounding when `e` approaches 1. The result is wrapped to `[0, 2π)` with `mod(ν, 2pi)`. The function is `@inline`, takes and returns `Float64`, and allocates nothing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `E` | Float64 | n/a | yes | Positional argument `E`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_eccentric_to_true_anomaly`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:92-92`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Valid only for `0 <= e < 1`; callers gate on that but the function itself accepts any `e` and silently returns a value for hyperbolic eccentricities. For `e` extremely close to 1 the `sqrt` term collapses to zero and `ν` snaps toward 0 or π regardless of `E`. NaN inputs propagate. `mod(ν, 2pi)` can return exactly `2pi` for tiny negative inputs due to floating rounding.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 67.
