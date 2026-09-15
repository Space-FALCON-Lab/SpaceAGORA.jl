---
id: simulation.targeting__gram_periapsis_target
label: _gram_periapsis_target
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_periapsis_target
  lines:
  - 112
  - 112
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: include_j2
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `include_j2` (default `true`).
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_gram_periapsis_target`.
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

# _gram_periapsis_target

## Purpose
Predicts the time to the next periapsis passage and the geodetic position at that passage, so the GRAM track cache can be split at periapsis where density gradients are steepest.

## Theory & Math
$M = (E - e\sin E) \bmod 2\pi$, $\Delta t_{peri} = \dfrac{2\pi - M}{n}$ with $n = \sqrt{\mu/a^3}$.

## Design & Implementation
Inside a `try`, extracts `a, e, ν` from `rvtoorbitalelement`, rejecting non-finite values, `a <= 0`, or `e` outside `[0, 1)`. Computes `n = sqrt(planet.μ / a^3)`, eccentric anomaly `E` via `_true_to_eccentric_anomaly`, mean anomaly `M = mod(E - e sin E, 2π)`, and `dt_peri = (2π - M) / n`, which is the time remaining to `M = 2π`. Rejects `dt_peri <= 1e-6`. Delegates to `_gram_kepler_target(pos, vel, planet, dt_peri; include_j2)` and returns `(dt_peri, alt_peri, lat_peri, lon_peri)` or `nothing`. Exceptions return `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `include_j2` | Bool | n/a | no | Keyword argument `include_j2` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_gram_periapsis_target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:93-93`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:120-120`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:119-119`
- `callees` → [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:135-135`
- `callees` → [[simulation.targeting__true_to_eccentric_anomaly|_true_to_eccentric_anomaly]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:129-129`
<!-- vulcan:connections:end -->

## Limitations
Exactly at periapsis (`M = 0`) the function returns a full period rather than zero, because `2π - 0 = 2π`; just past periapsis it correctly returns nearly a period. Drag shortens the true time to periapsis during aerobraking, so the target lags reality. Two-body assumptions apply as in `_gram_kepler_target`. Elements are converted twice per call.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 112.
