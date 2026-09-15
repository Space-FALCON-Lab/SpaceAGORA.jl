---
id: simulation.targeting__gram_kepler_target
label: _gram_kepler_target
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_kepler_target
  lines:
  - 72
  - 72
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
- id: dt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt`.
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
  description: Return value of `_gram_kepler_target`.
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

# _gram_kepler_target

## Purpose
Predicts the geodetic altitude, latitude, and longitude the spacecraft will occupy `dt` seconds ahead under two-body Keplerian motion, optionally with J2 secular drift of RAAN and argument of periapsis, to aim the GRAM track cache prefetch along the future ground track.

## Design & Implementation
Inside a `try`, converts `(pos, vel)` to classical elements with `rvtoorbitalelement`, extracting `a, e, i, Ω, ω, ν`. Returns `nothing` when `dt` is non-finite or `<= 1e-6`, or when `a <= 0`, `e < 0`, `e >= 1`, or any element is non-finite. Mean motion `n = sqrt(planet.μ / a^3)`; `E0` from `_true_to_eccentric_anomaly`, `M0 = E0 - e sin E0`, `E1 = _solve_kepler_elliptic(M0 + n dt, e)`, `ν1` from `_eccentric_to_true_anomaly`. With `include_j2=true`, `j2_secular_rates(a, e, i, planet)` supplies `Ωdot, ωdot` applied linearly over `dt`. The target elements `SVector{7}(a, e, i, Ω1, ω1, ν1, 0.0)` are converted back with `orbitalelemtorv`, rotated to planet-fixed by `r_intor_p!`, and reduced to `(alt, lat, lon)` by `rtolatlong`. Any exception returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `dt` | Float64 | n/a | yes | Positional argument `dt`. |
| in | `include_j2` | Bool | n/a | no | Keyword argument `include_j2` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_gram_kepler_target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- [[simulation.refresh__gram_kepler_or_linear_target|_gram_kepler_or_linear_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:53-53`
- [[simulation.targeting__gram_orbit_period_target|_gram_orbit_period_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:164-164`
- [[simulation.targeting__gram_periapsis_target|_gram_periapsis_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:135-135`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:81-81`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:103-103`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:104-104`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:105-105`
- `callees` → [[environment.gravity_models_j2_secular_rates|j2_secular_rates]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:97-97`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:80-80`
- `callees` → [[simulation.targeting__eccentric_to_true_anomaly|_eccentric_to_true_anomaly]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:92-92`
- `callees` → [[simulation.targeting__solve_kepler_elliptic|_solve_kepler_elliptic]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:91-91`
- `callees` → [[simulation.targeting__true_to_eccentric_anomaly|_true_to_eccentric_anomaly]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:89-89`
<!-- vulcan:connections:end -->

## Limitations
Elliptic orbits only; parabolic and hyperbolic cases are rejected. Drag, higher-order harmonics, and third-body effects are ignored, so during aerobraking passes the prediction diverges from truth at low altitude (the Allen-Eggers predictor is used there instead). The bare `catch` hides errors in element conversion. J2 rates are applied as constant secular drift, valid only for short `dt` relative to the nodal period.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 72.
