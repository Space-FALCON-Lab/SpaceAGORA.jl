---
id: environment.gravity_models_j2_secular_rates
label: j2_secular_rates
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: j2_secular_rates
  lines:
  - 120
  - 120
inputs:
- id: a
  type: Float64
  units: n/a
  required: true
  description: Positional argument `a`.
- id: e
  type: Float64
  units: n/a
  required: true
  description: Positional argument `e`.
- id: i
  type: Float64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  description: Return value of `j2_secular_rates`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# j2_secular_rates

## Purpose
Returns the first-order secular drift rates of the right ascension of the ascending node and the argument of periapsis due to J2, used by analytic checks and verification tooling to compare against integrated J2 propagation.

## Theory & Math
$$\dot\Omega = -\frac{3}{2} n J_2 \left(\frac{R_e}{p}\right)^2 \cos i, \qquad \dot\omega = \frac{3}{4} n J_2 \left(\frac{R_e}{p}\right)^2 \left(5\cos^2 i - 1\right)$$ with $n = \sqrt{\mu/a^3}$ and $p = a(1 - e^2)$.

## Design & Implementation
Marked `@inline`; arguments are semi-major axis `a` (m), eccentricity `e`, inclination `i` (rad) and a `planet` with `μ`, `J2` and `Rp_e`. Returns `(0.0, 0.0)` if `J2` is zero or non-finite, if `a <= 0`, `e < 0` or `e >= 1`, if the semi-latus rectum `p = a(1 - e^2)` is not positive, or if the mean motion `n = sqrt(μ/a^3)` is not positive. Otherwise `scale = J2 (Rp_e/p)^2`, `Ωdot = -1.5 n scale cos(i)` and `ωdot = 0.75 n scale (5 cos^2(i) - 1)`, both rad/s, following Vallado and Montenbruck & Gill.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | yes | Positional argument `a`. |
| in | `e` | Float64 | n/a | yes | Positional argument `e`. |
| in | `i` | Float64 | n/a | yes | Positional argument `i`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `j2_secular_rates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:112-112`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:97-97`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only first-order J2 secular terms are included; short-period, long-period and J2^2 terms are ignored, so the rates are inaccurate for very low or highly eccentric orbits. Hyperbolic and parabolic orbits return zero rather than throwing, which can mask misuse. Inclination is not range-checked, so values outside `[0, π]` are accepted.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 120.
