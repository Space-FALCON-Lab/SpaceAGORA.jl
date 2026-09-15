---
id: environment.simple_ephemerides__earth_gmst_iau82_rad
label: _earth_gmst_iau82_rad
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _earth_gmst_iau82_rad
  lines:
  - 99
  - 99
inputs:
- id: ut1_seconds_past_j2000
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ut1_seconds_past_j2000`.
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
  description: Return value of `_earth_gmst_iau82_rad`.
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

# _earth_gmst_iau82_rad

## Purpose
Computes Greenwich Mean Sidereal Time in radians from the IAU-82 polynomial, giving the rotation angle of Earth's prime meridian relative to the J2000 vernal equinox so the `SimpleEphemeridesModel` places geographic (lat/lon-keyed) atmosphere evaluations at the correct longitude.

## Theory & Math
With $T_u = t_{UT1} / (36525 \cdot 86400)$ (Julian centuries from J2000): $\mathrm{GMST}_s = 67310.54841 + 3.164400184812866\times10^{9}\,T_u + 0.093104\,T_u^2 - 6.2\times10^{-6}\,T_u^3$ seconds, and $\theta = \left(\mathrm{GMST}_s \bmod 86400\right)\cdot\dfrac{2\pi}{86400}$ wrapped to $[0, 2\pi)$.

## Design & Implementation
The argument `ut1_seconds_past_j2000` is converted to Julian centuries `Tu = t / (36525 * 86400)`. `@evalpoly(Tu, 67310.54841, 3.164400184812866e9, 0.093104, -6.2e-6)` evaluates the IAU-82 GMST polynomial in seconds (the linear coefficient is `876600 h * 3600 + 8640184.812866 s`). The result is reduced with `rem(gmst_s, 86400.0)`, scaled by `2π / 86400` to radians, and wrapped into `[0, 2π)` by `mod2pi`. The function is `@inline`, allocation-free, and pure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ut1_seconds_past_j2000` | Float64 | n/a | yes | Positional argument `ut1_seconds_past_j2000`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_earth_gmst_iau82_rad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`

**Downstream**

- `callees` → [[environment.simple_ephemerides__rotation_about_spin_axis|_rotation_about_spin_axis]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:119-119`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:105-105`
<!-- vulcan:connections:end -->

## Limitations
The simple model's timeline is leap-second-free UTC seconds treated directly as UT1; the comment bounds the error at |UT1 - UTC| < 0.9 s, or under 4e-3 degrees of Earth rotation. The polynomial loses precision far from J2000 because `Tu` is small while the linear coefficient is ~3e9; `rem` on a value of order 1e11 seconds after a century retains only ~1e-5 s resolution. Applies only to Earth; other planets use `ω[3] * elapsed`.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 99.
