---
id: environment.simple_ephemerides__rotation_about_spin_axis
label: _rotation_about_spin_axis
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _rotation_about_spin_axis
  lines:
  - 65
  - 65
inputs:
- id: theta
  type: Float64
  units: n/a
  required: true
  description: Positional argument `θ`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_rotation_about_spin_axis`.
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

# _rotation_about_spin_axis

## Purpose
Builds the 3x3 direction-cosine matrix that rotates J2000 inertial vectors into a planet-fixed frame whose z-axis coincides with the J2000 z-axis, given the planet's rotation angle `θ` about that axis. It is the sole frame transform used by `SimpleEphemeridesModel`.

## Theory & Math
$R(\theta) = \begin{pmatrix} \cos\theta & \sin\theta & 0 \\ -\sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix}$, so that $\mathbf{r}_{PCPF} = R(\theta)\,\mathbf{r}_{J2000}$.

## Design & Implementation
Computes `c = cos(θ)` and `s = sin(θ)` once and returns a static `@SMatrix [c s 0; -s c 0; 0 0 1]` of type `SMatrix{3,3,Float64}`. The sign convention (`+s` in row 1, `-s` in row 2) is the passive rotation that maps inertial coordinates to body-fixed coordinates for a positive spin angle, matching the direction of `pxform("J2000", body_frame, et)` in the SPICE path. Being `@inline` and static, it allocates nothing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Float64 | n/a | yes | Positional argument `θ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_rotation_about_spin_axis`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.simple_ephemerides__earth_gmst_iau82_rad|_earth_gmst_iau82_rad]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:119-119`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The transform ignores pole obliquity, precession, nutation, and polar motion: the planet spin axis is assumed to be exactly aligned with J2000 z. For Earth this misplaces the pole by roughly 23.4 degrees of ecliptic obliquity only if callers mistakenly use an ecliptic frame; within the J2000 equatorial frame the error is the neglected precession/nutation (tens of arcseconds per decade). No argument validation is done for non-finite `θ`.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 65.
