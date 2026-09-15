---
id: gnc.thruster_guidance_functions__oblate_surface_radius
label: _oblate_surface_radius
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _oblate_surface_radius
  lines:
  - 62
  - 62
inputs:
- id: u_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `u_pp`.
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
  type: Float64
  units: n/a
  description: Return value of `_oblate_surface_radius`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _oblate_surface_radius

## Purpose
Returns the geocentric distance from the planet centre to the reference ellipsoid surface along a given planet-fixed unit direction, serving as the lower bracket for the altitude-to-radius bisection in `_radius_for_oblate_altitude`.

## Theory & Math
$r_s(\hat{u}) = \left(\frac{u_x^2+u_y^2}{R_e^2}+\frac{u_z^2}{R_p^2}\right)^{-1/2}$ where $\hat{u}$ is the planet-fixed unit direction, $R_e$ the equatorial radius (m) and $R_p$ the polar radius (m).

## Design & Implementation
An `@inline` pure function of `u_pp::SVector{3,Float64}` and `planet` (fields `Rp_e`, `Rp_p` in m). It evaluates the ellipsoid equation for the ray `r*u_pp`: `inv(sqrt((u1^2+u2^2)/Rp_e^2 + u3^2/Rp_p^2))`, which is the radius `r` at which `x^2/Rp_e^2 + y^2/Rp_e^2 + z^2/Rp_p^2 = 1`. Returns a `Float64` in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_oblate_surface_radius`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions__radius_for_oblate_altitude|_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:68-68`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Assumes `u_pp` is unit length; a non-unit vector scales the result by `1/|u_pp|`. A zero vector divides by zero and returns `Inf`. No check is made that `Rp_e` and `Rp_p` are positive.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 62.
