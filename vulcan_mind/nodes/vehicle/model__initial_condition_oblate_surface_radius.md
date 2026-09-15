---
id: vehicle.model__initial_condition_oblate_surface_radius
label: _initial_condition_oblate_surface_radius
kind: function
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: _initial_condition_oblate_surface_radius
  lines:
  - 142
  - 142
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
  description: Return value of `_initial_condition_oblate_surface_radius`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# _initial_condition_oblate_surface_radius

## Purpose
Returns the distance from the planet centre to the reference ellipsoid along a planet-fixed unit direction, providing the lower bracket for the altitude-to-radius bisection in the oblate `InitialCondition` constructor.

## Theory & Math
$r_s(\hat u)=\left(\frac{u_x^2+u_y^2}{R_e^2}+\frac{u_z^2}{R_p^2}\right)^{-1/2}$, $R_e$ equatorial and $R_p$ polar radius (m).

## Design & Implementation
`@inline` pure function of `u_pp::SVector{3,Float64}` and `planet` returning `inv(sqrt((u1^2+u2^2)/Rp_e^2 + u3^2/Rp_p^2))`, the scalar `r` solving `(r u)` on the ellipsoid `x^2/Rp_e^2 + y^2/Rp_e^2 + z^2/Rp_p^2 = 1`. Units are metres. It is identical in form to `_oblate_surface_radius` in the thruster guidance module.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_initial_condition_oblate_surface_radius`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:154-154`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Assumes `u_pp` is unit length; a non-unit input scales the result. A zero vector produces `Inf`. Negative or zero planet radii are not checked. The duplication with the guidance module means the two must be kept consistent by hand.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 142.
