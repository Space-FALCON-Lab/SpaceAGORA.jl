---
id: vehicle.model__initial_condition_oblate_altitude
label: _initial_condition_oblate_altitude
kind: function
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: _initial_condition_oblate_altitude
  lines:
  - 126
  - 126
inputs:
- id: radius
  type: Float64
  units: n/a
  required: true
  description: Positional argument `radius`.
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
  description: Return value of `_initial_condition_oblate_altitude`.
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

# _initial_condition_oblate_altitude

## Purpose
Converts a geocentric radius along a planet-fixed unit direction into altitude above the oblate reference ellipsoid, used to bracket and bisect for the radius that yields a requested apsis altitude in the oblate `InitialCondition` constructor.

## Theory & Math
Bowring single-pass geodetic conversion with $p=\sqrt{x^2+y^2}$, $e^2=1-(1-f)^2$, $e'^2=e^2/(1-e^2)$, $\theta=\operatorname{atan2}(zR_e, pR_p)$, $\phi=\operatorname{atan2}(z+e'^2R_p\sin^3\theta,\ p-e^2R_e\cos^3\theta)$, $N=R_e/\sqrt{1-e^2\sin^2\phi}$. Standard altitude is $h=p\cos\phi+(z+e^2N\sin\phi)\sin\phi-N$; the code evaluates $h_{code}=p\cos\phi+(z+e^2N\sin^2\phi)\sin\phi-N$.

## Design & Implementation
`@inline` function of `radius::Float64` (m), `u_pp::SVector{3,Float64}` (unit vector, planet-fixed) and `planet` (fields `Rp_e`, `Rp_p` in m). It scales the direction to a point `(x, y, z)`, derives flattening `f = (Rp_e-Rp_p)/Rp_e`, `e2 = 1-(1-f)^2`, `ep2 = e2/(1-e2)`, `p_xy = sqrt(x^2+y^2)`, parametric latitude `θ = atan(z*Rp_e, p_xy*Rp_p)`, geodetic latitude `lat = atan(z + ep2*Rp_p*sin(θ)^3, p_xy - e2*Rp_e*cos(θ)^3)` and prime-vertical radius `N`. It then returns `p_xy*cos(lat) + (z + e2*N*sin(lat)^2)*sin(lat) - N`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `radius` | Float64 | n/a | yes | Positional argument `radius`. |
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_initial_condition_oblate_altitude`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:156-156`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The final expression uses `sin(lat)^2` inside the parenthesis whereas the standard Bowring altitude formula (and the sibling `_oblate_altitude_from_radius` in thruster guidance) uses `sin(lat)`; this introduces an altitude error that grows with latitude and eccentricity, so results differ from the guidance-side conversion for inclined orbits. A prolate planet (`Rp_e < Rp_p`) yields negative `e2` and NaN. `u_pp` is assumed unit length.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 126.
