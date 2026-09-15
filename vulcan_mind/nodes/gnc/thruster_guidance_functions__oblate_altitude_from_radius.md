---
id: gnc.thruster_guidance_functions__oblate_altitude_from_radius
label: _oblate_altitude_from_radius
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _oblate_altitude_from_radius
  lines:
  - 47
  - 47
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
  description: Return value of `_oblate_altitude_from_radius`.
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

# _oblate_altitude_from_radius

## Purpose
Converts a geocentric radius along a fixed planet-fixed unit direction into geodetic altitude above the oblate reference ellipsoid, used by periapsis-raise guidance to express the target periapsis in altitude rather than radius.

## Theory & Math
With $p=\sqrt{x^2+y^2}$, $f=(R_e-R_p)/R_e$, $e^2=1-(1-f)^2$, $e'^2=e^2/(1-e^2)$ and $\theta=\operatorname{atan2}(zR_e,\,pR_p)$: $\phi=\operatorname{atan2}\big(z+e'^2R_p\sin^3\theta,\; p-e^2R_e\cos^3\theta\big)$, $N=R_e/\sqrt{1-e^2\sin^2\phi}$, $h=p\cos\phi+(z+e^2N\sin\phi)\sin\phi-N$, where $R_e,R_p$ are the equatorial and polar radii (m), $\phi$ geodetic latitude (rad) and $h$ altitude (m).

## Design & Implementation
Inputs are `radius::Float64` (m), `u_pp::SVector{3,Float64}` (unit vector in the planet-fixed frame) and a `planet` providing equatorial radius `Rp_e` and polar radius `Rp_p` (m). The point `(x,y,z) = radius*u_pp` is converted with a Bowring-style closed-form geodetic latitude: flattening `f = (Rp_e-Rp_p)/Rp_e`, first eccentricity squared `e2 = 1-(1-f)^2`, second eccentricity squared `ep2 = e2/(1-e2)`, parametric angle `θ = atan(z*Rp_e, p_xy*Rp_p)`, then `lat = atan(z + ep2*Rp_p*sin(θ)^3, p_xy - e2*Rp_e*cos(θ)^3)`. With prime-vertical radius `N = Rp_e/sqrt(1-e2*sin(lat)^2)` the altitude is `p_xy*cos(lat) + (z + e2*N*sin(lat))*sin(lat) - N`. It is `@inline` and pure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `radius` | Float64 | n/a | yes | Positional argument `radius`. |
| in | `u_pp` | SVector{3, Float64} | n/a | yes | Positional argument `u_pp`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_oblate_altitude_from_radius`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions__radius_for_oblate_altitude|_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:70-70`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The single-pass Bowring formula is accurate to sub-millimetre for Earth-like flattening but is not iterated, so accuracy degrades for strongly oblate bodies. It assumes `Rp_e >= Rp_p`; a prolate planet yields negative `e2` and a NaN from `sqrt`. A spherical planet (`Rp_e == Rp_p`) gives `e2 = 0` and reduces cleanly to `radius - Rp_e`. `u_pp` is assumed to be unit length; it is not normalised here.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 47.
