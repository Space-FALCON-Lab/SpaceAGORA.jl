---
id: gncz.thruster_guidance_functions__osculating_elements_and_periapsis_direction
label: _osculating_elements_and_periapsis_direction
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _osculating_elements_and_periapsis_direction
  lines:
  - 19
  - 45
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC guidance namespace providing the planet gravitational parameter
    and the propulsive guidance models that consume these elements.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: elements
  type: NamedTuple
  units: m, dimensionless, rad
  description: Semi-major axis, eccentricity, true anomaly, and unit periapsis direction,
    or nothing when the state is not elliptical.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# _osculating_elements_and_periapsis_direction

## Purpose
`_osculating_elements_and_periapsis_direction` converts an inertial position and velocity pair into the instantaneous Keplerian description that apoapsis-triggered thruster guidance needs. Both propulsive guidance models in this file use it to locate the vehicle within its orbit before deciding whether a burn command is due.

## Theory & Math
Specific angular momentum is $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ and specific energy is $\varepsilon = \tfrac{1}{2}v^2 - \mu/r$. A bound orbit has $\varepsilon < 0$ and semi-major axis $a = -\mu/(2\varepsilon)$. The eccentricity vector is $\mathbf{e} = (\mathbf{v} \times \mathbf{h})/\mu - \mathbf{r}/r$, whose magnitude is the eccentricity and whose direction points at periapsis. True anomaly follows from $\cos\nu = (\mathbf{e} \cdot \mathbf{r})/(er)$, resolved into the second half of the orbit when the radial rate $\mathbf{r} \cdot \mathbf{v}$ is negative.

## Model & Assumptions
Only the point-mass gravitational parameter of the planet is used, so the elements are osculating rather than mean and ignore oblateness and drag. The routine returns nothing whenever the radius, speed, or angular momentum is non-finite or zero, whenever the energy is non-negative, and whenever the resulting semi-major axis or eccentricity falls outside the elliptical range, letting callers fall back to a neutral scaling.

## Design & Implementation
For near-circular orbits, where eccentricity falls below a small threshold, true anomaly is set to zero and the periapsis direction is taken as the inward radial direction, which keeps the returned frame well defined instead of dividing by a vanishing eccentricity. Static three-vectors are used throughout so the computation stays allocation free inside the right-hand side.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC guidance namespace providing the planet gravitational parameter and the propulsive guidance models that consume these elements. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `elements` | NamedTuple | m, dimensionless, rad | — | Semi-major axis, eccentricity, true anomaly, and unit periapsis direction, or nothing when the state is not elliptical. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions__flight_apoapsis_ratio_scale|_flight_apoapsis_ratio_scale]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:105-105`
- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:161-161`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Hyperbolic and parabolic states are rejected rather than described, and no attempt is made to smooth the discontinuity in the periapsis direction as eccentricity crosses the threshold. Osculating elements computed inside a drag passage oscillate strongly and should not be treated as a pass-averaged orbit.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:19-45`.
