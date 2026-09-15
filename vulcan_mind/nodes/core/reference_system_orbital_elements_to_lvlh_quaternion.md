---
id: core.reference_system_orbital_elements_to_lvlh_quaternion
label: orbital_elements_to_lvlh_quaternion
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: orbital_elements_to_lvlh_quaternion
  lines:
  - 470
  - 470
inputs:
- id: raan
  type: Float64
  units: n/a
  required: true
  description: Positional argument `raan`.
- id: inclination
  type: Float64
  units: n/a
  required: true
  description: Positional argument `inclination`.
- id: arg_of_perigee
  type: Float64
  units: n/a
  required: true
  description: Positional argument `arg_of_perigee`.
- id: true_anomaly
  type: Float64
  units: n/a
  required: true
  description: Positional argument `true_anomaly`.
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
  type: SVector{4,
  units: n/a
  description: Return value of `orbital_elements_to_lvlh_quaternion`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# orbital_elements_to_lvlh_quaternion

## Purpose
Computes the quaternion rotating inertial vectors into the LVLH frame — nadir z, anti-angular-momentum y — from four orbital angles.

## Theory & Math
$$
\hat{r} = \begin{bmatrix} \cos\Omega\cos u - \sin\Omega\sin u\cos i \\ \sin\Omega\cos u + \cos\Omega\sin u\cos i \\ \sin u \sin i \end{bmatrix},\quad \hat{h} = \begin{bmatrix} \sin\Omega\sin i \\ -\cos\Omega\sin i \\ \cos i \end{bmatrix},\quad u = \omega + \nu
$$

## Design & Implementation
Forms the argument of latitude, builds the inertial radial and angular-momentum unit vectors from the angles, sets `z = -r̂`, `y = -ĥ` and `x = y × z` normalised, and stacks them as rows of a direction-cosine matrix. It then converts the DCM to a quaternion with the trace-based branch selection that picks the largest diagonal term for numerical stability, and normalises the result. Returns scalar-last `[qx, qy, qz, qw]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raan` | Float64 | n/a | yes | Positional argument `raan`. |
| in | `inclination` | Float64 | n/a | yes | Positional argument `inclination`. |
| in | `arg_of_perigee` | Float64 | n/a | yes | Positional argument `arg_of_perigee`. |
| in | `true_anomaly` | Float64 | n/a | yes | Positional argument `true_anomaly`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `orbital_elements_to_lvlh_quaternion`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_latlongtoned|latlongtoNED]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:449-449`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system_rotate_vector_by_quaternion|rotate_vector_by_quaternion]] · `callers` · call · `src/core/interfaces/reference_system.jl:555-555`
<!-- vulcan:connections:end -->

## Limitations
The LVLH x-axis is only along velocity for circular orbits; for eccentric orbits it is horizontal but not aligned with velocity. The intermediate vectors are constructed via array literals, which allocate before conversion to static.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 470.
