---
id: dynx.coupled_perturbations_calculate_magnetic_torque
label: calculate_magnetic_torque
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: calculate_magnetic_torque
  lines:
  - 1828
  - 1838
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace under which magnetic effectors resolve this
    helper.
- id: dipole_moment
  type: AbstractVector
  units: A*m^2
  required: true
  description: Commanded or residual magnetic dipole moment of the spacecraft, expressed
    in the same frame as the field.
- id: field_vector
  type: AbstractVector
  units: T
  required: true
  description: Local geomagnetic flux density from the tilted-dipole or IGRF field
    model.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: torque
  type: SVector{3,Float64}
  units: N*m
  description: Magnetic torque acting on the spacecraft, in the common frame of the
    dipole and field inputs.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynx
origin: agent
---

# calculate_magnetic_torque

## Purpose
`calculate_magnetic_torque` evaluates the torque a magnetic dipole experiences in an ambient field. It is the shared kernel behind the magnetic torque rod effector, the eddy-current damping model and any residual-dipole disturbance, so the sign convention and frame handling are defined exactly once.

## Theory & Math
The governing relation is the classical dipole torque law

$$\vec{\tau} = \vec{m} \times \vec{B}$$

with $\vec{m}$ the magnetic dipole moment in A*m^2 and $\vec{B}$ the magnetic flux density in T, giving $\vec{\tau}$ in N*m. Its magnitude is $\lVert \vec{\tau} \rVert = m B \sin\theta$, where $\theta$ is the angle in rad between the dipole and the field, so torque vanishes when the dipole is aligned or anti-aligned with the field. The associated potential energy is $U = -\vec{m}\cdot\vec{B}$, and the torque is its negative gradient with respect to orientation. Because $\vec{\tau} \perp \vec{B}$ always, a magnetic actuator is instantaneously rank-two: no torque can be produced about the local field line, which is the fundamental underactuation of magnetorquer-only attitude control. Typical low-Earth-orbit values are $B \approx 2$--$5\times10^{-5}$ T and $m \approx 0.1$--$10$ A*m^2, giving torques of order $10^{-6}$ to $10^{-4}$ N*m.

## Model & Assumptions
Both inputs must be expressed in the same frame; the function performs no rotation, so the caller is responsible for having already converted the field from planet-fixed to body axes. It assumes the dipole is a point moment, valid when the spacecraft dimension is small compared with the field gradient scale, which holds comfortably in orbit. The field is treated as uniform across the vehicle, so no force term $\nabla(\vec{m}\cdot\vec{B})$ is returned.

## Design & Implementation
The implementation converts both arguments to `SVector{3,Float64}` before taking the cross product, which keeps the operation allocation-free and type-stable inside the integrator right-hand side even when callers pass plain `Vector` arguments. Returning a static vector lets the caller accumulate the result into a body-torque sum without heap traffic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace under which magnetic effectors resolve this helper. |
| in | `dipole_moment` | AbstractVector | A*m^2 | yes | Commanded or residual magnetic dipole moment of the spacecraft, expressed in the same frame as the field. |
| in | `field_vector` | AbstractVector | T | yes | Local geomagnetic flux density from the tilted-dipole or IGRF field model. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `torque` | SVector{3,Float64} | N*m | — | Magnetic torque acting on the spacecraft, in the common frame of the dipole and field inputs. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2064-2064`
- [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1808-1808`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No frame checking is performed, so a mismatch between a planet-fixed field and a body-frame dipole silently produces a wrong torque. The uniform-field assumption means the small magnetic force on a dipole in a field gradient is not returned. Saturation, hysteresis and coil dynamics of a physical torque rod are handled by the caller, not here.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl:1828-1838`, used by the effector methods at lines 1892 and 2033 of the same file.
