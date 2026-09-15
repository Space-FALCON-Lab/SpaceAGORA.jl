---
id: dynx.rotational_attitude_kinematics_quaternion_derivative
label: quaternion_derivative
kind: function
source:
  file: src/dynamics/rotational/attitude_kinematics.jl
  symbol: quaternion_derivative
  lines:
  - 16
  - 32
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which the rotational kinematics
    kernel is reached.
- id: omega_body
  type: SVector{3,Float64}
  units: rad/s
  required: true
  description: Angular velocity of the body frame relative to inertial, expressed
    in body axes.
- id: quaternion
  type: AbstractVector
  units: n/a
  required: true
  description: Scalar-last unit quaternion rotating body axes into inertial axes.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: q_dot
  type: SVector{4,Float64}
  units: 1/s
  description: Time derivative of the scalar-last attitude quaternion.
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

# quaternion_derivative

## Purpose
`quaternion_derivative` is the attitude kinematics kernel for the rotational dynamics package. Given a body-frame angular velocity and the current scalar-last quaternion it returns the quaternion rate, which the integrator advances alongside the Euler equations. It is the single definition of attitude propagation shared by every rotational model in the package.

## Theory & Math
For a body rate expressed in body axes, the quaternion kinematics are

$$\dot{q} = \tfrac{1}{2}\, q \otimes \begin{bmatrix} \vec{\omega}_b \\ 0 \end{bmatrix}$$

which, with the scalar-last convention $q = [\vec{q}_v, q_s]$, expands to

$$\dot{\vec{q}}_v = \tfrac{1}{2}\left( q_s \vec{\omega}_b - \vec{\omega}_b \times \vec{q}_v \right), \qquad \dot{q}_s = -\tfrac{1}{2}\, \vec{\omega}_b \cdot \vec{q}_v$$

exactly as implemented. Here $\vec{\omega}_b$ is in rad/s, $\vec{q}_v$ and $q_s$ are dimensionless, and $\dot{q}$ has units of 1/s. The sign of the cross-product term encodes the frame: the body-rate composition uses $-\vec{\omega}_b \times \vec{q}_v$, while an inertial-rate composition would use $+\vec{\omega} \times \vec{q}_v$. Because $q^{\mathsf{T}}\dot{q} = 0$ identically, the unit-norm constraint $\lVert q \rVert = 1$ is preserved to integrator accuracy rather than by projection.

## Model & Assumptions
The kernel assumes a unit quaternion on input and a body-frame rate consistent with the Euler equations solved by `angular_acceleration`, which uses body-frame torques and the body-frame gyroscopic term. Mixing frames here is not a cosmetic error: the docstring records that the pre-2026-07 inertial-rate form drifted inertial angular momentum by 68 percent over 2000 s on a torque-free asymmetric tumble with I = diag(1.5, 1.0, 2.0) and |omega| about 0.04 rad/s, while the body-rate form conserves it to integrator accuracy. A probe test pins that invariant.

## Design & Implementation
The function is marked `@inline` and takes an `SVector{3,Float64}` rate with a loosely typed quaternion, converting components to `Float64` explicitly so it accepts state views and plain vectors without allocating. It returns an `SVector{4,Float64}` in scalar-last order, matching the state layout used across the translational and rotational packages, and performs no normalisation of its own.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which the rotational kinematics kernel is reached. |
| in | `omega_body` | SVector{3,Float64} | rad/s | yes | Angular velocity of the body frame relative to inertial, expressed in body axes. |
| in | `quaternion` | AbstractVector | n/a | yes | Scalar-last unit quaternion rotating body axes into inertial axes. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `q_dot` | SVector{4,Float64} | 1/s | — | Time derivative of the scalar-last attitude quaternion. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1598-1598`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/rotational/attitude_kinematics.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
No unit-norm enforcement is applied, so quaternion drift must be corrected by the caller or absorbed by a sufficiently tight integrator tolerance. Passing an inertial-frame angular velocity produces a silently wrong but plausible-looking trajectory. The scalar-last ordering is not self-describing, so a scalar-first quaternion supplied by external code yields an incorrect rate with no error raised.

## Provenance
Mapped from `src/dynamics/rotational/attitude_kinematics.jl:16-32`, exported by `DynamicsRotational` at line 10 of `src/dynamics/rotational/rotational_models.jl`.
