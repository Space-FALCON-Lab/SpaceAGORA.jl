---
id: dynx.rotational_rigid_body_dynamics_angular_acceleration
label: angular_acceleration
kind: function
source:
  file: src/dynamics/rotational/rigid_body_dynamics.jl
  symbol: angular_acceleration
  lines:
  - 1
  - 12
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which the rigid-body rotational
    kernel is reached.
- id: omega_body
  type: SVector{3,Float64}
  units: rad/s
  required: true
  description: Body-frame angular velocity of the spacecraft relative to inertial.
- id: torque_body
  type: SVector{3,Float64}
  units: N*m
  required: true
  description: Net external torque about the mass centre, expressed in body axes.
- id: inertia_tensor
  type: AbstractMatrix
  units: kg*m^2
  required: true
  description: Body-frame inertia tensor about the mass centre.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: omega_dot
  type: SVector{3,Float64}
  units: rad/s^2
  description: Body-frame angular acceleration satisfying Euler's rotational equations.
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

# angular_acceleration

## Purpose
`angular_acceleration` solves Euler's rotational equations for the body-frame angular acceleration of a rigid spacecraft, optionally including stored reaction-wheel momentum. It is the rotational counterpart of `acceleration_from_force` and, together with `quaternion_derivative`, forms the complete attitude right-hand side used by the rotational dynamics package.

## Theory & Math
The governing equation for a rigid body with internal momentum storage is

$$\mathbf{J}\,\dot{\vec{\omega}}_b = \vec{\tau}_b - \vec{\omega}_b \times \left( \mathbf{J}\vec{\omega}_b + \vec{h}_w \right)$$

where $\mathbf{J}$ is the body-frame inertia tensor in kg*m^2, $\vec{\omega}_b$ the body angular velocity in rad/s, $\vec{\tau}_b$ the net external torque in N*m and $\vec{h}_w$ the reaction-wheel angular momentum in the body frame in N*m*s. The term $\vec{\omega}_b \times \mathbf{J}\vec{\omega}_b$ is the gyroscopic coupling that makes a torque-free asymmetric body tumble: for a principal-axis inertia $\mathbf{J} = \mathrm{diag}(J_1, J_2, J_3)$ it expands to

$$J_1\dot{\omega}_1 = (J_2 - J_3)\omega_2\omega_3, \quad J_2\dot{\omega}_2 = (J_3 - J_1)\omega_3\omega_1, \quad J_3\dot{\omega}_3 = (J_1 - J_2)\omega_1\omega_2$$

These admit stable spin about the major and minor axes and instability about the intermediate axis. Setting `include_gyroscopic=false` reduces the equation to $\dot{\vec{\omega}}_b = \mathbf{J}^{-1}\vec{\tau}_b$, which is the small-rate linearisation valid when $\lVert\vec{\omega}\rVert$ is small enough that the quadratic term is negligible.

## Model & Assumptions
The kernel assumes a rigid body with a constant inertia tensor about the mass centre and body-frame torque input, and that any wheel momentum passed in has already been rotated into body axes. Sloshing propellant, flexible appendages and time-varying inertia from deployments are outside its scope. The gyroscopic term is enabled by default because disabling it breaks angular-momentum conservation for any asymmetric body.

## Design & Implementation
The implementation uses a left division `inertia_tensor \ rhs` rather than forming an explicit inverse, which is both more accurate and cheaper for a 3x3 static matrix and works for a full non-diagonal inertia. The function is `@inline` with static-array arguments so it allocates nothing inside the integrator, and the gyroscopic branch is selected by a compile-time-visible `Bool` keyword rather than a runtime lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which the rigid-body rotational kernel is reached. |
| in | `omega_body` | SVector{3,Float64} | rad/s | yes | Body-frame angular velocity of the spacecraft relative to inertial. |
| in | `torque_body` | SVector{3,Float64} | N*m | yes | Net external torque about the mass centre, expressed in body axes. |
| in | `inertia_tensor` | AbstractMatrix | kg*m^2 | yes | Body-frame inertia tensor about the mass centre. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `omega_dot` | SVector{3,Float64} | rad/s^2 | — | Body-frame angular acceleration satisfying Euler's rotational equations. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/rotational/rigid_body_dynamics.jl`
- [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1609-1609`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A singular or non-positive-definite inertia tensor is not detected, so a malformed configuration produces a non-finite acceleration instead of an error. The wheel-momentum term accounts for the reaction torque only through $\vec{h}_w$; wheel acceleration torque must be supplied separately in $\vec{\tau}_b$. Structural flexibility and inertia variation during a burn are not represented.

## Provenance
Mapped from `src/dynamics/rotational/rigid_body_dynamics.jl:1-12`, exported by `DynamicsRotational` at line 11 of `src/dynamics/rotational/rotational_models.jl`.
