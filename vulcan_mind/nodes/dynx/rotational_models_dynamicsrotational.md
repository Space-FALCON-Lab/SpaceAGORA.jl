---
id: dynx.rotational_models_dynamicsrotational
label: DynamicsRotational
kind: struct
source:
  file: src/dynamics/rotational/rotational_models.jl
  symbol: DynamicsRotational
  lines:
  - 1
  - 16
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace that consumes the rotational kernels aggregated
    by this module.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rotational_api
  type: Module
  units: n/a
  description: 'Exported surface: quaternion_derivative, angular_acceleration, body_torque,
    body_angular_velocity and combine_torques.'
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

# DynamicsRotational

## Purpose
`DynamicsRotational` is the aggregating module for the attitude side of the dynamics package. It includes the three rotational kernel files and re-exports the five functions that together form the attitude right-hand side, giving callers one namespace for quaternion kinematics, Euler dynamics and torque marshalling.

## Theory & Math
The module's exported surface implements the coupled attitude system

$$\dot{q} = \tfrac{1}{2}\, q \otimes \begin{bmatrix}\vec{\omega}_b\\0\end{bmatrix}, \qquad \mathbf{J}\dot{\vec{\omega}}_b = \vec{\tau}_b - \vec{\omega}_b\times\left(\mathbf{J}\vec{\omega}_b + \vec{h}_w\right)$$

with $q$ a scalar-last unit quaternion, $\vec{\omega}_b$ in rad/s, $\mathbf{J}$ in kg*m^2, $\vec{\tau}_b$ in N*m and $\vec{h}_w$ in N*m*s. The net torque fed to the second equation is assembled as $\vec{\tau}_b = \vec{\tau}_{dyn} + \vec{\tau}_{ctrl}$ by `combine_torques`, where $\vec{\tau}_{dyn}$ collects environmental contributions (aerodynamic, magnetic, gravity-gradient, solar) and $\vec{\tau}_{ctrl}$ the actuator command. Together the two equations conserve inertial angular momentum $\vec{h}_{ii} = \mathbf{R}(q)\,(\mathbf{J}\vec{\omega}_b + \vec{h}_w)$ when $\vec{\tau}_b = 0$, which is the invariant the package's probe suite checks.

## Model & Assumptions
The module assumes every kernel it aggregates uses body-frame rates and body-frame torques with a scalar-last quaternion, and that the inertia tensor is constant. It deliberately keeps no state: all five exported functions are pure, which makes them safe to call from multiple threads propagating different spacecraft concurrently.

## Design & Implementation
Includes are resolved with `joinpath(@__DIR__, ...)` so load order is deterministic and independent of the process working directory: kinematics first, then rigid-body dynamics, then the torque marshalling helpers. Only `LinearAlgebra` and `StaticArrays` are pulled in, keeping the module free of simulation-level dependencies so it can be unit tested in isolation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace that consumes the rotational kernels aggregated by this module. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rotational_api` | Module | n/a | — | Exported surface: quaternion_derivative, angular_acceleration, body_torque, body_angular_velocity and combine_torques. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/rotational/rotational_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The module exports no attitude representation other than the scalar-last quaternion, so callers using Euler angles or Modified Rodrigues Parameters must convert externally. It offers no quaternion renormalisation entry point, and no validation that the inertia tensor is symmetric positive definite. Flexible-body and variable-inertia dynamics are outside its scope.

## Provenance
Mapped from `src/dynamics/rotational/rotational_models.jl:1-16` and the three files it includes under `src/dynamics/rotational/`.
