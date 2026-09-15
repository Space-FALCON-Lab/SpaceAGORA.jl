---
id: dynx.multibody_cloth_compliant_multibody_dynamics
label: compliant_multibody_dynamics
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: compliant_multibody_dynamics
  lines:
  - 581
  - 624
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which the compliant multibody solver
    is reached from the coupled dynamics.
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Body list, joint list and actuator set defining the compliant topology.
- id: state
  type: AbstractVector
  units: m,-,m/s,rad/s
  required: true
  description: Flat 13n state vector holding position, quaternion, linear velocity
    and body angular velocity for each body.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: state_derivative
  type: Vector{Float64}
  units: m/s,1/s,m/s^2,rad/s^2
  description: Flat 13n time derivative of the compliant multibody state.
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

# compliant_multibody_dynamics

## Purpose
`compliant_multibody_dynamics` is the right-hand side of the cloth-style compliant multibody model. Instead of enforcing joints as hard algebraic constraints, every joint is a stiff six-degree-of-freedom spring-damper, so the system stays an ordinary differential equation with no constraint projection and no mass-matrix factorisation of a coupled tree.

## Theory & Math
Each body carries the Newton-Euler equations in its own frame:

$$\dot{\vec{r}}_i = \vec{v}_i, \qquad m_i \dot{\vec{v}}_i = \vec{F}_i, \qquad \mathbf{J}_i \dot{\vec{\omega}}_i = \vec{\tau}_i - \vec{\omega}_i \times \mathbf{J}_i \vec{\omega}_i$$

with $m_i$ in kg, $\mathbf{J}_i$ the body-frame inertia in kg*m^2, $\vec{F}_i$ the world-frame force in N and $\vec{\tau}_i$ the body-frame torque in N*m. Attitude propagates through the scalar-last quaternion kinematics

$$\dot{q}_i = \tfrac{1}{2}\, q_i \otimes \begin{bmatrix} \vec{\omega}_i \\ 0 \end{bmatrix}$$

Joint loads are Kelvin-Voigt in both translation and rotation:

$$\vec{f}_j = \mathbf{K}_x \left(\vec{p}_c - \vec{p}_p\right) + \mathbf{C}_x \left(\dot{\vec{p}}_c - \dot{\vec{p}}_p\right), \qquad \vec{\tau}_j = \mathbf{K}_r \vec{\theta}_{err} + \mathbf{C}_r \Delta\vec{\omega}$$

where $\vec{p}$ are the world-frame joint attachment points in m, $\mathbf{K}_x$ is in N/m, $\mathbf{C}_x$ in N*s/m, $\mathbf{K}_r$ in N*m/rad, $\mathbf{C}_r$ in N*m*s/rad, and $\vec{\theta}_{err}$ is the axis-angle error in rad extracted by `_axis_angle_error` from the relative quaternion against the rest pose. The attachment torque about each body mass centre is $\vec{r}_{att/cm} \times \vec{f}_j$ rotated into body axes by $\mathbf{R}^{\mathsf{T}}$.

## Model & Assumptions
The model assumes rigid bodies joined by massless compliant elements, small joint deflections so the linear stiffness law holds, and quaternion normalisation maintained by `_normalize_state_quaternions!` between steps rather than by a constraint. Stiffness sets the fastest timescale: with $\mathbf{K}_x = 5\times10^{3}$ N/m and a link mass near 1 kg the joint natural frequency is roughly $\sqrt{K/m} \approx 70$ rad/s, so explicit steps must stay well under about 10 ms.

## Design & Implementation
The state is a flat `13n` vector laid out as position (3), quaternion (4), linear velocity (3) and body angular velocity (3) per body, indexed with the base offset `13*(i-1)`; `compliant_state_parts` decodes one body's slice. The function first accumulates joint loads through `compliant_joint_loads`, distributing equal and opposite forces to parent and child and adding the attachment cross-product torques, then writes the derivative block for each body. Angular acceleration uses `J \ (torques_body[i] - cross(s.ω, J * s.ω))`, a linear solve rather than an explicit inverse.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which the compliant multibody solver is reached from the coupled dynamics. |
| in | `model` | CompliantMultibodyModel | n/a | yes | Body list, joint list and actuator set defining the compliant topology. |
| in | `state` | AbstractVector | m,-,m/s,rad/s | yes | Flat 13n state vector holding position, quaternion, linear velocity and body angular velocity for each body. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `state_derivative` | Vector{Float64} | m/s,1/s,m/s^2,rad/s^2 | — | Flat 13n time derivative of the compliant multibody state. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_residual|residual]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:662-662`
- [[dynamics.cloth_multibody_step_compliant_multibody_implicit_midpoint|step_compliant_multibody_implicit_midpoint]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:657-657`
- [[dynamics.cloth_multibody_step_compliant_multibody_rk4|step_compliant_multibody_rk4]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:636-636`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__parent_kinematics|_parent_kinematics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:598-598`
- `callees` → [[dynamics.cloth_multibody__quat_raw_mul|_quat_raw_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:618-618`
- `callees` → [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:595-595`
- `callees` → [[dynamics.cloth_multibody_compliant_state_parts|compliant_state_parts]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:615-615`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_raw_mul|_quat_raw_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:618-618`
<!-- vulcan:connections:end -->

## Limitations
Stiff joints make the system numerically stiff, so an explicit RK4 step can go unstable; `step_compliant_multibody_implicit_midpoint` exists for that reason. The quaternion derivative uses the raw product without renormalisation inside the step, so drift accumulates and must be corrected between steps. Joint compliance is a linearised approximation: large deflections, joint limits, backlash and friction are not represented, and the returned `dx` is heap-allocated per call.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl:581-624`, with joint loads at line 528 and the integrators at lines 627 and 646 of the same file.
