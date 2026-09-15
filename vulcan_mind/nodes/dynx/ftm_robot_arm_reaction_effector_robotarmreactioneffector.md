---
id: dynx.ftm_robot_arm_reaction_effector_robotarmreactioneffector
label: RobotArmReactionEffector
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl
  symbol: RobotArmReactionEffector
  lines:
  - 13
  - 24
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace under which this AbstractForceTorqueModel
    subtype is registered.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: effector_config
  type: RobotArmReactionEffector
  units: N/m,N*s/m,N*m/rad,N*m*s/rad
  description: Mutable effector carrying the arm plan, joint compliance gains and
    the actuator list used to size the base reaction wrench.
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

# RobotArmReactionEffector

## Purpose
`RobotArmReactionEffector` is the mutable effector that feeds robot-arm joint reaction loads back onto the spacecraft base. It holds the target spacecraft index, the current `RobotArmPlan`, the epoch at which that plan was refreshed, translational and rotational compliance gains, force and torque scale factors, and the compliant joint actuator list shared with the cloth multibody solver.

## Theory & Math
A manipulator mounted on a free-flying base exchanges momentum with it. For an arm of $n$ links the reaction wrench on the base is

$$\vec{F}_{base} = -\sum_{k=1}^{n} m_k \ddot{\vec{r}}_k, \qquad \vec{\tau}_{base} = -\sum_{k=1}^{n} \left( \mathbf{I}_k \dot{\vec{\omega}}_k + \vec{\omega}_k \times \mathbf{I}_k \vec{\omega}_k + \vec{r}_{k/base} \times m_k \ddot{\vec{r}}_k \right)$$

with $m_k$ the link mass in kg, $\ddot{\vec{r}}_k$ the link mass-centre acceleration in m/s^2, $\mathbf{I}_k$ the link inertia in kg*m^2, and $\vec{\omega}_k$ the link angular velocity in rad/s. The compliance gains stored on the struct realise a Kelvin-Voigt joint model

$$\vec{f}_j = \mathbf{K}_x \Delta\vec{x} + \mathbf{C}_x \Delta\dot{\vec{x}}, \qquad \vec{\tau}_j = \mathbf{K}_r \vec{\theta}_{err} + \mathbf{C}_r \Delta\vec{\omega}$$

with $\mathbf{K}_x$ in N/m (default 5.0e3), $\mathbf{C}_x$ in N*s/m (default 30.0), $\mathbf{K}_r$ in N*m/rad (default 15.0) and $\mathbf{C}_r$ in N*m*s/rad (default 0.5); $\vec{\theta}_{err}$ is the axis-angle joint orientation error in rad.

## Model & Assumptions
The effector assumes a single mounting spacecraft, selected by `spacecraft_idx`, and that the arm plan is refreshed externally with `updated_at_s` tracking staleness. Compliance gains are isotropic scalars by default but accept full 3x3 matrices, so anisotropic joints are representable. The `force_scale` and `torque_scale` fields default to zero, which disables the reaction path until a scenario opts in.

## Design & Implementation
It is declared with `Base.@kwdef mutable struct ... <: AbstractForceTorqueModel` so scenarios construct it by keyword and the plan can be mutated in place between right-hand-side evaluations without reallocating the effector. The companion `calcForceTorque` method at lines 27-33 short-circuits when the satellite index does not match or the plan is unset, and currently returns a zero wrench because reaction estimation is owned by the cloth/RNEA dynamics path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace under which this AbstractForceTorqueModel subtype is registered. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `effector_config` | RobotArmReactionEffector | N/m,N*s/m,N*m/rad,N*m*s/rad | — | Mutable effector carrying the arm plan, joint compliance gains and the actuator list used to size the base reaction wrench. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The wrench computation is a placeholder: the current `calcForceTorque` method returns zero force and zero torque even with a valid plan, so momentum exchange reaches the base only through the coupled cloth robot-arm right-hand side. Gains are not checked for positive definiteness, and a negative damping value would inject energy. Only one spacecraft per effector instance is supported.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl:13-24`, with the dispatch method at lines 27-33 of the same file.
