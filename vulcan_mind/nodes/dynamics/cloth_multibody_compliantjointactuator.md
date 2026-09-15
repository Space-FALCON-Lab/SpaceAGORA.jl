---
id: dynamics.cloth_multibody_compliantjointactuator
label: CompliantJointActuator
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantJointActuator
  lines:
  - 88
  - 88
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: joint
  type: Int
  units: n/a
  required: true
  description: Field `joint`.
- id: torque_limit_n_m
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `torque_limit_n_m`.
- id: kp_n_m_rad
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `kp_n_m_rad`.
- id: kd_n_m_s_rad
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `kd_n_m_s_rad`.
- id: feedforward_torque_child_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `feedforward_torque_child_body`.
- id: efficiency
  type: Float64
  units: n/a
  required: true
  description: Field `efficiency`.
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
  type: CompliantJointActuator
  units: n/a
  description: Constructed `CompliantJointActuator`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# CompliantJointActuator

## Purpose
A PD torque actuator acting across one joint, with feedforward, saturation and efficiency, used to drive joints toward their rest orientation or a commanded one.

## Design & Implementation
Immutable with the joint index, per-axis `torque_limit_n_m`, three-by-three `kp` and `kd` gain matrices, a feedforward torque in the child body frame and an `efficiency` scalar. The keyword constructor accepts a scalar limit broadcast to three axes or a vector taken by absolute value, requires non-negative limits and efficiency, and defaults the limit to infinite and the gains to zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `joint` | Int | n/a | yes | Field `joint`. |
| in | `torque_limit_n_m` | SVector{3, Float64} | n/a | yes | Field `torque_limit_n_m`. |
| in | `kp_n_m_rad` | SMatrix{3, 3, Float64} | n/a | yes | Field `kp_n_m_rad`. |
| in | `kd_n_m_s_rad` | SMatrix{3, 3, Float64} | n/a | yes | Field `kd_n_m_s_rad`. |
| in | `feedforward_torque_child_body` | SVector{3, Float64} | n/a | yes | Field `feedforward_torque_child_body`. |
| in | `efficiency` | Float64 | n/a | yes | Field `efficiency`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantJointActuator | n/a | — | Constructed `CompliantJointActuator`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:249-249`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_actuators|cloth_robot_arm_actuators]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:283-283`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Gains act on the world-frame axis-angle error, but the limit and feedforward apply in the child body frame, so the saturation box rotates with the child; the interplay is correct but easy to misconfigure.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 88.
