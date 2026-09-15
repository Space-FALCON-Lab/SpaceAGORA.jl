---
id: dynamics.cloth_multibody__actuator_torque_child_world
label: _actuator_torque_child_world
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _actuator_torque_child_world
  lines:
  - 515
  - 515
inputs:
- id: actuator
  type: CompliantJointActuator
  units: n/a
  required: true
  description: Positional argument `actuator`.
- id: child
  type: Any
  units: n/a
  required: true
  description: Positional argument `child`.
- id: phi_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `ϕ_world`.
- id: Deltaomega_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `Δω_world`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_actuator_torque_child_world`.
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

# _actuator_torque_child_world

## Purpose
Evaluates one actuator's torque on the child body in world coordinates, applying PD gains, feedforward, saturation and efficiency.

## Theory & Math
$$
\tau_{pd} = K_p \vec{\phi} - K_d \Delta\vec{\omega},\qquad \tau_c = \eta\, R\, \operatorname{clamp}\left(R^\top \tau_{pd} + \tau_{ff},\; \pm\tau_{\max}\right)
$$

## Design & Implementation
Forms the PD torque in world frame from the axis-angle error and relative angular rate, rotates it into the child body frame, adds the feedforward, clamps each axis to the torque limit, scales by efficiency, and rotates back to world. Returns the world-frame torque.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `actuator` | CompliantJointActuator | n/a | yes | Positional argument `actuator`. |
| in | `child` | Any | n/a | yes | Positional argument `child`. |
| in | `phi_world` | SVector{3, Float64} | n/a | yes | Positional argument `ϕ_world`. |
| in | `Deltaomega_world` | SVector{3, Float64} | n/a | yes | Positional argument `Δω_world`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_actuator_torque_child_world`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:555-555`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Saturation is per axis in the child frame, so a diagonal command can exceed the scalar limit by up to a factor of the square root of three.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 515.
