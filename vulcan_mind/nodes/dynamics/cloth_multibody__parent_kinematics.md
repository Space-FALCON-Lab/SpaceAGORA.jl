---
id: dynamics.cloth_multibody__parent_kinematics
label: _parent_kinematics
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _parent_kinematics
  lines:
  - 479
  - 479
inputs:
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
- id: parent
  type: Int
  units: n/a
  required: true
  description: Positional argument `parent`.
- id: p_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `p_body`.
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
  type: Any
  units: n/a
  description: Return value of `_parent_kinematics`. Returns `(`.
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

# _parent_kinematics

## Purpose
Gathers everything the joint force law needs about one side of a joint: the body's pose and rates, its rotation matrix, and the attachment point's world position and velocity.

## Design & Implementation
For `parent == 0` it returns the model's fixed base with zero rates. Otherwise it unpacks the body's state, forms `R`, and computes the attachment point as `r + R p_body` and its velocity through `_body_offset_velocity`. Returns a named tuple. Despite the name it is used for both the parent and child sides.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `parent` | Int | n/a | yes | Positional argument `parent`. |
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_parent_kinematics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:538-538`
- [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:598-598`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__body_offset_velocity|_body_offset_velocity]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:505-505`
- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:485-485`
- `callees` → [[dynamics.cloth_multibody_compliant_state_parts|compliant_state_parts]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:496-496`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:485-485`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:485-485`
<!-- vulcan:connections:end -->

## Limitations
Called twice per joint in `compliant_joint_loads` and twice again in `compliant_multibody_dynamics` for the same state, so each body's kinematics are rebuilt up to four times per derivative.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 479.
