---
id: dynamics.cloth_multibody__joint_rest
label: _joint_rest
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _joint_rest
  lines:
  - 510
  - 510
inputs:
- id: joint
  type: CompliantJoint
  units: n/a
  required: true
  description: Positional argument `joint`.
- id: jidx
  type: Int
  units: n/a
  required: true
  description: Positional argument `jidx`.
- id: joint_rest_quaternions
  type: Any
  units: n/a
  required: true
  description: Positional argument `joint_rest_quaternions`.
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
  description: 'Return value of `_joint_rest`. Returns `joint_rest_quaternions ===
    nothing ? joint.rest_child_parent_quat : _unit_quat(j`.'
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

# _joint_rest

## Purpose
Selects the rest orientation for a joint, from the joint itself or from an optional per-joint override vector.

## Design & Implementation
Returns `joint.rest_child_parent_quat` when `joint_rest_quaternions` is `nothing`, else the normalised entry at `jidx`. `@inline`. The override is how a controller commands a new joint angle without rebuilding the model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `joint` | CompliantJoint | n/a | yes | Positional argument `joint`. |
| in | `jidx` | Int | n/a | yes | Positional argument `jidx`. |
| in | `joint_rest_quaternions` | Any | n/a | yes | Positional argument `joint_rest_quaternions`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_joint_rest`. Returns `joint_rest_quaternions === nothing ? joint.rest_child_parent_quat : _unit_quat(j`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:546-546`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:511-511`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:511-511`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:511-511`
<!-- vulcan:connections:end -->

## Limitations
The override vector must be indexed identically to `model.joints`; a shorter vector fails with a bounds error rather than a message.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 510.
