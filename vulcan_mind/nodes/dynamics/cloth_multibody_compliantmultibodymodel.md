---
id: dynamics.cloth_multibody_compliantmultibodymodel
label: CompliantMultibodyModel
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantMultibodyModel
  lines:
  - 41
  - 41
inputs:
- id: bodies
  type: Vector{CompliantBody}
  units: n/a
  required: true
  description: Field `bodies`.
- id: joints
  type: Vector{CompliantJoint}
  units: n/a
  required: true
  description: Field `joints`.
- id: base_position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `base_position`.
- id: base_quaternion
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `base_quaternion`.
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
  type: CompliantMultibodyModel
  units: n/a
  description: Constructed `CompliantMultibodyModel`.
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

# CompliantMultibodyModel

## Purpose
The complete compliant system: bodies, joints and the pose of the fixed base that anchor joints attach to.

## Design & Implementation
Immutable, holding `Vector{CompliantBody}`, `Vector{CompliantJoint}`, and the base position and unit quaternion. Body index zero in a joint refers to this base, which has zero velocity and angular rate in `_parent_kinematics`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bodies` | Vector{CompliantBody} | n/a | yes | Field `bodies`. |
| in | `joints` | Vector{CompliantJoint} | n/a | yes | Field `joints`. |
| in | `base_position` | SVector{3, Float64} | n/a | yes | Field `base_position`. |
| in | `base_quaternion` | SVector{4, Float64} | n/a | yes | Field `base_quaternion`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantMultibodyModel | n/a | — | Constructed `CompliantMultibodyModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:323-323`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:270-270`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The base is rigidly fixed in the world frame; a model attached to a moving spacecraft must transform loads externally, which `cloth_robot_arm_dynamics.jl` does.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 41.
