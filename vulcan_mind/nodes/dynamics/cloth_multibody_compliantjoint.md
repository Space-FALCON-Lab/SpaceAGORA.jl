---
id: dynamics.cloth_multibody_compliantjoint
label: CompliantJoint
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantJoint
  lines:
  - 27
  - 27
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: parent
  type: Int
  units: n/a
  required: true
  description: Field `parent`.
- id: child
  type: Int
  units: n/a
  required: true
  description: Field `child`.
- id: parent_point_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `parent_point_body`.
- id: child_point_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `child_point_body`.
- id: k_translation_n_m
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `k_translation_n_m`.
- id: c_translation_n_s_m
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `c_translation_n_s_m`.
- id: k_rotation_n_m_rad
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `k_rotation_n_m_rad`.
- id: c_rotation_n_m_s_rad
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `c_rotation_n_m_s_rad`.
- id: rest_child_parent_quat
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `rest_child_parent_quat`.
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
  type: CompliantJoint
  units: n/a
  description: Constructed `CompliantJoint`.
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

# CompliantJoint

## Purpose
A six-axis spring-damper connection between a parent body (or the fixed base) and a child body, defined by attachment points and stiffness and damping matrices.

## Design & Implementation
Immutable with `parent` and `child` indices (parent zero meaning the base), attachment points in each body's frame, three-by-three translational and rotational stiffness and damping matrices, and `rest_child_parent_quat`, the relative orientation at which the rotational spring is unloaded. Matrices rather than scalars allow anisotropic compliance.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `parent` | Int | n/a | yes | Field `parent`. |
| in | `child` | Int | n/a | yes | Field `child`. |
| in | `parent_point_body` | SVector{3, Float64} | n/a | yes | Field `parent_point_body`. |
| in | `child_point_body` | SVector{3, Float64} | n/a | yes | Field `child_point_body`. |
| in | `k_translation_n_m` | SMatrix{3, 3, Float64} | n/a | yes | Field `k_translation_n_m`. |
| in | `c_translation_n_s_m` | SMatrix{3, 3, Float64} | n/a | yes | Field `c_translation_n_s_m`. |
| in | `k_rotation_n_m_rad` | SMatrix{3, 3, Float64} | n/a | yes | Field `k_rotation_n_m_rad`. |
| in | `c_rotation_n_m_s_rad` | SMatrix{3, 3, Float64} | n/a | yes | Field `c_rotation_n_m_s_rad`. |
| in | `rest_child_parent_quat` | SVector{4, Float64} | n/a | yes | Field `rest_child_parent_quat`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantJoint | n/a | — | Constructed `CompliantJoint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:303-303`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:257-257`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Rest orientation is a single quaternion, so the rotational spring is linear in axis-angle error and only meaningful for moderate deflections; there is no representation of a rest translational offset other than coincident attachment points.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 27.
