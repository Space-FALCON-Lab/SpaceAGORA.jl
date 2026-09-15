---
id: vehx.kinematics_rotate_to_inertial
label: rotate_to_inertial
kind: function
source:
  file: src/vehicle/kinematics/kinematics.jl
  symbol: rotate_to_inertial
  lines:
  - 10
  - 20
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Assembled vehicle whose root attitude anchors every child link frame.
- id: body
  type: Link
  units: n/a
  required: true
  description: Link whose stored quaternion is being resolved into inertial axes.
- id: root_index
  type: Int
  units: n/a
  required: true
  description: Assembly index selecting which root chain the link belongs to.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: dcm
  type: SMatrix{3,3,Float64}
  units: n/a
  description: Direction cosine matrix mapping link-frame vectors into the inertial
    frame.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- kinematics
charts:
- vehx
origin: agent
---

# rotate_to_inertial

## Purpose
`rotate_to_inertial` produces the rotation that carries a vector expressed in a link frame into the inertial frame. Geometry, aerodynamic normals, solar panel areas and thruster directions are all stored in link axes, so every environment and dynamics computation that needs an inertial quantity passes through this function. It is the counterpart of `rotate_to_body`, and the two together define the frame convention for the whole vehicle package.

## Theory & Math
For a unit quaternion $\mathbf{q}$ the attitude matrix $R(\mathbf{q})$ maps inertial components to frame components, so the inverse map is $R^{\mathsf T}$. For a child link the composition is $C^{\mathcal I \leftarrow \mathcal L} = R(\mathbf{q}_{root})^{\mathsf T} R(\mathbf{q}_{link})^{\mathsf T}$, which is orthonormal with $\det = +1$ whenever both quaternions satisfy $\|\mathbf{q}\| = 1$.

## Model & Assumptions
Attitude is represented by unit quaternions. The root link stores its attitude directly with respect to the inertial frame, while every child link stores its attitude relative to the root body frame. The function encodes that convention explicitly: a root body returns the transpose of its own direction cosine matrix, and a child body returns the root rotation composed with the child rotation. The quaternion is assumed normalised; no renormalisation is performed here, so integration drift must be corrected upstream.

## Design & Implementation
The branch on `body.root` at line 13 is the whole decision. `rot` comes from the quaternion utilities included at the top of the module from the core numerics directory, and it returns the inertial-to-frame matrix, so the transpose is taken to obtain the frame-to-inertial direction. Composition for a child is written as `rot(model.root.q)' * rot(body.q)'`, applying the child rotation first and then the root rotation. Static matrices are used throughout so the product stays allocation free and can be inlined into tight dynamics loops.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `model` | SpacecraftModel | n/a | yes | Assembled vehicle whose root attitude anchors every child link frame. |
| in | `body` | Link | n/a | yes | Link whose stored quaternion is being resolved into inertial axes. |
| in | `root_index` | Int | n/a | yes | Assembly index selecting which root chain the link belongs to. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `dcm` | SMatrix{3,3,Float64} | n/a | — | Direction cosine matrix mapping link-frame vectors into the inertial frame. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:623-623`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:782-782`
- [[vehicle.geometry_properties_get_normal_vector|get_normal_vector]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:114-114`
- [[vehicle.geometry_properties_get_tangent_vector|get_tangent_vector]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:129-129`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/vehicle/kinematics/kinematics.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
Only a single level of nesting is supported: a link attached to another non-root link is still resolved through the root attitude, so deep kinematic chains are not composed joint by joint. The `root_index` argument is accepted for signature compatibility with the multi-assembly structure helpers but is not used in the body. No check is made that the quaternion is unit norm or that the link actually belongs to the supplied model.

## Provenance
Mapped from `src/vehicle/kinematics/kinematics.jl:10-20`; `rot` is defined in the quaternion utilities under `src/core/numerics/`.
