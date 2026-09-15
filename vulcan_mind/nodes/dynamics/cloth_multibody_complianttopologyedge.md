---
id: dynamics.cloth_multibody_complianttopologyedge
label: CompliantTopologyEdge
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantTopologyEdge
  lines:
  - 66
  - 66
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
  type: Union{Nothing, SVector{4, Float64}}
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
  type: CompliantTopologyEdge
  units: n/a
  description: Constructed `CompliantTopologyEdge`.
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

# CompliantTopologyEdge

## Purpose
Specification of one joint in a topology before it is built, with optional explicit rest orientation.

## Design & Implementation
Immutable mirroring `CompliantJoint` but with `rest_child_parent_quat` as `Union{Nothing, SVector{4}}`. The keyword outer constructor accepts scalar or matrix stiffness and damping via `_mat3`, defaulting all four to zero, and normalises a supplied rest quaternion. When rest is `nothing`, the builder computes it from the nodes' initial orientations.

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
| in | `rest_child_parent_quat` | Union{Nothing, SVector{4, Float64}} | n/a | yes | Field `rest_child_parent_quat`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantTopologyEdge | n/a | — | Constructed `CompliantTopologyEdge`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:388-388`
- [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:223-223`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
All-zero default stiffness means an edge constructed with only attachment points is a free hinge in every axis, which is rarely intended and produces no warning.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 66.
