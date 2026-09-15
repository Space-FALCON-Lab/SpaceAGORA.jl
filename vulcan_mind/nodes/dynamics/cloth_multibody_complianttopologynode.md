---
id: dynamics.cloth_multibody_complianttopologynode
label: CompliantTopologyNode
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantTopologyNode
  lines:
  - 55
  - 55
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `mass_kg`.
- id: inertia_body_kg_m2
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `inertia_body_kg_m2`.
- id: position
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `position`.
- id: quaternion
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `quaternion`.
- id: velocity
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `velocity`.
- id: angular_velocity
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `angular_velocity`.
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
  type: CompliantTopologyNode
  units: n/a
  description: Constructed `CompliantTopologyNode`.
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

# CompliantTopologyNode

## Purpose
Specification of one body in a topology before it is built: mass, inertia and initial pose and rates.

## Design & Implementation
Immutable with `name`, `mass_kg`, inertia, position, unit quaternion, velocity and angular velocity. A keyword outer constructor validates positive mass, accepts scalar or matrix inertia through `_mat3`, converts tuples to static vectors, and defaults quaternion to identity and rates to zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `mass_kg` | Float64 | n/a | yes | Field `mass_kg`. |
| in | `inertia_body_kg_m2` | SMatrix{3, 3, Float64} | n/a | yes | Field `inertia_body_kg_m2`. |
| in | `position` | SVector{3, Float64} | n/a | yes | Field `position`. |
| in | `quaternion` | SVector{4, Float64} | n/a | yes | Field `quaternion`. |
| in | `velocity` | SVector{3, Float64} | n/a | yes | Field `velocity`. |
| in | `angular_velocity` | SVector{3, Float64} | n/a | yes | Field `angular_velocity`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantTopologyNode | n/a | — | Constructed `CompliantTopologyNode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:375-375`
- [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:201-201`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only mass positivity is checked; a zero inertia passes and later makes the angular dynamics singular.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 55.
