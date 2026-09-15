---
id: dynamics.cloth_multibody_compliantjointload
label: CompliantJointLoad
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantJointLoad
  lines:
  - 99
  - 99
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
- id: translation_force_parent_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `translation_force_parent_world`.
- id: translation_force_child_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `translation_force_child_world`.
- id: compliance_torque_parent_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `compliance_torque_parent_world`.
- id: compliance_torque_child_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `compliance_torque_child_world`.
- id: actuator_torque_parent_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `actuator_torque_parent_world`.
- id: actuator_torque_child_world
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `actuator_torque_child_world`.
- id: compliance_torque_child_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `compliance_torque_child_body`.
- id: actuator_torque_child_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `actuator_torque_child_body`.
- id: total_torque_child_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `total_torque_child_body`.
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
  type: CompliantJointLoad
  units: n/a
  description: Constructed `CompliantJointLoad`.
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

# CompliantJointLoad

## Purpose
Per-joint force and torque diagnostics for one state, separating compliance and actuator contributions and reporting both world and child-body frames.

## Design & Implementation
Immutable with the joint name and indices, translational forces on parent and child in world frame, compliance and actuator torques on each side in world frame, and the compliance, actuator and total torque on the child in its body frame. Parent-side quantities are the negatives of child-side ones by construction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `parent` | Int | n/a | yes | Field `parent`. |
| in | `child` | Int | n/a | yes | Field `child`. |
| in | `translation_force_parent_world` | SVector{3, Float64} | n/a | yes | Field `translation_force_parent_world`. |
| in | `translation_force_child_world` | SVector{3, Float64} | n/a | yes | Field `translation_force_child_world`. |
| in | `compliance_torque_parent_world` | SVector{3, Float64} | n/a | yes | Field `compliance_torque_parent_world`. |
| in | `compliance_torque_child_world` | SVector{3, Float64} | n/a | yes | Field `compliance_torque_child_world`. |
| in | `actuator_torque_parent_world` | SVector{3, Float64} | n/a | yes | Field `actuator_torque_parent_world`. |
| in | `actuator_torque_child_world` | SVector{3, Float64} | n/a | yes | Field `actuator_torque_child_world`. |
| in | `compliance_torque_child_body` | SVector{3, Float64} | n/a | yes | Field `compliance_torque_child_body`. |
| in | `actuator_torque_child_body` | SVector{3, Float64} | n/a | yes | Field `actuator_torque_child_body`. |
| in | `total_torque_child_body` | SVector{3, Float64} | n/a | yes | Field `total_torque_child_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantJointLoad | n/a | — | Constructed `CompliantJointLoad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:562-562`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Twelve fields of redundant information per joint make this relatively heavy to construct for a large grid on every derivative evaluation, which `compliant_multibody_dynamics` does.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 99.
