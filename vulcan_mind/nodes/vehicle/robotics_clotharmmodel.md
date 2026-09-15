---
id: vehicle.robotics_clotharmmodel
label: ClothArmModel
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmModel
  lines:
  - 46
  - 46
inputs:
- id: links
  type: Vector{ClothArmLink}
  units: n/a
  required: true
  description: Field `links`.
- id: joints
  type: Vector{ClothArmJoint}
  units: n/a
  required: true
  description: Field `joints`.
- id: mount_offset_body
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `mount_offset_body`.
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
  type: ClothArmModel
  units: n/a
  description: Constructed `ClothArmModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# ClothArmModel

## Purpose
The serial chain: ordered links and joints plus where the chain attaches to the spacecraft body.

## Design & Implementation
Immutable, holding a `Vector{ClothArmLink}`, a `Vector{ClothArmJoint}` of the same length, and `mount_offset_body`, the vector from the base pose position to the first joint origin expressed in the base frame. Joint `i` precedes link `i`, so the chain alternates joint, link, joint, link from the mount outward.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `links` | Vector{ClothArmLink} | n/a | yes | Field `links`. |
| in | `joints` | Vector{ClothArmJoint} | n/a | yes | Field `joints`. |
| in | `mount_offset_body` | SVector{3, Float64} | n/a | yes | Field `mount_offset_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmModel | n/a | — | Constructed `ClothArmModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:145-145`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The equal-length invariant between links and joints is assumed by every consumer but not checked by the struct itself, only by `default_cloth_arm_model`; a hand-built model with mismatched vectors fails inside `cloth_fk` with a bounds error.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 46.
