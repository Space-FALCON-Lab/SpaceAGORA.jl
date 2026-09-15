---
id: vehicle.robotics_clotharmlink
label: ClothArmLink
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmLink
  lines:
  - 29
  - 29
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: vector_parent
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `vector_parent`.
- id: com_offset_parent
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `com_offset_parent`.
- id: radius_m
  type: Float64
  units: n/a
  required: true
  description: Field `radius_m`.
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `mass_kg`.
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
  type: ClothArmLink
  units: n/a
  description: Constructed `ClothArmLink`.
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

# ClothArmLink

## Purpose
Geometry and mass of one rigid link, expressed in the parent joint's frame.

## Design & Implementation
Immutable with a `name`, the tip vector `vector_parent` from the joint origin to the next joint, `com_offset_parent` from the joint origin to the centre of mass, a cylindrical `radius_m` used for contact and rendering, and `mass_kg`. The default model builder places the tip along local x and the centre of mass at half length.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `vector_parent` | SVector{3, Float64} | n/a | yes | Field `vector_parent`. |
| in | `com_offset_parent` | SVector{3, Float64} | n/a | yes | Field `com_offset_parent`. |
| in | `radius_m` | Float64 | n/a | yes | Field `radius_m`. |
| in | `mass_kg` | Float64 | n/a | yes | Field `mass_kg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmLink | n/a | — | Constructed `ClothArmLink`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:129-129`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No inertia tensor is stored, so the dynamics layer must reconstruct one from mass, length and radius under a uniform-cylinder assumption; the link cannot describe an asymmetric body.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 29.
