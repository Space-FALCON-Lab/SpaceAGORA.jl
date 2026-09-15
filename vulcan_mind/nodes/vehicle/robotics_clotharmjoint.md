---
id: vehicle.robotics_clotharmjoint
label: ClothArmJoint
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmJoint
  lines:
  - 38
  - 38
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: axis_parent
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `axis_parent`.
- id: lower_rad
  type: Float64
  units: n/a
  required: true
  description: Field `lower_rad`.
- id: upper_rad
  type: Float64
  units: n/a
  required: true
  description: Field `upper_rad`.
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
  type: ClothArmJoint
  units: n/a
  description: Constructed `ClothArmJoint`.
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

# ClothArmJoint

## Purpose
One revolute joint: its rotation axis in the parent frame and its angular travel limits.

## Design & Implementation
Immutable with a `name`, a unit `axis_parent`, and `lower_rad` and `upper_rad`. The default builder normalises the axis through `_normalize_axis` and sets symmetric limits of ±175 degrees. Forward kinematics rotates about `axis_parent` expressed in the accumulated parent orientation, so a joint's axis is fixed relative to the link before it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `axis_parent` | SVector{3, Float64} | n/a | yes | Field `axis_parent`. |
| in | `lower_rad` | Float64 | n/a | yes | Field `lower_rad`. |
| in | `upper_rad` | Float64 | n/a | yes | Field `upper_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmJoint | n/a | — | Constructed `ClothArmJoint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:138-138`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only revolute joints are representable; there is no prismatic variant. The limits are stored but enforced only by `cloth_ik`'s clamp, not by `cloth_fk`, which happily evaluates out-of-range angles.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 38.
