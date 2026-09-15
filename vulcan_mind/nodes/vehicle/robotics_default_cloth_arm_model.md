---
id: vehicle.robotics_default_cloth_arm_model
label: default_cloth_arm_model
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: default_cloth_arm_model
  lines:
  - 116
  - 116
inputs:
- id: link_lengths_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `link_lengths_m` (default `(0.18, 0.16, 0.12)`).
- id: link_radii_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `link_radii_m` (default `(0.018, 0.016, 0.014)`).
- id: link_masses_kg
  type: Any
  units: n/a
  required: false
  description: Keyword argument `link_masses_kg` (default `(0.15, 0.10, 0.05)`).
- id: joint_axes
  type: Any
  units: n/a
  required: false
  description: Keyword argument `joint_axes` (default `((0.0, 0.0, 1.0), (0.0, 1.0,
    0.0), (0.0, 1.0, 0.0))`).
- id: joint_limit_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `joint_limit_rad` (default `deg2rad(175.0)`).
- id: mount_offset_body
  type: Any
  units: n/a
  required: false
  description: Keyword argument `mount_offset_body` (default `(0.0, 0.0, 0.0)`).
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
  description: Return value of `default_cloth_arm_model`. Returns `ClothArmModel(links,
    joints, SVector{3, Float64}(mount_offset_body))`.
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

# default_cloth_arm_model

## Purpose
Constructs the three-link reference arm used by examples and tests, with every dimension overridable by keyword.

## Design & Implementation
Defaults to link lengths of 0.18, 0.16 and 0.12 m, radii of 18, 16 and 14 mm, masses of 150, 100 and 50 g, joint axes z then y then y, symmetric limits of 175 degrees and a zero mount offset. It validates that the radius, mass and axis tuples match the length tuple, raising `ArgumentError` naming the offending keyword, then builds links along local x with the centre of mass at half length and joints with normalised axes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link_lengths_m` | Any | n/a | no | Keyword argument `link_lengths_m` (default `(0.18, 0.16, 0.12)`). |
| in | `link_radii_m` | Any | n/a | no | Keyword argument `link_radii_m` (default `(0.018, 0.016, 0.014)`). |
| in | `link_masses_kg` | Any | n/a | no | Keyword argument `link_masses_kg` (default `(0.15, 0.10, 0.05)`). |
| in | `joint_axes` | Any | n/a | no | Keyword argument `joint_axes` (default `((0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0))`). |
| in | `joint_limit_rad` | Any | n/a | no | Keyword argument `joint_limit_rad` (default `deg2rad(175.0)`). |
| in | `mount_offset_body` | Any | n/a | no | Keyword argument `mount_offset_body` (default `(0.0, 0.0, 0.0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmModel | n/a | — | Return value of `default_cloth_arm_model`. Returns `ClothArmModel(links, joints, SVector{3, Float64}(mount_offset_body))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/vehicle/robotics/robotics.jl:131-131`
- `callees` → [[vehicle.robotics__normalize_axis|_normalize_axis]] · `callers` · call · `src/vehicle/robotics/robotics.jl:140-140`
- `callees` → [[vehicle.robotics_clotharmjoint|ClothArmJoint]] · `callers` · call · `src/vehicle/robotics/robotics.jl:138-138`
- `callees` → [[vehicle.robotics_clotharmlink|ClothArmLink]] · `callers` · call · `src/vehicle/robotics/robotics.jl:129-129`
- `callees` → [[vehicle.robotics_clotharmmodel|ClothArmModel]] · `callers` · call · `src/vehicle/robotics/robotics.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
All joints share one limit magnitude, so a model needing per-joint limits must construct `ClothArmJoint` values directly; the centre-of-mass placement at half length is hard-coded.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 116.
