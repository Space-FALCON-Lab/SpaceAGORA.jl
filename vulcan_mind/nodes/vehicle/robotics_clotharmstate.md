---
id: vehicle.robotics_clotharmstate
label: ClothArmState
kind: struct
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: ClothArmState
  lines:
  - 66
  - 66
inputs:
- id: pose
  type: ClothArmPose
  units: n/a
  required: true
  description: Field `pose`.
- id: q
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `q`.
- id: dq
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `dq`.
- id: link_linear_velocities
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_linear_velocities`.
- id: link_angular_velocities
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_angular_velocities`.
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
  type: ClothArmState
  units: n/a
  description: Constructed `ClothArmState`.
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

# ClothArmState

## Purpose
A pose bundled with the joint coordinates and rates that produced it, plus per-link velocity slots for the dynamics layer to fill.

## Design & Implementation
Immutable, holding the `ClothArmPose`, joint vector `q`, joint rate vector `dq`, and `link_linear_velocities` and `link_angular_velocities` as vectors of static three-vectors. `cloth_fk_state` constructs it with both velocity vectors zero-filled; the multibody dynamics populate them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pose` | ClothArmPose | n/a | yes | Field `pose`. |
| in | `q` | Vector{Float64} | n/a | yes | Field `q`. |
| in | `dq` | Vector{Float64} | n/a | yes | Field `dq`. |
| in | `link_linear_velocities` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_linear_velocities`. |
| in | `link_angular_velocities` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_angular_velocities`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothArmState | n/a | — | Constructed `ClothArmState`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_cloth_fk_state|cloth_fk_state]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:221-221`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the struct is immutable and the velocity vectors are zero at construction, a consumer reading them from a freshly built state gets zeros that look like a valid resting arm rather than a not-yet-computed marker.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 66.
