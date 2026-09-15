---
id: vehicle.kinematics_kinematics
label: Kinematics
kind: module
source:
  file: src/vehicle/kinematics/kinematics.jl
  symbol: Kinematics
  lines:
  - 1
  - 1
inputs:
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# Kinematics

## Purpose

`Kinematics` is the module that owns frame-conversion helpers for the articulated `SpacecraftModel` link tree. It exports `rotate_to_inertial`, `rotate_to_body` and `rotate_link`, giving callers the direction cosine matrices needed to move vectors between a link frame, the root body frame and the inertial frame, and a mutating setter for a link's attitude quaternion.

## Design & Implementation

The module pulls in `StaticArrays` and `LinearAlgebra`, imports `..SpacecraftModels` for the `SpacecraftModel` and `Link` types, and `include`s `core/numerics/quaternion_utils.jl` for `rot`, `project_unit_quaternion` and `dcm_to_quaternion`. The convention it enforces is that a root `Link` stores its attitude relative to inertial while every child stores attitude relative to the root body frame, so `rotate_to_inertial` composes `rot(model.root.q)' * rot(body.q)'` for children and returns `rot(body.q)'` alone for a root.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/kinematics/kinematics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The include of `quaternion_utils.jl` is path-relative and brings those helpers into `Kinematics` rather than sharing a single loaded copy, so the definitions are duplicated wherever else that file is included. Nothing in the module validates that `model.root` is the actual parent of a given child link; the two-level root/child convention is assumed, and deeper chains are not composed.

## Provenance
Mapped from `src/vehicle/kinematics/kinematics.jl` line 1.
