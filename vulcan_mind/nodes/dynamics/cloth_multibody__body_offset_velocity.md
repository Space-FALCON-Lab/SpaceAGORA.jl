---
id: dynamics.cloth_multibody__body_offset_velocity
label: _body_offset_velocity
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _body_offset_velocity
  lines:
  - 178
  - 178
inputs:
- id: v_world
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_world`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
- id: omega_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `ω_body`.
- id: p_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `p_body`.
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
  type: SVector
  units: n/a
  description: Return value of `_body_offset_velocity`. Returns `SVector{3, Float64}(v_world)
    + _rot(q) * cross(SVector{3, Float64}(ω_body), p_bo`.
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

# _body_offset_velocity

## Purpose
Computes the world-frame velocity of a point fixed in a body, given the body's velocity, orientation and body-frame angular rate.

## Theory & Math
$$
\vec{v}_p = \vec{v} + R(q)\,(\vec{\omega}_b \times \vec{p}_b)
$$

## Design & Implementation
Returns `v_world + R(q) (ω_body × p_body)`. `@inline`. Used for joint attachment points so damping acts on the relative velocity of the actual attachment, not the body centres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v_world` | Any | n/a | yes | Positional argument `v_world`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `omega_body` | Any | n/a | yes | Positional argument `ω_body`. |
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_body_offset_velocity`. Returns `SVector{3, Float64}(v_world) + _rot(q) * cross(SVector{3, Float64}(ω_body), p_bo`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__parent_kinematics|_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:505-505`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:179-179`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__diag3|_diag3]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:179-179`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:179-179`
<!-- vulcan:connections:end -->

## Limitations
Assumes `p_body` is a body-frame offset from the centre of mass, consistent with how attachment points are stored.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 178.
