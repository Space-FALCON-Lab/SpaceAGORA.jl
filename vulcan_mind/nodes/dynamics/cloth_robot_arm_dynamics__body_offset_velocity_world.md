---
id: dynamics.cloth_robot_arm_dynamics__body_offset_velocity_world
label: _body_offset_velocity_world
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _body_offset_velocity_world
  lines:
  - 130
  - 130
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
  description: Return value of `_body_offset_velocity_world`. Returns `SVector{3,
    Float64}(v_world) + _rot(q) * cross(SVector{3, Float64}(ω_body), SVec`.
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

# _body_offset_velocity_world

## Purpose
Computes the world-frame velocity of a point fixed in a body at body-frame offset `p_body`, given the body's world velocity, attitude, and body angular velocity.

## Theory & Math
$\mathbf{v}_P = \mathbf{v}_{cm} + R(q)\,(\boldsymbol\omega_b \times \mathbf{p}_b)$, with $\mathbf{v}_{cm}$ the body centre-of-mass velocity (world), $\boldsymbol\omega_b$ the body-frame angular velocity, and $\mathbf{p}_b$ the body-frame offset.

## Design & Implementation
Returns `v_world + R(q) * (ω_body × p_body)` using `_rot(q)` and `cross` on `SVector{3,Float64}` conversions of the inputs. Used by `_coupled_parent_kinematics` to obtain joint attachment-point velocities for the damping term.

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
| out | `result` | SVector | n/a | — | Return value of `_body_offset_velocity_world`. Returns `SVector{3, Float64}(v_world) + _rot(q) * cross(SVector{3, Float64}(ω_body), SVec`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics__coupled_parent_kinematics|_coupled_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:320-320`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:131-131`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:131-131`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:131-131`
<!-- vulcan:connections:end -->

## Limitations
Calls `_rot` (with its normalisation) on every evaluation. Assumes `ω_body` is expressed in the body frame; passing a world-frame rate gives wrong velocities without any diagnostic.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 130.
