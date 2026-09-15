---
id: dynamics.cloth_robot_arm_dynamics__coupled_parent_kinematics
label: _coupled_parent_kinematics
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _coupled_parent_kinematics
  lines:
  - 306
  - 306
inputs:
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: parent
  type: Int
  units: n/a
  required: true
  description: Positional argument `parent`.
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
  type: Any
  units: n/a
  description: Return value of `_coupled_parent_kinematics`. Returns `(`.
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

# _coupled_parent_kinematics

## Purpose
Returns the pose, rotation matrix, angular rates, and world position/velocity of a joint attachment point `p_body` on either the spacecraft base (`parent == 0`) or arm link `parent`.

## Design & Implementation
For the base: `q` from `sc_view.q` (or identity), `r = sc_view.pos`, `v = sc_view.vel`, `ω = sc_view.ω` (or zero). For a link: `_coupled_body_state(sc_view, parent)`. Both branches compute `R = _rot(q)`, `ω_world = R * ω`, `point = r + R * p_body`, and `point_velocity = _body_offset_velocity_world(v, q, ω, p_body)`. The NamedTuple fields are `r, q, R, ω, ω_world, point, point_velocity`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `parent` | Int | n/a | yes | Positional argument `parent`. |
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_coupled_parent_kinematics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:374-374`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:312-312`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:308-308`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__body_offset_velocity_world|_body_offset_velocity_world]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:320-320`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__coupled_body_state|_coupled_body_state]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:323-323`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:312-312`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:308-308`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:312-312`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:308-308`
<!-- vulcan:connections:end -->

## Limitations
`sc_view.vel` is accessed without a `hasproperty` guard on the base branch, unlike `q` and `ω`, so a view lacking `vel` throws. `_rot` is evaluated twice per attachment (once here and again inside `_body_offset_velocity_world`). The two branches return NamedTuples of identical field names but the same types only if `sc_view` fields are `Float64`.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 306.
