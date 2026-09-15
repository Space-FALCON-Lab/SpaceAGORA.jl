---
id: simulation.dynamics_rhs__assign_orientation_rhs_bang
label: _assign_orientation_rhs!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _assign_orientation_rhs!
  lines:
  - 1584
  - 1584
inputs:
- id: du_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `du_view`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: inertia_tensor
  type: AbstractMatrix{<:Real}
  units: n/a
  required: true
  description: Positional argument `inertia_tensor`.
- id: torques
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `torques`.
- id: propagate_quaternion
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `propagate_quaternion`.
- id: include_gyroscopic
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `include_gyroscopic`.
- id: rw_assembly
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rw_assembly` (default `nothing`).
- id: rw_torque_body
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Keyword argument `rw_torque_body` (default `SVector{3, Float64}(0.0,
    0.0, 0.0)`).
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
  type: Nothing
  units: n/a
  description: Return value of `_assign_orientation_rhs!`; mutates `du_view` in place.
    Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _assign_orientation_rhs!

## Purpose
Writes the attitude derivatives â€” quaternion rate, wheel momentum rate and angular acceleration â€” for one satellite.

## Design & Implementation
Computes the quaternion derivative from body rate if `propagate_quaternion`, otherwise zero; if a reaction wheel assembly exists, forms the wheel momentum in body frame and sets the wheel momentum rate from the pseudo-inverse of the commanded torque; then computes angular acceleration with optional gyroscopic term. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `inertia_tensor` | AbstractMatrix{<:Real} | n/a | yes | Positional argument `inertia_tensor`. |
| in | `torques` | AbstractVector{<:Real} | n/a | yes | Positional argument `torques`. |
| in | `propagate_quaternion` | Bool | n/a | yes | Keyword argument `propagate_quaternion`. |
| in | `include_gyroscopic` | Bool | n/a | yes | Keyword argument `include_gyroscopic`. |
| in | `rw_assembly` | Any | n/a | no | Keyword argument `rw_assembly` (default `nothing`). |
| in | `rw_torque_body` | SVector{3, Float64} | n/a | no | Keyword argument `rw_torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_assign_orientation_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1359-1359`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2135-2135`
- [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2232-2232`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2025-2025`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1877-1877`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1766-1766`

**Downstream**

- `callees` → [[dynamics.torque_models_body_angular_velocity|body_angular_velocity]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1594-1594`
- `callees` → [[dynamics.torque_models_body_torque|body_torque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1595-1595`
- `callees` → [[dynx.rotational_attitude_kinematics_quaternion_derivative|quaternion_derivative]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1598-1598`
- `callees` → [[dynx.rotational_rigid_body_dynamics_angular_acceleration|angular_acceleration]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1609-1609`
<!-- vulcan:connections:end -->

## Limitations
The wheel momentum vector is built with a runtime-sized `SVector`, which specialises per wheel count.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1584.
