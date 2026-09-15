---
id: dynamics.cloth_multibody_compliant_state_vector
label: compliant_state_vector
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: compliant_state_vector
  lines:
  - 436
  - 436
inputs:
- id: positions
  type: Any
  units: n/a
  required: true
  description: Positional argument `positions`.
- id: quaternions
  type: Any
  units: n/a
  required: true
  description: Positional argument `quaternions`.
- id: velocities
  type: Any
  units: n/a
  required: false
  description: Keyword argument `velocities` (default `fill(SVector{3, Float64}(0.0,
    0.0, 0.0), length(positions))`).
- id: angular_velocities
  type: Any
  units: n/a
  required: false
  description: Keyword argument `angular_velocities` (default `fill(SVector{3, Float64}(0.0,
    0.0, 0.0), length(positions))`).
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
  description: Return value of `compliant_state_vector`. Returns `x`.
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

# compliant_state_vector

## Purpose
Packs per-body positions, quaternions, velocities and angular rates into the flat thirteen-per-body state vector the dynamics operate on.

## Design & Implementation
Validates that all four collections have the same length, allocates `13n` zeros, and writes each body's block in the order position, normalised quaternion, velocity, angular rate. Velocities and angular rates default to zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `positions` | Any | n/a | yes | Positional argument `positions`. |
| in | `quaternions` | Any | n/a | yes | Positional argument `quaternions`. |
| in | `velocities` | Any | n/a | no | Keyword argument `velocities` (default `fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(positions))`). |
| in | `angular_velocities` | Any | n/a | no | Keyword argument `angular_velocities` (default `fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(positions))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compliant_state_vector`. Returns `x`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:316-316`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state|cloth_robot_arm_initial_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:182-182`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:450-450`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:450-450`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:450-450`
<!-- vulcan:connections:end -->

## Limitations
The thirteen-slot layout is implicit shared knowledge between this function, `compliant_state_parts` and the dynamics; there is no named accessor for the layout constants.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 436.
