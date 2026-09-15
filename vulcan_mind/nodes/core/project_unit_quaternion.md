---
id: core.project_unit_quaternion
label: project_unit_quaternion
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: project_unit_quaternion
  lines:
  - 19
  - 27
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SimulationModel namespace providing shared quaternion numerics to dynamics
    and control.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: quaternion
  type: SVector{4,Float64}
  units: n/a
  description: Scalar-last unit quaternion, or the identity quaternion when the input
    norm is invalid.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
- numerics
charts:
- core
origin: agent
---

# project_unit_quaternion

## Purpose
`project_unit_quaternion` converts a real quaternion-like vector into a finite scalar-last unit quaternion for attitude propagation and control. It centralizes normalization and supplies a deterministic identity fallback when the input norm cannot define a valid orientation.

## Theory & Math
For `q ∈ R⁴`, the projection is `q̂ = q / sqrt(qᵀq)` when `qᵀq` is finite and greater than machine epsilon. The function returns `q_id = (0,0,0,1)` otherwise. The scalar-last convention matches `quat_mult`, rotational kinematics, and the attitude state layout.

## Model & Assumptions
The input has at least four indexable values representing one quaternion in scalar-last order. Normalization changes magnitude but not orientation for finite nonzero inputs. Replacing invalid input with identity is an explicit recovery policy and assumes identity is safer than propagating NaNs.

## Design & Implementation
The function constructs an `SVector{4,Float64}`, computes `dot(q_unit,q_unit)`, checks finiteness and the epsilon threshold, then divides by the square root. Returning a static vector avoids heap allocation in the RHS and control paths. The identity constant is defined beside the function in `quaternion_utils.jl`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SimulationModel namespace providing shared quaternion numerics to dynamics and control. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `quaternion` | SVector{4,Float64} | n/a | — | Scalar-last unit quaternion, or the identity quaternion when the input norm is invalid. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control__robot_arm_control_axis_angle_about|_robot_arm_control_axis_angle_about]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:140-140`
- [[gnc.robot_arm_control__robot_arm_control_quat_conj|_robot_arm_control_quat_conj]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:134-134`
- [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:153-153`
- [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2338-2338`
- [[vehicle.kinematics_rotate_link|rotate_link]] · `callees` → `callers` · call · `src/vehicle/kinematics/kinematics.jl:38-38`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/core/numerics/quaternion_utils.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
The identity fallback hides the original invalid quaternion from the caller, so diagnostics must be added at the state-validation boundary when invalid attitude data is a fault. The function does not enforce continuity across the quaternion double cover or choose a sign relative to a previous attitude. Inputs with fewer than four elements fail by indexing.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl:19-27`.
