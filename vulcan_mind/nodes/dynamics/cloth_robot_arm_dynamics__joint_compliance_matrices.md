---
id: dynamics.cloth_robot_arm_dynamics__joint_compliance_matrices
label: _joint_compliance_matrices
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _joint_compliance_matrices
  lines:
  - 119
  - 119
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `name`.
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
  type: Vector{SMatrix{3,
  units: n/a
  description: Return value of `_joint_compliance_matrices`.
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

# _joint_compliance_matrices

## Purpose
Expands a stiffness or damping keyword into a `Vector{SMatrix{3,3,Float64}}` with one matrix per arm joint, supporting a shared value or a per-joint list.

## Design & Implementation
If `value` is a `Real`, `AbstractMatrix`, or `Tuple`, returns `fill(_compliance_matrix(value, name), n)`. If it is an `AbstractVector` of length `n`, maps `_compliance_matrix` over its elements. Otherwise throws `ArgumentError`. Called four times per model build and four times per RHS evaluation for Kx, Cx, Kr, Cr.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{SMatrix{3, | n/a | — | Return value of `_joint_compliance_matrices`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:238-238`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:353-353`

**Downstream**

- `callees` → [[dynamics.cloth_robot_arm_dynamics__compliance_matrix|_compliance_matrix]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
A length-3 vector is interpreted as per-joint only when `n == 3`, otherwise it falls through to the else branch and throws, so the diagonal-vector convenience is unavailable for 3-joint arms. The vectors are allocated fresh on every call, including inside `assign_coupled_cloth_robot_arm_rhs!`, producing per-step garbage.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 119.
