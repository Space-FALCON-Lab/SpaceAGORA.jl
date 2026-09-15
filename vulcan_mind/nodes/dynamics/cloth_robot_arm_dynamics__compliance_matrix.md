---
id: dynamics.cloth_robot_arm_dynamics__compliance_matrix
label: _compliance_matrix
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _compliance_matrix
  lines:
  - 101
  - 101
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_compliance_matrix`.
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

# _compliance_matrix

## Purpose
Normalises one joint's stiffness or damping specification (scalar, 3-vector, 3-tuple, or 3 x 3 matrix) into an `SMatrix{3,3,Float64}`, throwing a descriptive `ArgumentError` on bad shapes.

## Design & Implementation
Dispatches on runtime type: `Real` goes to `_diag3`; `AbstractMatrix` must be `size == (3,3)`; `AbstractVector` and `Tuple` must have length 3 and become `Diagonal(Float64.(value))`. Any other type throws `ArgumentError("$(name) must be a scalar, a 3x3 matrix, or a per-joint vector of those values.")`. The `name::Symbol` argument is only used in error text.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_compliance_matrix`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics__joint_compliance_matrices|_joint_compliance_matrices]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:121-121`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_robot_arm_dynamics__diag3|_diag3]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:103-103`
<!-- vulcan:connections:end -->

## Limitations
Non-symmetric or negative-definite matrices are accepted, allowing physically inconsistent joints. A `Vector` of length 3 is ambiguous with a per-joint list for a 3-link arm; `_joint_compliance_matrices` resolves that by checking `length == n` first, so a 3-link arm cannot pass a diagonal vector shared across joints. Element types must convert to `Float64` or a `MethodError` escapes.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 101.
