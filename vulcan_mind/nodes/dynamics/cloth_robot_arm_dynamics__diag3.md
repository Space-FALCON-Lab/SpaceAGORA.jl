---
id: dynamics.cloth_robot_arm_dynamics__diag3
label: _diag3
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _diag3
  lines:
  - 98
  - 98
inputs:
- id: v
  type: Real
  units: n/a
  required: true
  description: Positional argument `v`.
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
  type: SMatrix
  units: n/a
  description: Return value of `_diag3`. Returns `SMatrix{3, 3, Float64}(Float64(v)
    * I)`.
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

# _diag3

## Purpose
Builds a 3 x 3 static diagonal matrix `v * I` from a scalar, used to expand scalar stiffness or damping into an isotropic compliance matrix.

## Design & Implementation
One-line `@inline` definition `SMatrix{3,3,Float64}(Float64(v) * I)` where `I` is `LinearAlgebra.UniformScaling`. Called by `_compliance_matrix` for `Real` inputs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Real | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix | n/a | — | Return value of `_diag3`. Returns `SMatrix{3, 3, Float64}(Float64(v) * I)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__body_offset_velocity|_body_offset_velocity]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- [[dynamics.cloth_robot_arm_dynamics__compliance_matrix|_compliance_matrix]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:103-103`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:98-98`
<!-- vulcan:connections:end -->

## Limitations
Accepts negative or zero values without complaint, which yields a non-restoring or absent spring. Non-finite inputs propagate NaN or Inf into the joint model silently.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 98.
