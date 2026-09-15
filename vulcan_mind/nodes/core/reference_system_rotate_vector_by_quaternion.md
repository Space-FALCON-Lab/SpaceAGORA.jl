---
id: core.reference_system_rotate_vector_by_quaternion
label: rotate_vector_by_quaternion
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rotate_vector_by_quaternion
  lines:
  - 560
  - 560
inputs:
- id: v
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `v`.
- id: q
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `q`.
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
  description: Return value of `rotate_vector_by_quaternion`. Returns `v_rotated`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# rotate_vector_by_quaternion

## Purpose
Rotates a three-vector by a scalar-last quaternion using the efficient two-cross-product form.

## Theory & Math
$$
t = 2\, \vec{q} \times \vec{v},\qquad \vec{v}' = \vec{v} + q_w\, t + \vec{q} \times t
$$

## Design & Implementation
Splits `q` into its vector part and scalar, computes `t = 2 (q_v × v)`, and returns `v + q_w t + q_v × t`. This form uses two cross products and no matrix construction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Vector{Float64} | n/a | yes | Positional argument `v`. |
| in | `q` | Vector{Float64} | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rotate_vector_by_quaternion`. Returns `v_rotated`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_orbital_elements_to_lvlh_quaternion|orbital_elements_to_lvlh_quaternion]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:555-555`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system_rtn_dcm_from_inertial|rtn_dcm_from_inertial]] · `callers` · call · `src/core/interfaces/reference_system.jl:617-617`
<!-- vulcan:connections:end -->

## Limitations
Typed on `Vector{Float64}` rather than `AbstractVector` or `SVector`, so static-vector callers must convert and pay allocations; the quaternion is not normalised, so a non-unit input scales the result.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 560.
