---
id: core.quaternion_utils_error_quaternion
label: error_quaternion
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: error_quaternion
  lines:
  - 61
  - 61
inputs:
- id: current
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Positional argument `current`.
- id: target
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Positional argument `target`.
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
  description: Return value of `error_quaternion`. Returns `result / norm(result)`.
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

# error_quaternion

## Purpose

Produces the attitude-error quaternion between a current attitude and a commanded target attitude, which attitude controllers use as the rotational error signal driving corrective torque. The result is the rotation that still has to be applied to reach the target.

## Design & Implementation

`error_quaternion(current, target)` takes two `SVector{4,Float64}` scalar-last quaternions. It unpacks `current` into `a, b, c, s`, assembles the 4x4 composition matrix `q_matrix` directly in column-major form (equivalent to `[s*I - hat([a,b,c]) [a,b,c]; -[a,b,c]' s]`), forms the conjugate of the target as `qd = [-target[1], -target[2], -target[3], target[4]]`, multiplies, and returns `result / norm(result)`. Building the matrix inline avoids the slice allocations a naive block assembly would cause, and the explicit normalisation guarantees a unit output.

## Theory & Math

The error quaternion is the composition of the current attitude with the inverse of the target,

$$\delta q = q_{\text{cur}} \otimes q_{\text{tgt}}^{*}, \qquad q^{*} = (-q_1, -q_2, -q_3, q_4)$$

where $q_{\text{tgt}}^{*}$ is the conjugate, which equals the inverse for a unit quaternion. The final division by $\|\delta q\|$ projects the product back onto the unit sphere $S^3$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `current` | SVector{4, Float64} | n/a | yes | Positional argument `current`. |
| in | `target` | SVector{4, Float64} | n/a | yes | Positional argument `target`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `error_quaternion`. Returns `result / norm(result)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/numerics/quaternion_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

If `norm(result)` is zero or non-finite the division yields `NaN` or `Inf` with no guard, which can happen when either input is the zero vector. The sign ambiguity of the double cover is not resolved: no check forces the scalar part positive, so the controller may see an error representing the long way around unless it handles the sign itself. Inputs must already be unit quaternions in the scalar-last order, and the concrete `SVector{4,Float64}` signature excludes other vector or number types.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 61.
