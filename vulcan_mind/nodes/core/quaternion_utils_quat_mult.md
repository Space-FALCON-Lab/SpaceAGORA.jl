---
id: core.quaternion_utils_quat_mult
label: quat_mult
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: quat_mult
  lines:
  - 9
  - 9
inputs:
- id: q
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `q`.
- id: p
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `quat_mult`. Returns `SVector{4, Float64}(v1, v2, v3,
    s)`.
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

# quat_mult

## Purpose

Composes two attitude quaternions, returning the product `q * p` in the scalar-last convention. Attitude propagation, reference-frame chaining and error-quaternion construction all rely on this as the basic composition operation for successive rotations.

## Design & Implementation

`quat_mult(q, p)` unpacks both `AbstractVector` arguments into scalar locals (`q0 = q[4]` for the scalar part, `q1..q3` for the vector part, and likewise for `p`), evaluates the Hamilton product in four explicit expressions, and returns a freshly built `SVector{4,Float64}` ordered `[v1, v2, v3, s]` so the scalar part stays last. Returning a static vector keeps the operation heap-free, and `@inline` lets it fuse into propagation loops.

## Theory & Math

Writing $q = q_s + \mathbf{q}_v$ and $p = p_s + \mathbf{p}_v$ with scalar parts $q_s = q[4]$, $p_s = p[4]$ and vector parts $\mathbf{q}_v$, $\mathbf{p}_v$, the Hamilton product implemented here is

$$q \otimes p = \left( q_s p_s - \mathbf{q}_v \cdot \mathbf{p}_v,\; q_s \mathbf{p}_v + p_s \mathbf{q}_v + \mathbf{q}_v \times \mathbf{p}_v \right)$$

with the scalar component stored in slot 4 and the vector component in slots 1 to 3 of the result.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | AbstractVector | n/a | yes | Positional argument `q`. |
| in | `p` | AbstractVector | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `quat_mult`. Returns `SVector{4, Float64}(v1, v2, v3, s)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:159-159`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/numerics/quaternion_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The result is not renormalised, so repeated composition accumulates norm drift that callers must correct. The output element type is fixed to `Float64` regardless of the inputs, which converts away higher-precision or differentiable number types. Indices 1 through 4 are read without a length check, so a short vector raises a raw bounds error, and the scalar-last convention is assumed rather than validated.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 9.
