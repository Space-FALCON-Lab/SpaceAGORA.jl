---
id: core.quaternion_utils_hat
label: hat
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: hat
  lines:
  - 29
  - 29
inputs:
- id: omega
  type: AbstractVector{<:Float64}
  units: n/a
  required: true
  description: Positional argument `ω`.
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
  description: Return value of `hat`.
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

# hat

## Purpose

Builds the 3x3 skew-symmetric matrix that is the matrix equivalent of the cross product, so that `hat(w) * v` equals `w x v`. Attitude kinematics and rigid-body dynamics code in the ADCS stack uses it to express the angular-velocity cross term and the Euler equation gyroscopic torque as plain matrix products.

## Design & Implementation

`hat(w)` takes any `AbstractVector{<:Float64}`, reads the first three entries into locals `w1`, `w2`, `w3`, and returns an `SMatrix{3,3,Float64}` built directly from its column-major argument list (column 1 `[0, w3, -w2]`, column 2 `[-w3, 0, w1]`, column 3 `[w2, -w1, 0]`). Constructing the static matrix in one call avoids the temporary array a literal matrix would allocate, and the function is marked `@inline` so it disappears into the caller's stack-allocated arithmetic.

## Theory & Math

For a vector $\omega = [\omega_1, \omega_2, \omega_3]^\top \in \mathbb{R}^3$ the hat map produces

$$[\omega]_\times = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix}$$

where $\omega_i$ are the components of $\omega$ in the frame of interest (rad/s when $\omega$ is an angular rate). It satisfies $[\omega]_\times v = \omega \times v$ for every $v \in \mathbb{R}^3$, and $[\omega]_\times^\top = -[\omega]_\times$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `omega` | AbstractVector{<:Float64} | n/a | yes | Positional argument `ω`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `hat`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- [[vehicle.mass_properties_update_inertia_tensor|update_inertia_tensor]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:61-61`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The signature is restricted to `Float64` element types, so dual numbers or `Float32` states will not dispatch here and forward-mode automatic differentiation through this call fails. Only the first three entries are read; a longer vector is silently truncated and a shorter one raises a bounds error rather than a descriptive message. No finiteness check is performed, so `NaN` rates propagate into the returned matrix.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 29.
