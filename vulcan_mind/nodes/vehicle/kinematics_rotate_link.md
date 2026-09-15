---
id: vehicle.kinematics_rotate_link
label: rotate_link
kind: function
source:
  file: src/vehicle/kinematics/kinematics.jl
  symbol: rotate_link
  lines:
  - 32
  - 32
inputs:
- id: body
  type: Link
  units: n/a
  required: true
  description: Positional argument `body`.
- id: q
  type: SVector{4, Float64}
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
  description: Return value of `rotate_link`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# rotate_link

## Purpose

`rotate_link` overwrites the attitude quaternion of a non-root `Link` in place. Three methods accept the new orientation in different forms: a unit quaternion `q::SVector{4,Float64}`, a direction cosine matrix `dcm::SMatrix{3,3,Float64}`, or an axis/angle pair `axis::SVector{3,Float64}` with angle `θ` in radians.

## Design & Implementation

Every method opens with `@assert !body.root "Cannot rotate a root body directly"`, then broadcasts into `body.q` with `.=` so the link object is mutated rather than replaced. The quaternion method passes through `project_unit_quaternion(q)` to renormalise; the DCM method converts with `dcm_to_quaternion(dcm)`. The axis/angle method substitutes the default axis `(0, 1, 0)` when `norm(axis) <= 1e-6`, normalises the axis, and builds the scalar-last quaternion `[axis * sin(θ/2); cos(θ/2)]`. All three return the mutated `body.q`.

## Theory & Math

For a unit axis $\hat{a}$ and rotation angle $\theta$ the scalar-last quaternion assembled by the axis/angle method is

$$q = \begin{bmatrix} \hat{a}\sin(\theta/2) \\ \cos(\theta/2) \end{bmatrix}$$

where $\hat{a} = \text{axis}/\lVert \text{axis} \rVert$ and $\theta$ is in radians. The construction is unit-norm by identity, since $\lVert\hat{a}\rVert^2\sin^2(\theta/2) + \cos^2(\theta/2) = 1$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `q` | SVector{4, Float64} | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rotate_link`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__apply_solar_panel_aoa_bang|_apply_solar_panel_aoa!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:240-240`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/kinematics/kinematics.jl`

**Downstream**

- `callees` → [[core.project_unit_quaternion|project_unit_quaternion]] · `callers` · call · `src/vehicle/kinematics/kinematics.jl:38-38`
- `callees` → [[core.quaternion_utils_dcm_to_quaternion|dcm_to_quaternion]] · `callers` · call · `src/vehicle/kinematics/kinematics.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations

The `@assert` guard is compiled out when Julia runs with `--check-bounds=no`-style optimisation of assertions disabled, so a root link could be silently corrupted. Setting orientation absolutely means the call discards the previous attitude rather than composing with it, which surprises callers expecting an incremental rotation. The axis fallback silently rotates about the y-axis instead of raising on a degenerate axis, and the DCM method trusts the caller to supply an orthonormal matrix.

## Provenance
Mapped from `src/vehicle/kinematics/kinematics.jl` line 32.
