---
id: core.quaternion_utils_qtoeulerangles
label: qToEulerAngles
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: qToEulerAngles
  lines:
  - 75
  - 75
inputs:
- id: q
  type: Any
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
  type: AbstractArray
  units: n/a
  description: Return value of `qToEulerAngles`. Returns `[roll, pitch, yaw]`.
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

# qToEulerAngles

## Purpose

Converts a quaternion into three Euler angles in a 3-2-1 (yaw-pitch-roll) sequence, giving telemetry, plotting and reporting code a human-readable attitude representation instead of raw quaternion components.

## Design & Implementation

`qToEulerAngles(q)` first calls `normalize(q)` and destructures the result as `qx, qy, qz, qw`, again scalar-last. Roll uses `atan(cosr_cosp, sinr_cosp)` with `sinr_cosp = 2(qw*qx + qy*qz)` and `cosr_cosp = 1 - 2(qx^2 + qy^2)`; pitch is formed from `sqrt(1 + 2(qw*qy - qx*qz))` and `sqrt(1 - 2(qw*qy - qx*qz))` combined as `2*atan(cosp, sinp) - pi/2`; yaw mirrors the roll form using the z components. It returns a plain three-element `Vector{Float64}` of angles in radians.

## Theory & Math

For a unit quaternion $(q_x, q_y, q_z, q_w)$ the standard 3-2-1 extraction is

$$\phi = \operatorname{atan2}\!\left(2(q_w q_x + q_y q_z),\; 1 - 2(q_x^2 + q_y^2)\right)$$
$$\theta = \arcsin\!\left(2(q_w q_y - q_x q_z)\right)$$
$$\psi = \operatorname{atan2}\!\left(2(q_w q_z + q_x q_y),\; 1 - 2(q_y^2 + q_z^2)\right)$$

with $\phi$ roll, $\theta$ pitch and $\psi$ yaw in radians. The pitch branch here uses the numerically stable half-angle form $\theta = 2\operatorname{atan2}(\sqrt{1-2u}, \sqrt{1+2u}) - \pi/2$ with $u = q_w q_y - q_x q_z$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `qToEulerAngles`. Returns `[roll, pitch, yaw]`. |
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

The roll and yaw calls pass their arguments to `atan` in the order (cosine, sine) rather than the conventional (sine, cosine), so those two angles come out as the complement of the usual definition; consumers must know which convention this project treats as canonical. The 3-2-1 sequence is singular at pitch of plus or minus 90 degrees, where roll and yaw are no longer separable. The function allocates a heap `Vector` on every call, which is costly inside a per-step telemetry loop, and it returns radians with no wrapping applied.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 75.
