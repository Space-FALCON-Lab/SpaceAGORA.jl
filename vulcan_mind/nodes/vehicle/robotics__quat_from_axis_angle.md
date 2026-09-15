---
id: vehicle.robotics__quat_from_axis_angle
label: _quat_from_axis_angle
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _quat_from_axis_angle
  lines:
  - 94
  - 94
inputs:
- id: axis
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `axis`.
- id: theta
  type: Float64
  units: n/a
  required: true
  description: Positional argument `θ`.
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
  type: SVector{4,
  units: n/a
  description: Return value of `_quat_from_axis_angle`.
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

# _quat_from_axis_angle

## Purpose
Builds the rotation quaternion for a joint angle about its axis.

## Theory & Math
$$
q = \left( \hat{u} \sin\tfrac{\theta}{2},\; \cos\tfrac{\theta}{2} \right),\qquad \hat{u} = \frac{a}{\|a\|}
$$

## Design & Implementation
Normalises the axis, returning identity if it is zero or non-finite, then forms the quaternion with vector part `sin(θ/2)` times the unit axis and scalar part `cos(θ/2)`. `@inline` with static inputs and output.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `axis` | SVector{3, Float64} | n/a | yes | Positional argument `axis`. |
| in | `theta` | Float64 | n/a | yes | Positional argument `θ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `_quat_from_axis_angle`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:185-185`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The identity fallback for a degenerate axis means a misconfigured joint silently does not rotate rather than failing; `_normalize_axis` in the model builder is what actually rejects such axes.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 94.
