---
id: vehicle.robotics__normalize_axis
label: _normalize_axis
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _normalize_axis
  lines:
  - 149
  - 149
inputs:
- id: axis
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `axis`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_normalize_axis`.
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

# _normalize_axis

## Purpose
Produces a unit joint axis and refuses a zero-length one at model construction time, where the error is actionable.

## Design & Implementation
Takes the norm and throws `ArgumentError` if it is not above machine epsilon, otherwise returns the axis divided by its norm. Contrast with `_quat_from_axis_angle`, which tolerates a bad axis at evaluation time; the strictness lives here so bad models are rejected once rather than silently degraded on every FK call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `axis` | SVector{3, Float64} | n/a | yes | Positional argument `axis`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_normalize_axis`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`
- [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:140-140`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A non-finite axis passes the comparison as false and is rejected, which is correct, but the message does not distinguish NaN from zero.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 149.
