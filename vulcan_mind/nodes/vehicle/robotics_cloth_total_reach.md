---
id: vehicle.robotics_cloth_total_reach
label: cloth_total_reach
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: cloth_total_reach
  lines:
  - 156
  - 156
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: Float64
  units: n/a
  description: Return value of `cloth_total_reach`.
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

# cloth_total_reach

## Purpose
Reports the arm's maximum extension, the sum of link lengths, used to reject targets that cannot possibly be reached before running inverse kinematics.

## Design & Implementation
Sums `norm(link.vector_parent)` over `model.links`. Using the tip vector norm rather than a stored length means a link defined with a non-axial tip vector is still measured correctly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `cloth_total_reach`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It ignores joint limits and the mount offset, so it overstates the reachable workspace for any real configuration; it is an upper bound, not a reachability test.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 156.
