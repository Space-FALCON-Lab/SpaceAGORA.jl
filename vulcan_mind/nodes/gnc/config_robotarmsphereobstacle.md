---
id: gnc.config_robotarmsphereobstacle
label: RobotArmSphereObstacle
kind: struct
source:
  file: src/gnc/robotics/robot_arm_hypr/config.jl
  symbol: RobotArmSphereObstacle
  lines:
  - 2
  - 2
inputs:
- id: center
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `center`.
- id: radius_m
  type: Float64
  units: n/a
  required: true
  description: Field `radius_m`.
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
  type: RobotArmSphereObstacle
  units: n/a
  description: Constructed `RobotArmSphereObstacle`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# RobotArmSphereObstacle

## Purpose
`RobotArmSphereObstacle` models a workspace keep-out volume as a sphere for the robot-arm HYPR planner's clearance checks. Cost evaluation and collision sampling test candidate arm configurations against a collection of these spheres, making them the obstacle primitive of the whole robot-arm planning path.

## Design & Implementation
An immutable struct with two fields, `center::SVector{3, Float64}` in metres and `radius_m::Float64`. A convenience outer constructor `RobotArmSphereObstacle(center, radius_m::Real)` converts any 3-element container via `SVector{3, Float64}(center)` and the radius via `Float64(radius_m)`, so callers may pass a plain `Vector` or a tuple. The `SVector` choice keeps obstacles stack-allocated and isbits, which matters because the inner collision loop evaluates distances to every obstacle at every one of the configured `n_samples` path samples per particle per iteration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `center` | SVector{3, Float64} | n/a | yes | Field `center`. |
| in | `radius_m` | Float64 | n/a | yes | Field `radius_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmSphereObstacle | n/a | — | Constructed `RobotArmSphereObstacle`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/config.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/config.jl:8-8`
<!-- vulcan:connections:end -->

## Limitations
No validation is performed: a negative or `NaN` `radius_m` is accepted and will make every clearance test behave unpredictably, and a `center` of the wrong length throws a `StaticArrays` conversion error rather than an explanatory one. Only spheres are representable, so a long or flat obstacle must be approximated by a conservative bounding sphere or by many overlapping spheres, which inflates the per-sample cost linearly.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/config.jl` line 2.
