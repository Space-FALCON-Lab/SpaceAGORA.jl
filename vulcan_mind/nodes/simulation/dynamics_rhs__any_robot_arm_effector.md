---
id: simulation.dynamics_rhs__any_robot_arm_effector
label: _any_robot_arm_effector
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _any_robot_arm_effector
  lines:
  - 1651
  - 1651
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Bool
  units: n/a
  description: Return value of `_any_robot_arm_effector`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _any_robot_arm_effector

## Purpose
Determines once per run whether any effector carries a `RobotArmPlan`, so the per-satellite coupling scan can be skipped entirely in the common no-arm case.

## Design & Implementation
Checks the control effector tuple and then the dynamic effector tuple, guarded by `hasproperty` on the configuration, returning true on the first effector whose `plan` property is a `RobotArmPlan`. The result is cached in `shared_buffers.robot_arm_present` by `_robot_arm_present`. Returns `Bool`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_any_robot_arm_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__robot_arm_present|_robot_arm_present]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1673-1673`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Duck-typed on a property named `plan`, so an unrelated effector with a `plan` field of another type is inspected but correctly rejected; the scan itself is cheap but runs `hasproperty` dynamically.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1651.
