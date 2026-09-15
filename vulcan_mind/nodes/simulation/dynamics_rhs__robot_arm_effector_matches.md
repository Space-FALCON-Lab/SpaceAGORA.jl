---
id: simulation.dynamics_rhs__robot_arm_effector_matches
label: _robot_arm_effector_matches
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _robot_arm_effector_matches
  lines:
  - 1620
  - 1620
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  description: Return value of `_robot_arm_effector_matches`.
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

# _robot_arm_effector_matches

## Purpose
Tests whether an effector is the robot-arm plan holder for a specific satellite, so the coupling lookup can find the one effector among the control and dynamic tuples that drives that satellite's arm.

## Design & Implementation
Requires the effector to expose both `plan` and `spacecraft_idx` properties, the index to equal `sat_idx`, and the plan to be a `RobotArmPlan`. All four conditions are combined with short-circuit `&&` so the type test runs only after the cheaper property checks. Declared `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_robot_arm_effector_matches`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__robot_arm_coupling|_robot_arm_coupling]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1682-1682`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1629-1629`
<!-- vulcan:connections:end -->

## Limitations
Duck-typed on property names, so an unrelated effector that happens to expose `plan` and `spacecraft_idx` is inspected on every call even though the final type test rejects it.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1620.
