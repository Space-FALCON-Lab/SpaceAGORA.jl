---
id: simulation.dynamics_rhs__robot_arm_coupling_from_effector
label: _robot_arm_coupling_from_effector
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _robot_arm_coupling_from_effector
  lines:
  - 1636
  - 1636
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_robot_arm_coupling_from_effector`. Returns `(`.
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

# _robot_arm_coupling_from_effector

## Purpose
Extracts the plan, elapsed plan time, stiffness and damping parameters and joint actuators from a robot-arm effector, supplying defaults for any it does not declare.

## Design & Implementation
Returns a named tuple with `plan`, `t_s` as `t - updated_at_s` floored at zero, translational and rotational stiffness and damping defaulting to 5e3 N/m, 30 N·s/m, 15 N·m/rad and 0.5 N·m·s/rad, and `joint_actuators` defaulting to an empty vector. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_coupling_from_effector`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__robot_arm_coupling|_robot_arm_coupling]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1682-1682`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__effector_value|_effector_value]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1640-1640`
<!-- vulcan:connections:end -->

## Limitations
Duck-typed property access with defaults that mirror the grid builder's, so a mismatch between the two default sets would be invisible.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1636.
