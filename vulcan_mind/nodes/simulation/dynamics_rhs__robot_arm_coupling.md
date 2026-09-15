---
id: simulation.dynamics_rhs__robot_arm_coupling
label: _robot_arm_coupling
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _robot_arm_coupling
  lines:
  - 1679
  - 1679
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: Nothing
  units: n/a
  description: Return value of `_robot_arm_coupling`. Returns `nothing`.
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

# _robot_arm_coupling

## Purpose
Finds the robot-arm coupling parameters for a satellite by locating the effector that carries its plan.

## Design & Implementation
Scans the control effector tuple and then the dynamic effector tuple for an effector matching `_robot_arm_effector_matches` with this `sat_idx`, returning `_robot_arm_coupling_from_effector` for the first hit or `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_robot_arm_coupling`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1698-1698`
- [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2315-2315`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__robot_arm_coupling_from_effector|_robot_arm_coupling_from_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1682-1682`
- `callees` → [[simulation.dynamics_rhs__robot_arm_effector_matches|_robot_arm_effector_matches]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1682-1682`
<!-- vulcan:connections:end -->

## Limitations
A linear scan on every RHS call for every satellite when any arm is present; the presence cache only avoids it in the no-arm case.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1679.
