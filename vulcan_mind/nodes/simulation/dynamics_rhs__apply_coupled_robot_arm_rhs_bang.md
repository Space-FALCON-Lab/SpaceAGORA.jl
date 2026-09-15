---
id: simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang
label: _apply_coupled_robot_arm_rhs!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _apply_coupled_robot_arm_rhs!
  lines:
  - 1694
  - 1694
inputs:
- id: du_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `du_view`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
- id: forces
  type: Any
  units: n/a
  required: true
  description: Positional argument `forces`.
- id: torques
  type: Any
  units: n/a
  required: true
  description: Positional argument `torques`.
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
  description: Return value of `_apply_coupled_robot_arm_rhs!`; mutates `du_view`
    in place. Returns `nothing`.
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

# _apply_coupled_robot_arm_rhs!

## Purpose
Adds the coupled cloth-robot-arm dynamics to a satellite's RHS when that satellite has an active arm plan.

## Design & Implementation
Returns immediately unless a robot arm is present in the run and a coupling for this satellite is found, then calls `assign_coupled_cloth_robot_arm_rhs!` with the plan, elapsed time, forces, torques and stiffness parameters. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `forces` | Any | n/a | yes | Positional argument `forces`. |
| in | `torques` | Any | n/a | yes | Positional argument `torques`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_apply_coupled_robot_arm_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2117-2117`
- [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2221-2221`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1860-1860`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1748-1748`

**Downstream**

- `callees` → [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1700-1700`
- `callees` → [[simulation.dynamics_rhs__robot_arm_coupling|_robot_arm_coupling]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1698-1698`
- `callees` → [[simulation.dynamics_rhs__robot_arm_present|_robot_arm_present]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1697-1697`
<!-- vulcan:connections:end -->

## Limitations
The presence check is cached, but the per-satellite coupling lookup scans effector tuples on every call when an arm exists.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1694.
