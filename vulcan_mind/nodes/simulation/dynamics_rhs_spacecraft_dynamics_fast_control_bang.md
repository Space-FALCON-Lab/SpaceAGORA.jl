---
id: simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang
label: spacecraft_dynamics_fast_control!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics_fast_control!
  lines:
  - 2200
  - 2200
inputs:
- id: du
  type: ComponentVector
  units: n/a
  required: true
  description: Positional argument `du`.
- id: u
  type: ComponentVector
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `spacecraft_dynamics_fast_control!`; mutates `du` in
    place. Returns `merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.p`
    or `state`.
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

# spacecraft_dynamics_fast_control!

## Purpose
The fast-rate RHS evaluating only control effectors and robot-arm coupling, for solver configurations that split control from slow dynamics.

## Design & Implementation
Batched or serial loop over satellites accumulating control effectors, applying coupling, and assigning derivatives with zero dynamic forces.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du` | ComponentVector | n/a | yes | Positional argument `du`. |
| in | `u` | ComponentVector | n/a | yes | Positional argument `u`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `spacecraft_dynamics_fast_control!`; mutates `du` in place. Returns `merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.p` or `state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang|assign_control_only_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2223-2223`
- `callees` → [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2220-2220`
- `callees` → [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2221-2221`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2232-2232`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2206-2206`
<!-- vulcan:connections:end -->

## Limitations
Assumes the slow RHS supplies gravity and aerodynamics.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 2200.
