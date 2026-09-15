---
id: simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang
label: spacecraft_dynamics_implicit_atmosphere!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics_implicit_atmosphere!
  lines:
  - 1970
  - 1970
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
  type: Nothing
  units: n/a
  description: Return value of `spacecraft_dynamics_implicit_atmosphere!`; mutates
    `du` in place. Returns `nothing` or `_spacecraft_dynamics_flat_constellation_effector_queue!(`
    or `merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.p`
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

# spacecraft_dynamics_implicit_atmosphere!

## Purpose
The implicit-partition RHS: atmosphere-dependent effectors evaluated implicitly, skipped entirely when all satellites are above a vanishing atmosphere.

## Design & Implementation
Short-circuits to zero derivatives when outside the atmosphere, obtains the plan, routes to the flat queue with `partition=:implicit` if chosen, otherwise accumulates implicit-partition effectors per satellite.

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
| out | `result` | Nothing | n/a | — | Return value of `spacecraft_dynamics_implicit_atmosphere!`; mutates `du` in place. Returns `nothing` or `_spacecraft_dynamics_flat_constellation_effector_queue!(` or `merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.p` or `state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang|assign_force_only_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2017-2017`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2015-2015`
- `callees` → [[simulation.dynamics_rhs__all_active_spacecraft_outside_atmosphere|_all_active_spacecraft_outside_atmosphere]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1982-1982`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2025-2025`
- `callees` → [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2001-2001`
- `callees` → [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1988-1988`
- `callees` → [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2006-2006`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1999-1999`
- `callees` → [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1986-1986`
<!-- vulcan:connections:end -->

## Limitations
The early exit is only available for `NoAtmosphereModel`.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1970.
