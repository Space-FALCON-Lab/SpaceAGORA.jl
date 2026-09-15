---
id: simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang
label: spacecraft_dynamics_explicit_remainder!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics_explicit_remainder!
  lines:
  - 2077
  - 2077
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
  description: Return value of `spacecraft_dynamics_explicit_remainder!`; mutates
    `du` in place. Returns `_spacecraft_dynamics_flat_constellation_effector_queue!(`
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

# spacecraft_dynamics_explicit_remainder!

## Purpose
The explicit-partition RHS for the split solver: everything not in the implicit backbone, including explicit effectors, control and coupling.

## Design & Implementation
Obtains the execution plan, routes to the flat queue with `partition=:explicit` if chosen, otherwise runs batched or serial per-satellite accumulation of explicit-partition effectors plus control, and assigns derivatives.

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
| out | `result` | Any | n/a | — | Return value of `spacecraft_dynamics_explicit_remainder!`; mutates `du` in place. Returns `_spacecraft_dynamics_flat_constellation_effector_queue!(` or `merge(base_shape, SimulationModel.coupled_cloth_robot_arm_state_shape(coupling.p` or `state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2126-2126`
- `callees` → [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2116-2116`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2114-2114`
- `callees` → [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2117-2117`
- `callees` → [[simulation.dynamics_rhs__assign_heat_rate_derivative_bang|_assign_heat_rate_derivative!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2147-2147`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2135-2135`
- `callees` → [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2100-2100`
- `callees` → [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2087-2087`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2098-2098`
- `callees` → [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2085-2085`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2118-2118`
<!-- vulcan:connections:end -->

## Limitations
Duplicates the per-satellite tail of the other RHS variants.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 2077.
