---
id: simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang
label: _spacecraft_dynamics_flat_constellation_effector_queue!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _spacecraft_dynamics_flat_constellation_effector_queue!
  lines:
  - 1295
  - 1295
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
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: rhs_kind
  type: Symbol
  units: n/a
  required: true
  description: Keyword argument `rhs_kind`.
- id: partition
  type: Union{Nothing, Symbol}
  units: n/a
  required: false
  description: Keyword argument `partition` (default `nothing`).
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
  description: Return value of `_spacecraft_dynamics_flat_constellation_effector_queue!`;
    mutates `du` in place. Returns `nothing`.
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

# _spacecraft_dynamics_flat_constellation_effector_queue!

## Purpose
The complete RHS under the flat constellation plan: shared samples are prefilled once, dynamic effectors are evaluated across the worker pool, then each satellite's control effectors, robot-arm coupling and derivatives are assigned.

## Design & Implementation
Sets the current time, prefills shared body samples, prefills planet frame and atmosphere when any wrench effector needs them, calls `_accumulate_dynamic_effectors_flat_batch!`, then loops satellites reading force and torque from the totals, accumulating control effectors, applying robot-arm coupling, and assigning translational, mass, heat-load and orientation derivatives according to `rhs_kind` (`:slow`, `:explicit` or `:implicit`) and the partition.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du` | ComponentVector | n/a | yes | Positional argument `du`. |
| in | `u` | ComponentVector | n/a | yes | Positional argument `u`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `rhs_kind` | Symbol | n/a | yes | Keyword argument `rhs_kind`. |
| in | `partition` | Union{Nothing, Symbol} | n/a | no | Keyword argument `partition` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_spacecraft_dynamics_flat_constellation_effector_queue!`; mutates `du` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2087-2087`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1988-1988`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1840-1840`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1726-1726`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1318-1318`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1318-1318`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1318-1318`
- `callees` → [[dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang|assign_force_only_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1352-1352`
- `callees` → [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1410-1410`
- `callees` → [[dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang|assign_slow_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1378-1378`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1318-1318`
- `callees` → [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1341-1341`
- `callees` → [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1399-1399`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1332-1332`
- `callees` → [[simulation.dynamics_rhs__assign_heat_rate_derivative_bang|_assign_heat_rate_derivative!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1395-1395`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1359-1359`
- `callees` → [[simulation.dynamics_rhs__flat_totals_force_torque|_flat_totals_force_torque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1350-1350`
- `callees` → [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1326-1326`
- `callees` → [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1312-1312`
- `callees` → [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1343-1343`
- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1317-1317`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1371-1371`
<!-- vulcan:connections:end -->

## Limitations
Roughly 350 lines; the per-satellite tail duplicates the batch and slow RHS variants, and the `rhs_kind` branching means three different derivative layouts share one function.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1295.
