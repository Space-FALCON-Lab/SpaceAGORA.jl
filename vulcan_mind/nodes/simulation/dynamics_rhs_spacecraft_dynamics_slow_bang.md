---
id: simulation.dynamics_rhs_spacecraft_dynamics_slow_bang
label: spacecraft_dynamics_slow!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics_slow!
  lines:
  - 1831
  - 1831
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
  description: Return value of `spacecraft_dynamics_slow!`; mutates `du` in place.
    Returns `_spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, t,
    plan; rhs_k` or `sat_idx <= length(times) && times[sat_idx] == t`.
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

# spacecraft_dynamics_slow!

## Purpose
The main monolithic RHS: all dynamic and control effectors for all satellites, routed through the execution plan.

## Design & Implementation
Obtains the plan, routes to the flat queue if chosen, otherwise prefills shared body samples when batching and loops satellites batched or serially accumulating dynamic and control effectors, coupling, and derivatives.

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
| out | `result` | Any | n/a | — | Return value of `spacecraft_dynamics_slow!`; mutates `du` in place. Returns `_spacecraft_dynamics_flat_constellation_effector_queue!(du, u, p, t, plan; rhs_k` or `sat_idx <= length(times) && times[sat_idx] == t`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang|assign_slow_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1869-1869`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1859-1859`
- `callees` → [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1860-1860`
- `callees` → [[simulation.dynamics_rhs__assign_heat_rate_derivative_bang|_assign_heat_rate_derivative!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1888-1888`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1877-1877`
- `callees` → [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1845-1845`
- `callees` → [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1840-1840`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1843-1843`
- `callees` → [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1838-1838`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1861-1861`
<!-- vulcan:connections:end -->

## Limitations
Around 360 lines with the per-satellite tail duplicated across variants.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1831.
