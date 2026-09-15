---
id: simx.engine_dynamics_rhs_spacecraft_dynamics_bang
label: spacecraft_dynamics!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics!
  lines:
  - 1716
  - 1829
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: state_u
  type: ComponentVector
  units: m,m/s,kg,J
  required: true
  description: Full stacked constellation state; `u.sc[i]` is spacecraft i's translational,
    attitude, mass and heat-load block.
- id: params_p
  type: ODEParams
  units: n/a
  required: true
  description: Solver parameter object carrying args, shared_buffers, is_active flags
    and every preallocated workspace.
- id: time_t
  type: Float64
  units: s
  required: true
  description: Mission-elapsed time at which the derivative is requested; also stamped
    into shared_buffers.current_time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: derivative_du
  type: ComponentVector
  units: m/s,m/s^2,kg/s,W
  description: In-place state derivative written per spacecraft; inactive spacecraft
    blocks are zeroed.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# spacecraft_dynamics!

## Purpose
`spacecraft_dynamics!` is the full-fidelity right-hand side handed to the ODE solver. For every active spacecraft it accumulates dynamic-effector and control-effector wrenches, resolves robot-arm coupling, evaluates stage heat rates, and writes translational, rotational and thermal derivatives into the supplied `du` in place.

## Theory & Math
Per spacecraft the assembled derivative is the SI-native set

$$\dot{\vec r} = \vec v, \qquad \dot{\vec v} = \frac{1}{m}\sum_k \vec F_k, \qquad \dot m = \sum_k \dot m_k$$

for translation, and, when `orientation_sim` is enabled, Euler's equation with reaction-wheel exchange

$$\mathbf{J}\,\dot{\vec\omega} = \sum_k \vec\tau_k - \vec\tau_{rw} - \vec\omega \times \left(\mathbf{J}\vec\omega + \vec h_{rw}\right)$$

together with the quaternion kinematics

$$\dot{q} = \tfrac{1}{2}\, q \otimes \begin{bmatrix} 0 \\ \vec\omega \end{bmatrix}$$

where $\mathbf{J}$ is `spacecraft[i].inertia_tensor` in kg m^2, $\vec h_{rw}$ is the wheel angular momentum from the `rw_assembly`, and the gyroscopic term is retained because `include_gyroscopic=true` on both branches. Heat loads integrate $\dot{Q}_j$ in W supplied by `_compute_stage_heat_rates!`.

## Model & Assumptions
Forces and torques are accumulated into `MVector{3,Float64}` scratch so no heap allocation occurs per spacecraft per call. Inactive spacecraft, flagged by `p.is_active[i]`, have their derivative block set to zero rather than being skipped, which keeps the stacked state vector shape constant for the solver. Attitude propagation is conditional on `p.args.mission_configuration.orientation_sim`; translation-only runs never touch the quaternion or angular-rate slots.

## Design & Implementation
The function first stamps `shared_buffers.current_time[]` and asks `_rhs_execution_plan` for the execution plan. If the plan mode is `:flat_constellation_effector_queue` it delegates immediately to `_spacecraft_dynamics_flat_constellation_effector_queue!` with `rhs_kind=:full` and returns. Otherwise it chooses between a threaded and a serial satellite loop: threading is enabled when the plan mode is not `:serial` and `_rhs_batch_parallel_enabled` agrees for the current spacecraft count. The threaded branch first calls `_prefill_shared_body_samples!` so all shared ephemeris and environment samples exist before workers start, then runs Polyester's `@batch` with `minbatch = max(1, ceil(n_sats / Polyester.num_cores()))`. Both branches execute an identical body under `@views`, and the serial branch adds `@inbounds`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `state_u` | ComponentVector | m,m/s,kg,J | yes | Full stacked constellation state; `u.sc[i]` is spacecraft i's translational, attitude, mass and heat-load block. |
| in | `params_p` | ODEParams | n/a | yes | Solver parameter object carrying args, shared_buffers, is_active flags and every preallocated workspace. |
| in | `time_t` | Float64 | s | yes | Mission-elapsed time at which the derivative is requested; also stamped into shared_buffers.current_time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `derivative_du` | ComponentVector | m/s,m/s^2,kg/s,W | — | In-place state derivative written per spacecraft; inactive spacecraft blocks are zeroed. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:276-276`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1757-1757`
- `callees` → [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1747-1747`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1745-1745`
- `callees` → [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1748-1748`
- `callees` → [[simulation.dynamics_rhs__assign_heat_rate_derivative_bang|_assign_heat_rate_derivative!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1778-1778`
- `callees` → [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1766-1766`
- `callees` → [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1731-1731`
- `callees` → [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1726-1726`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1729-1729`
- `callees` → [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1724-1724`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1749-1749`
<!-- vulcan:connections:end -->

## Limitations
The threaded and serial bodies are duplicated verbatim, so a change to one must be mirrored in the other or the two paths silently diverge. Prefilling shared samples only happens on the threaded branch, which means any effector reaching for a shared sample must tolerate lazy population when running serially. `@batch` requires the per-spacecraft work to be free of data races on `shared_buffers`, so effectors that write unbuffered shared scratch must be routed to the serial plan by the calibration layer.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl:1716-1829`, with the flat constellation queue at line 1295 and the shared-sample prefill at line 1212 of the same file.
