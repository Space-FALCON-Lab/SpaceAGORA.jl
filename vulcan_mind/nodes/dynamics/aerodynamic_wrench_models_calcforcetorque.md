---
id: dynamics.aerodynamic_wrench_models_calcforcetorque
label: calcForceTorque
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: calcForceTorque
  lines:
  - 588
  - 588
inputs:
- id: model
  type: AerodynamicCoefficientConstant
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: AbstractVector{Float64}
  units: n/a
  required: true
  description: Positional argument `x`.
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `calcForceTorque`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# calcForceTorque

## Purpose
Legacy `ODEParams`-based force/torque entry point for the three aerodynamic models; only the `AerodynamicCoefficientfM` method is functional and computes a threaded per-link free-molecular force sum.

## Theory & Math
Per link: $\mathbf{F}_D = q\,C_D A\,\hat{\mathbf{d}}$, $\mathbf{F}_L = q\cos\phi\,C_L A\,\hat{\mathbf{l}}$, $\mathbf{F}_S = q\,C_S A\,\hat{\mathbf{s}}$ with $q = \tfrac{1}{2}\rho\,|\mathbf{v}_{rw}|^2$, molecular speed ratio $S = \sqrt{\gamma/2}\,M$, $M = |\mathbf{v}_{rw}|/\sqrt{\gamma R T}$, bank angle $\phi = 0$.

## Design & Implementation
The fM method stages the environment via `SimulationCallbacks._stage_environment_state(x, param, i, current_time; write_buffers=true)`, validates the incidence symbol, returns zeros early for vacuum (`rho <= eps`), computes speed of sound `sqrt(γ R T)`, wind-relative velocity `vel_pp - wind_pp` (with NED wind from `latlongtoNED`), Mach, molecular speed ratio `S = sqrt(γ/2) * Mach`, dynamic pressure `q = 0.5 ρ v²`, and drag/lift/cross unit vectors from the planet-relative angular momentum. It obtains a thread decision, a per-satellite `AeroScratchWorkspace`, fills per-link slots through `compute_link_wrench!` (threaded via `ParallelPolicy.threaded_foreach_worker_persistent(:rhs_aero, ...)` or serially), sums slots in link order for determinism, records a policy observation, normalises CL/CD by total area, stores caches, and returns `(force_ii, zero torque)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerodynamicCoefficientConstant | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector{Float64} | n/a | yes | Positional argument `x`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `calcForceTorque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:106-106`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1655-1655`
- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:916-916`
- [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1500-1500`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:265-265`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:654-654`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:748-748`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:603-603`
- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:773-773`
- `callees` → [[dynamics.aerodynamic_wrench_models__store_aero_caches_bang|_store_aero_caches!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:721-721`
- `callees` → [[dynamics.aerodynamic_wrench_models__validate_fm_incidence|_validate_fm_incidence]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:717-717`
- `callees` → [[dynamics.aerodynamic_wrench_models_collect_and_reset_link_wrenches_bang|collect_and_reset_link_wrenches!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:694-694`
- `callees` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:775-775`
- `callees` → [[simulation.runtime__stage_environment_state|_stage_environment_state]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:706-706`
- `callees` → [[vehx.kinematics_rotate_to_inertial|rotate_to_inertial]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:623-623`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:593-593`
<!-- vulcan:connections:end -->

## Limitations
The `AerodynamicCoefficientConstant` and `NoBallisticFlight` methods are dead code: they reference undefined names (`vel_pp_rw_hat`, `q`, `args`, `CS_body`, `rot_body_to_inertial`, `drag_ii`, `drag_pp`, `m.body`) and the latter contains a `println`, so both throw `UndefVarError` if dispatched. The fM method always returns zero torque even with `orientation_sim`, unlike `_aero_pure_wrench`. Several computed quantities (`h_ii`, `pos_pp_mag`, the first `mach`/`S`) are unused.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 588.
