---
id: dynamics.perturbations_calcforcetorque
label: calcForceTorque
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: calcForceTorque
  lines:
  - 918
  - 918
inputs:
- id: model
  type: NBodyGravityModel
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
The legacy force-torque entry point for the N-body effector, computing the third-body perturbation with cache, memo and optional threading over bodies.

## Design & Implementation
Reads position and mass, resolves the epoch from shared buffers, asks the multibody thread policy whether to parallelise over bodies, obtains the satellite's workspace, and for each body fetches its position from the ephemeris cache or SPICE. It then evaluates `_nbody_body_force_ii` per body — threaded when decided — sums in body order, and returns force plus zero torque.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | NBodyGravityModel | n/a | yes | Positional argument `model`. |
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
- [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `callees` → `callers` · feedback · `src/dynamics/coupled/perturbations.jl:1500-1500`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2245-2245`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2033-2033`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1892-1892`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:265-265`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:920-920`
- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:925-925`
- `callees` → [[dynamics.perturbations__nbody_body_force_ii|_nbody_body_force_ii]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:951-951`
- `callees` → [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:936-936`
- `callees` → [[dynamics.perturbations__nbody_body_position_from_spice_j2000_m|_nbody_body_position_from_spice_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:937-937`
- `callees` → [[dynamics.perturbations__nbody_workspace_for_sat_bang|_nbody_workspace_for_sat!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:928-928`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:921-921`
- `callees` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:927-927`
- `callees` → [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:950-950`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:961-961`
<!-- vulcan:connections:end -->

## Limitations
Body-level threading only helps with many bodies; the position fetch loop is serial and takes the SPICE lock on cache misses.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 918.
