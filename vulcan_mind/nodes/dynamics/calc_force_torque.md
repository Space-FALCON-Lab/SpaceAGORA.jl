---
id: dynamics.calc_force_torque
label: calcForceTorque
kind: function
source:
  file: src/dynamics/coupled/force_torque_models.jl
  symbol: calcForceTorque
  lines:
  - 5
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace registering concrete force and torque dispatch
    methods.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: wrench
  type: Tuple{SVector{3,Float64},SVector{3,Float64}}
  units: N,Nm
  description: Inertial force and body torque returned by a concrete effector method.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
- effectors
charts:
- dynamics
origin: agent
---

# calcForceTorque

## Purpose
`calcForceTorque` is the legacy generic dispatch surface for coupled dynamic effectors. The declaration in `force_torque_models.jl` establishes one function name; concrete gravity, aerodynamic, perturbation, thruster, guidance, and robot-arm models add methods in the included files. The RHS calls the method appropriate for each effector instance.

## Theory & Math
Each method returns a wrench pair `(F_ii, τ_b)`, with force in inertial axes and torque in body axes. The dynamics RHS sums compatible force and torque contributions before applying `m a = F` and the rigid-body rotational equation. The generic declaration itself performs no computation.

## Model & Assumptions
Concrete methods assume the state vector, ODE parameters, and satellite index use the same layout and frame conventions. Effector contributions are additive, and models must not silently return a body-frame force or inertial-frame torque. Environment requirements are exposed separately through the sampling API.

## Design & Implementation
Line 5 declares `function calcForceTorque end` before the six effector includes. Gravity and perturbation files define model-specific methods, while `dynamics_rhs.jl` invokes `SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)`. The sampled `wrench` path is a newer parallel interface, but this generic remains part of the public API and compatibility surface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace registering concrete force and torque dispatch methods. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `wrench` | Tuple{SVector{3,Float64},SVector{3,Float64}} | N,Nm | — | Inertial force and body torque returned by a concrete effector method. |
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
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:265-265`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Calling the generic without a concrete method raises a dispatch error. The contract does not itself check units, frame labels, finite values, or mass validity; individual methods and the RHS provide those guards inconsistently by model. Legacy methods can allocate or access shared state differently from the sampled wrench path.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models.jl:5` and its concrete methods under `src/dynamics/coupled/` and `src/environment/gravity/`.
