---
id: environment.gravity_models_calcforcetorque
label: calcForceTorque
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: calcForceTorque
  lines:
  - 184
  - 184
inputs:
- id: model
  type: ConstantGravityModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: ComponentVector
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
- environment
charts:
- environment
origin: agent
---

# calcForceTorque

## Purpose
Legacy force/torque entry point (line 184) for `InverseSquaredGravityModel` on the `ComponentVector` state layout, used by older simulation paths that call `calcForceTorque(model, x, param, i)` for each effector. Sibling methods with the same shape exist for `ConstantGravityModel` and `InverseSquaredJ2GravityModel`.

## Design & Implementation
Extracts `pos_ii = SVector(x[1], x[2], x[3])` (m) and `mass = Float64(x[7])` (kg) by fixed index from the component vector, evaluates `_inverse_squared_gravity_accel(pos_ii, param.args.environment_model.planet)`, multiplies by mass, and computes `torque_body` with the legacy `_gravity_gradient_torque_body(model, pos_ii, x, param, i)`. Returns `(force_ii, torque_body)` as a tuple of two `SVector{3,Float64}` in N and N m.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ConstantGravityModel | n/a | yes | Positional argument `model`. |
| in | `x` | ComponentVector | n/a | yes | Positional argument `x`. |
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
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:265-265`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/gravity/gravity_models.jl:186-186`
- `callees` → [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callers` · call · `src/environment/gravity/gravity_models.jl:189-189`
- `callees` → [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:187-187`
<!-- vulcan:connections:end -->

## Limitations
The hard-coded indices 1:3 and 7 assume a state layout of position, velocity, mass; any reordering of the `ComponentVector` silently reads wrong values. The J2 sibling passes an inertial position to a routine that expects body-fixed coordinates. This path does not participate in the gravity-backbone or IMEX partitioning; those require the `wrench` interface. `i` is only used for the gravity-gradient inertia lookup.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 184.
