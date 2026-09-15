---
id: environment.gravity_models_environment_requirements
label: environment_requirements
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: environment_requirements
  lines:
  - 257
  - 257
inputs:
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
  type: EffectorEnvironmentRequirements
  units: n/a
  description: Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements(planet_frame=true)`.
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

# environment_requirements

## Purpose
Declares (line 257) that `InverseSquaredJ2GravityModel` needs the planet-fixed frame sample so its `wrench` and `gravity_backbone_acceleration_ii` methods can evaluate J2 in body-fixed coordinates and rotate the result back to J2000. The other two gravity models rely on the default (no requirements).

## Design & Implementation
`@inline environment_requirements(::InverseSquaredJ2GravityModel) = EffectorEnvironmentRequirements(planet_frame=true)`. The engine ORs this with every other effector's requirements when building the per-stage `EnvironmentSample`, causing `PlanetFrameSample` (with `l_pi`, `pos_pp`, `vel_pp`, altitude, latitude, longitude) to be computed once per stage.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EffectorEnvironmentRequirements | n/a | — | Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements(planet_frame=true)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:145-145`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:15-15`
- [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1500-1500`
- [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1517-1517`
- [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1318-1318`
- [[simulation.setup__any_effector_consumes_atmosphere|_any_effector_consumes_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:75-75`
- [[simulation.setup__warn_density_without_atmospheric_effector|_warn_density_without_atmospheric_effector]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:100-100`
- [[simulation.solver_policy__auto_stiff_smooth_gravity_reject_reason|_auto_stiff_smooth_gravity_reject_reason]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:131-131`
- [[simulation.solver_policy__gravity_backbone_reject_reason|_gravity_backbone_reject_reason]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:199-199`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:286-286`

**Downstream**

- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · feedback · `src/environment/gravity/gravity_models.jl:257-257`
- `callees` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/environment/gravity/gravity_models.jl:274-274`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:259-259`
- `callees` → [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/environment/gravity/gravity_models.jl:276-276`
- `callees` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/environment/gravity/gravity_models.jl:274-274`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:259-259`
- `callees` → [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callers` · call · `src/environment/gravity/gravity_models.jl:270-270`
- `callees` → [[environment.gravity_models__inverse_squared_j2_gravity_accel|_inverse_squared_j2_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:267-267`
- `callees` → [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/environment/gravity/gravity_models.jl:276-276`
- `callees` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/environment/gravity/gravity_models.jl:274-274`
- `callees` → [[environment.gravity_models_gravity_gradient|gravity_gradient]] · `callers` · call · `src/environment/gravity/gravity_models.jl:289-289`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:259-259`
<!-- vulcan:connections:end -->

## Limitations
Requesting the planet frame triggers a SPICE `pxform` evaluation and geodetic conversion on every ODE stage even though only `l_pi` and `pos_pp` are used here; the altitude and latitude computation is wasted work for this model. The requirement is static, so it cannot be relaxed when the model is inactive. If the sampler cannot furnish the frame (missing kernels), the failure surfaces as an `ArgumentError` inside `wrench` at the first evaluation, not at configuration time.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 257.
