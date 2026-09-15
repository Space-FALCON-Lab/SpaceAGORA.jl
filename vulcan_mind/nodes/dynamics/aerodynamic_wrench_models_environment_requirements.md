---
id: dynamics.aerodynamic_wrench_models_environment_requirements
label: environment_requirements
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: environment_requirements
  lines:
  - 224
  - 224
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
  description: Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements(planet_frame=true,
    atmosphere=true)`.
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

# environment_requirements

## Purpose
Declares that all three aerodynamic effector types need a planet-fixed frame and an atmosphere sample from the environment stage.

## Design & Implementation
Three `@inline` methods, one per model type, each returning `EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)`. The solver policy consults these to reject smooth-gravity `Tsit5` routing and gravity-backbone splitting whenever an aerodynamic effector is active.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EffectorEnvironmentRequirements | n/a | — | Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)`. |
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

- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · feedback · `src/dynamics/coupled/aerodynamic_wrench_models.jl:224-224`
<!-- vulcan:connections:end -->

## Limitations
Solar and third-body requirements default to false, which is correct, but the requirement is static; per-link atmosphere sampling (which needs the live density model through `ODEParams`) is not expressed here.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 224.
