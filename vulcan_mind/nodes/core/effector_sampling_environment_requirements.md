---
id: core.effector_sampling_environment_requirements
label: environment_requirements
kind: function
source:
  file: src/core/types/effector_sampling.jl
  symbol: environment_requirements
  lines:
  - 151
  - 151
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
  description: Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements()`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# environment_requirements

## Purpose
Additive declaration hook that an effector model overloads to state which sampled environment fields its `wrench` needs. The default method makes every model opt out, so adding a new effector type never breaks sampling of existing ones.

## Design & Implementation
Defined as `@inline environment_requirements(::Any) = EffectorEnvironmentRequirements()`, returning the all-false request with an empty `third_body_names` tuple. Effector packages add methods such as `environment_requirements(::DragModel) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)`. The engine calls this once per effector when assembling the environment sampling plan, not per ODE evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EffectorEnvironmentRequirements | n/a | — | Return value of `environment_requirements`. Returns `EffectorEnvironmentRequirements()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:103-103`
- [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callees` → `callers` · feedback · `src/core/types/effector_sampling.jl:145-145`
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

- `callees` → [[core.effector_sampling_effectorenvironmentrequirements|EffectorEnvironmentRequirements]] · `callers` · call · `src/core/types/effector_sampling.jl:151-151`
- `callees` → [[core.effector_sampling_wrench_caching_bang|wrench_caching!]] · `callers` · call · `src/core/types/effector_sampling.jl:165-165`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:154-154`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `callers` · call · `src/core/types/effector_sampling.jl:165-165`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:154-154`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:154-154`
<!-- vulcan:connections:end -->

## Limitations
Since the fallback matches `::Any`, a typo in an effector's method signature (for example the wrong type name) silently selects the default and the effector's `wrench` then receives `nothing` for the fields it needs, failing later with a `MethodError` on `nothing`. The requirement cannot depend on runtime state or configuration flags because it is evaluated from the model value alone.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 151.
