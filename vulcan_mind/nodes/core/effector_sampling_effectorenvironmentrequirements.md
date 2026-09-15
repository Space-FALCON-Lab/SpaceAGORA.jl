---
id: core.effector_sampling_effectorenvironmentrequirements
label: EffectorEnvironmentRequirements
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: EffectorEnvironmentRequirements
  lines:
  - 130
  - 130
inputs:
- id: planet_frame
  type: Bool
  units: n/a
  required: true
  description: Field `planet_frame`.
- id: atmosphere
  type: Bool
  units: n/a
  required: true
  description: Field `atmosphere`.
- id: solar
  type: Bool
  units: n/a
  required: true
  description: Field `solar`.
- id: third_body_names
  type: TB
  units: n/a
  required: true
  description: Field `third_body_names`.
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
  description: Constructed `EffectorEnvironmentRequirements`.
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

# EffectorEnvironmentRequirements

## Purpose
Capability request that an effector returns from `environment_requirements(model)` to tell the sampling layer which optional `EnvironmentSample` fields must be built for its `wrench` evaluation, so the engine only pays for planet-frame, atmosphere, solar or third-body sampling when some effector needs it.

## Design & Implementation
`struct EffectorEnvironmentRequirements{TB <: Tuple{Vararg{String}}}` with `planet_frame::Bool`, `atmosphere::Bool`, `solar::Bool` and `third_body_names::TB`. A keyword constructor defaults every flag to `false` and `third_body_names` to the empty tuple `()`, which is the value returned by the `::Any` fallback of `environment_requirements`. The engine merges requirements across all effectors with logical OR on the flags and a union of body names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_frame` | Bool | n/a | yes | Field `planet_frame`. |
| in | `atmosphere` | Bool | n/a | yes | Field `atmosphere`. |
| in | `solar` | Bool | n/a | yes | Field `solar`. |
| in | `third_body_names` | TB | n/a | yes | Field `third_body_names`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EffectorEnvironmentRequirements | n/a | — | Constructed `EffectorEnvironmentRequirements`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:104-104`
- [[core.effector_sampling_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:151-151`
- [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callees` → `callers` · feedback · `src/dynamics/coupled/aerodynamic_wrench_models.jl:224-224`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1661-1661`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2261-2261`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2068-2068`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1923-1923`
- [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callees` → `callers` · feedback · `src/dynamics/coupled/perturbations.jl:972-972`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · feedback · `src/environment/gravity/gravity_models.jl:257-257`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `callers` · feedback · `src/dynamics/coupled/aerodynamic_wrench_models.jl:224-224`
- [[grp.src_environment_gravity|environment/gravity/]] · `members_out` → `callers` · feedback · `src/environment/gravity/gravity_models.jl:257-257`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · feedback · `src/core/types/effector_sampling.jl:145-145`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/core/types/effector_sampling.jl:145-145`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/core/types/effector_sampling.jl:145-145`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/core/types/effector_sampling.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
The type parameter `TB` makes the struct's type depend on the exact tuple of body names, so merging requirements from heterogeneous effectors yields a new concrete type and cannot be stored in a homogeneous typed collection without abstraction. The constructor accepts any `Tuple`, so a tuple of non-strings passes the keyword layer and only fails at the inner constructor's type bound. Requirements are static per model; an effector cannot request atmosphere only below a given altitude.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 130.
