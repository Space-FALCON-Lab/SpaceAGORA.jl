---
id: environment.gravity_models_inversesquaredgravitymodel
label: InverseSquaredGravityModel
kind: struct
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: InverseSquaredGravityModel
  lines:
  - 13
  - 13
inputs:
- id: gravity_gradient
  type: Bool
  units: n/a
  required: false
  description: Field `gravity_gradient` (default `false`).
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
  type: InverseSquaredGravityModel
  units: n/a
  description: Constructed `InverseSquaredGravityModel` (keyword constructor via @kwdef).
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

# InverseSquaredGravityModel

## Purpose
Effector model type for point-mass (two-body) gravity with optional gravity-gradient torque. It is the canonical central-field model used when no oblateness or harmonics are needed, and it participates in the gravity-backbone fast path.

## Design & Implementation
`@kwdef struct InverseSquaredGravityModel <: AbstractForceTorqueModel` with one field `gravity_gradient::Bool = false`. The force is `mass * _inverse_squared_gravity_accel(pos_ii, planet)` with `planet.μ` in m^3/s^2 taken from the environment; torque comes from `_gravity_gradient_torque_body`. It declares `gravity_backbone_structure(::InverseSquaredGravityModel) = :position_only_static_gravity` and implements `gravity_backbone_acceleration_ii` returning the same inverse-square acceleration, and provides both `calcForceTorque` and `wrench` methods.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gravity_gradient` | Bool | n/a | no | Field `gravity_gradient` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | InverseSquaredGravityModel | n/a | — | Constructed `InverseSquaredGravityModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__base_gravity_effector|_base_gravity_effector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:22-22`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:27-27`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No `μ` override field exists, so per-model gravitational parameters (for example a modified value for testing) are impossible without changing the planet. The model does not include the indirect third-body term or any relativistic correction. The commented-out defaults in the source reference Earth values that are not used, which can mislead readers into thinking the model is Earth-specific.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 13.
