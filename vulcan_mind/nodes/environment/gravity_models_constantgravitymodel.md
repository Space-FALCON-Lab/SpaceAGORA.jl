---
id: environment.gravity_models_constantgravitymodel
label: ConstantGravityModel
kind: struct
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: ConstantGravityModel
  lines:
  - 8
  - 8
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
  type: ConstantGravityModel
  units: n/a
  description: Constructed `ConstantGravityModel` (keyword constructor via @kwdef).
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

# ConstantGravityModel

## Purpose
Effector model type selecting central inverse-square gravity for a spacecraft. Despite the name, it evaluates `-μ/r^2` at the current position rather than a constant field; the name survives from the legacy configuration vocabulary. It is one of three gravity model types in this file and is a subtype of `AbstractForceTorqueModel`.

## Design & Implementation
Declared with `@kwdef struct ConstantGravityModel <: AbstractForceTorqueModel` and a single field `gravity_gradient::Bool = false` that enables the body-frame gravity-gradient torque. The gravitational parameter is not stored here (the `μ` field is commented out) but read from `planet.μ` at evaluation time. Methods dispatched on this type in the same file are `calcForceTorque` (legacy `ComponentVector` path), `wrench` (sampled `StateSample` path), `gravity_backbone_structure` returning `:position_only_static_gravity`, and `gravity_backbone_acceleration_ii`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gravity_gradient` | Bool | n/a | no | Field `gravity_gradient` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstantGravityModel | n/a | — | Constructed `ConstantGravityModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`ConstantGravityModel` and `InverseSquaredGravityModel` are functionally identical, differing only in type; configuration files that expect a truly constant acceleration will get an inverse-square field. Because `μ` comes from the planet object, the model cannot be used with a body whose `μ` differs from the environment planet. `environment_requirements` is not overridden, so no planet-frame sampling is requested; the gravity-gradient torque therefore depends on `x.q_ib` and `x.spacecraft` being populated by the sampler.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 8.
