---
id: dynamics.perturbations_gravity_backbone_kick_acceleration_ii
label: gravity_backbone_kick_acceleration_ii
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: gravity_backbone_kick_acceleration_ii
  lines:
  - 1005
  - 1005
inputs:
- id: model
  type: NBodyGravityModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: StateSample
  units: n/a
  required: true
  description: Positional argument `x`.
- id: env
  type: EnvironmentSample
  units: n/a
  required: true
  description: Positional argument `env`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `gravity_backbone_kick_acceleration_ii`.
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

# gravity_backbone_kick_acceleration_ii

## Purpose
Exposes the N-body acceleration as an explicit velocity kick for the split solver, which treats third-body forces as a perturbation to the static backbone.

## Design & Implementation
Requires `env.third_bodies` and calls `_nbody_acceleration_ii`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | NBodyGravityModel | n/a | yes | Positional argument `model`. |
| in | `x` | StateSample | n/a | yes | Positional argument `x`. |
| in | `env` | EnvironmentSample | n/a | yes | Positional argument `env`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `gravity_backbone_kick_acceleration_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:221-221`
- [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1387-1387`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1519-1519`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:342-342`

**Downstream**

- `callees` → [[dynamics.perturbations__nbody_acceleration_ii|_nbody_acceleration_ii]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1013-1013`
<!-- vulcan:connections:end -->

## Limitations
Raises if the environment sample lacks third bodies, which happens when the effector's requirements were not honoured.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1005.
