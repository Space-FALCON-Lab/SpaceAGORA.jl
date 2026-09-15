---
id: core.effector_sampling_effectorsampling
label: EffectorSampling
kind: module
source:
  file: src/core/types/effector_sampling.jl
  symbol: EffectorSampling
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# EffectorSampling

## Purpose
Small leaf module (`module EffectorSampling`) that defines the typed sample structs and the additive extension hooks through which dynamic effectors (aerodynamics, SRP, N-body, thrusters) contribute forces and torques to the simulation without touching integrator internals. It depends only on `StaticArrays` so that effector packages can extend it without pulling in the engine.

## Design & Implementation
Exports the sample types `StateSample`, `PlanetFrameSample`, `AtmosphereSample`, `SolarEphemerisSample`, `ThirdBodyEphemerisSample`, `EnvironmentSample`, the capability request `EffectorEnvironmentRequirements`, and the hook generics `wrench`, `wrench_caching!`, `environment_requirements`, `solver_partition`, `gravity_backbone_structure`, `gravity_backbone_acceleration_ii`, `gravity_backbone_kick_structure`, `gravity_backbone_kick_acceleration_ii`. `wrench`, `gravity_backbone_acceleration_ii` and `gravity_backbone_kick_acceleration_ii` are declared with `function ... end` and have no default method, so an effector that opts in must define them; the declaration hooks (`environment_requirements`, `solver_partition`, the two `_structure` functions) have `@inline` fallbacks on `::Any` that opt out. Every sample field is an `SVector`/`SMatrix` of `Float64` for allocation-free evaluation in the ODE right-hand side.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The module has no runtime validation: an effector that declares `:position_only_static_gravity` but does not implement `gravity_backbone_acceleration_ii` fails only at first call with a `MethodError`. Third-body names are `String`s in an `NTuple`, so each distinct set is its own type and can trigger compilation per constellation configuration. Units are fixed to SI by convention in docstrings only; nothing enforces metres, kilograms or radians.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 1.
