---
id: parcore.abstract_types_abstracttypes
label: AbstractTypes
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractTypes
  lines:
  - 1
  - 67
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Loaded first in the SimulationModel submodule block; depends on nothing
    but Base.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: types
  type: Module
  units: n/a
  description: The eight abstract supertypes that anchor the package's dispatch hierarchy
    for force models, planets, atmospheres, thermal models, thrusters, effectors,
    ephemerides and guidance.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# AbstractTypes

## Purpose
`AbstractTypes` declares the root of the package's type hierarchy. Every physical model in the simulation subtypes one of the eight abstract types defined here, which is what lets configuration records constrain their type parameters and lets model evaluation dispatch on the concrete subtype without a runtime branch.

## Model & Assumptions
The hierarchy is deliberately flat: `AbstractForceTorqueModel`, `AbstractControlEffectorModel`, `AbstractPlanet`, `AbstractDensityModel`, `AbstractThermalModel`, `AbstractEphemeridesModel`, `AbstractThrusterModel` and `AbstractGuidanceModel` are all direct children of `Any` with no shared intermediate supertype. Each therefore expresses a capability contract enforced by convention and by the methods other modules define on it, not by any structural constraint in this file.

## Design & Implementation
The module contains one export line and eight `abstract type` declarations, each preceded by a docstring naming the role and the expected interface. It is the first file included in the submodule block of `simulation_model.jl`, ahead of `effector_sampling.jl`, because the sampling contracts and the configuration records both annotate against these names. Keeping the declarations in a dedicated module with no dependencies means the whole package can be type-annotated without creating an include cycle.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Loaded first in the SimulationModel submodule block; depends on nothing but Base. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `types` | Module | n/a | — | The eight abstract supertypes that anchor the package's dispatch hierarchy for force models, planets, atmospheres, thermal models, thrusters, effectors, ephemerides and guidance. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/abstract_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the types are empty and unrelated, nothing prevents a model from subtyping the wrong root, and the mistake surfaces only when a configuration field rejects it. The interface each abstract type implies is documented in prose rather than encoded, so a subtype that omits a required method compiles fine and raises a MethodError during propagation. Adding a ninth root here forces a recompile of essentially the entire package.

## Provenance
Mapped from `src/core/types/abstract_types.jl:1-67`.
