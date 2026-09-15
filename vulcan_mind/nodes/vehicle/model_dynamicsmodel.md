---
id: vehicle.model_dynamicsmodel
label: DynamicsModel
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: DynamicsModel
  lines:
  - 417
  - 417
inputs:
- id: spacecraft
  type: Vector{SpacecraftModel}
  units: n/a
  required: true
  description: Field `spacecraft`.
- id: dynamic_effectors
  type: T_Effectors
  units: n/a
  required: true
  description: Field `dynamic_effectors`.
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
  type: DynamicsModel
  units: n/a
  description: Constructed `DynamicsModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# DynamicsModel

## Purpose
Top-level container pairing the list of spacecraft with the tuple of dynamic effectors (gravity, drag, SRP, third-body, thrusters) that the RHS applies to each of them.

## Design & Implementation
`struct DynamicsModel{T_Effectors<:Tuple}` with `spacecraft::Vector{SpacecraftModel}` and `dynamic_effectors::T_Effectors`. The single inner constructor `DynamicsModel(roots::Vector{SpacecraftModel}, dynamic_effectors::T_Effectors)` stores them without validation. The docstring notes that the spacecraft vector is the axis over which threads are parallelised, each thread handling one root and its sub-links. The tuple type parameter lets `calcForceTorque` dispatch statically per effector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Vector{SpacecraftModel} | n/a | yes | Field `spacecraft`. |
| in | `dynamic_effectors` | T_Effectors | n/a | yes | Field `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DynamicsModel | n/a | — | Constructed `DynamicsModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:160-160`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[simulation.constellation_ensemble__ensemble_member_configuration|_ensemble_member_configuration]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:42-42`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No check that `spacecraft` is non-empty or that `id` fields are unique, and no check that every effector's `environment_requirements` are satisfiable. Because `dynamic_effectors` is a `Tuple`, adding or removing an effector changes the type and forces recompilation of the RHS.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 417.
