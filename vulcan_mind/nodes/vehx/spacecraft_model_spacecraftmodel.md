---
id: vehx.spacecraft_model_spacecraftmodel
label: SpacecraftModel
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: SpacecraftModel
  lines:
  - 374
  - 386
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: initial_condition
  type: AbstractInitialCondition
  units: n/a
  required: true
  description: Orbit and attitude initialisation, either Keplerian or Cartesian.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model
  type: SpacecraftModel
  units: n/a
  description: Assembled vehicle definition consumed by dynamics, guidance, navigation
    and control effectors.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- spacecraft
charts:
- vehx
origin: agent
---

# SpacecraftModel

## Purpose
`SpacecraftModel` is the central vehicle object of the simulator. It aggregates the multi-body topology, the mass and inertia budget, the actuator counts and the initial condition into one mutable record that is handed to the dynamics right-hand side, the structure property helpers, the environment samplers and the plotting layer. Every other type in the vehicle package is either a component of it or a function over it.

## Model & Assumptions
The vehicle is a rooted collection of rigid links joined by joints. Mass is split into a dry part and a propellant part so that a burn can deplete the latter while the former stays fixed. The inertia tensor is stored as a single static three-by-three matrix in body axes, representing the composite assembly about the vehicle centre of mass. Reaction wheel and thruster counts are cached so that state vector layout can be sized without walking the link list. The `instant_actuation` flag selects whether commanded articulation such as a solar panel angle is applied immediately or driven through actuator dynamics.

## Design & Implementation
The struct is `mutable` because assembly proceeds incrementally: `add_body!` and `add_joint!` push into the link and joint vectors and update the cached counts after the object already exists. Fields are concretely typed, with `SMatrix{3,3,Float64}` for the inertia tensor and `Int64` counters, so the compiler can specialise the dynamics loop. `initial_condition` is typed on the abstract supertype, which lets a Keplerian `InitialCondition` and a `CartesianInitialCondition` be used interchangeably at the cost of a dynamic dispatch at setup time. The `id` field distinguishes vehicles in a constellation run. A keyword constructor immediately below the struct supplies defaults for every field.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `initial_condition` | AbstractInitialCondition | n/a | yes | Orbit and attitude initialisation, either Keplerian or Cartesian. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model` | SpacecraftModel | n/a | — | Assembled vehicle definition consumed by dynamics, guidance, navigation and control effectors. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__link_q|_link_q]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:104-104`
- [[core.runtime_types_model|Model]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:159-159`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only one root is representable, so a genuinely disconnected multi-assembly vehicle cannot be held in a single model even though several structure helpers accept a root index. `links` and `joints` are abstractly parameterised vectors, which costs type stability when links of different reaction wheel counts are mixed. The cached counters can drift from the actual hardware if a link is mutated after registration, and the inertia tensor is not recomputed automatically when mass properties change.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl:374-386`, within the `SpacecraftModels` module that also defines `Link`, `Joint`, the initial condition types and the guidance, navigation and control model wrappers.
