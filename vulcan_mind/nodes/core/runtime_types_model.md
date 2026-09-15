---
id: core.runtime_types_model
label: Model
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Model
  lines:
  - 158
  - 158
inputs:
- id: body
  type: SpacecraftModel
  units: n/a
  required: false
  description: Field `body` (default `SpacecraftModel()`).
- id: aerodynamics
  type: Aerodynamics
  units: n/a
  required: false
  description: Field `aerodynamics` (default `Aerodynamics()`).
- id: engines
  type: Engines
  units: n/a
  required: false
  description: Field `engines` (default `Engines()`).
- id: initial_condition
  type: Initial_condition
  units: n/a
  required: false
  description: Field `initial_condition` (default `Initial_condition()`).
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
  type: Model
  units: n/a
  description: Constructed `Model` (keyword constructor via @kwdef).
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

# Model

## Purpose
Legacy top-level bundle tying a `SpacecraftModel` together with scalar aerodynamic, engine, and initial-condition records for the Python-port simulation entry points.

## Design & Implementation
`@kwdef struct Model` with four fields, each defaulting to its type's default constructor: `body::SpacecraftModel`, `aerodynamics::Aerodynamics`, `engines::Engines`, and `initial_condition::Initial_condition`. Being immutable, replacing any component requires reconstructing the `Model`, although `Aerodynamics` and `Engines` are themselves mutable and can be edited in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | SpacecraftModel | n/a | no | Field `body` (default `SpacecraftModel()`). |
| in | `aerodynamics` | Aerodynamics | n/a | no | Field `aerodynamics` (default `Aerodynamics()`). |
| in | `engines` | Engines | n/a | no | Field `engines` (default `Engines()`). |
| in | `initial_condition` | Initial_condition | n/a | no | Field `initial_condition` (default `Initial_condition()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Model | n/a | — | Constructed `Model` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:105-105`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- `callees` → [[core.runtime_types_aerodynamics|Aerodynamics]] · `callers` · call · `src/core/types/runtime_types.jl:160-160`
- `callees` → [[core.runtime_types_engines|Engines]] · `callers` · call · `src/core/types/runtime_types.jl:161-161`
- `callees` → [[core.runtime_types_initial_condition|Initial_condition]] · `callers` · call · `src/core/types/runtime_types.jl:162-162`
- `callees` → [[vehx.spacecraft_model_spacecraftmodel|SpacecraftModel]] · `callers` · call · `src/core/types/runtime_types.jl:159-159`
<!-- vulcan:connections:end -->

## Limitations
Constructing a default `Model()` builds a full `SpacecraftModel()` which may be expensive. The mix of immutable outer and mutable inner records makes aliasing easy: two `Model` values built from the same `Aerodynamics` instance share mutable state.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 158.
