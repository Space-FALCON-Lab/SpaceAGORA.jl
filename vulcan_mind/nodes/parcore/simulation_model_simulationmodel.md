---
id: parcore.simulation_model_simulationmodel
label: SimulationModel
kind: struct
source:
  file: src/core/simulation_model.jl
  symbol: SimulationModel
  lines:
  - 2
  - 141
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Physics, vehicle, GNC and IO submodule namespaces pulled in by the
    include and @reexport sequence.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: namespace
  type: Module
  units: n/a
  description: Single flattened namespace re-exporting spacecraft, environment, effector,
    sampling and configuration symbols to the rest of the package.
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

# SimulationModel

## Purpose
`SimulationModel` is the canonical aggregator module. Its own comment states that it owns no behavior: it exists to define the load order of the simulation's constituent modules and to re-export their public APIs under one name so that campaign, CLI and analysis code depends on a single namespace instead of thirty include paths.

## Model & Assumptions
The module assumes each included file defines its own module and that the include order respects the dependency order, because Julia resolves `using .Submodule` at include time. It also assumes symbol names across submodules do not collide, since `@reexport` flattens all of them into one namespace; a collision surfaces as an ambiguity error at load, not at call time.

## Design & Implementation
Loading starts with LinearAlgebra, StaticArrays, CSV, DataFrames and Reexport, then conditionally includes `runtime_services.jl` into the parent module using an `isdefined` guard so repeated loads do not redefine it. Utility files for quaternions, reference-system records and planet shapes are included without re-export. The submodule block then interleaves `include` and `@reexport using` in dependency order: abstract types and effector sampling first, then command types, robotics and arm planning, then cloth multibody and rotational and translational dynamics, then environment, IO, mission and vehicle layers, closing with the navigation, guidance and control hook modules.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Physics, vehicle, GNC and IO submodule namespaces pulled in by the include and @reexport sequence. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `namespace` | Module | n/a | — | Single flattened namespace re-exporting spacecraft, environment, effector, sampling and configuration symbols to the rest of the package. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/simulation_model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because every symbol is flattened, the module gives no import-level indication of which submodule owns a name, which makes ownership questions a grep exercise. The include-order coupling is implicit: moving one include earlier can break a downstream type annotation with an obscure UndefVarError. Load time is dominated by this file, since it forces compilation of the whole physics stack even for callers that need only one subsystem.

## Provenance
Mapped from `src/core/simulation_model.jl:2-141`.
