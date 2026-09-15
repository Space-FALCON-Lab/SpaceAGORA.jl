---
id: simx.engine_config_solver_config_solverconfig
label: SolverConfig
kind: struct
source:
  file: src/simulation/engine/config/solver_config.jl
  symbol: SolverConfig
  lines:
  - 1
  - 2
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: solver_config_binding
  type: SolverConfig
  units: n/a
  description: The SolverConfig name as it becomes visible in SimulationEngine scope,
    re-exported from SimulationModel.SimConfig.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# SolverConfig

## Purpose
This file is a placeholder in the engine's configuration include chain. It declares no code; its two lines record that `SolverConfig` is defined in `SimulationModel.SimConfig` and re-exported through `SimulationModel`, so it is already visible in `SimulationEngine` scope via `using ..SimulationModel`.

## Model & Assumptions
The engine module includes five configuration files in a fixed order, and `SolverConfig` occupies the second slot. Keeping the slot occupied by a documented note rather than deleting the include means the configuration chain reads as a complete list of the four configuration types plus the aggregate, and a reader looking for the solver record is told exactly where it lives instead of finding a gap.

## Design & Implementation
The include at `src/simulation/engine/simulation_engine.jl:12` evaluates this file for its comments only. `SolverConfig` itself is a shared type because it is referenced from both the model layer and the engine layer: `SimulationConfiguration` carries an optional `solver_config` field, and `run_simulation` falls back to `_solver_config_from_env()` when that field is nothing. Defining the type in the model package avoids a circular dependency between the type that describes a solve and the engine that performs it. Downstream, the name is exported again from `SimulationEngine` so callers see one solver configuration type regardless of entry point.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `solver_config_binding` | SolverConfig | n/a | — | The SolverConfig name as it becomes visible in SimulationEngine scope, re-exported from SimulationModel.SimConfig. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callees` → `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:10-10`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the file contains no declaration, tooling that locates a type by its defining file will not find `SolverConfig` here; the real definition and its field defaults live in the `SimConfig` submodule of `SimulationModel`. Renaming the type there leaves this note stale with no compile-time signal.

## Provenance
Mapped from `src/simulation/engine/config/solver_config.jl:1-2`, an include-manifest slot; the include site is `src/simulation/engine/simulation_engine.jl:12`.
