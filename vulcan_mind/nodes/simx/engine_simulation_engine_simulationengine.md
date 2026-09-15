---
id: simx.engine_simulation_engine_simulationengine
label: SimulationEngine
kind: struct
source:
  file: src/simulation/engine/simulation_engine.jl
  symbol: SimulationEngine
  lines:
  - 2
  - 42
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
- id: simulation_model
  type: Module
  units: n/a
  required: true
  description: SimulationModel namespace supplying the configuration types, effector
    interfaces and IO submodules the engine builds on.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: engine_api
  type: Module
  units: n/a
  description: Namespace exporting the four configuration structs, SimulationEngineConfig,
    simulation_engine_config_from_env, run_simulation, SolverIntegratorCache and the
    two n-body ephemeris cache entry points.
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
# SimulationEngine

## Purpose
`SimulationEngine` is the canonical aggregator for the propagation engine. Its own comment states the rule it follows: no behaviour ownership. The file declares the module, brings in the dependencies the included files assume, includes twelve implementation files in dependency order, and exports nine names.

## Model & Assumptions
Include order is the module's real contract, because Julia evaluates includes into a single flat module scope. Configuration comes first so the types exist, then the environment adapter that constructs them, then sampling and state accessors, then setup, then the solver policy, then the right-hand side, calibration, persistence, checkpointing, reporting, execution and finally the public API. `execution.jl` sits near the end because `run_simulation` calls into nearly everything above it.

## Design & Implementation
Dependencies are declared once at the top: `using ..SimulationModel` for the model layer, `import ..RuntimeServices` so the shared SPICE lock is reachable by qualified name, `import DiffEqBase` and `import ADTypes: AutoFiniteDiff` for the solver interface, `using SPICE` for the raw ephemeris calls and `using SparseArrays` for the block-diagonal Jacobian prototype. The five configuration includes cover parallel, solver, runtime policy, artifact and the aggregate `simulation_engine_config.jl`; note that the solver slot is a comment-only file because `SolverConfig` is re-exported from `SimulationModel`. The export list is deliberately narrow: everything prefixed with an underscore stays internal.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `simulation_model` | Module | n/a | yes | SimulationModel namespace supplying the configuration types, effector interfaces and IO submodules the engine builds on. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `engine_api` | Module | n/a | — | Namespace exporting the four configuration structs, SimulationEngineConfig, simulation_engine_config_from_env, run_simulation, SolverIntegratorCache and the two n-body ephemeris cache entry points. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/simulation_engine.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Flat include scope means every internal helper across twelve files shares one namespace, so a duplicated private name in two files silently redefines a method rather than colliding. The aggregator carries no conditional includes, so extension-only code paths must be handled inside the included files. Reordering the includes is a breaking change that surfaces as an undefined-symbol error at precompile time.

## Provenance
Mapped from `src/simulation/engine/simulation_engine.jl:2-42`; the twelve included files live under `src/simulation/engine/` and `src/simulation/engine/config/`.
