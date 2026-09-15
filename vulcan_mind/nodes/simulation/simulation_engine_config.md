---
id: simulation.simulation_engine_config
label: SimulationEngineConfig
kind: function
source:
  file: src/simulation/engine/config/simulation_engine_config.jl
  symbol: SimulationEngineConfig
  lines:
  - 8
  - 14
outputs:
- id: config
  type: SimulationEngineConfig
  units: n/a
  description: Engine configuration record consumed by run_simulation and callback
    setup.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
- configuration
charts:
- simulation
origin: agent
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
---

# SimulationEngineConfig

## Purpose
`SimulationEngineConfig` is the typed configuration boundary for a single simulation run. It groups solver, timestep, callback, output, and execution choices so `run_simulation` can construct the integration problem without parsing user-facing options inside the stepping loop.

## Theory & Math
The configuration determines the numerical integration of `du/dt = f(u,t,p)`. Solver choice and tolerances control the approximation error, while callback settings determine event handling and persistence. The record itself performs no integration and carries no mutable integrator state.

## Model & Assumptions
The configuration is assumed to contain compatible state dimensions, initial conditions, solver options, and callback plans. Tolerances are interpreted by the chosen DifferentialEquations solver. Output and checkpoint paths are expected to be derived through `IOConfig` so parallel runs do not overwrite one another unexpectedly.

## Design & Implementation
The configuration file declares the record and its defaults. `run_simulation` receives it as a required input, builds the problem and callback set, and delegates stepping to the solver. Runtime services and environment models are referenced through the module aggregation layer rather than copied into the configuration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config` | SimulationEngineConfig | n/a | — | Engine configuration record consumed by run_simulation and callback setup. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:191-191`

**Downstream**

- `callees` → [[core.simulation_configuration_solverconfig|SolverConfig]] · `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:10-10`
- `callees` → [[simx.engine_config_artifact_config_artifactconfig|ArtifactConfig]] · `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:12-12`
- `callees` → [[simx.engine_config_parallel_config_parallelconfig|ParallelConfig]] · `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:9-9`
- `callees` → [[simx.engine_config_runtime_policy_config_runtimepolicyconfig|RuntimePolicyConfig]] · `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:11-11`
- `callees` → [[simx.engine_config_solver_config_solverconfig|SolverConfig]] · `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
The record cannot validate every combination of solver and callback settings before integration. An invalid ODE parameter object, missing output directory, or incompatible state layout fails during problem construction or stepping. Defaults are intended for package examples and should be reviewed for production campaigns.

## Provenance
Mapped from `src/simulation/engine/config/simulation_engine_config.jl:1-14`.
