---
id: simulation.run_simulation
label: run_simulation
kind: function
source:
  file: src/simulation/engine/public_api.jl
  symbol: run_simulation
  lines:
  - 19
  - 32
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: config
  type: SimulationEngineConfig
  units: n/a
  required: true
  description: Validated single-run engine configuration.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: SimulationResult
  units: n/a
  description: Solver result and associated simulation metadata returned to callers
    or campaign sampling.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
- public-api
charts:
- simulation
origin: agent
---

# run_simulation

## Purpose
`run_simulation` is the public single-scenario execution entrypoint. It turns `SimulationEngineConfig` into a configured ODE problem, installs callbacks and persistence policies, advances the solver, and returns the result consumed by examples, verification, and Monte Carlo campaign code.

## Theory & Math
The engine integrates the configured state derivative `du/dt = f(u,t,p)` from the initial time to the terminal time. Error control is delegated to the selected solver and its absolute and relative tolerances. Callbacks can mutate or observe state at events such as burns, checkpoints, termination conditions, or callback-defined mission boundaries.

## Model & Assumptions
The configuration must provide a consistent state and parameter object, a valid RHS, and compatible callback and output settings. Dynamics and environment models are expected to obey the frame and unit conventions established by `SimulationModel`. The public entrypoint assumes native SPICE/GRAM access follows `RuntimeServices` locking rules.

## Design & Implementation
The function in `public_api.jl` is intentionally small: it delegates construction to engine configuration and problem-building helpers, invokes the integrator, and packages the result. `dynamics_rhs.jl` supplies the derivative path, while persistence and callback modules observe the run. `run_monte_carlo` can use this entrypoint as its sample function.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `config` | SimulationEngineConfig | n/a | yes | Validated single-run engine configuration. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationResult | n/a | — | Solver result and associated simulation metadata returned to callers or campaign sampling. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support_run_and_report|run_and_report]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:179-179`
- [[analysis.runner__run_once|_run_once]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:47-47`
- [[misc.precompile_workload_run_spaceagora_precompile_workload|_run_spaceagora_precompile_workload]] · `callees` → `callers` · call · `src/precompile_workload.jl:47-47`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:160-160`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:492-492`

**Downstream**

- `callees` → [[simulation.from_env__with_engine_env_overrides|_with_engine_env_overrides]] · `callers` · call · `src/simulation/engine/public_api.jl:20-20`
- `callees` → [[simulation.public_api__depwarn_untyped_run_simulation|_depwarn_untyped_run_simulation]] · `callers` · call · `src/simulation/engine/public_api.jl:24-24`
- `callees` → [[simulation.public_api__require_simulation_configuration|_require_simulation_configuration]] · `callers` · call · `src/simulation/engine/public_api.jl:25-25`
<!-- vulcan:connections:end -->

## Limitations
Solver failure, callback exceptions, invalid initial state, and native-library errors propagate through the run result or exception path. Output files may be partially written when a run stops after a checkpoint or callback side effect. Determinism depends on solver settings, callback order, and any randomized model or campaign input.

## Provenance
Mapped from `src/simulation/engine/public_api.jl:18-32`.
