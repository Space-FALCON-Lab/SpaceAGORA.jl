---
id: spaceagora.run_simulation
label: run_simulation
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: run_simulation
  lines:
  - 516
  - 517
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORA public namespace forwarding to SimulationEngine.run_simulation.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: SimulationEngine.run_simulation
  units: n/a
  description: Return value of `run_simulation`. Returns `SimulationEngine.run_simulation(args...;
    kwargs...)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- public-api
- simulation
charts:
- spaceagora
origin: agent
---

# run_simulation

## Purpose
The top-level `run_simulation` methods in `src/SpaceAGORA.jl` forward public package calls into `SimulationEngine.run_simulation`. They keep the package namespace ergonomic while preserving the engine’s configuration and keyword arguments. Two methods cover the ordinary argument path and the explicit `SimulationEngineConfig` path.

## Theory & Math
The wrapper does not integrate a state itself. It preserves the engine contract `du/dt = f(u,t,p)` and returns the engine’s solver result unchanged. Any solver tolerances, callback semantics, and output policies are therefore defined by the downstream `SimulationEngine` implementation rather than duplicated in this public namespace.

## Model & Assumptions
Callers must supply arguments accepted by the engine overload and must keep the package’s state, frame, and unit conventions. Keyword arguments such as `isolate_state` and `return_solution` are forwarded to the engine. The wrapper assumes `SimulationEngine` has been included by the package load order.

## Design & Implementation
The first method at line 516 forwards variadic positional and keyword arguments; the second immediately below forwards an explicit `SimulationEngineConfig`. This is a deliberate thin boundary: `src/SpaceAGORA.jl` owns exports and documentation, while `src/simulation/engine/public_api.jl` owns problem construction, integration, callbacks, and persistence.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | yes | SpaceAGORA public namespace forwarding to SimulationEngine.run_simulation. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationEngine.run_simulation | n/a | — | Return value of `run_simulation`. Returns `SimulationEngine.run_simulation(args...; kwargs...)`. |
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

- *none*
<!-- vulcan:connections:end -->

## Limitations
Errors from configuration, solver setup, callbacks, persistence, or native environment calls propagate from the engine. The wrapper adds no validation and no retry behavior. Because it forwards keywords, a change in the engine public API can affect package callers even when this file remains unchanged.

## Provenance
Mapped from `src/SpaceAGORA.jl:516-517`.
