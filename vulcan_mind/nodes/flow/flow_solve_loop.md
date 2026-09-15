---
id: flow.solve_loop
label: Solve loop
kind: group
inputs:
- id: initial_state
  type: ComponentVector
  units: n/a
  description: Initial conditions.
- id: ode_params
  type: ODEParams
  units: n/a
  description: Prepared parameters.
- id: sample_config
  type: SimulationConfiguration
  units: n/a
  description: One campaign sample's configuration.
  required: false
- id: scenario_config
  type: SimulationConfiguration
  units: n/a
  description: One verification scenario's configuration.
  required: false
outputs:
- id: callback_hooks
  type: CallbackSet
  units: n/a
  description: Per-step and event callbacks the integrator fires.
- id: rhs_calls
  type: ODEFunction
  units: n/a
  description: The derivative the integrator evaluates.
- id: saved_values
  type: SavedValues / DataFrame
  units: n/a
  description: Time series collected by the save callback.
- id: solution
  type: ODESolution + solver metadata
  units: n/a
  description: The in-memory result when requested.
tags:
- master-flow
charts:
- master
origin: agent
opens: run-pipeline
---

# Solve loop

## Purpose
The heart of a run: `run_simulation` builds the ODE problem from the initial state and parameters, chooses a solver and tolerances under the configured policy, integrates the mission — in checkpointed segments when enabled — and hands the collected time series to the writers.

## Design & Implementation
`execution.jl` validates the configuration, runs setup, assembles the `CallbackSet`, resolves the solver mode (default, split IMEX, multirate or gravity-backbone split), builds tolerances and an optional block-diagonal Jacobian prototype, and calls `_solve_with_solver_policy`. With checkpointing the mission is solved in segments and a checkpoint written after each; on failure partial results are still saved. The solution or solver trace is returned when the caller asks for it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_state` | ComponentVector | n/a | — | Initial conditions. |
| in | `ode_params` | ODEParams | n/a | — | Prepared parameters. |
| in | `sample_config` | SimulationConfiguration | n/a | no | One campaign sample's configuration. |
| in | `scenario_config` | SimulationConfiguration | n/a | no | One verification scenario's configuration. |
| out | `callback_hooks` | CallbackSet | n/a | — | Per-step and event callbacks the integrator fires. |
| out | `rhs_calls` | ODEFunction | n/a | — | The derivative the integrator evaluates. |
| out | `saved_values` | SavedValues / DataFrame | n/a | — | Time series collected by the save callback. |
| out | `solution` | ODESolution + solver metadata | n/a | — | The in-memory result when requested. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.campaigns|Campaigns]] · `sample_config` → `sample_config` · dataflow · `src/simulation/campaigns/monte_carlo.jl`
- [[flow.setup_run|Set up the run]] · `initial_state` → `initial_state` · dataflow · `src/simulation/engine/execution.jl`
- [[flow.setup_run|Set up the run]] · `ode_params` → `ode_params` · dataflow · `src/simulation/engine/execution.jl`
- [[flow.verification|Telemetry verification study]] · `scenario_config` → `scenario_config` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callback_hooks` → [[flow.callbacks|Integration callbacks]] · `callback_hooks` · dataflow · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `rhs_calls` → [[flow.rhs|Dynamics right-hand side]] · `rhs_calls` · dataflow · `src/simulation/engine/dynamics_rhs.jl`
- `saved_values` → [[flow.write_results|Write results & checkpoints]] · `saved_values` · dataflow · `src/io/outputs/io_outputs.jl`
- `solution` → [[output.solution_object|Solution / solver metadata (in memory)]] · `solution` · dataflow · `src/simulation/engine/execution.jl`
<!-- vulcan:connections:end -->

## Limitations
A solve holds the whole `SavedValues` history in memory until it ends, and the checkpoint segment loop re-creates the ODE problem per segment; only the `gravity_backbone_split` mode can resume from a backbone checkpoint.
