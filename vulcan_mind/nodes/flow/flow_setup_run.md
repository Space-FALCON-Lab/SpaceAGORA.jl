---
id: flow.setup_run
label: Set up the run
kind: group
inputs:
- id: run_config
  type: SimulationConfiguration
  units: n/a
  description: The configuration to prepare.
outputs:
- id: initial_state
  type: ComponentVector
  units: n/a
  description: Packed per-satellite initial conditions.
- id: ode_params
  type: ODEParams
  units: n/a
  description: Shared buffers, caches and environment snapshots.
- id: cache_file
  type: serialized NBodyEphemerisCache
  units: n/a
  description: Optional prewarmed ephemeris table written to disk.
tags:
- master-flow
charts:
- master
origin: agent
opens: simulation-simulation-engine-setup-jl
---

# Set up the run

## Purpose
Prepares everything a solve needs before the first derivative is evaluated: initial conditions for every satellite, the shared buffers and scratch workspaces, per-satellite density model instances, the run-scoped environment snapshots, and the ephemeris tables that let the solve avoid SPICE on the hot path.

## Design & Implementation
`build_initial_conditions` packs position, velocity, mass, heat loads and optionally attitude into a `ComponentVector`; `ODEParams` allocates `SharedBuffers` sized by satellite count; the `_initialize_*` family in `setup.jl` resets caches and captures `PolicyDecisionEnvConfig`, `RhsPlanEnvConfig` and `CallbackEnvConfig`; and the three ephemeris cache builders tabulate third-body, Sun and planet-orientation data over the mission, reusing tables across runs when the parameters match. A prewarmed N-body table can be serialised for later campaigns.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `run_config` | SimulationConfiguration | n/a | — | The configuration to prepare. |
| out | `initial_state` | ComponentVector | n/a | — | Packed per-satellite initial conditions. |
| out | `ode_params` | ODEParams | n/a | — | Shared buffers, caches and environment snapshots. |
| out | `cache_file` | serialized NBodyEphemerisCache | n/a | — | Optional prewarmed ephemeris table written to disk. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.configure|Configure a run]] · `run_config` → `run_config` · dataflow · `src/simulation/engine/execution.jl`

**Downstream**

- `cache_file` → [[output.ephemeris_cache_file|Prewarmed N-body ephemeris cache]] · `cache_file` · dataflow · `src/simulation/engine/setup.jl`
- `initial_state` → [[flow.solve_loop|Solve loop]] · `initial_state` · dataflow · `src/simulation/engine/execution.jl`
- `ode_params` → [[flow.solve_loop|Solve loop]] · `ode_params` · dataflow · `src/simulation/engine/execution.jl`
<!-- vulcan:connections:end -->

## Limitations
Environment snapshots freeze at setup, so changing a variable mid-run has no effect; ephemeris tables are keyed on exact epoch and duration, so campaigns that vary either do not share them.
