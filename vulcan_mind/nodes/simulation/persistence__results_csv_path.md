---
id: simulation.persistence__results_csv_path
label: _results_csv_path
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _results_csv_path
  lines:
  - 2
  - 2
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: SimulationModel.IOConfig._results_csv_path
  units: n/a
  description: Return value of `_results_csv_path`. Returns `SimulationModel.IOConfig._results_csv_path(args)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _results_csv_path

## Purpose
Yields the canonical CSV output path for a simulation run, the file the results DataFrame is written to at the end of propagation.

## Design & Implementation
An `@inline` forwarder to `SimulationModel.IOConfig._results_csv_path(args)`, which returns `joinpath(args.simulation_settings.results_directory, "simulation_results.csv")`. Keeping the name in the IO configuration module means the engine, the bundle writer and the collision-avoidance fallback all agree on one location without passing strings around. It is a pure function of `args` and performs no filesystem access.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOConfig._results_csv_path | n/a | — | Return value of `_results_csv_path`. Returns `SimulationModel.IOConfig._results_csv_path(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:97-97`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The filename is a constant, so repeated runs into the same results directory silently overwrite previous output, and parallel runs sharing a directory race on the same handle, which is precisely why `_collision_results_csv_path` exists as an escape hatch. No check is made that the path is writable or that the parent directory exists.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 2.
