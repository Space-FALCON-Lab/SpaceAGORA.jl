---
id: simulation.persistence__collision_results_csv_path
label: _collision_results_csv_path
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _collision_results_csv_path
  lines:
  - 3
  - 3
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
  type: SimulationModel.IOConfig._collision_results_csv_path
  units: n/a
  description: Return value of `_collision_results_csv_path`. Returns `SimulationModel.IOConfig._collision_results_csv_path(args)`.
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

# _collision_results_csv_path

## Purpose
Generates a collision-resistant alternative CSV path, used when the canonical `simulation_results.csv` cannot be written because another process holds or owns it.

## Design & Implementation
Forwards to `SimulationModel.IOConfig._collision_results_csv_path(args)`, which formats the current UTC instant with `dateformat"yyyymmddTHHMMSSsss"`, concatenates it with `getpid()` and a `rand(UInt)` draw into a token, and returns `joinpath(results_directory, "simulation_results.$token.csv")`. The three-part token combines millisecond time, process identity and a random draw so that same-millisecond writers within one process still differ.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOConfig._collision_results_csv_path | n/a | — | Return value of `_collision_results_csv_path`. Returns `SimulationModel.IOConfig._collision_results_csv_path(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:111-111`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`

**Downstream**

- `callees` → [[io.io_serialization__atomic_write_file|_atomic_write_file]] · `callers` · call · `src/simulation/engine/persistence.jl:5-5`
<!-- vulcan:connections:end -->

## Limitations
Uniqueness is probabilistic rather than guaranteed: the random component comes from the default task-local RNG, so a run that seeds that RNG for reproducibility can reproduce the same token across two same-millisecond calls in the same process. The function does not test whether the produced path already exists, and the timestamped files accumulate in the results directory with no cleanup, so repeated fallbacks leave orphaned outputs behind.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 3.
