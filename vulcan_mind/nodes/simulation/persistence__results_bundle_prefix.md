---
id: simulation.persistence__results_bundle_prefix
label: _results_bundle_prefix
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _results_bundle_prefix
  lines:
  - 1
  - 1
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
  type: SimulationModel.IOConfig._results_bundle_prefix
  units: n/a
  description: Return value of `_results_bundle_prefix`. Returns `SimulationModel.IOConfig._results_bundle_prefix(args)`.
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

# _results_bundle_prefix

## Purpose
Computes the common path prefix used to name every file in a results bundle, so that the CSV, Arrow/Feather data file and TOML manifest of one run share a stem.

## Design & Implementation
The engine-level definition is an `@inline` forwarder to `SimulationModel.IOConfig._results_bundle_prefix(args)`, keeping the simulation engine free of direct filesystem policy. The implementation returns `joinpath(args.simulation_settings.results_directory, "simulation_results")`, so the stem is fixed and only the configured results directory varies between runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOConfig._results_bundle_prefix | n/a | — | Return value of `_results_bundle_prefix`. Returns `SimulationModel.IOConfig._results_bundle_prefix(args)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:127-127`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The stem is a hard-coded literal, so two concurrent runs pointed at the same `results_directory` collide and the later one overwrites the earlier bundle. The function neither creates nor validates the directory, so a nonexistent or unwritable `results_directory` is only discovered when a writer later opens the path. No extension is appended; each caller must supply its own suffix consistently.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 1.
