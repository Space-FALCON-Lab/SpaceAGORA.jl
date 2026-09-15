---
id: io.io_config__results_csv_path
label: _results_csv_path
kind: function
source:
  file: src/io/config/io_config.jl
  symbol: _results_csv_path
  lines:
  - 10
  - 10
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
  type: String
  units: n/a
  description: Return value of `_results_csv_path`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# _results_csv_path

## Purpose

Gives the canonical path of the aggregate simulation results CSV for a run: `simulation_results.csv` inside the configured results directory.

## Design & Implementation

An `@inline` one-liner returning `joinpath(args.simulation_settings.results_directory, "simulation_results.csv")` as a `String`. The filename is a literal, so every writer and every downstream reader agrees on it without passing the name around. It shares its stem with `_results_bundle_prefix`, which returns the same directory plus `simulation_results` for callers that append their own extension.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_results_csv_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:97-97`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Because the name is fixed, two simulations sharing a results directory will collide on this file — that is exactly the case `_collision_results_csv_path` exists to handle. The function performs no existence or permission check and does not create the parent directory.

## Provenance
Mapped from `src/io/config/io_config.jl` line 10.
