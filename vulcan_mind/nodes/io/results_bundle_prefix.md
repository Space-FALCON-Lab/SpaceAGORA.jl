---
id: io.results_bundle_prefix
label: _results_bundle_prefix
kind: function
source:
  file: src/io/config/io_config.jl
  symbol: _results_bundle_prefix
  lines:
  - 6
  - 8
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: IOConfig namespace supplying result-directory configuration and path
    conventions.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: prefix
  type: String
  units: n/a
  description: Path prefix used by simulation output writers for the current results
    directory.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
- paths
charts:
- io
origin: agent
---

# _results_bundle_prefix

## Purpose
`_results_bundle_prefix` derives the stable base path for a simulation result bundle from the active simulation settings. Output writers use this helper so CSV, metadata, and related artifacts agree on the configured result directory and common filename stem.

## Theory & Math
The function is a deterministic path mapping: `prefix = joinpath(results_directory, "simulation_results")`. It does not open a file, create a directory, or inspect existing results. The mapping is intentionally independent of timestamps so a single run can refer to a coherent bundle.

## Model & Assumptions
The input object must expose `simulation_settings.results_directory`, and that value must be a valid path for the host filesystem. Callers decide whether a relative path is resolved against the current working directory or normalized earlier in configuration construction.

## Design & Implementation
`io_config.jl` defines the helper as an inline function and exports it for the output module. `io_outputs.jl` calls `IOConfig._results_bundle_prefix(args)` before writing simulation artifacts. The sibling `_results_csv_path` and `_checkpoint_paths` functions derive specialized paths from the same configuration object.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | IOConfig namespace supplying result-directory configuration and path conventions. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `prefix` | String | n/a | — | Path prefix used by simulation output writers for the current results directory. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:127-127`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The helper does not guarantee that the directory exists or is writable and does not prevent two runs from using the same stable bundle prefix. Concurrent writers should use the collision-safe CSV helper or an external run-isolation policy. Missing settings fields fail at access rather than returning a sentinel path.

## Provenance
Mapped from `src/io/config/io_config.jl:6-8`.
