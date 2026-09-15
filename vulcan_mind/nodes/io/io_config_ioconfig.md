---
id: io.io_config_ioconfig
label: IOConfig
kind: module
source:
  file: src/io/config/io_config.jl
  symbol: IOConfig
  lines:
  - 1
  - 1
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
outputs:
- id: result
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# IOConfig

## Purpose

Small module that centralises the filesystem layout of a simulation run: where the results CSV is written, where a per-process collision-free results file goes, and where checkpoint data and its manifest live. Every path is derived from `args.simulation_settings` so no other module hard-codes an output filename.

## Design & Implementation

The module pulls in `Dates` and `Random` purely for the unique-token construction in `_collision_results_csv_path`, and exports the five path builders `_results_bundle_prefix`, `_results_csv_path`, `_collision_results_csv_path`, `_checkpoint_directory` and `_checkpoint_paths`. All of them are `@inline`, take the run-wide `args` object and return `String` paths (or a `NamedTuple` of them) built with `joinpath`, so the same logic works on Windows and POSIX separators.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The functions only compute path strings; they never create directories, check writability, or verify that `args.simulation_settings.results_directory` exists, so callers must `mkpath` before opening a stream. The `args` argument is untyped, so a settings object missing `results_directory` or `checkpoint_directory` fails at runtime rather than at method dispatch.

## Provenance
Mapped from `src/io/config/io_config.jl` line 1.
