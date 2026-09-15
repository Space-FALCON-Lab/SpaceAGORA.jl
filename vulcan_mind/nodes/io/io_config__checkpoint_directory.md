---
id: io.io_config__checkpoint_directory
label: _checkpoint_directory
kind: function
source:
  file: src/io/config/io_config.jl
  symbol: _checkpoint_directory
  lines:
  - 20
  - 20
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
  description: Return value of `_checkpoint_directory`.
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

# _checkpoint_directory

## Purpose

Resolves the directory that checkpoint files are written to, honouring an explicit override in the settings and otherwise falling back to a `checkpoints` subdirectory of the results directory.

## Design & Implementation

If `args.simulation_settings.checkpoint_directory` is an empty string the function returns `joinpath(args.simulation_settings.results_directory, "checkpoints")`; otherwise it returns the configured value verbatim. Emptiness, not `nothing`, is the sentinel for `unset`, which keeps the settings field a plain `String` and avoids a `Union` type in the configuration struct.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_checkpoint_directory`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `callees` → `callers` · call · `src/io/config/io_config.jl:28-28`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

An explicitly configured directory is passed through without normalisation, so a relative path is interpreted against whatever the process working directory happens to be at write time, and a path containing whitespace only is treated as configured rather than unset. The directory is not created or tested for writability here.

## Provenance
Mapped from `src/io/config/io_config.jl` line 20.
