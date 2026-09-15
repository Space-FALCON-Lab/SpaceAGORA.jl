---
id: io.io_config__checkpoint_paths
label: _checkpoint_paths
kind: function
source:
  file: src/io/config/io_config.jl
  symbol: _checkpoint_paths
  lines:
  - 27
  - 27
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
  type: Any
  units: n/a
  description: Return value of `_checkpoint_paths`. Returns `(`.
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

# _checkpoint_paths

## Purpose

Bundles the two files that make up a checkpoint — the binary state blob and its TOML manifest — into one `NamedTuple` so callers save and restore them as a pair.

## Design & Implementation

Calls `_checkpoint_directory(args)` once, then returns `(data=joinpath(ckpt_dir, "simulation_checkpoint.bin"), manifest=joinpath(ckpt_dir, "simulation_checkpoint.manifest.toml"))`. Returning a `NamedTuple` rather than a tuple means call sites read `paths.data` and `paths.manifest`, which prevents the two from being swapped, and both names are fixed literals so a checkpoint written by one run is discoverable by any other pointed at the same directory.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_checkpoint_paths`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_serialization__clear_checkpoint_bang|_clear_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:84-84`
- [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:63-63`
- [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:35-35`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/config/io_config.jl`

**Downstream**

- `callees` → [[io.io_config__checkpoint_directory|_checkpoint_directory]] · `callers` · call · `src/io/config/io_config.jl:28-28`
- `callees` → [[simulation.resume_checkpoint__checkpoint_directory|_checkpoint_directory]] · `callers` · call · `src/io/config/io_config.jl:28-28`
<!-- vulcan:connections:end -->

## Limitations

Only one checkpoint per directory is representable: a new save overwrites the previous `.bin` and manifest, and there is no generation counter or rotation. The two writes are not atomic with respect to each other, so a crash between them leaves a manifest that disagrees with the blob. Neither file nor the containing directory is created here.

## Provenance
Mapped from `src/io/config/io_config.jl` line 27.
