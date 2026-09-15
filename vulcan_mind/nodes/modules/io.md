---
id: module.io
label: IOConfig
kind: module
source:
  file: src/io/config/io_config.jl
  symbol: IOConfig
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Result-directory, collision-safe CSV, and checkpoint path helpers exported
    by IOConfig.
tags:
- module
charts:
- master
origin: agent
---

# IOConfig

## Purpose
`IOConfig` owns the path conventions used by simulation output and checkpoint persistence. It does not serialize state itself; instead, it derives stable paths from `args.simulation_settings` so the serialization and output modules share one configuration surface. The functions are marked `@inline` because each performs a small field lookup or path join on a hot configuration path.

## Theory & Math
There is no numerical model. The path contract is a deterministic mapping from simulation settings to filesystem locations. The collision-safe CSV path adds a UTC millisecond timestamp, the current process identifier, and a random unsigned token, making concurrent result writers unlikely to select the same filename.

## Model & Assumptions
`args` must expose `simulation_settings` with `results_directory` and `checkpoint_directory` fields. An empty checkpoint directory means the checkpoint files belong in a `checkpoints` child of the results directory. Path joining follows the host operating system through `joinpath`, so callers should not embed platform-specific separators.

## Design & Implementation
`_results_bundle_prefix` returns the base name for a simulation result bundle, and `_results_csv_path` appends the stable CSV filename. `_collision_results_csv_path` formats the current UTC time and combines it with process and random identifiers. `_checkpoint_directory` selects the explicit directory or derives the default. `_checkpoint_paths` returns a named tuple containing the binary checkpoint and TOML manifest paths. `SimulationModel` re-exports this module at the IO owner section of `simulation_model.jl`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Result-directory, collision-safe CSV, and checkpoint path helpers exported by IOConfig. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[io.io_config__checkpoint_directory|_checkpoint_directory]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[io.io_config__checkpoint_paths|_checkpoint_paths]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[io.io_config__collision_results_csv_path|_collision_results_csv_path]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[io.io_config__results_csv_path|_results_csv_path]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[io.io_config_ioconfig|IOConfig]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[io.io_outputs__append_save_field_columns_bang|_append_save_field_columns!]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__append_saved_segment_bang|_append_saved_segment!]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__append_series_columns_bang|_append_series_columns!]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__build_results_dataframe|_build_results_dataframe]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__find_sample_value|_find_sample_value]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__is_flat_scalar|_is_flat_scalar]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_outputs_iooutputs|IOOutputs]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[io.io_serialization__clear_checkpoint_bang|_clear_checkpoint!]] · `module_api` · call · `src/io/serialization/io_serialization.jl`
- `api` → [[io.io_serialization_ioserialization|IOSerialization]] · `module_api` · call · `src/io/serialization/io_serialization.jl`
- `api` → [[io.results_bundle_prefix|_results_bundle_prefix]] · `module_api` · call · `src/io/config/io_config.jl`
- `api` → [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `module_api` · call · `src/io/outputs/io_outputs.jl`
- `api` → [[module.core|SimulationModel]] · `io` · call · `src/core/simulation_model.jl:109-113`
<!-- vulcan:connections:end -->

## Limitations
The helpers do not create directories, verify writability, or reserve filenames. The collision-resistant name is probabilistic because it includes a random token, and the timestamp uses the wall clock. Invalid or missing settings fields fail at field access. Consumers must still coordinate writes when a result bundle contains multiple related files.

## Provenance
Mapped from `src/io/config/io_config.jl` and its re-export site in `src/core/simulation_model.jl`.
