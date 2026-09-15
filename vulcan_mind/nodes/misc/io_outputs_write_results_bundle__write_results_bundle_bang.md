---
id: misc.io_outputs_write_results_bundle__write_results_bundle_bang
label: _write_results_bundle!
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _write_results_bundle!
  lines:
  - 120
  - 169
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: IOConfig path helpers (_results_bundle_prefix) and IOSerialization
    atomic-write and hashing helpers used to emit the bundle.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: bundle_files
  type: Feather+TOML
  units: n/a
  description: Arrow feather table, optional CSV twin, and a TOML manifest recording
    schema version, sizes and SHA-256 digests.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
- arrow
- manifest
charts:
- misc
origin: agent
---

# _write_results_bundle!

## Purpose
`_write_results_bundle!` persists a completed simulation's results table as a self-describing bundle rather than a bare data file. The bundle is what downstream analysis, regression comparison and archiving consume, so it pairs the binary table with a manifest that states the schema version it was written under and the cryptographic digest of every file it names. That makes an incomplete or corrupted write detectable after the fact instead of silently feeding a study.

## Model & Assumptions
The bundle base path comes from `IOConfig._results_bundle_prefix(args)`; the function appends `.feather` and `.manifest.toml` to it. The results `DataFrame` is assumed already flattened into scalar columns by `_build_results_dataframe` and its recursive column expander, because Arrow cannot store the nested named tuples and dictionaries that the solver's save callback produces. The `args` object must expose `simulation_settings.save_csv`, `mission_configuration.mission_time`, `mission_configuration.orientation_sim` and `dynamics_model.spacecraft`.

## Design & Implementation
Every file goes through `IOSerialization._atomic_write_file`, which writes to a uniquely named temporary sibling and then moves it into place, so a reader never observes a half-written table. The feather table is written first with `Arrow.write`. A `files` dictionary then records, per artifact, its path, `filesize`, and `IOSerialization._sha256_hex` digest. When `save_csv` is set the function reuses an already-written CSV if `csv_path` was supplied by `_write_results_csv!`, and only writes its own `<prefix>.csv` when none was passed — this avoids duplicating the collision-handling logic that the CSV writer implements. The manifest finally records `schema_version`, `created_utc` from `now(UTC)`, `mission_time_s`, the row count as `steps`, `spacecraft_count`, `orientation_sim`, and the `files` table, printed with `TOML.print`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | IOConfig path helpers (_results_bundle_prefix) and IOSerialization atomic-write and hashing helpers used to emit the bundle. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `bundle_files` | Feather+TOML | n/a | — | Arrow feather table, optional CSV twin, and a TOML manifest recording schema version, sizes and SHA-256 digests. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`
- [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:129-129`

**Downstream**

- `callees` → [[io.io_serialization__atomic_write_file|_atomic_write_file]] · `callers` · call · `src/io/outputs/io_outputs.jl:131-131`
- `callees` → [[io.io_serialization__sha256_hex|_sha256_hex]] · `callers` · call · `src/io/outputs/io_outputs.jl:137-137`
- `callees` → [[io.results_bundle_prefix|_results_bundle_prefix]] · `callers` · call · `src/io/outputs/io_outputs.jl:127-127`
- `callees` → [[simulation.persistence__results_bundle_prefix|_results_bundle_prefix]] · `callers` · call · `src/io/outputs/io_outputs.jl:127-127`
- `callees` → [[simulation.persistence__sha256_hex|_sha256_hex]] · `callers` · call · `src/io/outputs/io_outputs.jl:137-137`
<!-- vulcan:connections:end -->

## Limitations
The manifest is written after the data files, so a crash between the two leaves an unmanifested feather table on disk. Digests are computed by re-reading each file, which doubles IO for large results. The bundle prefix is stable rather than run-unique, so two concurrent runs sharing a results directory overwrite each other's feather and manifest without the collision fallback that protects the CSV path.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl:120-169`, with the atomic writer and hasher at `src/io/serialization/io_serialization.jl:10-32`.
