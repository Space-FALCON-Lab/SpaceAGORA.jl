---
id: io.io_outputs_iooutputs
label: IOOutputs
kind: module
source:
  file: src/io/outputs/io_outputs.jl
  symbol: IOOutputs
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

# IOOutputs

## Purpose
`IOOutputs` is the module that turns the simulation's saved-callback snapshots into tabular outputs: it concatenates per-segment `SavedValues` into flat time and data vectors, flattens nested snapshot values into a `DataFrame`, and writes that frame to CSV or to a Feather+manifest results bundle. It sits between the solver's save callbacks and the on-disk artifacts read by analysis tooling.

## Design & Implementation
The module depends on `Arrow`, `CSV`, `DataFrames`, `Dates`, `TOML` and the sibling `IOConfig` (path resolution) and `IOSerialization` (`_atomic_write_file`, `_sha256_hex`). Functions are layered: `_append_saved_segment!` accumulates segments, `_build_results_dataframe` drives `_append_save_field_columns!` which recurses through `_append_series_columns!` using `_find_sample_value` and `_is_flat_scalar` to decide column splitting, and `_write_results_csv!` / `_write_results_bundle!` perform the atomic writes. `_write_results_bundle!` writes `<prefix>.feather`, optionally `<prefix>.csv`, and a `<prefix>.manifest.toml` containing schema version, UTC creation time, `mission_time_s`, step count, spacecraft count, `orientation_sim` and per-file size and SHA-256. The five underscore-prefixed functions are explicitly exported for use by the simulation engine.

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

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Column naming follows the recursive `prefix_key` scheme with no collision detection, so two save fields whose prefixes and keys concatenate to the same string overwrite each other in the `DataFrame`. All snapshots are materialised into Julia vectors before writing, so memory scales with steps times fields. Feather output is always written even when only CSV was requested. The manifest's `created_utc` uses `now(UTC)` so bundles are not byte-reproducible across runs.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 1.
