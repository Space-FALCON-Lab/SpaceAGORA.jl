---
id: simx.engine_persistence_write_results_bundle__write_results_bundle_bang
label: _write_results_bundle!
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _write_results_bundle!
  lines:
  - 24
  - 37
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: results_df
  type: DataFrame
  units: n/a
  required: true
  description: Assembled results table whose columns were expanded from the saved
    state series by _build_results_dataframe.
- id: times
  type: Vector{Float64}
  units: s
  required: true
  description: Saved sample times aligned row-for-row with the results table.
- id: csv_path
  type: Union{Nothing,String}
  units: n/a
  required: true
  description: Optional already-written CSV path to reference instead of writing the
    table a second time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: bundle_manifest
  type: NamedTuple
  units: n/a
  description: Bundle description returned by the IO layer, including the written
    paths and the schema version stamped into the manifest.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _write_results_bundle!

## Purpose
`_write_results_bundle!` is the engine-side entry point for emitting a run's results bundle. It forwards the assembled table, the sample times and the run configuration to `SimulationModel.IOOutputs._write_results_bundle!`, adding the one piece of state the IO layer cannot know: the engine's `RESULTS_BUNDLE_SCHEMA_VERSION` constant.

## Model & Assumptions
The whole file is a thin forwarding layer. Every other definition in it is a one-line `@inline` alias onto `IOConfig`, `IOSerialization` or `IOOutputs`, which keeps the engine's call sites short while leaving the actual format ownership in the IO modules. Because the aliases are `@inline`, the indirection has no runtime cost and the engine can be read without chasing module prefixes on every persistence call.

## Design & Implementation
The forwarded call passes `results_df`, `times`, `args` and `RESULTS_BUNDLE_SCHEMA_VERSION` positionally and `csv_path` as a keyword. Passing an existing `csv_path` is the normal case when the CSV was written first by `_write_results_csv!`, so the bundle references that file rather than re-serialising the table. Sibling aliases in this file cover path construction (`_results_bundle_prefix`, `_results_csv_path`, `_collision_results_csv_path`), integrity (`_sha256_hex`), durability (`_atomic_write_file`, which writes through a temporary and renames), and series assembly (`_append_saved_segment!`, `_append_series_columns!`, `_find_sample_value`).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `results_df` | DataFrame | n/a | yes | Assembled results table whose columns were expanded from the saved state series by _build_results_dataframe. |
| in | `times` | Vector{Float64} | s | yes | Saved sample times aligned row-for-row with the results table. |
| in | `csv_path` | Union{Nothing,String} | n/a | yes | Optional already-written CSV path to reference instead of writing the table a second time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `bundle_manifest` | NamedTuple | n/a | — | Bundle description returned by the IO layer, including the written paths and the schema version stamped into the manifest. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`
- [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:129-129`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The schema version is a compile-time constant string, so writing a bundle at an older schema requires editing the constant rather than passing a parameter. The function returns whatever the IO layer returns and performs no verification that the bundle landed, so a partially written bundle is only detected by the atomic-write helper's rename semantics and the recorded SHA-256 digests.

## Provenance
Mapped from `src/simulation/engine/persistence.jl:24-37`; `RESULTS_BUNDLE_SCHEMA_VERSION` is defined at `src/simulation/engine/setup.jl:18`.
