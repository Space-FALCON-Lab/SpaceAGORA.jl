---
id: io.io_outputs__append_save_field_columns_bang
label: _append_save_field_columns!
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _append_save_field_columns!
  lines:
  - 74
  - 74
inputs:
- id: results_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `results_df`.
- id: field
  type: Any
  units: n/a
  required: true
  description: Positional argument `field`.
- id: saved_data
  type: Vector
  units: n/a
  required: true
  description: Positional argument `saved_data`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: Nothing
  units: n/a
  description: Return value of `_append_save_field_columns!`; mutates `results_df`
    in place. Returns `nothing`.
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

# _append_save_field_columns!

## Purpose
`_append_save_field_columns!` adds the columns for one declared save field to the results `DataFrame`, expanding per-satellite fields into one column group per spacecraft. `_build_results_dataframe` calls it once per entry in `save_fields`.

## Design & Implementation
Signature `(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)`. It extracts `field_series = [snapshot[field.name] for snapshot in saved_data]`, indexing each snapshot (a dictionary-like object) by the field's `name`. If `field.per_satellite` is true it loops `sat_idx in 1:num_sats`, slices `value[sat_idx]` from every element, and calls `_append_series_columns!` with the prefix `"sc$(sat_idx)_$(field.column_prefix)"`; otherwise it calls it once with `field.column_prefix`. The `DataFrame` is mutated in place; the return value is `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results_df` | DataFrame | n/a | yes | Positional argument `results_df`. |
| in | `field` | Any | n/a | yes | Positional argument `field`. |
| in | `saved_data` | Vector | n/a | yes | Positional argument `saved_data`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_save_field_columns!`; mutates `results_df` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__build_results_dataframe|_build_results_dataframe]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:91-91`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`

**Downstream**

- `callees` → [[io.io_outputs__append_series_columns_bang|_append_series_columns!]] · `callers` · call · `src/io/outputs/io_outputs.jl:79-79`
- `callees` → [[simulation.persistence__append_series_columns_bang|_append_series_columns!]] · `callers` · call · `src/io/outputs/io_outputs.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
For per-satellite fields every snapshot value must be indexable up to `num_sats`; a snapshot whose vector is shorter (for example after a satellite is deactivated) raises `BoundsError`. `nothing` values in `field_series` are not tolerated in the per-satellite path because `value[sat_idx]` is applied unconditionally, unlike the recursive flattener which guards against `nothing`. A `field.name` missing from any snapshot throws `KeyError`. Column prefixes are not checked for uniqueness across fields.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 74.
