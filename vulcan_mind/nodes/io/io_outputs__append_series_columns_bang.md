---
id: io.io_outputs__append_series_columns_bang
label: _append_series_columns!
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _append_series_columns!
  lines:
  - 39
  - 39
inputs:
- id: results_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `results_df`.
- id: prefix
  type: String
  units: n/a
  required: true
  description: Positional argument `prefix`.
- id: series
  type: Any
  units: n/a
  required: true
  description: Positional argument `series`.
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
  description: Return value of `_append_series_columns!`; mutates `results_df` in
    place. Returns `nothing`.
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

# _append_series_columns!

## Purpose
`_append_series_columns!` recursively flattens a time series of possibly nested values into one or more columns of a `DataFrame`, naming child columns by joining the parent `prefix` with the field key or index. It is how per-step snapshots containing `NamedTuple`s, dictionaries, tuples and arrays become wide tabular CSV/Feather output.

## Design & Implementation
Signature `(results_df::DataFrame, prefix::String, series)`. It takes a sample with `_find_sample_value`; when the sample is `nothing` or `_is_flat_scalar`, the whole `collect(series)` is assigned to `results_df[!, prefix]`. For a `NamedTuple` sample it recurses per key with `getproperty`; for an `AbstractDict` it recurses over keys sorted by `string` for deterministic ordering; for a `Tuple` or `AbstractArray` it recurses over `eachindex(sample)` with `getindex`. In each branch the child series maps `nothing` elements through unchanged so gaps propagate. Any other sample type falls through to a raw column. Column names are `prefix * "_" * key` (or index). The function mutates `results_df` in place and returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results_df` | DataFrame | n/a | yes | Positional argument `results_df`. |
| in | `prefix` | String | n/a | yes | Positional argument `prefix`. |
| in | `series` | Any | n/a | yes | Positional argument `series`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_series_columns!`; mutates `results_df` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__append_save_field_columns_bang|_append_save_field_columns!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:79-79`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`

**Downstream**

- `callees` → [[io.io_outputs__find_sample_value|_find_sample_value]] · `callers` · call · `src/io/outputs/io_outputs.jl:40-40`
- `callees` → [[io.io_outputs__is_flat_scalar|_is_flat_scalar]] · `callers` · call · `src/io/outputs/io_outputs.jl:41-41`
- `callees` → [[simulation.persistence__find_sample_value|_find_sample_value]] · `callers` · call · `src/io/outputs/io_outputs.jl:40-40`
<!-- vulcan:connections:end -->

## Limitations
Structure is inferred from a single sample; if later elements have different keys, lengths, or types, the recursion throws `KeyError`, `BoundsError` or `MethodError` mid-flattening, leaving the frame partially populated. Multi-dimensional arrays are indexed by `CartesianIndex`, producing column names like `prefix_CartesianIndex(1, 2)`. Dict keys that stringify identically collide. Each recursion level allocates a full new vector, so deeply nested or wide fields cost O(steps × leaves) allocations. No column-name collision check exists.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 39.
