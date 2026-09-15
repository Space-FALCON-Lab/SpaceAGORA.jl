---
id: simulation.persistence__append_series_columns_bang
label: _append_series_columns!
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _append_series_columns!
  lines:
  - 13
  - 13
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
  type: Any
  units: n/a
  description: Return value of `_append_series_columns!`; mutates `results_df` in
    place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _append_series_columns!

## Purpose
Flattens an arbitrarily nested saved series into flat DataFrame columns, giving each leaf a name derived from the dotted path it occupies inside the saved value.

## Design & Implementation
Forwards to `SimulationModel.IOOutputs._append_series_columns!(results_df, prefix, series)`, which mutates `results_df` in place. It takes a template from `_find_sample_value`, and if that is `nothing` or a flat scalar, meaning `missing`, `nothing`, a `Number`, an `AbstractString`, a `Symbol` or a `Bool`, it assigns `collect(series)` directly at `prefix`. A `NamedTuple` recurses over `keys(sample)`, an `AbstractDict` recurses over keys sorted by their string form for deterministic column order, and a `Tuple` or `AbstractArray` recurses over `eachindex(sample)`, each child series mapping `nothing` through unchanged and the name growing as `prefix_key`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results_df` | DataFrame | n/a | yes | Positional argument `results_df`. |
| in | `prefix` | String | n/a | yes | Positional argument `prefix`. |
| in | `series` | Any | n/a | yes | Positional argument `series`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_append_series_columns!`; mutates `results_df` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__append_save_field_columns_bang|_append_save_field_columns!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:79-79`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Recursion depth is unbounded, so a self-referential or very deep saved structure overflows the stack. Column names are built by string concatenation with an underscore, so a NamedTuple field already containing an underscore can collide with a nested path and one column silently overwrites another. Ragged arrays are expanded against the first sample's `eachindex`, so longer later elements lose their tail and shorter ones raise. Any composite type outside the recognised cases falls through to being stored whole in one column.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 13.
