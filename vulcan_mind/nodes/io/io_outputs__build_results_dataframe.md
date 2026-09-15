---
id: io.io_outputs__build_results_dataframe
label: _build_results_dataframe
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _build_results_dataframe
  lines:
  - 87
  - 87
inputs:
- id: times
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times`.
- id: saved_data
  type: Vector
  units: n/a
  required: true
  description: Positional argument `saved_data`.
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
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
  type: DataFrame
  units: n/a
  description: Return value of `_build_results_dataframe`.
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

# _build_results_dataframe

## Purpose
`_build_results_dataframe` assembles the final tabular results for a simulation run: a `DataFrame` whose first column is `time` and whose remaining columns are the flattened save fields, one group per spacecraft where applicable. It is the single entry point the engine uses before writing CSV or Feather.

## Design & Implementation
Signature `(times::Vector{Float64}, saved_data::Vector, save_fields, args)::DataFrame`. It creates `DataFrame(time=times)`, reads `num_sats = length(args.dynamics_model.spacecraft)`, and iterates `save_fields` in order, delegating each to `_append_save_field_columns!`. Column order therefore follows the order of `save_fields` and, within per-satellite fields, spacecraft index then nested key order. The `times` vector is stored by reference as the `time` column (no copy). Returns the new frame; `args` and `saved_data` are not mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `times` | Vector{Float64} | n/a | yes | Positional argument `times`. |
| in | `saved_data` | Vector | n/a | yes | Positional argument `saved_data`. |
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DataFrame | n/a | — | Return value of `_build_results_dataframe`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`
- [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:126-126`
- [[simulation.persistence__find_sample_value|_find_sample_value]] · `callees` → `callers` · feedback · `src/simulation/engine/persistence.jl:18-18`

**Downstream**

- `callees` → [[io.io_outputs__append_save_field_columns_bang|_append_save_field_columns!]] · `callers` · call · `src/io/outputs/io_outputs.jl:91-91`
<!-- vulcan:connections:end -->

## Limitations
`length(times)` must equal `length(saved_data)` or the `DataFrame` constructor / column assignment throws a `DimensionMismatch` after partially building. Because `times` is inserted without copying, later mutation of the caller's vector changes the frame. Fields are processed sequentially with no error isolation, so one malformed field aborts the whole frame. `num_sats` is taken from the model definition rather than from the data, so the frame silently assumes all spacecraft were active for all saved steps.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 87.
