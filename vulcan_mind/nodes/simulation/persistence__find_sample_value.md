---
id: simulation.persistence__find_sample_value
label: _find_sample_value
kind: function
source:
  file: src/simulation/engine/persistence.jl
  symbol: _find_sample_value
  lines:
  - 16
  - 16
inputs:
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
  type: SimulationModel.IOOutputs._find_sample_value
  units: n/a
  description: Return value of `_find_sample_value`. Returns `SimulationModel.IOOutputs._find_sample_value(series)`.
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

# _find_sample_value

## Purpose
Scans a saved series and returns the first element that is not `nothing`, giving the column-expansion logic a concrete exemplar whose structure determines how the series is flattened into DataFrame columns.

## Design & Implementation
Forwards to `SimulationModel.IOOutputs._find_sample_value(series)`, a linear loop that returns the first `value !== nothing` and returns `nothing` when the whole series is empty or entirely missing. The identity comparison is used rather than equality so that it works on any element type without invoking user-defined `==`. Its result drives the branch selection in `_append_series_columns!` between scalar, NamedTuple, dictionary and array handling.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `series` | Any | n/a | yes | Positional argument `series`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.IOOutputs._find_sample_value | n/a | — | Return value of `_find_sample_value`. Returns `SimulationModel.IOOutputs._find_sample_value(series)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__append_series_columns_bang|_append_series_columns!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:40-40`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/persistence.jl`

**Downstream**

- `callees` → [[io.io_outputs__build_results_dataframe|_build_results_dataframe]] · `callers` · feedback · `src/simulation/engine/persistence.jl:18-18`
- `callees` → [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callers` · call · `src/simulation/engine/persistence.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
Only the first present value is inspected, so a heterogeneous series whose later elements have a different shape, for instance a NamedTuple that gains a field partway through a run, is expanded against the wrong template and the mismatched entries fail at property access. Scanning is worst case linear in the series length, paid once per column, and an all-`nothing` series degrades to a single undifferentiated column.

## Provenance
Mapped from `src/simulation/engine/persistence.jl` line 16.
