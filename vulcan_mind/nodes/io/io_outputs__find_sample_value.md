---
id: io.io_outputs__find_sample_value
label: _find_sample_value
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _find_sample_value
  lines:
  - 30
  - 30
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
  type: Nothing
  units: n/a
  description: Return value of `_find_sample_value`. Returns `value` or `nothing`.
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

# _find_sample_value

## Purpose
`_find_sample_value` returns the first non-`nothing` element of a series so that `_append_series_columns!` can inspect a representative value and decide how to flatten the whole column. Series can contain `nothing` placeholders for time steps where a per-satellite or optional quantity was absent.

## Design & Implementation
It iterates `for value in series` and returns the first `value` for which `value !== nothing`; if every element is `nothing` or the series is empty it returns `nothing`. The comparison is identity-based (`!==`) so `missing` and `NaN` count as valid samples. The function is O(n) in the worst case but typically returns on the first element. It does not mutate its argument and works with any iterable, not only `Vector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `series` | Any | n/a | yes | Positional argument `series`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_find_sample_value`. Returns `value` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__append_series_columns_bang|_append_series_columns!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:40-40`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the first non-`nothing` element is examined, so a series whose elements change structure over time (for example a `NamedTuple` with different keys after a mode switch) is flattened according to that first sample and later elements may throw on `getproperty`/`getindex` in the caller. `missing` is treated as a sample, which makes the caller store the column raw even if later values are composite. An all-`nothing` series yields `nothing`, which the caller handles by storing the raw vector.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 30.
