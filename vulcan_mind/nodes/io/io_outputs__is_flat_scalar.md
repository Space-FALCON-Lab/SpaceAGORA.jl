---
id: io.io_outputs__is_flat_scalar
label: _is_flat_scalar
kind: function
source:
  file: src/io/outputs/io_outputs.jl
  symbol: _is_flat_scalar
  lines:
  - 26
  - 26
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
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
  type: Bool
  units: n/a
  description: Return value of `_is_flat_scalar`.
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

# _is_flat_scalar

## Purpose
`_is_flat_scalar` is the predicate `_append_series_columns!` uses to decide that a sampled value can be stored directly as a single `DataFrame` column rather than being decomposed into child columns. It marks the leaf case of the recursive flattening.

## Design & Implementation
Declared `@inline` with return type `Bool`. A value is flat when it is `missing`, `nothing`, a `Number`, an `AbstractString`, a `Symbol`, or a `Bool` (the last is redundant with `Number` but kept explicit). The checks use `===` for the singletons and `isa` for the abstract types, so `Float64`, `Int`, `Complex`, `String`, `SubString` and `Symbol` are all leaves. Anything else, including `NamedTuple`, `AbstractDict`, `Tuple`, `AbstractArray` and arbitrary structs, is treated as composite by the caller. Pure function, no allocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_flat_scalar`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[io.io_outputs__append_series_columns_bang|_append_series_columns!]] · `callees` → `callers` · call · `src/io/outputs/io_outputs.jl:41-41`
- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/outputs/io_outputs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A `StaticArrays.SVector` or any `AbstractArray` is not flat, so vector-valued fields always expand into indexed columns even when the caller wants a single object column. Custom struct types that are not `Number` are also non-flat but not decomposable, so the caller falls through to storing them raw; this predicate cannot distinguish that case. Because `Bool <: Number`, the explicit `Bool` branch never changes the result.

## Provenance
Mapped from `src/io/outputs/io_outputs.jl` line 26.
