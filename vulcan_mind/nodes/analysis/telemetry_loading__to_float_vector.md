---
id: analysis.telemetry_loading__to_float_vector
label: _to_float_vector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _to_float_vector
  lines:
  - 4
  - 4
inputs:
- id: values
  type: Any
  units: n/a
  required: true
  description: Positional argument `values`.
- id: context
  type: String
  units: n/a
  required: true
  description: Positional argument `context`.
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
  type: Vector{Float64}
  units: n/a
  description: Return value of `_to_float_vector`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _to_float_vector

## Purpose
Converts an arbitrary indexable column of numbers (typically a `DataFrame` column read from Arrow, possibly with `Union{Missing, T}` element type) into a dense `Vector{Float64}`, failing loudly on any `missing` so downstream numerical code never sees non-numeric entries.

## Design & Implementation
Marked `@inline`; allocates `Vector{Float64}(undef, length(values))` and loops `@inbounds for i in eachindex(values)`. Each element is compared with `v === missing`; a hit throws `ArgumentError("Missing value in $context at index $i")` where `context` is a caller-supplied label such as `"telemetry-altitude"`. Otherwise `Float64(v)` is stored, which accepts any `Real`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Any | n/a | yes | Positional argument `values`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Float64} | n/a | — | Return value of `_to_float_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:139-139`
- [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:79-79`
- [[analysis.telemetry_loading__load_telemetry_curve|_load_telemetry_curve]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:208-208`
- [[analysis.telemetry_loading__require_column|_require_column]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:66-66`
- [[envana.ana_calibration_estimate_event_biases|_estimate_event_biases]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:43-43`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:9-9`
<!-- vulcan:connections:end -->

## Limitations
`eachindex(values)` may not be `1:length(values)` for offset arrays, in which case `out[i]` under `@inbounds` writes out of range without an error. `Float64(v)` throws a `MethodError` (not the labelled `ArgumentError`) for strings, so a text column produces a less helpful message. `NaN` and `Inf` pass through unchanged and must be handled by callers. The context string is interpolated even on the success path only lazily, so cost is one allocation for the output vector.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 4.
