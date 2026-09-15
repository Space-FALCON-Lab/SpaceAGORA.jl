---
id: analysis.telemetry_loading__require_column
label: _require_column
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _require_column
  lines:
  - 63
  - 63
inputs:
- id: df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `df`.
- id: candidates
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `candidates`.
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
  description: Return value of `_require_column`.
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

# _require_column

## Purpose
Fetches the first available column from a list of candidate names in a `DataFrame` and returns it as `Vector{Float64}`, so telemetry loaders can accept both legacy and current column naming (for example `sc1_pos_1` versus `sc1_position_1`) while still failing clearly when neither is present.

## Design & Implementation
Marked `@inline`; iterates `candidates::Vector{String}` in order and tests `col in names(df)`. On the first hit it returns `_to_float_vector(df[!, col], "$context:$col")`, using the non-copying `df[!, col]` accessor. If no candidate matches it throws `ArgumentError("Missing required column for $context. Tried: ...")` listing every candidate joined with commas.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `df` | DataFrame | n/a | yes | Positional argument `df`. |
| in | `candidates` | Vector{String} | n/a | yes | Positional argument `candidates`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Float64} | n/a | — | Return value of `_require_column`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:140-140`
- [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:73-73`
- [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:216-216`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:66-66`
<!-- vulcan:connections:end -->

## Limitations
`names(df)` allocates a fresh `Vector{String}` on every candidate check, which is wasteful for wide frames, though loaders call it only a handful of times. Candidate order defines priority silently, so if both a legacy and a new column exist with different data the first listed wins with no warning. `df[!, col]` returns the underlying column, so the subsequent conversion copies but the original is not protected from aliasing elsewhere.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 63.
