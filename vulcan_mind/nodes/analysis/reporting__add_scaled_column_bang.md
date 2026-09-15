---
id: analysis.reporting__add_scaled_column_bang
label: _add_scaled_column!
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _add_scaled_column!
  lines:
  - 54
  - 54
inputs:
- id: src
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `src`.
- id: dst
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `dst`.
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
  description: Return value of `_add_scaled_column!`; mutates `src` in place. Returns
    `nothing`.
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

# _add_scaled_column!

## Purpose
Local closure defined inside both `_append_display_metric_columns!` and `_append_display_error_columns!` that copies a numeric DataFrame column into a new column after multiplying each element by that row's display scale factor.

## Design & Implementation
Signature `_add_scaled_column!(src::Symbol, dst::Symbol)`. If the captured DataFrame lacks `src` it returns `nothing` without creating `dst`. Otherwise it reads `df[!, src]`, allocates `Vector{Float64}(undef, n)`, and in an `@inbounds` loop sets `out[i] = Float64(src_values[i]) * scales[i]`, where `scales` is the closure-captured per-row vector computed by the enclosing function. The result is assigned with `df[!, dst] = out`, which replaces any existing column of that name. The metric variant is applied to `mae_km`, `rmse_km`, `max_abs_km`, `p95_abs_km`, `bias_km`, `limit_max_abs_km`, `limit_max_rmse_km`; the error variant to `telemetry_value_km`, `sim_interp_value_km`, `error_km`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `src` | Symbol | n/a | yes | Positional argument `src`. |
| in | `dst` | Symbol | n/a | yes | Positional argument `dst`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_add_scaled_column!`; mutates `src` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:100-100`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:59-59`
<!-- vulcan:connections:end -->

## Limitations
`Float64(src_values[i])` throws on `missing` or non-numeric entries; there is no guard. The `@inbounds` loop assumes `length(scales) == nrow(df)`, which holds only because the enclosing function built `scales` from the same frame immediately before. Being a closure, it is recompiled per enclosing call and cannot be tested in isolation.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 54.
