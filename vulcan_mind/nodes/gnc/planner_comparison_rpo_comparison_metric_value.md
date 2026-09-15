---
id: gnc.planner_comparison_rpo_comparison_metric_value
label: rpo_comparison_metric_value
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_metric_value
  lines:
  - 705
  - 705
inputs:
- id: row
  type: Any
  units: n/a
  required: true
  description: Positional argument `row`.
- id: metric
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `metric`.
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
  type: Float64
  units: n/a
  description: Return value of `rpo_comparison_metric_value`. Returns `Float64(value)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_comparison_metric_value

## Purpose
Reads one comparison metric from a flattened result row as a `Float64`, synthesising the `:success_pct` metric from the Boolean `success` field and coercing other Booleans to 0/1.

## Design & Implementation
`rpo_comparison_metric_value(row, metric::Symbol)` returns `100.0` or `0.0` for `:success_pct` based on `row.success`; otherwise fetches `getproperty(row, metric)`, returns `1.0`/`0.0` for a `Bool`, and `Float64(value)` for anything else. Used by the summary plot to build y-values.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `row` | Any | n/a | yes | Positional argument `row`. |
| in | `metric` | Symbol | n/a | yes | Positional argument `metric`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `rpo_comparison_metric_value`. Returns `Float64(value)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:763-763`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:709-709`
<!-- vulcan:connections:end -->

## Limitations
`Float64(value)` throws `MethodError` for non-numeric fields such as `planner_label` or `case_label`, so callers must only pass metric symbols from `rpo_comparison_metric_specs`. Missing fields raise `ErrorException` from `getproperty` rather than returning `NaN`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 705.
