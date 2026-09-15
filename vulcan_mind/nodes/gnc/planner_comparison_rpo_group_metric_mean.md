---
id: gnc.planner_comparison_rpo_group_metric_mean
label: rpo_group_metric_mean
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_group_metric_mean
  lines:
  - 677
  - 677
inputs:
- id: results
  type: Any
  units: n/a
  required: true
  description: Positional argument `results`.
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
  type: Any
  units: n/a
  description: Return value of `rpo_group_metric_mean`. Returns `sum(values) / length(values)`.
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

# rpo_group_metric_mean

## Purpose
Computes the arithmetic mean of one numeric metric across a collection of result rows, coercing Booleans to 0/1 and ignoring non-finite values, for summary tables of planner comparison batches.

## Design & Implementation
`rpo_group_metric_mean(results, metric::Symbol)` loops over rows, reads `getproperty(row, metric)`, maps `true`/`false` to `1.0`/`0.0`, converts to `Float64`, and pushes finite values into a `Float64[]`. Returns `NaN` when no finite values remain, otherwise `sum / length`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results` | Any | n/a | yes | Positional argument `results`. |
| in | `metric` | Symbol | n/a | yes | Positional argument `metric`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_group_metric_mean`. Returns `sum(values) / length(values)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:682-682`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:683-683`
<!-- vulcan:connections:end -->

## Limitations
Silently excludes `NaN` and `Inf` entries (for example `fuel_used_pct` when propellant mass is zero, or `min_clearance = Inf` for zero-step tracking), so the reported mean can be over a smaller sample than the row count without indication. It does not apply `rpo_comparison_metric_includes_result`, so failed cases are averaged in unless the caller pre-filters.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 677.
