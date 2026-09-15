---
id: gnc.planner_comparison_rpo_comparison_trace_axis_name
label: rpo_comparison_trace_axis_name
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_trace_axis_name
  lines:
  - 722
  - 722
inputs:
- id: prefix
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `prefix`.
- id: idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `rpo_comparison_trace_axis_name`. Returns `string(prefix,
    idx)`.
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

# rpo_comparison_trace_axis_name

## Purpose
Produces the string a Plotly trace uses to reference its axis (`"x"`, `"x2"`, `"y3"`), complementing `rpo_comparison_axis_name` which produces the layout-side key.

## Design & Implementation
`rpo_comparison_trace_axis_name(prefix::AbstractString, idx::Integer)` returns `prefix` unchanged for `idx == 1` and `string(prefix, idx)` otherwise. Called twice per subplot in `rpo_comparison_metric_summary_plot` with prefixes `"x"` and `"y"`, and the results are also passed as `anchor` values in the axis layout attributes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prefix` | AbstractString | n/a | yes | Positional argument `prefix`. |
| in | `idx` | Integer | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_comparison_trace_axis_name`. Returns `string(prefix, idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:753-753`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returns a `String` where the layout helper returns a `Symbol`; the distinction is required by PlotlyJS but easy to confuse. No range check on `idx`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 722.
