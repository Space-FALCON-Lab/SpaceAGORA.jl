---
id: gnc.planner_comparison_rpo_comparison_axis_name
label: rpo_comparison_axis_name
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_axis_name
  lines:
  - 716
  - 716
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
  type: Symbol
  units: n/a
  description: Return value of `rpo_comparison_axis_name`. Returns `Symbol(prefix,
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

# rpo_comparison_axis_name

## Purpose
Produces the Plotly layout key (`:xaxis`, `:xaxis2`, ...) for the axis belonging to subplot `idx`, matching Plotly's convention that the first axis carries no numeric suffix.

## Design & Implementation
`rpo_comparison_axis_name(prefix::AbstractString, idx::Integer)` returns `Symbol(prefix)` for `idx == 1` and `Symbol(prefix, idx)` otherwise. The summary plot stores these as keys in `layout_kwargs` so each subplot's domain and titles are attached to its own axis.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prefix` | AbstractString | n/a | yes | Positional argument `prefix`. |
| in | `idx` | Integer | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `rpo_comparison_axis_name`. Returns `Symbol(prefix, idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:795-795`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Does not validate `idx >= 1`; `idx = 0` yields `:xaxis0`, which Plotly ignores. The paired `rpo_comparison_trace_axis_name` must be used for trace references, since Plotly expects `"x2"` on traces but `xaxis2` in layout, and mixing the two helpers produces a silently mis-anchored subplot.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 716.
