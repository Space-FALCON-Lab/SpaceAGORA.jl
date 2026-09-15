---
id: gnc.planner_comparison_rpo_comparison_metric_summary_plot
label: rpo_comparison_metric_summary_plot
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_metric_summary_plot
  lines:
  - 728
  - 728
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
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
  type: PlotlyJS.Plot
  units: n/a
  description: Return value of `rpo_comparison_metric_summary_plot`. Returns `PlotlyJS.Plot(traces,
    PlotlyJS.Layout(; layout_kwargs...))`.
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

# rpo_comparison_metric_summary_plot

## Purpose
Builds the 3 x 3 grid of box-plus-scatter subplots summarising every metric in `rpo_comparison_metric_specs` across planners, the primary visual output of a comparison batch.

## Design & Implementation
`rpo_comparison_metric_summary_plot(batch)` computes subplot domains from `xgap = 0.055`, `ygap = 0.095`, `xspan = (1 - xgap * 2) / 3`, `yspan = (1 - ygap * 2) / 3`. For each metric it gathers `(planner_label, value, "case <label>")` triples over rows passing `rpo_comparison_metric_includes_result`, then pushes a `PlotlyJS.box` (no points, `boxmean = true`, blue fill) and an overlaid `PlotlyJS.scatter` of orange case markers with hover text, both bound to axes from `rpo_comparison_trace_axis_name`. Layout entries per axis set `domain`, `anchor`, and titles (metric label on x, units on y). Returns `PlotlyJS.Plot(traces, PlotlyJS.Layout(; layout_kwargs...))` at 1200 x 920 px with the legend hidden.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.Plot | n/a | — | Return value of `rpo_comparison_metric_summary_plot`. Returns `PlotlyJS.Plot(traces, PlotlyJS.Layout(; layout_kwargs...))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1167-1167`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_axis_name|rpo_comparison_axis_name]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:795-795`
- `callees` → [[gnc.planner_comparison_rpo_comparison_metric_includes_result|rpo_comparison_metric_includes_result]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:762-762`
- `callees` → [[gnc.planner_comparison_rpo_comparison_metric_specs|rpo_comparison_metric_specs]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:729-729`
- `callees` → [[gnc.planner_comparison_rpo_comparison_metric_value|rpo_comparison_metric_value]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:763-763`
- `callees` → [[gnc.planner_comparison_rpo_comparison_trace_axis_name|rpo_comparison_trace_axis_name]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:753-753`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:764-764`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:760-760`
<!-- vulcan:connections:end -->

## Limitations
Grid size and pixel dimensions are constants, so more than nine metrics overflow the layout. Colours and gaps are hard-coded RGB strings. Box statistics over a handful of cases per planner are of limited significance, and planners with no successful cases produce empty categories for most subplots. Every planner label is normalised via `rpo_comparison_planner_label`, which throws for unknown symbols.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 728.
