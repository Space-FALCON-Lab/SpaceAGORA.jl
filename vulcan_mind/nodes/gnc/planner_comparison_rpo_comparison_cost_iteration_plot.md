---
id: gnc.planner_comparison_rpo_comparison_cost_iteration_plot
label: rpo_comparison_cost_iteration_plot
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_cost_iteration_plot
  lines:
  - 1072
  - 1072
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: planner
  type: Any
  units: n/a
  required: true
  description: Keyword argument `planner`.
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
  description: Return value of `rpo_comparison_cost_iteration_plot`. Returns `PlotlyJS.Plot(`.
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

# rpo_comparison_cost_iteration_plot

## Purpose
Creates a 2-D Plotly figure of planner cost versus iteration for every case of one planner, giving a convergence view for HYPR, CHOMP, STOMP, and any other planner that records `cost_history`.

## Design & Implementation
`rpo_comparison_cost_iteration_plot(batch; planner)` normalises the planner and throws `ArgumentError("No planner results found for ...")` if `batch.plans_by_planner` lacks it. For each plan with a `cost_history` property it filters to finite `Float64` costs, skips empty histories, computes x values via `rpo_comparison_cost_iteration_xvalues`, and pushes a `scatter` trace in `"markers"` mode for a single sample or `"lines+markers"` otherwise, named by case label with a hover template showing iteration and cost. If no traces result, a placeholder empty trace named "no finite cost history" is added. Layout is 920 x 620 px with y `rangemode = "tozero"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `planner` | Any | n/a | yes | Keyword argument `planner`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.Plot | n/a | — | Return value of `rpo_comparison_cost_iteration_plot`. Returns `PlotlyJS.Plot(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1180-1180`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1081-1081`
- `callees` → [[gnc.planner_comparison_rpo_comparison_cost_iteration_xvalues|rpo_comparison_cost_iteration_xvalues]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1083-1083`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1084-1084`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1073-1073`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1077-1077`
<!-- vulcan:connections:end -->

## Limitations
Planners such as RRT-Connect that have no `cost_history` produce only the placeholder trace, yet `rpo_write_planner_comparison_outputs` still writes an HTML file for them. Filtering non-finite costs before computing x values can misalign iteration numbers. Dimensions and colours are fixed.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 1072.
