---
id: gnc.planner_comparison_rpo_comparison_single_path_plot
label: rpo_comparison_single_path_plot
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_single_path_plot
  lines:
  - 924
  - 924
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
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
  description: Return value of `rpo_comparison_single_path_plot`. Returns `PlotlyJS.Plot(`.
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

# rpo_comparison_single_path_plot

## Purpose
Renders one planner/case combination in 3-D with fixed colours: the desired reference (blue dotted), the tracked trajectory (orange), the raw planner path before refinement (grey dashed with markers, when present), and start/goal markers, intended for inspecting individual failed cases.

## Design & Implementation
`rpo_comparison_single_path_plot(batch, plan; planner)` normalises the planner symbol, reads `plan.case.label`, `plan.tracking.r_ref_rtn`, and `plan.tracking.x_hist[1:3, :]`, starts from the station trace, and pushes `scatter3d` traces with widths 5 (desired, tracked) and 2 (raw). The raw path is included only if `hasproperty(plan, :raw_path)`. The title is fixed to "<label> Failed Case <case> Path" with the same RTN scene settings as the family plot.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `planner` | Any | n/a | yes | Keyword argument `planner`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.Plot | n/a | — | Return value of `rpo_comparison_single_path_plot`. Returns `PlotlyJS.Plot(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_station_trace|rpo_comparison_station_trace]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:930-930`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:931-931`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:925-925`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:926-926`
<!-- vulcan:connections:end -->

## Limitations
The title always says "Failed Case" even when called on a successful plan. Not referenced by `rpo_write_planner_comparison_outputs`, which uses the multi-case `rpo_comparison_failed_paths_plot` instead; it exists for interactive use. `Matrix{Float64}(plan.raw_path)` throws if `raw_path` is not numeric.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 924.
