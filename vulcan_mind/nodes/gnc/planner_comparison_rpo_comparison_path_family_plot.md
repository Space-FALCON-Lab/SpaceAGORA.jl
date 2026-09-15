---
id: gnc.planner_comparison_rpo_comparison_path_family_plot
label: rpo_comparison_path_family_plot
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_path_family_plot
  lines:
  - 864
  - 864
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: planner
  type: Any
  units: n/a
  required: false
  description: Keyword argument `planner` (default `nothing`).
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
  description: Return value of `rpo_comparison_path_family_plot`. Returns `PlotlyJS.Plot(`.
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

# rpo_comparison_path_family_plot

## Purpose
Produces a 3-D Plotly figure overlaying the station geometry with the desired (retimed reference) and tracked (LQ-MPC flown) paths plus start and goal markers for every case, optionally restricted to one planner.

## Design & Implementation
`rpo_comparison_path_family_plot(batch; planner = nothing)` seeds `traces` with `rpo_comparison_station_trace(batch)`, then for each planner in `batch.planner_types` (or the single normalised `planner`) and each plan in `batch.plans_by_planner[planner_type]`, pushes four `scatter3d` traces: a dotted `plan.tracking.r_ref_rtn` line named "<label> case <case> desired", a solid `plan.tracking.x_hist[1:3, :]` line "tracked", a green circle at `plan.case.start_rtn`, and a red diamond at `plan.case.goal_rtn`. The layout uses `aspectmode = "data"` with axis titles radial, along-track, and cross-track in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `planner` | Any | n/a | no | Keyword argument `planner` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.Plot | n/a | — | Return value of `rpo_comparison_path_family_plot`. Returns `PlotlyJS.Plot(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1170-1170`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_station_trace|rpo_comparison_station_trace]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:865-865`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:873-873`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:866-866`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:868-868`
<!-- vulcan:connections:end -->

## Limitations
Four traces per plan means a seven-planner, ten-case batch yields 281 traces and a large HTML file. Colours for the path lines are left to Plotly's default cycle, so desired and tracked lines of the same case may not share a hue. Raw planner paths are not shown here; see `rpo_comparison_single_path_plot`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 864.
