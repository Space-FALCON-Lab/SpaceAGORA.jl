---
id: gnc.planner_comparison_rpo_comparison_failed_paths_plot
label: rpo_comparison_failed_paths_plot
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_failed_paths_plot
  lines:
  - 990
  - 990
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
- id: failed_indices
  type: Any
  units: n/a
  required: true
  description: Keyword argument `failed_indices`.
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
  description: Return value of `rpo_comparison_failed_paths_plot`. Returns `PlotlyJS.Plot(`.
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

# rpo_comparison_failed_paths_plot

## Purpose
Builds a single 3-D figure containing the desired, tracked, raw, start, and goal traces for every failed case of one planner, so all of a planner's failures can be inspected together against the station geometry.

## Design & Implementation
`rpo_comparison_failed_paths_plot(batch; planner, failed_indices)` normalises the planner, fetches `plans = batch.plans_by_planner[planner_type]`, and for each `idx in failed_indices` pushes up to five `scatter3d` traces per plan (dotted desired, solid tracked, dashed raw when `hasproperty(plan, :raw_path)`, green start circle, red goal diamond), each named with the case label. The title reports the planner label and `length(failed_indices)`. Axes are RTN metres with `aspectmode = "data"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `planner` | Any | n/a | yes | Keyword argument `planner`. |
| in | `failed_indices` | Any | n/a | yes | Keyword argument `failed_indices`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlotlyJS.Plot | n/a | — | Return value of `rpo_comparison_failed_paths_plot`. Returns `PlotlyJS.Plot(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_failed_path_outputs|rpo_write_failed_path_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1142-1142`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_station_trace|rpo_comparison_station_trace]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:993-993`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1000-1000`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:991-991`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:992-992`
<!-- vulcan:connections:end -->

## Limitations
Indices are trusted; an out-of-range index raises `BoundsError`. Per-case line colours come from Plotly's default cycle, so with many failures colours repeat. The function does no filtering itself; `rpo_write_failed_path_outputs` computes `failed_indices` from `results[idx].success`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 990.
