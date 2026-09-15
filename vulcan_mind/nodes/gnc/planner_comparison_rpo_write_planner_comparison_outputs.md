---
id: gnc.planner_comparison_rpo_write_planner_comparison_outputs
label: rpo_write_planner_comparison_outputs
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_write_planner_comparison_outputs
  lines:
  - 1149
  - 1149
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: output_dir
  type: AbstractString
  units: n/a
  required: false
  description: Keyword argument `output_dir` (default `batch.config.output_dir`).
- id: write_plotly_outputs
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `write_plotly_outputs` (default `batch.config.write_plotly_outputs`).
- id: write_failed_path_outputs
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `write_failed_path_outputs` (default `batch.config.write_failed_path_outputs`).
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
  description: Return value of `rpo_write_planner_comparison_outputs`. Returns `(`.
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

# rpo_write_planner_comparison_outputs

## Purpose
Top-level exporter for a comparison batch: writes the flattened results CSV and, when enabled, the metric summary plot, the all-planner path family plot, per-planner path plots, per-planner cost-versus-iteration plots, and per-planner failed-path plots, returning the paths of everything written.

## Design & Implementation
Signature `rpo_write_planner_comparison_outputs(batch; output_dir = batch.config.output_dir, write_plotly_outputs = batch.config.write_plotly_outputs, write_failed_path_outputs = batch.config.write_failed_path_outputs)`. It `mkpath`s the directory, writes `rpo_planner_comparison_results.csv` via `rpo_write_namedtuple_csv(rpo_flatten_planner_results(batch))`, then if Plotly output is on saves `rpo_planner_comparison_metrics.html`, `rpo_planner_comparison_paths.html`, `rpo_planner_paths_<planner>.html`, and `rpo_planner_cost_vs_iteration_<planner>.html` for each planner. Failed-path plots are delegated to `rpo_write_failed_path_outputs`. Returns a NamedTuple with `csv`, `metrics_plot`, `path_plot`, `method_path_plots`, `cost_iteration_plots`, `failed_path_plots`, and the two Boolean flags.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `output_dir` | AbstractString | n/a | no | Keyword argument `output_dir` (default `batch.config.output_dir`). |
| in | `write_plotly_outputs` | Bool | n/a | no | Keyword argument `write_plotly_outputs` (default `batch.config.write_plotly_outputs`). |
| in | `write_failed_path_outputs` | Bool | n/a | no | Keyword argument `write_failed_path_outputs` (default `batch.config.write_failed_path_outputs`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_write_planner_comparison_outputs`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1180-1180`
- `callees` → [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1167-1167`
- `callees` → [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1170-1170`
- `callees` → [[gnc.planner_comparison_rpo_flatten_planner_results|rpo_flatten_planner_results]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1156-1156`
- `callees` → [[gnc.planner_comparison_rpo_write_failed_path_outputs|rpo_write_failed_path_outputs]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1185-1185`
- `callees` → [[gnc.planner_comparison_rpo_write_namedtuple_csv|rpo_write_namedtuple_csv]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1157-1157`
<!-- vulcan:connections:end -->

## Limitations
With Plotly output disabled `metrics_plot` and `path_plot` are `nothing` and the Dicts are empty. Files are overwritten silently and file names embed the raw planner symbol rather than `rpo_comparison_artifact_slug`. Cost-iteration HTML files are written even for planners without cost histories, producing placeholder plots. All plotting is serial, and `PlotlyJS.savefig` failures propagate as exceptions after the CSV has already been written.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 1149.
