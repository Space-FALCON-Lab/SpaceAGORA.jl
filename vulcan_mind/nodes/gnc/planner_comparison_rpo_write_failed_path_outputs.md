---
id: gnc.planner_comparison_rpo_write_failed_path_outputs
label: rpo_write_failed_path_outputs
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_write_failed_path_outputs
  lines:
  - 1133
  - 1133
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
- id: output_dir
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `output_dir`.
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
  description: Return value of `rpo_write_failed_path_outputs`. Returns `failed_path_plots`.
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

# rpo_write_failed_path_outputs

## Purpose
Writes one HTML figure per planner containing all of that planner's failed cases under `<output_dir>/rpo_failed_path_plots/`, returning a map from planner symbol to file path.

## Design & Implementation
`rpo_write_failed_path_outputs(batch, output_dir::AbstractString)` iterates `batch.planner_types`, computes `failed_indices = [idx for idx in eachindex(results) if !results[idx].success]` from `batch.results_by_planner[planner]`, skips planners with no failures, lazily `mkpath`s the failed directory, and calls `PlotlyJS.savefig(rpo_comparison_failed_paths_plot(...), path)` with `path = joinpath(failed_dir, "rpo_failed_paths_<planner>.html")`. Returns a `Dict{Symbol, String}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `output_dir` | AbstractString | n/a | yes | Positional argument `output_dir`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_write_failed_path_outputs`. Returns `failed_path_plots`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1185-1185`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1142-1142`
<!-- vulcan:connections:end -->

## Limitations
Existing files are overwritten without warning. If every planner succeeds, the directory is never created and the returned Dict is empty, which callers must handle. `PlotlyJS.savefig` to HTML requires the PlotlyJS rendering stack and can be slow for batches with many failures.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 1133.
